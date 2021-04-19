#include "FindHB.h"

using namespace std;
using namespace OpenBabel;

FindHB::FindHB()
{
    OBMessageHandler messageHandler;
    messageHandler.SetOutputLevel(obError);
    OpenBabel::obErrorLog = messageHandler;
}

int FindHB::find_atom(string atomname, Mol2* Rec, int astart, int aend){
    int index=-1;
    for (int i=astart; i<=aend; i++){
        if (Rec->atomnames[i] == atomname){
            index=i;
        }
    }
    return index;
}

double FindHB::distance(vector<double> xyz1, vector<double> xyz2){
    double d = ( (xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0])) + ((xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1])) + ((xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]));
    d = sqrt(d);
    return d;
}

double FindHB::distance_squared(vector<double> xyz1, vector<double> xyz2){
    double d = ( (xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0])) + ((xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1])) + ((xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]));
    return d;
}

double FindHB::angle(vector<double> xyz1, vector<double> xyz2, vector<double> xyz3){
    //
    // Law of cosines: c** = a** + b** - 2ab cos(C)
    //
    double ab = distance(xyz1, xyz2);
    double ac = distance(xyz1, xyz3);
    double bc = distance(xyz2, xyz3);
    double angle = acos(((ab*ab)+(bc*bc)-(ac*ac))/(2*ab*bc));
    angle = angle * 180.0 / PI;
    return (angle);
}

void FindHB::find_ligandHB(string molfile, Mol2* Lig){
    OBMol mol;
    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
        printf("Could not find input format for file\n");
    }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading. Maybe an OpenBabel issue?\n", molfile.c_str());
    }
    if (!conv.Read(&mol, &ifs)) {
        printf("Could not read molecule from file. Maybe an OpenBabel issue?\n");
    }

    OBAtom *atom;
    vector<int> vtemp(2);

    FOR_ATOMS_OF_MOL(atom, mol){
        if (atom->IsHbondAcceptor()){
            Lig->HBacceptors.push_back(int(atom->GetIdx()-1));
        }
        else if (atom->IsHbondDonorH()){
            FOR_NBORS_OF_ATOM(nbr, &*atom){
                vtemp[0] = int(nbr->GetIdx()-1);
                vtemp[1] = int(atom->GetIdx()-1);
                Lig->HBdonors.push_back(vtemp);
            }
        }
    }
}

void FindHB::parse_residue(int atom_start, int atom_end, string resname, Mol2* Rec, Mol2* Lig, double dist_cutoff){
    vector<int> vtemp(2);
    int res = this->find_atom(string("CA"), Rec, atom_start, atom_end);
    if (res > 0){
        COORD_MC* Coord = new COORD_MC;
        vector<double> com = Coord->compute_com(Lig);
        double d = distance(Rec->xyz[res], com);        // distance from residue alpha-carbon and ligand center of mass
        delete Coord;

        if (d < dist_cutoff){

            map<string, int> aa =  {
                {"ALA", 0},
                {"ARG", 1},
                {"ASN", 2},
                {"ASP", 3},
                {"CYS", 4},
                {"GLN", 5},
                {"GLU", 6},
                {"GLY", 7},
                {"HIS", 8},
                {"ILE", 9},
                {"LEU", 10},
                {"LYS", 11},
                {"MET", 12},
                {"PHE", 13},
                {"PRO", 14},
                {"SER", 15},
                {"THR", 16},
                {"TRP", 17},
                {"TYR", 18},
                {"VAL", 19},
                {"ASH", 20},
                {"CYX", 21},
                {"GLH", 22},
                {"HIE", 23},
                {"HID", 24},
                {"HIP", 25}
            };
            res = aa[resname];
            bool hasHD1 = false;
            bool hasHE2 = false;

            switch (res) {
            case 1:             // ARG
                vtemp[0] = find_atom(string("NH1"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HH11"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                vtemp[0] = find_atom(string("NH1"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HH12"),Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                vtemp[0] = find_atom(string("NH2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HH21"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                vtemp[0] = find_atom(string("NH2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HH22"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
                break;

            case 2:             //ASN
                Rec->HBacceptors.push_back(find_atom(string("OD1"), Rec, atom_start, atom_end));

                vtemp[0] = find_atom(string("ND2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HD21"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                vtemp[0] = find_atom(string("ND2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HD22"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
                break;

            case 3:             // ASP
                Rec->HBacceptors.push_back(find_atom(string("OD1"), Rec, atom_start, atom_end));
                Rec->HBacceptors.push_back(find_atom(string("OD2"), Rec, atom_start, atom_end));
                break;

            case 5:             // GLN
                Rec->HBacceptors.push_back(find_atom(string("OE1"), Rec, atom_start, atom_end));

                vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HE21"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HE22"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
                break;

            case 6:             // GLU
                Rec->HBacceptors.push_back(find_atom(string("OE1"), Rec, atom_start, atom_end));
                Rec->HBacceptors.push_back(find_atom(string("OE2"), Rec, atom_start, atom_end));
                break;

            case 8:             // HIS
                if (find_atom(string("HD1"), Rec, atom_start, atom_end) > 0){
                    hasHD1 = true;
                }
                if (find_atom(string("HE2"), Rec, atom_end, atom_end) > 0){
                    hasHE2 = true;
                }

                if (hasHD1 and hasHE2){     // HIP
                    vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
                    vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
                    Rec->HBdonors.push_back(vtemp);

                    vtemp[0] = find_atom(string("ND1"), Rec, atom_start, atom_end);
                    vtemp[1] = find_atom(string("HD1"), Rec, atom_start, atom_end);
                    Rec->HBdonors.push_back(vtemp);
                }
                else if (hasHD1){           // HID
                    vtemp[0] = find_atom(string("ND1"), Rec, atom_start, atom_end);
                    vtemp[1] = find_atom(string("HD1"), Rec, atom_start, atom_end);
                    Rec->HBdonors.push_back(vtemp);
                }
                else if (hasHE2){           // HIE
                    vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
                    vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
                    Rec->HBdonors.push_back(vtemp);
                }
                break;

            case 11:                        // LYS
                vtemp[0] = find_atom(string("NZ"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HZ1"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                vtemp[0] = find_atom(string("NZ"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HZ2"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                vtemp[0] = find_atom(string("NZ"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HZ3"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
                break;

            case 15:                        // SER
                vtemp[0] = find_atom(string("OG"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HG"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                Rec->HBacceptors.push_back(find_atom(string("OG"), Rec, atom_start, atom_end)); // ???
                break;

            case 16:                        // THR
                vtemp[0] = find_atom(string("OG1"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HG1"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                Rec->HBacceptors.push_back(find_atom(string("OG1"), Rec, atom_start, atom_end)); // ???
                break;

            case 18:                        // TYR
                vtemp[0] = find_atom(string("OH"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HH"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                Rec->HBacceptors.push_back(find_atom(string("OH"), Rec, atom_start, atom_end)); // ???
                break;

            case 20:                        // ASH
                Rec->HBacceptors.push_back(find_atom(string("OD1"), Rec, atom_start, atom_end));

                vtemp[0] = find_atom(string("OD2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HD2"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                Rec->HBacceptors.push_back(find_atom(string("OD2"), Rec, atom_start, atom_end));
                break;

            case 22:                        // GLH
                Rec->HBacceptors.push_back(find_atom(string("OE1"), Rec, atom_start, atom_end));

                vtemp[0] = find_atom(string("OE2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                Rec->HBacceptors.push_back(find_atom(string("OE2"), Rec, atom_start, atom_end));
                break;

            case 23:                        // HIE
                vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
                break;

            case 24:                        // HID
                vtemp[0] = find_atom(string("ND1"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HD1"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
                break;

            case 25:                        // HIP
                vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);

                vtemp[0] = find_atom(string("ND1"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HD1"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
                break;

            }

            // main chain atoms

            Rec->HBacceptors.push_back(find_atom(string("O"), Rec, atom_start, atom_end));

            if (resname != "PRO"){
                vtemp[0] = find_atom(string("N"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("H"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
            }
        }
    }
    else if (resname == "HOH"){
        int res = find_atom(string("O"), Rec, atom_start, atom_end);
        if (res > 0){
            COORD_MC* Coord = new COORD_MC;
            vector<double> com = Coord->compute_com(Lig);
            double d = distance(Rec->xyz[res], com);        // distance from residue alpha-carbon and ligand center of mass
            delete Coord;
            if (d < dist_cutoff){
                Rec->HBacceptors.push_back(find_atom(string("O"), Rec, atom_start, atom_end));
                vtemp[0] = find_atom(string("O"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("H1"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
                vtemp[0] = find_atom(string("O"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("H2"), Rec, atom_start, atom_end);
                Rec->HBdonors.push_back(vtemp);
            }
        }
    }
}

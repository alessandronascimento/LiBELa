#include <iostream>
#include <cstdio>
#include <map>
#include <openbabel/atom.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <openbabel/math/align.h>
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/PARSER.cpp"
#include "../LiBELa/COORD_MC.cpp"

using namespace std;
using namespace OpenBabel;

struct HBond {
    vector<int> acceptors;
    vector<vector<int> > donors;
};

int find_atom(string atomname, Mol2* Rec, int astart, int aend){
    int index=-1;
    for (int i=astart; i<=aend; i++){
        if (Rec->atomnames[i] == atomname){
            index=i;
        }
    }
    return index;
}


void usage(){
    printf("Usage: Find_HB -r <rec.mol> -l <lig.mol2> \n");
    printf("\t -l <Ligand.mol> \n");
    printf("\t -r <Receptor.mol2> \n");
}

double distance(vector<double> xyz1, vector<double> xyz2){
    double d = ( (xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0])) + ((xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1])) + ((xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]));
    d = sqrt(d);
    return d;
}


void find_ligandHB(string molfile, HBond* HB){
    OBMol mol;
    OBConversion* conv = new OBConversion;
    OBFormat *format = conv->FormatFromExt(molfile.c_str());
    if (!format || !conv->SetInFormat(format)) {
        printf("Could not find input format for file\n");
    }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading. Maybe an OpenBabel issue?\n", molfile.c_str());
    }
    // Read the molecule.
    if (!conv->Read(&mol, &ifs)) {
        printf("Could not read molecule from file. Maybe an OpenBabel issue?\n");
    }
    delete conv;

    OBAtom *atom;
    int index = -1;

    vector<int> vtemp (2);

    printf("PArsing atoms... \n");

    FOR_ATOMS_OF_MOL(atom, mol){
        if (atom->IsHbondAcceptor()){
            index = atom->GetIdx()-1;
            HB->acceptors.push_back(index);
        }
        else if (atom->IsHbondDonorH()){
            FOR_NBORS_OF_ATOM(nbr, &*atom){
                vtemp[0] = nbr->GetIdx()-1;
                vtemp[1] = atom->GetIdx()-1;
                HB->donors.push_back(vtemp);

            }
        }
    }
}


void parse_residue(int atom_start, int atom_end, string resname, Mol2* Rec, Mol2* Lig, HBond* HB, double dist_cutoff=10.){
    int res = find_atom(string("CA"), Rec, atom_start, atom_end);
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
        vector<int> vtemp(2);
        bool hasHD1 = false;
        bool hasHE2 = false;

        switch (res) {
        case 1:             // ARG
            vtemp[0] = find_atom(string("NH1"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HH11"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            vtemp[0] = find_atom(string("NH1"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HH12"),Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            vtemp[0] = find_atom(string("NH2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HH21"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            vtemp[0] = find_atom(string("NH2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HH22"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);
            break;

        case 2:             //ASN
            HB->acceptors.push_back(find_atom(string("OD1"), Rec, atom_start, atom_end));

            vtemp[0] = find_atom(string("ND2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HD21"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            vtemp[0] = find_atom(string("ND2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HD22"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);
            break;

        case 3:             // ASP
            HB->acceptors.push_back(find_atom(string("OD1"), Rec, atom_start, atom_end));
            HB->acceptors.push_back(find_atom(string("OD2"), Rec, atom_start, atom_end));
            break;

        case 5:             // GLN
            HB->acceptors.push_back(find_atom(string("OE1"), Rec, atom_start, atom_end));

            vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HE21"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HE22"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);
            break;

        case 6:             // GLU
            HB->acceptors.push_back(find_atom(string("OE1"), Rec, atom_start, atom_end));
            HB->acceptors.push_back(find_atom(string("OE2"), Rec, atom_start, atom_end));
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
                HB->donors.push_back(vtemp);

                vtemp[0] = find_atom(string("ND1"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HD1"), Rec, atom_start, atom_end);
                HB->donors.push_back(vtemp);
            }
            else if (hasHD1){           // HID
                vtemp[0] = find_atom(string("ND1"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HD1"), Rec, atom_start, atom_end);
                HB->donors.push_back(vtemp);
            }
            else if (hasHE2){           // HIE
                vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
                vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
                HB->donors.push_back(vtemp);
            }
            break;

        case 11:                        // LYS
            vtemp[0] = find_atom(string("NZ"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HZ1"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            vtemp[0] = find_atom(string("NZ"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HZ2"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            vtemp[0] = find_atom(string("NZ"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HZ3"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);
            break;

        case 15:                        // SER
            vtemp[0] = find_atom(string("OG"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HG"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            HB->acceptors.push_back(find_atom(string("OG"), Rec, atom_start, atom_end)); // ???
            break;

        case 16:                        // THR
            vtemp[0] = find_atom(string("OG1"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HG1"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            HB->acceptors.push_back(find_atom(string("OG1"), Rec, atom_start, atom_end)); // ???
            break;

        case 18:                        // TYR
            vtemp[0] = find_atom(string("OH"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HH"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            HB->acceptors.push_back(find_atom(string("OH"), Rec, atom_start, atom_end)); // ???
            break;

        case 20:                        // ASH
            HB->acceptors.push_back(find_atom(string("OD1"), Rec, atom_start, atom_end));

            vtemp[0] = find_atom(string("OD2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HD2"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            HB->acceptors.push_back(find_atom(string("OD2"), Rec, atom_start, atom_end));
            break;

        case 22:                        // GLH
            HB->acceptors.push_back(find_atom(string("OE1"), Rec, atom_start, atom_end));

            vtemp[0] = find_atom(string("OE2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            HB->acceptors.push_back(find_atom(string("OE2"), Rec, atom_start, atom_end));
            break;

        case 23:                        // HIE
            vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);
            break;

        case 24:                        // HID
            vtemp[0] = find_atom(string("ND1"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HD1"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);
            break;

        case 25:                        // HIP
            vtemp[0] = find_atom(string("NE2"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HE2"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);

            vtemp[0] = find_atom(string("ND1"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("HD1"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);
            break;

        }

        // main chain atoms

        HB->acceptors.push_back(find_atom(string("O"), Rec, atom_start, atom_end));

        if (resname != "PRO"){
            vtemp[0] = find_atom(string("N"), Rec, atom_start, atom_end);
            vtemp[1] = find_atom(string("H"), Rec, atom_start, atom_end);
            HB->donors.push_back(vtemp);
        }
    }
}

int main(int argc, char* argv[]){

    int c;
    char *ligfile;
    char *recfile;
    double cutoff=10.0;
    PARSER* Input = new PARSER;
    vector<int> hb_donors;
    vector<int> hb_acceptors;

    if (argc < 2){
        usage();
        exit(1);
    }

    while ((c = getopt (argc, argv, "r:l:d:")) != -1)
        switch (c)
        {
        case 'r':
            recfile = optarg;
            Input->rec_mol2 = string(recfile);
            break;
        case 'l':
            ligfile = optarg;
            Input->reflig_mol2 = string(ligfile);
            Input->lig_mol2 = string(ligfile);
            break;
        case 'd':
            cutoff = double(atof(optarg));
            break;
        case '?':
            if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt)){
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                usage();
            }
            else{
                fprintf (stderr, "Unknown option character `\\x%x'.\n",optopt);
                usage();
            }
            return 1;
            break;
        default:
            usage();
        }


    Mol2* Rec = new Mol2(Input, Input->rec_mol2);
    Mol2* Lig = new Mol2(Input, Input->lig_mol2);

    HBond* HB = new HBond;
    for (int i=0; i< Rec->residue_pointer.size()-1; i++){
        parse_residue(Rec->residue_pointer[i]-1, Rec->residue_pointer[i+1]-2, Rec->resnames[i], Rec, Lig, HB, cutoff);
    }

    printf("Number of acceptors found: %7lu\n", HB->acceptors.size());
    printf("Number of donors found: %7lu\n", HB->donors.size());

    HBond* HBLig = new HBond;
    find_ligandHB(Input->lig_mol2, HBLig);
    printf("Number of acceptors found: %7lu\n", HBLig->acceptors.size());
    printf("Number of donors found: %7lu\n", HBLig->donors.size());

    return  0;

}

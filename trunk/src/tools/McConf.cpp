#include <iostream>
#include <cstdlib>
#include <vector>
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/Conformer.cpp"
#include "../LiBELa/PARSER.cpp"
#include "../LiBELa/WRITER.cpp"
#include "../LiBELa/COORD_MC.cpp"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/rotor.h>
#include <openbabel/conformersearch.h>
#include <openbabel/shared_ptr.h>
#include <openbabel/forcefield.h>
#include <openbabel/math/vector3.h>


using namespace std;


double* copy_to_obmol(vector<vector<double> > vec_xyz){
    double *myxyz = new double[vec_xyz.size()*3];
    for (unsigned i=0; i<vec_xyz.size(); i++){
        myxyz[3*i] = vec_xyz[i][0];
        myxyz[(3*i)+1] = vec_xyz[i][1];
        myxyz[(3*i)+2] = vec_xyz[i][2];
    }
    return myxyz;
}

vector<vector<double> > copy_from_obmol(shared_ptr<OBMol> mymol){
    vector<vector<double > > vec_xyz;
    vector<double> tmp(3);
    double *myxyz = new double[mymol->NumAtoms()*3];
    myxyz = mymol->GetCoordinates();
    for (unsigned i=0; i < mymol->NumAtoms(); i++){
        tmp[0] = (myxyz[3*i]);
        tmp[1] = (myxyz[(3*i)+1]);
        tmp[2] = (myxyz[(3*i)+2]);
        vec_xyz.push_back(tmp);
    }

    tmp.clear();
//    delete[] myxyz;
    return vec_xyz;
}

shared_ptr<OBMol> GetMol(const std::string &molfile){
    shared_ptr<OBMol> mol(new OBMol);

    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
    std::cout << "Could not find input format for file " << molfile << endl;
    return mol;
  }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        std::cout << "Could not open " << molfile << " for reading." << endl;
        return mol;
    }

    if (!conv.Read(mol.get(), &ifs)) {
        std::cout << "Could not read molecule from file " << molfile << endl;
        return mol;
    }
    return mol;
}

int main(int argc, char* argv[]) {

    printf("****************************************************************************************************\n");
    printf("****************************************************************************************************\n");
    printf("*                            McConf - Conformer Generation for McLiBELa                            *\n");
    printf("*                                                                                                  *\n");
    printf("* University of Sao Paulo                                                                          *\n");
    printf("* More Info:                                                                                       *\n");
    printf("*      http://www.biotechmol.ifsc.usp.br/                                                          *\n");
    printf("*                                                                                                  *\n");
    printf("****************************************************************************************************\n");
    printf("****************************************************************************************************\n");


    if (argc < 2){
            printf("Usage: %s file.mol2 [input_file]\n", argv[0]);
            exit(1);
    }

    PARSER* Input = new PARSER;

    if (argc > 2){
        Input->set_parameters(argv[2]);
    }

    Input->write_mol2 = true;
    Input->dock_mode = true;

    Mol2* Mol = new Mol2(Input, string(argv[1]));

    WRITER* Writer = new WRITER("McConf", Input);

    shared_ptr<OBMol> mol;
    OBForceField* OBff;
    OBRotorList RotorList;
    OBRotorIterator RotorIterator;
    OBRotor *Rotor;
    vector<vector<int> > atoms_in_dihedrals;
    mol = GetMol(string(argv[1]));
    OBff = OBForceField::FindForceField("GAFF");


    if (!OBff){
        cout << "Could not find FF GAFF!" << endl;
        exit(1);
    }

    OBff->Setup(*mol);
    OBff->GetCoordinates(*mol);
    OBff->SteepestDescent(100, 1.0e-10);
    OBff->GetCoordinates(*mol);
    Mol->xyz = copy_from_obmol(mol);
    RotorList.Setup(*mol);
    Rotor = RotorList.BeginRotor(RotorIterator);
    mol->ToInertialFrame();


    vector<int> tmp(4);
    char info[98];
    sprintf(info, "Found %d rotatable bonds in ligand %s.", int(RotorList.Size()), Mol->molname.c_str());
    Writer->print_info(info);
    Writer->print_line();

    for (unsigned i = 0; i < RotorList.Size() ; ++i, Rotor = RotorList.NextRotor(RotorIterator)) {
        tmp = Rotor->GetDihedralAtoms();
        printf("%d %d %d %d\n", tmp[0], tmp[1], tmp[2], tmp[3]);
        atoms_in_dihedrals.push_back(tmp);
        tmp.clear();
    }

// End of file parsing

    double* myxyz = new double[mol->NumAtoms()*3];

    myxyz = copy_to_obmol(Mol->xyz);
    mol->SetCoordinates(myxyz);
    vector<vector<double> > new_xyz;

    double ene;

    double current_angle, new_angle;
    for (unsigned i=0; i< RotorList.Size(); i++){
        for (int j=0; j<360; j=j+10){
            current_angle = mol->GetTorsion(mol->GetAtom(atoms_in_dihedrals[i][0]), mol->GetAtom(atoms_in_dihedrals[i][1]), mol->GetAtom(atoms_in_dihedrals[i][2]), mol->GetAtom(atoms_in_dihedrals[i][3]));
            new_angle = current_angle + (j*1.0) ;
            mol->SetTorsion(mol->GetAtom(atoms_in_dihedrals[i][0]), mol->GetAtom(atoms_in_dihedrals[i][1]), mol->GetAtom(atoms_in_dihedrals[i][2]), mol->GetAtom(atoms_in_dihedrals[i][3]), new_angle*PI/180.);
            OBff->Setup(*mol);
            ene = OBff->Energy();
            new_xyz = copy_from_obmol(mol);
            Writer->writeMol2(Mol, new_xyz, ene, 0.0);
        }
    }

    delete Writer;
    delete Mol;
    delete Input;

    return 0;
}

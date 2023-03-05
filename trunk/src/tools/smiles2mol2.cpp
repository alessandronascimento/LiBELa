#include <iostream>
#include <string>
#include <ctime>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>
#include <openbabel/builder.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>

using namespace std;
using namespace OpenBabel;

int main(int argc,char **argv){
    OBMol mol;
    string smiles_input = string(argv[1]);
    OBConversion conv;
    if (!(conv.SetInFormat("smi") && conv.ReadString(&mol, smiles_input))){
        printf("Could not read SMILES from input. Please, check! \n");
    }

    int ti = clock();

    OBBuilder builder;
    builder.Build(mol);

    mol.AddHydrogens(false, true); // adding H atoms: false= not for polar atoms only. true=correct for pH 7.4.

    OBForceField* OBff = OBForceField::FindForceField("MMFF94");
    OBff->Setup(mol);
    OBff->SteepestDescent(100);
//    OBff->SteepestDescent(100, 1.0e-4);
//    OBff->WeightedRotorSearch(250, 50);
//    OBff->SteepestDescent(100, 1.0e-4);

    OBff->UpdateCoordinates(mol);

    int tf = clock()-ti;
    char aname[4];


    FOR_ATOMS_OF_MOL(atom, mol){
        sprintf(aname, "%s%d", OBElements::GetSymbol(atom->GetAtomicNum()), atom->GetIdx());
        printf("Atom name: %s\n", aname);
        printf("Atom type: %s\n", atom->GetType());
        printf("Atom charge: %f\n\n", atom->GetPartialCharge());
    }
/*
    vector<double> charges;
    vector<string> atomnames;
*/
    conv.SetOutFormat("mol2");
    conv.Write(&mol, &cout);

    printf("3D generation computations took %f second(s)", float((tf*1.0)/CLOCKS_PER_SEC));

    return 0;
}

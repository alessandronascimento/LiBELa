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
#include <openbabel/forcefield.h>
#include <openbabel/math/vector3.h>


using namespace std;

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
        printf("Usage: %s file.mol2 [input_file] [Ntrials=1000]\n", argv[0]);
        exit(1);
    }

    PARSER* Input = new PARSER;
    Mol2* Mol = new Mol2(Input, string(argv[1]));
    WRITER* Writer = new WRITER("McConf", Input);
    double confrmsd = 0.5;
    double confene = 50.0;
    int ntrial=1000;
    if (argc > 2){
        ntrial = atoi(argv[2]);
    }

    Input->write_mol2 = true;
    Input->dock_mode = true;


    bool file_read;
    OBMol* mol = new OBMol;
    OBConversion* conv = new OBConversion;
    OBFormat *format = conv->FormatFromExt(argv[1]);
    if (!format || !conv->SetInFormat(format)) {
        printf("Could not find input format for file\n");
        exit(1);
    }

    ifstream ifs(argv[1]);
    if (!ifs) {
        printf("Could not open %s for reading\n", argv[1]);
        exit(1);
    }

    if (!conv->Read(mol, &ifs)) {
        printf("Could not read molecule from file\n");
        exit(1);
    }

    OBMol* ref_mol = new OBMol;
    ifstream ifs2(argv[1]);
    conv->Read(ref_mol, &ifs2);

    delete conv;
    format->~OBFormat();

    OBForceField* OBff = OBForceField::FindForceField("GAFF");

    if (!OBff){
        printf("Could not find OpenBabel FF parameters!\n");
        exit(1);
    }

    if (Input->verbose){
        OBff->SetLogFile(&cout);
        OBff->SetLogLevel(OBFF_LOGLVL_LOW);
    }

    // Original conformation energy
    OBff->Setup(*mol);
    mol->SetTotalCharge(mol->GetTotalCharge());
    double energy = OBff->Energy();
    if (OBff->GetUnit() == "kJ/mol"){       // Converting to kcal/mol, if needed.
        energy = energy/4.18;
    }

    printf("Original energy for molecule %s: %10.4f\n", argv[1], energy);

    // Energy minimization
    OBff->GetCoordinates(*mol);
    OBff->SteepestDescent(1000);
    OBff->GetCoordinates(*mol);
    energy = OBff->Energy();
    if (OBff->GetUnit() == "kJ/mol"){       // Converting to kcal/mol, if needed.
        energy = energy/4.18;
    }
    printf("Energy for molecule after minimization %s: %10.4f\n", argv[1], energy);


    // Conformer Search
    OBff->DiverseConfGen(confrmsd, ntrial, confene, false);
    OBff->GetConformers(*mol);

    double rmsd;
    int generated_conformers=0;
    if (mol->NumConformers() > Input->lig_conformers){
        generated_conformers = Input->lig_conformers;
    }
    else {
        generated_conformers = mol->NumConformers();
    }

    if (mol->NumConformers() > 0){
        for (int i=0; i<generated_conformers; i++){
            double *xyz = new double [mol->NumAtoms()*3];
            vector<double> v3;
            vector<vector<double> > xyz_tmp;
            mol->SetConformer(i);
            OBff->Setup(*mol);
            OBff->GetCoordinates(*mol);
            energy = OBff->Energy();
            if (OBff->GetUnit() == "kJ/mol"){       // Converting to kcal/mol, if needed.
                energy = energy/4.18;
            }

            OBAlign* align = new OBAlign;
            align->SetRefMol(*ref_mol);
            align->SetTargetMol(*mol);
            align->Align();
            rmsd = align->GetRMSD();
            align->UpdateCoords(mol);

            printf("Energy/RMSD for conformer [%3d]: %10.4f kcal/mol / %10.4f Ang\n", i, energy, rmsd);

            delete align;

            xyz = mol->GetCoordinates();
            for (unsigned j=0; j<mol->NumAtoms(); j++){
                v3.push_back(xyz[3*j]);
                v3.push_back(xyz[(3*j)+1]);
                v3.push_back(xyz[(3*j)+2]);
                xyz_tmp.push_back(v3);
                v3.clear();
            }

            Writer->writeMol2(Mol,xyz_tmp, energy, rmsd, "McConf");
            xyz_tmp.clear();
            delete[] xyz;
        }
        file_read = true;
    }
    else {
        file_read = false;
    }

    printf("Number of conformers: %4d\n", generated_conformers);

    delete ref_mol;

    delete Writer;
    delete Mol;
    delete Input;

    return (file_read);
}

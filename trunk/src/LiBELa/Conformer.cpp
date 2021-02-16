/*
 * Conformer.cpp
 *
 *  Created on: 02/08/2012
 *      Author: Nascimento
 */

#include "Conformer.h"

using namespace std;
using namespace OpenBabel;

Conformer::Conformer() {

// Supressing OpenBabel warning messages
messageHandler = new OBMessageHandler;
messageHandler->SetOutputLevel(obError);
OpenBabel::obErrorLog = *messageHandler;
}

Conformer::~Conformer() {
    delete messageHandler;
}

OBMol Conformer::GetMol(const std::string &molfile){
    // Create the OBMol object.
    OBMol mol;

    // Create the OBConversion object.
    OBConversion* conv = new OBConversion;
    OBFormat *format = conv->FormatFromExt(molfile.c_str());
    if (!format || !conv->SetInFormat(format)) {
    printf("Could not find input format for file\n");
    return mol;
  }

    // Open the file.
    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading. Maybe an OpenBabel issue?\n", molfile.c_str());
        return mol;
    }
    // Read the molecule.
    if (!conv->Read(&mol, &ifs)) {
        printf("Could not read molecule from file. Maybe an OpenBabel issue?\n");
        return mol;
    }

    delete conv;

    return mol;
}

bool Conformer::generate_conformers_confab(PARSER* Input, Mol2* Lig, string molfile){
    bool file_read;
    OBMol* mol = new OBMol;
    OBConversion* conv = new OBConversion;

    OBFormat *format = conv->FormatFromExt(molfile.c_str());
    if (!format || !conv->SetInFormat(format)) {
    printf("Could not find input format for file\n");
    exit(1);
  }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading. Maybe an OpenBabel issue?\n", molfile.c_str());
        exit(1);
    }

    if (!conv->Read(mol, &ifs)) {
        printf("Could not read molecule from file. Maybe an OpenBabel issue?\n");
        exit(1);
    }

    OBMol* ref_mol = new OBMol;
    ifstream ifs2(molfile.c_str());
    conv->Read(ref_mol, &ifs2);

    delete conv;

    OBForceField* OBff;
    if (Input->ligand_energy_model == "GAFF" or Input->ligand_energy_model == "gaff"){
        OBff = OBForceField::FindForceField("GAFF");
    }
    else {
        OBff = OBForceField::FindForceField("MMFF94");
    }

#ifdef DEBUG
    OBff->SetLogFile(&cout);
    OBff->SetLogLevel(OBFF_LOGLVL_LOW);
#endif

    if (!OBff){
        printf("Could not find OpenBabel FF parameters!\n");
        exit(1);
    }

// Original conformation energy
    OBff->Setup(*mol);
    mol->SetTotalCharge(mol->GetTotalCharge());
    double energy = OBff->Energy();
    if (OBff->GetUnit() == "kJ/mol"){       // Converting to kcal/mol, if needed.
        energy = energy/4.18;
    }

// Conformer Search
    OBff->DiverseConfGen(0.5, Input->conf_search_trials, 50.0, false);
    OBff->GetConformers(*mol);

    int generated_conformers=0;
    if (mol->NumConformers() > Input->lig_conformers){
        generated_conformers = Input->lig_conformers;
    }
    else {
        generated_conformers = mol->NumConformers();
    }

    if (mol->NumConformers() > 0){
        for (int i=0; i<generated_conformers; i++){
            double* xyz = new double[mol->NumAtoms()*3];
            vector<double> v3;
            vector<vector<double> > xyz_tmp;
            mol->SetConformer(i);
            OBff->Setup(*mol);
            OBff->GetCoordinates(*mol);
            energy = OBff->Energy();
            if (OBff->GetUnit() == "kJ/mol"){       // Converting to kcal/mol, if needed.
                energy = energy/4.18;
            }
            Lig->conformer_energies.push_back(energy);

            OBAlign* align = new OBAlign;
            align->SetRefMol(*ref_mol);
            align->SetTargetMol(*mol);
            align->Align();
            align->UpdateCoords(mol);
            delete align;

            xyz = mol->GetCoordinates();
            for (unsigned j=0; j<mol->NumAtoms(); j++){
                v3.push_back(xyz[3*j]);
                v3.push_back(xyz[(3*j)+1]);
                v3.push_back(xyz[(3*j)+2]);
                xyz_tmp.push_back(v3);
                v3.clear();
            }
            Lig->mcoords.push_back(xyz_tmp);
            xyz_tmp.clear();

            delete[] xyz;
        }
        file_read = true;
    }
    else {
        file_read = false;
    }
    delete ref_mol;
    return file_read;
 }

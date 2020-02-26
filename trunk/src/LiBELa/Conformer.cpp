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
#ifndef DEBUG
    OBMessageHandler* messageHandler = new OBMessageHandler;
	messageHandler->SetOutputLevel(obError);
	OpenBabel::obErrorLog = *messageHandler;
#endif


}

Conformer::~Conformer() {
}

OBMol Conformer::GetMol(const std::string &molfile){
    // Create the OBMol object.
    OBMol mol;

    // Create the OBConversion object.
    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
    printf("Could not find input format for file\n");
    return mol;
  }

    // Open the file.
    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading\n", molfile.c_str());
        return mol;
    }
    // Read the molecule.
    if (!conv.Read(&mol, &ifs)) {
        printf("Could not read molecule from file\n");
        return mol;
    }
    return mol;
}

bool Conformer::generate_conformers_confab(PARSER* Input, Mol2* Lig, string molfile){
    bool file_read;
    OBMol* mol = new OBMol;
    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
    printf("Could not find input format for file\n");
    exit(1);
  }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading\n", molfile.c_str());
        exit(1);
    }

    if (!conv.Read(mol, &ifs)) {
        printf("Could not read molecule from file\n");
        exit(1);
    }

    OBMol* ref_mol = new OBMol;
    ref_mol = mol;

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
    Lig->conformer_energies.push_back(energy);
    Lig->mcoords.push_back(Lig->xyz);

// Conformer Search
    OBff->DiverseConfGen(0.5, Input->lig_conformers, 50.0, false);
    OBff->GetConformers(*mol);

    if (mol->NumConformers() > 0){
        for (int i=0; i<mol->NumConformers(); i++){
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
    delete mol;
    return file_read;
 }

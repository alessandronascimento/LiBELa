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
    OBMessageHandler* messageHandler;
	messageHandler = new OBMessageHandler;
	messageHandler->SetOutputLevel(obError);
	OpenBabel::obErrorLog = *messageHandler;
#endif


}

Conformer::~Conformer() {
}

shared_ptr<OBMol> Conformer::GetMol(const std::string &molfile){
    // Create the OBMol object.
    shared_ptr<OBMol> mol(new OBMol);

    // Create the OBConversion object.
    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
    std::cout << "Could not find input format for file " << molfile << endl;
    return mol;
  }

    // Open the file.
    ifstream ifs(molfile.c_str());
    if (!ifs) {
        std::cout << "Could not open " << molfile << " for reading." << endl;
        return mol;
    }
    // Read the molecule.
    if (!conv.Read(mol.get(), &ifs)) {
        std::cout << "Could not read molecule from file " << molfile << endl;
        return mol;
    }
    return mol;
}


bool Conformer::generate_conformers_GA(PARSER* Input, Mol2* Lig, string molfile){
    bool file_read;
    double charge;
    shared_ptr<OBMol> mol = this->GetMol(molfile);
    shared_ptr<OBMol> ref_mol = mol;
    OBConformerSearch cs;
    OBForceField* OBff;
    OBEnergyConformerScore* EneScore = new OBEnergyConformerScore;
    if (Input->ligand_energy_model == "GAFF"){
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
        cout << "Could not find OpenBabel FF parameters!" << endl;
        exit(1);
    }

// Original conformation energy
    OBff->Setup(*mol);
    mol->SetTotalCharge(mol->GetTotalCharge());
    Lig->conformer_energies.push_back(OBff->Energy());
    Lig->mcoords.push_back(Lig->xyz);

// Conformer Search
    cs.Setup(*mol.get(),
             Input->lig_conformers, // numConformers
             5, // numChildren=5
             5, // mutability=5
             25); // convergence=25
    cs.SetScore(EneScore);
    cs.Search();
    cs.GetConformers(*mol.get());

    if (mol->NumConformers() > 0){
        for (int i=0; i<mol->NumConformers(); i++){
            double* xyz = new double[mol->NumAtoms()*3];
            vector<double> v3;
            vector<vector<double> > xyz_tmp;
            mol->SetConformer(i);

            OBff->Setup(*mol);

            OBff->GetCoordinates(*mol);
            OBff->ConjugateGradients(Input->conformer_min_steps);
            OBff->GetCoordinates(*mol);
            Lig->conformer_energies.push_back(OBff->Energy());

            OBAlign* align = new OBAlign;
            align->SetRefMol(*ref_mol.get());
            align->SetTargetMol(*mol.get());
            align->Align();
            align->UpdateCoords(mol.get());
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
        }
        file_read = true;
    }
    else {
        file_read = false;
    }
    return file_read;
 }

bool Conformer::generate_conformers_WRS(PARSER* Input, Mol2* Lig, string molfile){
    OBMol* OBmol = new OBMol;
    OBConversion* OBconv = new OBConversion;
    OBForceField* OBff;
    OBFormat* format;

    if (Input->ligand_energy_model == "GAFF"){
        OBff = OBForceField::FindForceField("GAFF");
    }
    else {
        OBff = OBForceField::FindForceField("MMFF94");
    }
    if (!OBff){
        cout << "Could not find OpenBabel FF parameters!" << endl;
        exit(1);
    }

#ifdef DEBUG
    OBff->SetLogFile(&cout);
    OBff->SetLogLevel(OBFF_LOGLVL_LOW);
#endif

    format = OBconv->FormatFromExt(molfile.c_str());

    if (!format || !OBconv->SetInFormat(format)){
        printf("Conformer Class could not open file %s. Is this file format acceptable?\n", molfile.c_str());
        exit(1);
    }
    ifstream input_file(molfile.c_str());
    bool file_read = input_file.is_open();

    if (file_read){
        OBconv->Read(OBmol, &input_file);
        double total=0.00;
        FOR_ATOMS_OF_MOL(atom, OBmol){
            total+=atom->GetPartialCharge();
        }

        OBmol->SetTotalCharge(int(total));
        OBff->Setup(*OBmol);
        OBff->SteepestDescent(Input->conformer_min_steps, 1.0e-10); 		//ConjugateGradients(100);
        OBff->WeightedRotorSearch(Input->lig_conformers, Input->WRS_geom_steps);
        OBff->GetConformers(*OBmol);


// Copyng the conformer coordinates to Mol2::m_coords...
        if (OBmol->NumConformers() > 1){
            for (int i=0; i<OBmol->NumConformers(); i++){
                double* xyz = new double[OBmol->NumAtoms()*3];
                xyz = OBmol->GetConformer(i);
                OBmol->SetConformer(i);
                OBmol->SetConformer(i);
                OBff->Setup(*OBmol);
                Lig->conformer_energies.push_back(OBff->Energy());

                vector<double> v3;
                vector<vector<double> > xyz_tmp;

                for (unsigned j=0; j<OBmol->NumAtoms(); j++){
                    v3.push_back(xyz[3*j]);
                    v3.push_back(xyz[(3*j)+1]);
                    v3.push_back(xyz[(3*j)+2]);
                    xyz_tmp.push_back(v3);
                    v3.clear();
                }
                Lig->mcoords.push_back(xyz_tmp);
                xyz_tmp.clear();
            }
            file_read = true;
        }
        else {
            file_read = false;
        }
    }
    else {
        printf("Skipping conformer generation for file %s\n", molfile.c_str());
        Lig->mcoords.push_back(Lig->xyz);
        file_read = true;
    }
/*
    delete [] OBff;
    delete [] format;
    delete OBconv;
    delete OBmol;
*/
    return (file_read);
}

#include<iostream>
#include<zlib.h>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<string>
#include <vector>
#include "../LiBELa/iMcLiBELa.h"
#include "../LiBELa/main.h"
#include "../LiBELa/Mol2.h"
#include "../LiBELa/COORD_MC.h"
#include "../LiBELa/PARSER.h"
#include "../LiBELa/Optimizer.h"
#include "../LiBELa/McEntropy.h"
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/rotor.h>
#include <openbabel/conformersearch.h>
#include <openbabel/forcefield.h>
#include <openbabel/math/vector3.h>

#define k 0.0019872041

using namespace std;
using namespace OpenBabel;

OBMol GetMol(const std::string &molfile){
    OBMol mol;

    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
    printf("Could not find input format for file\n");
    return mol;
  }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading.\n", molfile.c_str());
        return mol;
    }

    if (!conv.Read(&mol, &ifs)) {
        printf("Could not read molecule from file\n");
        return mol;
    }
    return mol;
}

vector<vector<double> > copy_from_obmol(OBMol mymol){
    double* myxyz = new double[mymol.NumAtoms()*3];
    vector<vector<double > > vec_xyz;
    vector<double> tmp(3);
    myxyz = mymol.GetCoordinates();
    for (unsigned i=0; i < mymol.NumAtoms(); i++){
        tmp[0] = (myxyz[3*i]);
        tmp[1] = (myxyz[(3*i)+1]);
        tmp[2] = (myxyz[(3*i)+2]);
        vec_xyz.push_back(tmp);
    }

    tmp.clear();
    return vec_xyz;
}

int main(int argc, char* argv[]){

    string trajfile;
    string ref_mol2;
    int stride;
    int c;
    double dx, dy, dz, dalpha, dbeta, dgamma, angle;
    string ff = "GAFF";
    vector<vector<int> > atoms_in_dihedrals;

    if (argc < 2){
      printf("Usage %s -t <traj_file> -s <stride> -r <ref_mol2> [-h]\n", argv[0]);
      exit(1);
    }

    while ((c = getopt(argc, argv, "t:r:s:g:h")) != -1)
      switch (c){
      case 't':
          trajfile = string(optarg);
          break;
      case 'h':
          printf("Usage %s -t <traj_file> -s <stride> -r <ref_mol2>\n", argv[0]);
          break;
          exit(1);
      case '?':
          printf("Usage %s -t <traj_file> -s <stride> -r <ref_mol2>\n\n", argv[0]);
          break;
          exit(1);
      case 'r':
          ref_mol2 = string(optarg);
          break;
      case 's':
          stride = atoi(optarg);
          break;
      case 'f':
          ff = string(optarg);
      }

    printf("#*****************************************************************************************\n");
    printf("#                                                                                        *\n");
    printf("#         McEnergy - A program for Binding Energy computation within McLiBELa            *\n");
    printf("#                    Written by Alessandro S. Nascimento - 2018                          *\n");
    printf("#                      University of Sao Paulo - USP - Brazil                            *\n");
    printf("#                                                                                        *\n");
    printf("#                             asnascimento@ifsc.usp.br                                   *\n");
    printf("#                                                                                        *\n");
    printf("#*****************************************************************************************\n");
    printf("# Trajectory file:  %-74.74s\n", trajfile.c_str());
    printf("# Reference file: %-74.74s\n", ref_mol2.c_str());


    PARSER* Input = new PARSER;
    COORD_MC* Coord = new COORD_MC;

    Mol2* RefMol = new Mol2(Input, ref_mol2);                   // read the initial coordinates of the ligand

    OBMol mol = GetMol(Input->lig_mol2);
    double* myxyz = new double[mol.NumAtoms()*3];
    OBForceField* OBff;
    OBRotorList RotorList;
    OBRotorIterator RotorIterator;
    OBRotor *Rotor;

    if ( ff == "GAFF"){
        OBff = OBForceField::FindForceField("GAFF");
    }
    else if (ff == "MMFF94"){
        OBff = OBForceField::FindForceField("MMFF94");
    }
    else{
        OBff = OBForceField::FindForceField(ff);
    }

    if (!OBff){
        printf("Could not find OpenBabel FF parameters!\n");
        exit(1);
    }

    OBff->SetLogFile(&cout);
    OBff->SetLogLevel(OBFF_LOGLVL_LOW);

    OBff->GetCoordinates(mol);
    RotorList.Setup(mol);
    Rotor = RotorList.BeginRotor(RotorIterator);
    mol.ToInertialFrame();
    vector<int> tmp(4);
    for (unsigned i = 0; i < RotorList.Size(); ++i, Rotor = RotorList.NextRotor(RotorIterator)) {
        tmp = Rotor->GetDihedralAtoms();
        atoms_in_dihedrals.push_back(tmp);
        tmp.clear();
    }

    printf("Found %lu    rotatable bonds in ligand %s.\n", RotorList.Size(), RefMol->molname.c_str());

    vector<double> torsions(RotorList.Size());

    Mol2* TrajMol2 = new Mol2;
    TrajMol2->parse_gzipped_ensemble(Input, trajfile, stride);              // loads the trajectory at once;

    printf("#%8s %7.7s %7.7s %7.7s %7.7s %7.7s %7.7s ", "Frame", "DX", "DY", "DZ", "DALPHA", "DBETA", "DGAMMA");

    for (int j=0; j< RotorList.Size(); j++){
        printf("ROT[%3d] ", j);
    }

    printf("\n");


    for (unsigned i=0; i< TrajMol2->mcoords.size(); i++){
        Optimizer::align_t* align_data = new Optimizer::align_t;
        align_data->ref_xyz = RefMol->xyz;
        align_data->current_xyz = TrajMol2->mcoords[i];

        Optimizer::align_result_t* opt_result = new Optimizer::align_result_t;

        Optimizer* opt = new Optimizer(RefMol, RefMol, Input);
        opt->minimize_alignment_nlopt_simplex(align_data, opt_result);
        dx = opt_result->translation[0];
        dy = opt_result->translation[1];
        dz = opt_result->translation[2];
        dalpha = opt_result->rotation[0];
        dbeta = opt_result->rotation[1];
        dgamma = opt_result->rotation[2];

        for (unsigned j=0; j< RotorList.Size(); j++){
            torsions[j] = mol.GetTorsion(mol.GetAtom(atoms_in_dihedrals[j][0]), mol.GetAtom(atoms_in_dihedrals[j][1]), mol.GetAtom(atoms_in_dihedrals[j][2]), mol.GetAtom(atoms_in_dihedrals[j][3]));
        }

        printf("%8d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f ", i, dx, dy, dz, dalpha, dbeta, dgamma);

        for (int j=0; j< RotorList.Size(); j++){
            printf("%7.3f ", torsions[j]);
        }

        printf("\n");

        delete align_data;
        delete opt_result;
        delete opt;
    }

    return 0;
}


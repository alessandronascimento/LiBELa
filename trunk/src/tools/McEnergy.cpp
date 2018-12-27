#include<iostream>
#include<zlib.h>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<string>
#include <vector>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/rotor.h>
#include "../LiBELa/iMcLiBELa.h"
#include "../LiBELa/main.h"
#include "../LiBELa/Mol2.h"
#include "../LiBELa/COORD_MC.h"
#include "../LiBELa/PARSER.h"
#include "../LiBELa/Optimizer.h"
#include "../LiBELa/McEntropy.h"


#define k 0.0019872041

using namespace std;
using namespace OpenBabel;

OBMol GetMol(const string &molfile){
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
    double dx, dy, dz, dalpha, dbeta, dgamma, angle, rmsdi, rmsdf, T=300.0;
    double energy;
    vector<vector<int> > atoms_in_dihedrals;

    if (argc < 2){
      printf("Usage %s -t <traj_file> -s <stride> -r <ref_mol2> -T <temperature> [-h]\n", argv[0]);
      exit(1);
    }

    while ((c = getopt(argc, argv, "t:r:s:T:h")) != -1)
      switch (c){
      case 't':
          trajfile = string(optarg);
          break;
      case 'h':
          printf("Usage %s -t <traj_file> -s <stride> -r <ref_mol2> -f <OB Force Field> -T <temperature> [-h]\n", argv[0]);
          break;
          exit(1);
      case '?':
          printf("Usage %s -t <traj_file> -s <stride> -r <ref_mol2> -f <OB Force Field> -T <temperature> [-h]\n", argv[0]);
          break;
          exit(1);
      case 'r':
          ref_mol2 = string(optarg);
          break;
      case 's':
          stride = atoi(optarg);
          break;
      case 'T':
          T = double(atof(optarg));
          break;
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
    printf("#*****************************************************************************************\n");


    PARSER* Input = new PARSER;
    Input->temp = T;
    COORD_MC* Coord = new COORD_MC;

    Mol2* RefMol = new Mol2(Input, ref_mol2);                   // read the initial coordinates of the ligand

    OBMol mol = GetMol(ref_mol2);
    double* myxyz = new double[mol.NumAtoms()*3];
    OBRotorList RotorList;
    OBRotorIterator RotorIterator;
    OBRotor *Rotor;

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

    vector<double> com = Coord->compute_com(RefMol);

    McEntropy* Entropy = new McEntropy(Input, Coord, com, int(RotorList.Size()));

    printf("#%10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s ", "Frame", "DX", "DY", "DZ", "DALPHA", "DBETA", "DGAMMA", "RMSDi", "RMSDf");

    for (int j=0; j< RotorList.Size(); j++){
        printf(" ROT[%3d] ", j);
    }

    printf("\n");

    energy = 0.0;

    for (unsigned i=0; i< TrajMol2->mcoords.size(); i++){
        Optimizer::align_t* align_data = new Optimizer::align_t;
        align_data->ref_xyz = RefMol->xyz;
        align_data->current_xyz = TrajMol2->mcoords[i];
        rmsdi = Coord->compute_rmsd(RefMol->xyz, TrajMol2->mcoords[i], RefMol->N);

        Optimizer::align_result_t* opt_result = new Optimizer::align_result_t;

        Optimizer* opt = new Optimizer(RefMol, RefMol, Input);
        opt->minimize_alignment_nlopt_simplex(align_data, opt_result);
        dx = opt_result->translation[0];
        dy = opt_result->translation[1];
        dz = opt_result->translation[2];
        dalpha = opt_result->rotation[0];
        dbeta = opt_result->rotation[1];
        dgamma = opt_result->rotation[2];
        rmsdf = opt_result->rmsd;

        energy+= TrajMol2->ensemble_energies[i];

        for (unsigned j=0; j< RotorList.Size(); j++){
            torsions[j] = mol.GetTorsion(mol.GetAtom(atoms_in_dihedrals[j][0]), mol.GetAtom(atoms_in_dihedrals[j][1]), mol.GetAtom(atoms_in_dihedrals[j][2]), mol.GetAtom(atoms_in_dihedrals[j][3]));
        }

        Entropy->update(dx, dy, dz, dalpha, dbeta, dgamma, torsions);

        printf("%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f ", i, dx, dy, dz, dalpha, dbeta, dgamma, rmsdi, rmsdf);

        for (int j=0; j< RotorList.Size(); j++){
            printf("%10.4f ", torsions[j]);
        }
        printf("\n");

        delete align_data;
        delete opt_result;
        delete opt;
    }

    McEntropy::entropy_t* McEnt = new McEntropy::entropy_t;
    McEntropy::entropy_t* Max_Ent = new McEntropy::entropy_t;
    McEnt->Srot = 0.0; McEnt->Storsion = 0.0; McEnt->Strans = 0.0;

    Entropy->get_results(McEnt, Max_Ent, int(TrajMol2->mcoords.size()));

    energy = energy / TrajMol2->mcoords.size();

    printf("#*****************************************************************************************\n");

    printf("First-Order Approximation Translation Entropy (TS): %10.4g kcal/mol @ %7.2f K\n", McEnt->Strans*T, T);
    printf("First-Order Approximation Rotation Entropy (TS):    %10.4g kcal/mol @ %7.2f K\n", McEnt->Srot*T, T);
    printf("First-Order Approximation Torsion Entropy (TS):     %10.4g kcal/mol @ %7.2f K\n", McEnt->Storsion*T, T);
    printf("First-Order Approximation Total Entropy (S):        %10.4g kcal/(mol.K)@ %7.2f K\n", McEnt->S, T);
    printf("First-Order Approximation -TS (-TS):                %10.4g kcal/mol @ %7.2f K\n", -McEnt->TS, T);
    printf("First-Order Approximation -TS @ 300K:               %10.4g kcal/mol @ %7.2f K\n", -McEnt->S*300., 300.);

    printf("#*****************************************************************************************\n");

    printf("Maximal Entropies Computed for this System:\n");
    printf("First-Order Approximation Translation Entropy (TS): %10.4g kcal/mol @ %7.2f K\n", Max_Ent->Strans*T, T);
    printf("First-Order Approximation Rotation Entropy (TS):    %10.4g kcal/mol @ %7.2f K\n", Max_Ent->Srot*T, T);
    printf("First-Order Approximation Torsion Entropy (TS):     %10.4g kcal/mol @ %7.2f K\n", Max_Ent->Storsion*T, T);
    printf("First-Order Approximation Total Entropy (S):        %10.4g kcal/(mol.K)@ %7.2f K\n", Max_Ent->S, T);
    printf("First-Order Approximation -TS (-TS):                %10.4g kcal/mol @ %7.2f K\n", -Max_Ent->TS, T);
    printf("First-Order Approximation -TS @ 300K:               %10.4g kcal/mol @ %7.2f K\n", -Max_Ent->S*300., 300.);

    printf("#*****************************************************************************************\n");

    printf("Entropy loss (-TdS): %10.4g kcal/mol (%10.4f %s) @ %7.2f K\n", (-McEnt->TS - (-Max_Ent->TS)), ((-McEnt->TS/-Max_Ent->TS)*100), "%", T);

    printf("#*****************************************************************************************\n");

    printf("Ebind = <E> -TS\n");
    printf("Ebind = %8.4f - %8.4f = %8.4f kcal/mol\n", energy, McEnt->TS, (energy-McEnt->TS));

    printf("#*****************************************************************************************\n");

    delete [] myxyz;
    delete McEnt;
    delete Max_Ent;
    delete Entropy;
    delete RefMol;
    delete Coord;
    delete Input;
    delete TrajMol2;

    return 0;
}


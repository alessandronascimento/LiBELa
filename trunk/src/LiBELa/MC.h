#ifndef MC_H
#define MC_H

#include "PARSER.h"
#include "Grid.h"
#include "Mol2.h"
#include "COORD_MC.h"
#include "Energy2.h"
#include "Optimizer.h"
#include "WRITER.h"
#include "gsl/gsl_rng.h"
#include <zlib.h>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/rotor.h>
#include <openbabel/conformersearch.h>
#include <openbabel/shared_ptr.h>
#include <openbabel/forcefield.h>
#include <openbabel/math/align.h>
#include <openbabel/math/vector3.h>


class MC
{
public:
    double XSize, YSize, ZSize;
    vector <double> MaxMin;
    vector<vector<double> > xyz;
    double average_energy;
    double energy_standard_deviation;
    long double Boltzmann_weighted_average_energy;

    struct step_t{
        vector<vector<double> > xyz;
        double dx, dy, dz;
        double dalpha, dbeta, dgamma;
        int nconf;
        vector<double> torsion_angles;
    };

    gsl_rng * r;
    WRITER* Writer;
    char info[98];


    shared_ptr<OBMol> mol;
    OBForceField* OBff;
    OBRotorList RotorList;
    OBRotorIterator RotorIterator;
    OBRotor *Rotor;

    vector<vector<int> > atoms_in_dihedrals;

    MC(WRITER* _Writer);
    MC(Mol2* Lig, PARSER* Input);
    ~MC();
    void run(Mol2 *Rec, Mol2* Reflig , Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T);
    void run(Grid* Grids, Mol2* Reflig , Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T);
    double Boltzmman(double ene, double new_ene, double t, double b);
    void take_step(PARSER* Input, Mol2* Lig, step_t* step);
    void take_step_flex(PARSER* Input, Mol2* Lig, step_t* step);
    void write_conformers(Mol2* Lig);
    void MaxMinCM(double XCOM, double YCOM, double ZCOM, vector<double> Max);
    void ligand_run(Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T);

};

#endif // MC_H

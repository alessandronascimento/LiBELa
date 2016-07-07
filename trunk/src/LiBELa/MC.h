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
#include <openbabel/math/vector3.h>

using namespace OpenBabel;
//using namespace boost;

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
        double internal_energy;
    };

    gsl_rng * r;
    WRITER* Writer;
    char info[98];
    double* myxyz;


    shared_ptr<OBMol> mol;
    OBForceField* OBff;
    OBRotorList RotorList;
    OBRotorIterator RotorIterator;
    OBRotor *Rotor;

    vector<vector<int> > atoms_in_dihedrals;

    MC(WRITER* _Writer);
    MC(Mol2* Lig, PARSER* Input, WRITER *_Writer);
    ~MC();
    void run(Mol2 *Rec, Mol2* Reflig , Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T);
    void run(Grid* Grids, Mol2* Reflig , Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T);
    double Boltzmman(double ene, double new_ene, double t, double b);
    void take_step(PARSER* Input, Mol2* Lig, step_t* step);
    void take_step_flex(PARSER* Input, Mol2* Lig, step_t* step);
    void take_step_torsion(PARSER* Input, Mol2* Lig, step_t* step);
    void write_conformers(Mol2* Lig);
    void MaxMinCM(double XCOM, double YCOM, double ZCOM, vector<double> Max);
    void ligand_run(Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T);
    double* copy_to_obmol(vector<vector<double> > vec_xyz);
    vector<vector<double> > copy_from_obmol(shared_ptr<OBMol> mymol);
    shared_ptr<OBMol> GetMol(const std::string &molfile);
    double check_angle(double angle);

};

#endif // MC_H

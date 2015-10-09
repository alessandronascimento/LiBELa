#ifndef MCRECURSION_H
#define MCRECURSION_H


#include "../PARSER.h"
#include "../Grid.h"
#include "../Mol2.h"
#include "../COORD_MC.h"
#include "../Energy2.h"
#include "../Optimizer.h"
#include "../WRITER.h"
#include "gsl/gsl_rng.h"
#include <zlib.h>

class MC
{
public:
    long double XSize, YSize, ZSize;
    vector < long double > MaxMin;
    vector<vector<double> > xyz;

    struct step_t{
        vector<vector<double> > xyz;
        double dx, dy, dz;
        double dalpha, dbeta, dgamma;
        int nconf;
    };

    gsl_rng * r;
    WRITER* Writer;
    char info[98];
    MC(WRITER* _Writer);
    ~MC();
    void run(Mol2 *Rec, Mol2* Reflig , Mol2* Lig, vector<vector<double> > xyz, PARSER* Input);
    void run(Grid* Grids, Mol2* Reflig , Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double b, double T);
    double Boltzmman_recursive(double ene, double new_ene, double b, double t);
    void take_step(PARSER* Input, Mol2* Lig, step_t* step);
    void take_step_flex(PARSER* Input, Mol2* Lig, step_t* step);
    void MaxMinCM(long double XCOM, long double YCOM, long double ZCOM, vector<long double> Max);

};

#endif // MCRECURSION_H

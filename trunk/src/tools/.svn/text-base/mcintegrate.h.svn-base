#ifndef MCINTEGRATE_H
#define MCINTEGRATE_H

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
    void run(Grid* Grids,Mol2* Reflig ,Mol2* Lig, vector<vector<double> > xyz, PARSER* Input);
    void take_step(PARSER* Input, Mol2* Lig, step_t* step);
    void take_step_flex(PARSER* Input, Mol2* Lig, step_t* step);
    void write_conformers(Mol2* Lig);

};


#endif // MCINTEGRATE_H


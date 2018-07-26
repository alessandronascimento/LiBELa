/*
 * SA.h
 *
 *  Created on: Jan 23, 2013
 *      Author: asn
 */

#ifndef SA_H_
#define SA_H_
#include<gsl/gsl_rng.h>
#include<cmath>
#include <vector>
#include <time.h>
#include <cstdlib>
#include "Mol2.h"
#include "Grid.h"
#include "Energy2.h"
#include "COORD_MC.h"
#include "PARSER.h"

using namespace std;


class SA {
public:
	SA();
	virtual ~SA();
	vector<vector<double> > take_step(Mol2* Lig, PARSER* Input, gsl_rng* r);
	vector<vector<double> > take_step(Mol2* Lig, PARSER* Input, gsl_rng* r, vector<vector<double> > xyz);
	double evaluate_energy(Mol2* Lig2, vector<vector<double> > new_xyz, PARSER* Input, Grid* Grids);
    double evaluate_energy(Mol2* Lig2, vector<vector<double> > new_xyz, PARSER* Input, Mol2* Rec);
	double Boltzmman(double ene, double new_ene, double t);
	vector<vector<double> > optimize(Mol2* Lig, PARSER* Input, Grid* Grids, gsl_rng* r);
    vector<vector<double> > optimize(Mol2* Lig, PARSER* Input, Mol2* Rec, gsl_rng* r);
};

#endif /* SA_H_ */

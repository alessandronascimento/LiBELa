/*
 * SA.cpp
 *
 *  Created on: Jan 23, 2013
 *      Author: asn
 */

#include "SA.h"

SA::SA() {
}

SA::~SA() {
}

vector<vector<double> > SA::take_step(Mol2* Lig, PARSER* Input, gsl_rng* r){
	vector<vector<double> > new_xyz;
	double rnumber, x, y, z, a, b, g;
	COORD_MC* Coord = new COORD_MC;
	rnumber = gsl_rng_uniform(r);
	x = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
	rnumber = gsl_rng_uniform(r);
	y = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
	rnumber = gsl_rng_uniform(r);
	z = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

	rnumber = gsl_rng_uniform(r);
	a = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
	rnumber = gsl_rng_uniform(r);
	b = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
	rnumber = gsl_rng_uniform(r);
	g = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

	new_xyz = Coord->rototranslate(Lig->xyz, Lig, a, b,g, x, y, z);
	delete Coord;

	return(new_xyz);
}

vector<vector<double> > SA::take_step(Mol2* Lig, PARSER* Input, gsl_rng* r, vector<vector<double> >xyz){
	vector<vector<double> > new_xyz;
	double rnumber, x, y, z, a, b, g;
	COORD_MC* Coord = new COORD_MC;
	rnumber = gsl_rng_uniform(r);
	x = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
	rnumber = gsl_rng_uniform(r);
	y = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
	rnumber = gsl_rng_uniform(r);
	z = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

	rnumber = gsl_rng_uniform(r);
	a = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
	rnumber = gsl_rng_uniform(r);
	b = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
	rnumber = gsl_rng_uniform(r);
	g = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

	new_xyz = Coord->rototranslate(xyz, Lig, a, b,g, x, y, z);
	delete Coord;

	return(new_xyz);
}

double SA::evaluate_energy(Mol2* Lig2, vector<vector<double> > new_xyz, PARSER* Input, Grid* Grids){
    Energy2* Ene = new Energy2(Input);
	double energy = 0.00;
    energy = Ene->compute_ene(Grids, Lig2, new_xyz);
	delete Ene;
	return(energy);
}

double SA::evaluate_energy(Mol2* Lig2, vector<vector<double> > new_xyz, PARSER* Input, Mol2* Rec){
    Energy2* Ene = new Energy2(Input);
    double energy = 0.00;
    energy = Ene->compute_ene(Rec, Lig2, new_xyz);
    delete Ene;
    return(energy);
}

double SA::Boltzmman(double ene, double new_ene, double t){
    double x = - (new_ene - ene) / (0.0019858775203792202*t); // k=0.0019858775203792202
	return(exp(x));
}

vector<vector<double> > SA::optimize(Mol2* Lig, PARSER* Input, Grid* Grids, gsl_rng* r){
	double best_ene, ene, new_ene, p, rnumber;
	vector<vector<double> > xyz;
	vector<vector<double> > new_xyz;
	vector<vector<double> > best_xyz;
	double t = Input->sa_start_temp;

	ene = this->evaluate_energy(Lig, Lig->xyz, Input, Grids);
	best_ene = ene;
	xyz = Lig->xyz;
	while (t > Input->temp){
		for (int i=0; i<Input->sa_steps; i++){
			new_xyz = this->take_step(Lig, Input, r, xyz);
			new_ene = this->evaluate_energy(Lig, new_xyz, Input, Grids);
			if (new_ene <= ene){
				ene = new_ene;
				xyz = new_xyz;
				best_ene = new_ene;
				best_xyz = new_xyz;
			}
			else{
				p = this->Boltzmman(ene, new_ene, t);
                rnumber = gsl_rng_uniform(r);
				if (p > rnumber){
					ene = new_ene;
					xyz = new_xyz;
				}
			}
		}
		t = t/Input->sa_mu_t;

#ifdef DEBUG
		printf("Best Energy: %7.3f    @ temperature %7.3f\n", best_ene, t);
#endif

		ene = best_ene;
		xyz = best_xyz;
	}
	return(xyz);
}

vector<vector<double> > SA::optimize(Mol2* Lig, PARSER* Input, Mol2* Rec, gsl_rng* r){
    double best_ene, ene, new_ene, p, rnumber;
    vector<vector<double> > xyz;
    vector<vector<double> > new_xyz;
    vector<vector<double> > best_xyz;
    double t = Input->sa_start_temp;

    ene = this->evaluate_energy(Lig, Lig->xyz, Input, Rec);
    best_ene = ene;
    xyz = Lig->xyz;
    while (t > Input->temp){
        for (int i=0; i<Input->sa_steps; i++){
            new_xyz = this->take_step(Lig, Input, r, xyz);
            new_ene = this->evaluate_energy(Lig, new_xyz, Input, Rec);
            if (new_ene <= ene){
                ene = new_ene;
                xyz = new_xyz;
                best_ene = new_ene;
                best_xyz = new_xyz;
            }
            else{
                p = this->Boltzmman(ene, new_ene, t);
                rnumber = gsl_rng_uniform(r);
                if (p > rnumber){
                    ene = new_ene;
                    xyz = new_xyz;
                }
            }
        }
        t = t/Input->sa_mu_t;

#ifdef DEBUG
        printf("Best Energy: %7.3f    @ temperature %7.3f\n", best_ene, t);
#endif

        ene = best_ene;
        xyz = best_xyz;
    }
    return(xyz);
}

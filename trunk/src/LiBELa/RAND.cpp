/*
 * RAND.cpp
 *
 *  Created on: 07/10/2010
 *      Author: Nascimento
 */

#include "RAND.h"

using namespace std;

RAND::RAND(){
	srand(rand());
	r = gsl_rng_alloc (gsl_rng_ranlxs2);
	double seed = 30.0 * (rand()/(RAND_MAX + 1.0));

#ifdef DEBUG
	printf("*%-78s%20.5f*\n", "Random seed: ", seed);
#endif

	gsl_rng_set(r, seed);
}

void RAND::random(double cushion, double rotation_step){

    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	transx = -cushion + (1.0 * (rnumber*(2*cushion)));
    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	transy = -cushion + (1.0 * (rnumber*(2*cushion)));
    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	transz = -cushion + (1.0 * (rnumber*(2*cushion)));

    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	a = -rotation_step + (rnumber*(2*rotation_step));
    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	b = -rotation_step + (rnumber*(2*rotation_step));
    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	g = -rotation_step + (rnumber*(2*rotation_step));
	lign=0;

#ifdef DEBUG
	printf("a: %.1f b: %.1f c:%.1f x:%.1f y:%.1f z:%.1f\n", a, b, g, transx, transy, transz);
#endif
}

void RAND::random(double cushion, double rotation_step, Mol2 *CRec, Mol2 *CLig){


    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	transx = -cushion + (1.0 * (rnumber*(2*cushion)));
    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	transy = -cushion + (1.0 * (rnumber*(2*cushion)));
    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	transz = -cushion + (1.0 * (rnumber*(2*cushion)));

    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	a = -rotation_step + (rnumber*(2*rotation_step));
    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	b = -rotation_step + (rnumber*(2*rotation_step));
    rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
	g = -rotation_step + (rnumber*(2*rotation_step));

	if (CRec->mcoords.size() > 1){
		recn = rand() % (CRec->mcoords.size()-2);
	}
	if (CLig->mcoords.size() > 1){
		lign = rand() % (CLig->mcoords.size()-2);
	}

#ifdef DEBUG
	printf("a: %.1f b: %.1f c:%.1f x:%.1f y:%.1f z:%.1f LIGn:%d\n", a, b, g, transx, transy, transz, lign);
#endif
}

void RAND::print(){
	char info[98];
	sprintf(info,"Changing parameters: %.2f %.2f %.2f %.2f %.2f %.2f", transx, transy, transz, a, b, g);
	printf("*%-98s*\n", info);
};

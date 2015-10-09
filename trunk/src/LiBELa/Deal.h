/*
 * Deal.h
 *
 *  Created on: 01/06/2012
 *      Author: Nascimento
 */

#ifndef DEAL_H_
#define DEAL_H_

#include<QtGui>
#include<vector>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<time.h>
#include"McLiBELa.h"
#include <QList>
#include <QVector>
#include <QThread>
#include "PARSER.h"
#include "Mol2.h"
#include "WRITER.h"
#include "GUI/QtWriter.h"
#include <QtConcurrent>


/** @brief This class was written to use the Density Estination Algorithm to optimize binding energy after a rough
 * ligand docking based on ligand similarity and using a reference small molecule to 'guide' the docking. The density
 * estimation here uses gaussian estimates of rotation angles and translations.
 */

class Deal {
public:


/*!
 * Constructor. Here, the other methods are also called. So, calling the constructor, should be all that is necessary to
 * perform the optimization.
 */
	Deal(Mol2* _Rec, Mol2* _Lig, PARSER* Input);
//! Deconstructor. Deletes the gsl rng interface;
	virtual ~Deal();
//! This method sorts the energies array together with the angles/translations. It's important to select the top ranked rototranslations
//! and define the next estimative. Here, we used N**2 comparisong. This method can be further optimized.
	void sort_array(void);
//! This method iterates PARSER::deal_generation 'generations', each generation with PARSER::deal_pop_size sampling;
	void sort_and_evaluate(double *mu, int iteration);
//! Here, the estimates are computed (simple averages) for further use to sample angles and rotations in a gaussian distribution centered on mu;
	void update_mu(int top_ranked);
//! This method couples iterations and mu updates.
	void iterate(PARSER* Input);
//!
	double evaluate_energy(QVector<double> temp);
//!
	vector<int> define_chunck_sizes(void);
//!
	vector<vector<double> > rototranslate(Mol2* Lig, double alpha, double beta, double gamma, double transx, double transy, double transz);
//!
	vector<double> compute_com(Mol2 *Cmol);
//!
	double compute_ene_softcore_solvation(Mol2* Rec, Mol2* Lig, PARSER* Input, vector<vector<double> >coords);
//!
	double distance_squared (double x1, double x2, double y1, double y2, double z1, double z2);
//!
	double distance (double x1, double x2, double y1, double y2, double z1, double z2);




//! Size of population
	int population_size;
//! GSL type
	const gsl_rng_type *T;
//! GSL Random number generator interface
	gsl_rng *r;
//! Double array used to keep 6 variables: mean for alpha, beta, gamma, x, y and z;
	double *mu;
//! Double array used to keep sigmas for gaussian distribution generations. Initially defined as {60.0, 30.0,0 60.0, 1.0, 1.0, 1.0}.
	double *sigma;
//! C++ vector with the rotations and translations (angles and displacements). Important to define the next estimate.
	vector<vector<double> >rotrans;
//! C++ vector with the energies of each 'pose', defined as a set of 3 angles and translation in 3 directions.
	vector <double> energies;
//!
	Mol2* Rec;
//!
	Mol2* Lig;
//!
	PARSER* Input;
//!
	vector<int> slices;
//!
	int chunck_size;
//!
	clock_t time_i, time_f;
//!
	vector<vector<double> > final_coords;
};

#endif /* DEAL_H_ */

/*
 * ENERGY.h
 *
 *  Created on: 23/03/2010
 *      Author: Alessandro
 */

#ifndef ENERGY_H_
#define ENERGY_H_

#include<vector>
#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<stdio.h>
#include "PARSER.h"
#include "Mol2.h"
#include"RAND.h"
#include "McLiBELa.h"
#include "Grid.h"

using namespace std;

/*! @class ENERGY
 * The class ENERGY contains three methods dedicated to Binding Energy
 * computation. Here we implemented a method to compute a typical
 * VDW term (AMBER FF Style), typical electrostatic term (Coulomb potential),
 * and soft-core potentials for both VDW and electrostatic terms.
 */

class ENERGY {

public:

	// variables for energy computation

	//! van der Waals energy term
	double	vdw,
	//! Electrostatic energy term
	elec,
	//! Sum of radius of atoms i and j
	rij,
	//! Geometric mean of the well depths of for atoms i and j
	eij,
	//! Distance between atoms i and j
	dij,
	//! Squared distance (dij*dij) between atom i and j
	dij2,
	//! "Effective" distance given as (dij2*dij2*dij2)+deltaij_es6
	deff,
	//! Repulsive term of VDW potential, given as eij*(rij^12)
	acoef,
	//! Attractive term of VDW potential, given as 2*eij*(rij^6)
	bcoef,
	//! Receptor solvation energy (kcal/mol). Given as Sum over all receptor atoms (rec_solv_affinity*rec_solv_distf)
	rec_solv,
	//! Ligand solvation energy (kcal/mol). Given as Sum over all ligand atoms (lig_solv_affinity*lig_solv_distf)
	lig_solv,
	//! Phenomenologically determined affinity of each receptor atom for the solvent
	rec_solv_affinity,
	//! Phenomenologically determined affinity of each ligand atom for the solvent
	lig_solv_affinity,
	//! A gaussian weighted function of the atomic volumes displaced by ligand atoms when binding to receptor
	rec_solv_distf,
	//! A gaussian weighted function of the atomic volumes displaced by receptor atoms when binding to ligand
	lig_solv_distf;


	//
	// END OF VARIABLES
	//
	// METHODS DECLARATION
	//

	ENERGY();

	double distance(double x1, double x2, double y1, double y2, double z1, double z2);
	/*!
	 *
	 */
	double distance_squared (double x1, double x2, double y1, double y2, double z1, double z2);
	/*!
	 * This function computes the electrostatic term of the binding potential energy
	 */
	double compute_elec(vector<double> rec_charges, vector<double> lig_charges, vector<double> rec_crd, vector<double> lig_crd );

	/*!
	 * This function computes the van der Waals term of the binding potential energy
	 */
	double compute_vdw(vector<double> rec_epsilons, vector<double> lig_epsilons, vector<double> rec_radii, vector<double> lig_radii, vector<double> rec_crd, vector<double> lig_crd);

	/*!
	 * This function computes both, electrostatic and vdw term in a single method, i.e., looping just once over atoms i and j.
	 */
	double compute_ene(vector<double> rec_epsilons, vector<double> lig_epsilons, vector<double> rec_radii, vector<double> lig_radii, vector<double> rec_charges, vector<double> lig_charges, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd);

	/*!
	 * This function computes the VDW binding term using a Soft-core potential.
	 * Details about the potential can de found in Verkhivker GM, et al (1999).
	 * J. Mol. Recog. 12:371-389.
	 */
	double compute_vdw_softcore(vector<double> rec_epsilons, vector<double> lig_epsilons, vector<double> rec_radii, vector<double> lig_radii, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double deltaij6);

	/*!
	 * This function computes the electrostatic binding term using a Soft-core potential.
	 * Details about the potential can de found in Verkhivker GM, et al (1999).
	 * J. Mol. Recog. 12:371-389.
	 */
	double compute_elec_softcore(vector<double> rec_charges, vector<double> lig_charges, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double diel, double deltaij_es6 );

	/*! This function computed the energy penalty due to the rupture of interactions with
	 * solvent when a ligand binds to a receptor. It is derived from Verkhivker GM (1999)
	 * with AMBER radii and atomic charges. If necessary, other solvent modeling schemes
	 * such as GB, can be implemented in the future.
	 */
    double compute_solvation(vector<double>rec_charges, vector<double> lig_charges, vector<double>rec_radii, vector<double>lig_radii, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double deltaij_es6, double sigma);


	/*! This method computes the total system energy, i.e., the VDW term, in a softcore fashion;
	 * the electrostatic term, also in the softcore fashion; and the solvation correction as described
	 * in the method above. The main difference is that the atoms in the system are looped only once.
	 */
	double compute_ene_softcore_solvation(vector<double>rec_charges, vector<double> lig_charges, vector<double>rec_radii, vector<double>lig_radii,vector<double> rec_epsilons, vector<double> lig_epsilons,vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double diel, double deltaij_es6, double deltaij6, double sigma);
//! Method overloaded
	double compute_ene_softcore_solvation(Mol2* REC, Mol2* LIG,PARSER* Input, vector<vector<double> >coord);
//! Method overloaded
	double compute_ene_softcore_solvation(Mol2* REC, Mol2* LIG,PARSER* Input);

	/*!
	 *
	 * @param REC
	 * @param LIG
	 * @param Input
	 * @param Rand
	 * @return
	 */
	double compute_ene_softcore_solvation(Mol2* REC, Mol2* LIG,PARSER* Input, RAND *Rand);

	/*!
	 *
	 */
	double compute_ene_softcore_solvation(Mol2* REC, vector<vector<double> >coord, Mol2* LIG,PARSER* Input, RAND *Rand);

	double compute_ene_softcore_solvation(Mol2* REC, vector<vector<double> >coord, vector<vector<double> >lig_coords, Mol2* LIG,PARSER* Input);

	double compute_ene(Mol2* REC, vector<vector<double> >coords, Mol2* LIG, /* PARSER* Input,*/ RAND* Rand);

	double compute_ene(Mol2* REC, Mol2* LIG, /*PARSER* Input,*/ RAND *Rand);

	double compute_ene(Mol2* REC, Mol2* LIG,/*,PARSER* Input,*/ vector<vector<double> >coord);

	double compute_ene(Mol2* REC, vector<vector<double> >coord, vector<vector<double> >lig_coords, Mol2* LIG/*,PARSER* Input*/);

    double compute_ene_solvation(vector<double> rec_epsilons, vector<double> lig_epsilons, vector<double> rec_radii, vector<double> lig_radii, vector<double> rec_charges, vector<double> lig_charges, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double deltaij_es6, double sigma);

	double compute_ene_softcore(vector<double>rec_charges, vector<double> lig_charges, vector<double>rec_radii, vector<double>lig_radii,vector<double> rec_epsilons, vector<double> lig_epsilons,vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double diel, double deltaij_es6, double deltaij6);

	double compute_ene_from_grids(Grid* Grids, Mol2* Lig);
	double compute_ene_from_grids(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
	double compute_ene_from_grids_nosolv(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
	double compute_ene_from_grids_hardcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
	double compute_ene_from_grids_hardcore_nosolv(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);

};


#endif /* ENERGY_H_ */


/*
 * Evdw = [ Aij/(rij6 + dij6)^2 ] - [ Bij / (rij6 + dij6) ]
 * dEvdw/dr = [ (-12 Aij rij5) / (rij6+dij6)^3 ] - [ (-6 Bij rij5) / (rij6+dij6)^2 ]
 *
 * Eele = [qi*qj / 4*pi*e0*D*[(rij6+dij6)^(1/3)]
 * dEele/dr = -6*qi*qj*rij5 / 3*4*pi*e0*D*[(rij6+dij6)^(4/3)]
 *
 * Esol = SiXj + SjXi // Aqui a Fun��o S � invariante com a posi��o. A derivada de X �
 * dX/dr = [fj/(sigma^3)] * exp{- [ (rij6 + dij6)^(1/3)]/(2*sigma^2)} * (-1/{ 6*[ (rij6+dij6)^(2/3)]* sigma^2}) * 6rij^5
 */


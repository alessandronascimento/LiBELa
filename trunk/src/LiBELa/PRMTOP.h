/*
 * PRMTOP.h
 *
 *  Created on: 23/03/2010
 *      Author: Alessandro
 */

#ifndef PRMTOP_H_
#define PRMTOP_H_

#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<stdio.h>
#include<stdlib.h>

using namespace std;
/*! @class PRMTOP
 * The PRMTOP class has a number of methods dedicated to parse a AMBER prmtop
 * file and get all the necessary information for energy evaluation in a receptor-ligand
 * system. There is also a method to parse the atomic coordinates of a inpcrd file. The
 * parser works well for ligand and receptor prmtop/inpcrd files.
 */
class PRMTOP {

public:

	// variables

	//! Number of atoms
	int N;
	//! Number of residues
	int Nres;
	//! Number of AMBER atomtypes
	int Natomtypes;
	//! Counter variable
	int count;
	//! String to parse atom names
	string atm,
	//! String to parse comment lines
	line;
	//!
	vector<string> atomnames;
	//! C++ vector to store atomic charges
	vector<double>charges;
	//! C++ vector to store atomic masses
	vector<double>masses;
	//! C++ vector to store the residue labels
	vector<string> resnames;
	//! C++ vector to store the AMBERatoms
	vector<string>amberatoms;
	//! C++ vector with pointer to the number of the atoms in the residue
	vector<int> residue_pointer;
	//! C++ vector to store the atomtypes in the FF
	vector<string> a_atomtypes_prm;
	//! C++ vector to store FF atomic radii
	vector<double>a_radius;
	//! C++ vector to store FF well depth for Amberatoms
	vector<double>a_welldepth;
	//! C++ vector with N size to store the epsilons of the amberatoms;
	vector<double>epsilons;
	//! C++ vector with N size to store the radii of the amberatoms;
	vector<double>radii;
	//! C++ vector to store the atomic coordinates (3*N size)
	vector<double>coords;
	//! Atomic mass
	double mass,
	//! Atomic charge
	charge,
	//! Atomic coordinate;
	coord;
	//!
	vector<vector<double> >mcoords;
	//!
	vector<vector<double> >new_mcoords;
	PRMTOP *Cmol;

	//methods used to parse prmtop/inpcrd file
/*!
 * Class constructor. The some internal methods are called to parse the
 * number of atoms, the atomic charges, the atomic masses, the amberatom
 * names, and the coordinates. Also, the parameters of the force field are
 * taken from the file "vdw.param" (epsilons, radius)
 */
	PRMTOP(ifstream &prmtop, ifstream &inpcrd);

	/*!
	 * This method parses the prmtop file to get the number
	 * of atoms in the system
	 */
	void get_N(ifstream &prmtop);

	//!
	void get_atomnames(ifstream &prmtop, int N);

	/*!
	 * This method parses the prmtop file to get the atomic
	 * charges (electron charge multiplied by 18.2223).
	 */
	void get_charges(ifstream &prmtop, int N);

	/*!
	 * This method parses the PRMTOP file to get the residue
	 * label.
	 */
	void get_resnames(ifstream &prmtop);

	/*!
	 * This method parses the PRMTOP file to get the pointers to
	 * residues. That is, the number of the initial atom of each residue.
	 */
	void get_res_pointers(ifstream &prmtop);

	/*!
	 * This method parses the prmtop file to get the atomic
	 * masses of the atoms in the system and stores them in a
	 * C++ vector
	 */
	void get_masses(ifstream &prmtop, int N);

	/*!
	 * This method parses the atomic types of the system in
	 * Amber FF type.
	 */
	void get_amberatoms(ifstream &prmtop, int N);

	/*!
	 * This method parses the file "vdw.param" to get the parameters
	 * of the AMBER force field.
	 */
	void read_atomtypes_prm();

	/*!
	 * This method parses the amberatoms in the system and returns two
	 * C++ vectors with the radii and well depths for each atom type
	 * according to the AMBER FF (Cornell, 99).
	 */
	void get_vdw_parms(vector<string>amberatoms, vector<string>a_atomtypes_prm, vector<double> a_radius, vector<double>a_welldepth);

	/*!
	 * This method parses the inpcrd file to get the atomic coordinates of
	 * every atom in the system. The coordinates are them stored in a C++
	 * vector.
	 */
	void get_coordinates(ifstream &inpcrd);

	/*!
	 *
	 * @param mdcrd
	 */
	void get_trajectory(ifstream &mdcrd);
};

#endif /* PRMTOP_H_ */

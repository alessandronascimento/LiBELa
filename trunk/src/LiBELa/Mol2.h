/*
 * Mol2.h
 *
 *  Created on: 10/10/2011
 *      Author: Nascimento
 */

#ifndef MOL2_H_
#define MOL2_H_

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<cmath>
#include"PARSER.h"
#include <zlib.h>

using namespace std;
//using namespace OpenBabel;

class Mol2 {
public:


// variaveis

	//! Number of atoms
	int N ;
	//! Number of residues
	int Nres;
	//! Number of atomtypes
	int Natomtypes;
	//! Number of bonds;
	int Nbonds;
	//!
	string molname;
	//! Charges of the atoms. Must be divided by 18.2223 to get electron charges
	vector<double>charges;
	//! Atomic masses for the system.
	vector<double>masses;
	//! Atom names for the system according to AMBER FF.
	vector<string>amberatoms;
	//! Atomic parameters
	vector<string>atomtypes_prm;
	//! Atomic names
	vector<string> atomnames;
	//! Atomic radii
	vector<double>radius;
	//! Welldepth parameters to each atom according to AMBER FF
	vector<double>welldepth;
	//! Atomic coordinates
	vector<vector<double> > xyz;
	//! Atomic coordinates after rotation/translation
	vector<vector<double> > new_xyz;
	//! Atomic coordinates from the last step accepted.
	vector<vector<double> > old_xyz;
	//! Atomic coordinates from the optimized pose.
	vector<vector<double> > opt_overlay_xyz;
	//! Atomic coordinates from the optimized pose.
	vector<vector<double> > opt_energy_xyz;
	//! Trajectory vector
	vector<vector<vector<double> > > mcoords;
	//! Atomic parameter for LJ computation
	vector<double>epsilons;
	//! Squared root of atoms epsilons
	vector<double> epsilons_sqrt;
	//! Atomic radii
	vector<double>radii;
	//! C++ vector with the names of the residues
	vector<string> resnames;
	//! Pointers to the number of residues / number of atoms.
	vector<int> residue_pointer;
	//! Temporary string to parse prmtop file.
	string line;
	//!
	char str[80];
	//!
	char* elsa_dir_path;
	//! Keeps Vaa for RefMol/CompMol
	double self_obj_function;
	//!
	vector<vector<string> >bonds;
	//!
	vector<string> sybyl_atoms;

	vector<vector<vector<double> > > new_mcoords;
	string atm;
    //! Energies evaluated for the conformers generated. Uses GAFF.
    vector<double> conformer_energies;

	/*!
	 * Initializer. This class has, as arguments, a pointer to the class PARSER.
	 * The class uses some information given by the user there. The filename is also
	 * given as argument. This makes possible to use the same object to reference and
	 * comparing molecules.
	 */
	Mol2();
	Mol2(PARSER *Input, ifstream &mol2file);
	Mol2(PARSER *Input, string molfile);
    bool parse_gzipped_file(PARSER* Input, string molfile);
    bool parse_mol2file(PARSER* Input, string molfile);

	/*!
	 * This method is used to manually convert SYBYL atom types to AMBER (GAFF) atom
	 * types. The conversion does not need to be very accurate since the VDW parameters
	 * are the only parameters used here and only for alignment purposes.
	 */
	string convert2gaff(string atom);

	/*!
	 * This method parses the file "vdw.param" to get GAFF atomic VDW parameters
	 */
	void read_atomtypes_prm();

	/*!
	 *
	 */
	void get_epsilon(vector<string>atomtypes_prm, string amberatom, vector<double>welldepth);

	/*
	 *
	 */
	void get_radius(vector<string>atomtypes_prm, string amberatom, vector<double>radius);

	/*!
	 *
	 */
	void get_masses(string atomname);

	/*!
	 *
	 */
	void get_trajectory(ifstream &mdcrd);

	/*!
	 *
	 */
	void get_gaff_parameters();

	/*!
	 *
	 */
	~Mol2();

    bool parse_gzipped_ensemble(string molfile, int skipper);
    vector<double> ensemble_energies;
    vector<double> ensemble_rmsd;
};

#endif /* MOL2_H_ */

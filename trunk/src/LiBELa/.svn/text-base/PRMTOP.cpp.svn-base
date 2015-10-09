/*
 * PRMTOP.cpp
 *
 *  Created on: 23/03/2010
 *      Author: Alessandro
 */

#include "PRMTOP.h"



PRMTOP::PRMTOP(ifstream &prmtop, ifstream &inpcrd){
	PRMTOP::get_N(prmtop);
	PRMTOP::get_atomnames(prmtop, PRMTOP::N);
	PRMTOP::get_charges(prmtop, PRMTOP::N);
	PRMTOP::get_masses(prmtop, PRMTOP::N); // unused
	PRMTOP::get_resnames(prmtop);
	PRMTOP::get_res_pointers(prmtop);
	PRMTOP::get_amberatoms(prmtop, PRMTOP::N);
	PRMTOP::read_atomtypes_prm();
	PRMTOP::get_vdw_parms(PRMTOP::amberatoms, PRMTOP::a_atomtypes_prm, PRMTOP::a_radius, PRMTOP::a_welldepth);
	PRMTOP::get_coordinates(inpcrd);
}

void PRMTOP::get_N(ifstream &prmtop) {
	if (prmtop.is_open()){
		getline (prmtop, line);
		for (int i=1; i<=5; i++) {
			getline (prmtop, line); }
		prmtop >> PRMTOP::N >> PRMTOP::Natomtypes;
		for (int i=1; i<=9; i++){
			prmtop >> line;
		}
		prmtop >> PRMTOP::Nres;
		getline (prmtop, line);
	}
	else {
		cout << "Couldn't open prmtop file;" << endl;
		exit(1);
	}
}

void PRMTOP::get_atomnames(ifstream &prmtop, int N){
	string name, name2;
	getline (prmtop, line);
	while (line.size() < 6 or line.substr(6,9) != "ATOM_NAME")  {
		getline(prmtop, line);
	}
	getline(prmtop, line);
	int i=0;
	while (i<N){
		prmtop >> name;
		if (name.size() > 4){
			while (name.size()>4) {
				name2=name.substr(0,4);
				PRMTOP::atomnames.push_back(name2);
				i++;
				name=name.substr(4);
			}
			PRMTOP::atomnames.push_back(name);
			i++;
		}
		else {
			PRMTOP::atomnames.push_back(name);
			i++;
		}
	}
}

void PRMTOP::get_charges(ifstream &prmtop, int N) {
	getline (prmtop, line);
	while (line.size() < 6 or line.substr(6,6) != "CHARGE")  {
		getline (prmtop, line); }
	getline (prmtop, line);		// FORMAT
	for (int i=1; i<=N; i++) {
		prmtop >> charge;
		PRMTOP::charges.push_back(charge); }
}


void PRMTOP::get_masses(ifstream &prmtop, int N) {
	getline (prmtop, line);
	getline (prmtop, line);
	while (line.size() < 6 or line.substr(6,4) != "MASS")  {
		getline (prmtop, line); }
	getline (prmtop, line);		// FORMAT
	for (int i=1; i<=N; i++) {
		prmtop >> mass;
		PRMTOP::masses.push_back(mass); }
	getline (prmtop, line);
	getline (prmtop, line);
}

void PRMTOP::get_resnames(ifstream &prmtop){
	getline (prmtop, line);
	string resname;
	while (line.size() < 6 or line.substr(6,13) != "RESIDUE_LABEL")  {
		getline (prmtop, line);
	}
	getline (prmtop, line);		// FORMAT
	for (int i=1; i<= this->Nres; i++){
		prmtop >> resname;
		this->resnames.push_back(resname);
	}
	getline (prmtop, line);
}

void PRMTOP::get_res_pointers(ifstream &prmtop){
	getline(prmtop, line);
	int respointer;
	while (line.size() < 6 or line.substr(6,15)!= "RESIDUE_POINTER")  {
		getline (prmtop, line);
	}
	getline (prmtop, line);		// FORMAT
	for (int i=0; i < this->Nres; i++){
		prmtop >> respointer;
		PRMTOP::residue_pointer.push_back(respointer);
	}
	getline(prmtop, line);
}

void PRMTOP::get_amberatoms(ifstream &prmtop, int N) {
	getline (prmtop, line);
	while (! prmtop.eof()) {
		getline (prmtop, line);
		if (line.size() >= 15) {
			if (line.substr(6,15) == "AMBER_ATOM_TYPE")  {
				getline (prmtop, line);		// FORMAT
				for (int i=1; i<=N; i++) {
					prmtop >> atm;
					PRMTOP::amberatoms.push_back(atm); }
			}
		}
	}
}

void PRMTOP::read_atomtypes_prm() {
	double rad, well;
	ifstream vdwprm("vdw.param");
	if (!vdwprm.is_open()) {
		cout << "LiBELa found a problem when opening the file \"vdw.param\"." << endl;
		exit(1) ; }
	else {
		getline(vdwprm, line);
		while (!vdwprm.eof()) {
			vdwprm >> atm >> rad >> well;
			PRMTOP::a_atomtypes_prm.push_back(atm);
			PRMTOP::a_radius.push_back(rad);
			PRMTOP::a_welldepth.push_back(well);
		}
	}
}

void PRMTOP::get_vdw_parms(vector<string>amberatoms, vector<string>a_atomtypes_prm, vector<double> a_radius, vector<double>a_welldepth) {
	for (unsigned i=0; i< amberatoms.size(); i++){
		count = 0;
		for(unsigned j=0; j< a_atomtypes_prm.size(); j++){
			if (amberatoms[i] == a_atomtypes_prm[j]) {
				PRMTOP::epsilons.push_back(this->a_welldepth[j]);
				PRMTOP::radii.push_back(this->a_radius[j]);
				count=1;
			}
		}
		if (count < 1 ) {
			cout << "Missing VDW Parameter for atomtype: " << amberatoms[i] << endl;
			cout << "Program execution will be aborted !" << endl;
			exit(1);
		}
	}
}

void PRMTOP::get_coordinates(ifstream &inpcrd){
	if (inpcrd.is_open()){
		getline(inpcrd, atm); // comment
		getline(inpcrd, atm); // # of atoms
		for(int i =0;i<(3*N);i++){
			inpcrd >> coord;
			PRMTOP::coords.push_back(coord);
		}
	}
	else {
		printf("Could not open coordinate file. Please check!\n");
		exit(1);
	}
}

void PRMTOP::get_trajectory(ifstream &mdcrd){
	if (mdcrd.is_open()){
		getline(mdcrd, atm);
		vector<double> temp;
		while(!mdcrd.eof()){
			for (int i=0; i<(3*N); i++){
				mdcrd >> coord;
				temp.push_back(coord);
			}
			PRMTOP::mcoords.push_back(temp);
			temp.clear();
		}
	}
	else {
		printf("Could not open trajectory file. Please check!\n");
		exit(1);
	}
}

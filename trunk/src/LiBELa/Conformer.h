/*
 * Conformer.h
 *
 *  Created on: 02/08/2012
 *      Author: Nascimento
 */

#ifndef CONFORMER_H_
#define CONFORMER_H_

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/conformersearch.h>
#include <openbabel/shared_ptr.h>
#include <openbabel/forcefield.h>
#include <openbabel/math/align.h>
#include <vector>
#include "PARSER.h"
#include "Mol2.h"
#include <cstdlib>
#include <cstdio>
#include <memory>

using namespace OpenBabel;
using namespace std;

/**
 * @brief The Conformer class is dedicated to the conformer generation during the docking
 * procedure. It heavily relies on the OpenBabel API to such task.
 */

class Conformer {
public:
	Conformer();
	~Conformer();
    shared_ptr<OBMol> GetMol(const string &molfile);
    //! Generates molecule conformers using Genetic Algorithm, as implemented in OpenBabel
	bool generate_conformers_GA(PARSER* Input, Mol2* Lig, string molfile);
    //! Generates molecule conformers using the weighted rotor search, as implemented in OpenBabel
	bool generate_conformers_WRS(PARSER* Input, Mol2* Lig, string molfile);

};

#endif /* CONFORMER_H_ */

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

//#ifdef HAS_GUI

//#define shared_ptr boost::shared_ptr

//#endif


using namespace OpenBabel;
namespace std {

class Conformer {
public:
	Conformer();
	~Conformer();
    shared_ptr<OBMol> GetMol(const std::string &molfile);
	bool generate_conformers_GA(PARSER* Input, Mol2* Lig, string molfile);
	bool generate_conformers_WRS(PARSER* Input, Mol2* Lig, string molfile);

};

} /* namespace std */
#endif /* CONFORMER_H_ */

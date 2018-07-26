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
 * procedure.
 */

/**
 * @details The Conformer class is dedicated to the conformer generation in docking and MC simulations.
 * It is entirely based on the OpenBabel API and, as implemented in this library, it can use two
 * different functions to generate conformers. The first one is the weighted rotor search (WRS) and the
 * second one is the genetic algorithm (GA). These functions are implemented in this class.
 */

class Conformer {
public:
    //! Constructor
    Conformer();
    //! Destructor
	~Conformer();
    /**
     * @brief GetMol opens a MOL2 file as OBMol type
     * @param molfile a C++ string with the mol2 file to be opened
     * @return a OBMol type with the molecule characteristics.
     */
    shared_ptr<OBMol> GetMol(const string &molfile);
    //!
    /**
     * @brief generate_conformers_GA Generates molecule conformers using Genetic Algorithm, as implemented in OpenBabel
     * @param Input Pointer to the PARSER object
     * @param Lig Pointer to the Mol2 object
     * @param molfile C++ string with the MOL2 file to be processed.
     * @return true if generations successeds or false elsewhere
     */
    bool generate_conformers_GA(PARSER* Input, Mol2* Lig, string molfile);

    /**
     * @brief generate_conformers_WRS Generates molecule conformers using the weighted rotor search, as implemented in OpenBabel
     * @param Input Pointer to the PARSER object
     * @param Lig Pointer to the Mol2 object
     * @param molfile C++ string with the MOL2 file to be processed.
     * @return true if generations successeds or false elsewhere
     */
    bool generate_conformers_WRS(PARSER* Input, Mol2* Lig, string molfile);

};

#endif /* CONFORMER_H_ */

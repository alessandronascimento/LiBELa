/*
 * Grid.h
 *
 *  Created on: Jan 7, 2013
 *      Author: asn
 */

#ifndef GRID_H_
#define GRID_H_

#include <vector>
#include <cmath>
#include <sstream>
#include "Mol2.h"
#include "PARSER.h"

class Grid {
public:
	double grid_spacing;
	int npointsx, npointsy, npointsz;
	double xbegin, ybegin, zbegin, xend, yend, zend;
	vector<vector<vector<double> > > elec_grid;
    vector<vector<vector<double> > > pbsa_grid;
	vector<vector<vector<double> > > vdwA_grid;
	vector<vector<vector<double> > > vdwB_grid;
	vector<vector<vector<double> > > solv_gauss;
	vector<vector<vector<double> > > rec_solv_gauss;
	double rec_si;
	PARSER* Input;
    bool pbsa_loaded; // = false;

    /*!
     * \brief Grid This function initializes the class. It passes a copy of the PARSER
     * class to the Grid class and defines the grid spacing
     * \param _Input Pointer to the PARSER class.
     */
	Grid(PARSER* _Input);
    /*!
     * \brief Grid Initializer of the Grid class used through the code. This method copies
     * a pointer to the class PARSER to this class, defines the grid spacing and calls the
     * appropriate methods to start grid computation. Also, if the keyword "write_grids" is
     * set to "yes" in the PARSER, it calls the method to write the grids to a file.
     * \param _Input Pointer to the PARSER class.
     * \param Rec POINTER to a Mol2 class with Receptor information.
     * \param com C++ vector with the coordinates of the center of mass of the reference
     * ligand. It is used to define the center of the computation box.
     */
	Grid(PARSER* _Input, Mol2* Rec, vector<double> com);
	double distance(double x1, double x2, double y1, double y2, double z1, double z2) ;
	double distance_squared(double x1, double x2, double y1, double y2, double z1, double z2) ;
	void generate_points(vector<double> ref_com);
    void compute_grid_softcore(Mol2* Lig);
	void compute_grid_hardcore(Mol2* Lig);
	void write_grids_to_file(void);
	void load_grids_from_file(void);
    void load_Ambergrids_from_file(void);
	virtual ~Grid();
};

#endif /* GRID_H_ */

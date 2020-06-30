#ifndef RUNENGINE_H
#define RUNENGINE_H

#include<vector>
#include"PARSER.h"
#include"COORD_MC.h"
#include"ENERGY.h"
#include"RAND.h"
#include"WRITER.h"
#include"Docker.h"
#include"Mol2.h"
#include "MC.h"
#include"Conformer.h"
#include<iostream>
#include<fstream>
#include<cstdlib>
#include <string>
#include <sstream>
//#include<time.h>
#include<cstdio>
#include<ctime>
#include<unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include<gsl/gsl_rng.h>

#ifdef HAS_GUI
#include "GUI/QtWriter.h"
#include <QElapsedTimer>
#include <QProgressBar>
#include <QVector>
#include <QPlainTextEdit>
#include <QThread>
#include <QtGui>
#include <QtConcurrent>
#endif

#ifdef HAS_MPI
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
using namespace boost;
#endif

using namespace std;

class TEMP_SCHEME{

public:
	//! Atomic coordinate
	double crd;
	//! Electrostatic potential energy
	double e_elec;
	//! VDW potential energy
	double e_vdw;
	//! Total potential energy
	double e_tot;
	//! Solvation Energy
	double e_solv;
	//! Energy of the system in the initial state.
	double start_energy;
	//! Energy of the system after a change in the coordinates.
	double new_energy;
	//! Clock variables to time functions.
	clock_t start, end;
	//! Random number. Sorted using STDLIB.
	double rnumber;
	//! Alpha angle in Euler notation. Used for rotation.
	double a;
	//! Beta angle in Euler notation. Used for rotation.
	double b;
	//! Gamma angle in Euler notation. Used for rotation.
	double g;
	//! Shift in X axis (X translation)
	double transx;
	//! Shift in Y axis (Y translation)
	double transy;
	//!Shift in Z axis (Z translation)
	double transz;
	//! C++ vector used to store the modified coordinates for the ligand
	vector<vector<double> >new_coord;
	//! //! C++ vector used to store the modified coordinates for the ligand
	vector<vector<double> >old_coord;
	//! A three-element C++ vector used to store the Center Of Mass (COM) of the ligand coordinates
	vector<double> com;
	//! Root mean square deviation.
	double rmsd;
	//!
	double prob;
	//!
	char info[98];
	//!
	FILE *output;
	//!
	char rmsd_ene[100];
	//!
	FILE *rmsd_energy;
	//!
	int step;
	//!
	int step_accepted;
	//!
	vector<double> center;
	//!
	clock_t time0, time1;
	//!
	int time_elapsed;
	//!
	vector<double> com_tmp;
	//!
	double acceptance_rate;
	//!
	Grid* Grids;
    //!
    int argc;
    //!
    char** argv;


    ~TEMP_SCHEME();

#ifdef HAS_GUI
    QtWriter* QWriter;
    QProgressBar* progressbar;

    TEMP_SCHEME(PARSER* Input, QPlainTextEdit* Editor, QProgressBar* _progressbar );

	void sa_run(PARSER* Input);
#endif

/*!
 *
 * Input is defined as a PARSER class to parse the parameters
 * required to McLibela to run.
 *
 */
	PARSER* Input;


/*! REC is defined as a Mol2 class to include parameters and
 * coordinates of the receptor atoms.
 */
	Mol2 *REC;


/*!
 *  LIG is defined as a Mol2* class to include parameters and
 * coordinates of the ligand atoms.
 *
 */
	Mol2 *LIG;

/*!
 * \brief RefLig is defined as a Mol2* class to include parameters and
 * coordinates of the ligand atoms.
 */
	Mol2 *RefLig;

/*!
 * Ene is defined as a Energy2 class to include methods for binding
 * energy computation
 */
    Energy2* Ene;


/*!
 * COORD is defines as a COORD_MC class to include the methods for
 * changes in the coordinates (rotation, translation, roto-translation
 * and rmsd).
 *
 */
    COORD_MC COORD;

/*!
 * Rand is defined as a RAND class used in McLibela to generate random
 * numbers that define random walks (movements) of the ligand in the
 * binding site
 */
	RAND Rand;

/*! The Class Writer has the methods for writing output files, including
 * log files and coordinate files.
 */
	WRITER* Writer;

/*!
 * The constructor defines the Input object and parses the input file
 */
    TEMP_SCHEME(int ac, char* av[]);

/*!
 * The evaluation method evaluates the initial state of ligand and protein
 * and generates the first random movement of the ligand to start the simulation
 */
	void evaluation();


/*!
 * The sa_run method defines a simulated annealing run with initial temperature
 * defined by the user. The user also defines how many steps will be simulated
 * before the temperature is decreased. The temperature decrease is not linear
 * in McLibela. Rather than linearly, the temperature is decreased by 10*log(t)
 * in each SA step.
 */
	void sa_run();

/*!
 * \brief dock_run is a general method to start the docking engine within LiBELa.
 * From this general method, other methods are called to allow parallel docking,
 * and docking of multiple conformers.
 */

    void dock_run();
/*!
 * \brief dock_parallel_function is a general function to dock a molecule using a
 * parallel implementation. It is called from the current MPI implementation to
 * allow docking of a single molecule by a MPI job.
 * \param Lig2 is a Mol2 object containing the molecule to be docked
 * \return zero if docking is performed.
 */
	int dock_parallel_function(Mol2* Lig2);

/*!
 * \brief dock_concurrent_function Legacy method. Not used anymore
 * \param mylig
 * \return
 */
    int dock_concurrent_function(string mylig);

/*!
 * \brief dock_serial This method used to dock a list of ligands in a serial fashion.
 * \param ligand_list vector with the list of ligands to be docked
 * \param count number of ligands that will be docked
 * \param chunck_size Chunck size.
 * \return
 */

    int dock_serial(vector<string> ligand_list, int count, int chunck_size);

/*!
 * \brief dock_parallel Implementation of the parallel execution for docking.
 */
    void dock_parallel(void);

/*!
 * \brief dock_mpi Implementation of docking via MPI. Called by dock_parallel method.
 */
    void dock_mpi();


/*!
 * \brief eq_run General method for the Monte Carlo Simulation of the Ligand-Receptor
 * complex.
 */
    void eq_run();

/*!
 * \brief mcr_run Recursive Monte Carlo method as, described by Scheraga.
 */
    void mcr_run();


#ifdef HAS_GUI
/*!
 * \brief dock_concurrent Legacy method.
 */
	void dock_concurrent(void);
#endif

    /*!
     * \brief print_info Acessory method to print information using Writer class or Qwriter class.
     * \param info Array of char with info to be printed
     */
    void print_info(char info[98]);

/*!
 * \brief print_line Displays a line
 */
    void print_line();

/*!
 * \brief write_box Writes a box used for computing the grid
 * \param center C++ vector with the coordinates for the center of the box
 * \param min_x Lower limit for x-coordinate
 * \param min_y Lower limit for y-coordinate
 * \param min_z Lower limit for z-coordinate
 * \param max_x Upper limit for x-coordinate
 * \param max_y Upper limit for y-coordinate
 * \param max_z Upper limit for z-coordinate
 */
    void write_box(vector<double>center, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z);

};

#endif

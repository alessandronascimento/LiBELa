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
#include<time.h>
#include<cstdio>
#include<ctime>
#include<gsl/gsl_rng.h>

#ifdef HAS_GUI
#include "GUI/QtWriter.h"
#include <QProgressBar>
#include <QVector>
#include <QTextEdit>
#include <QThread>
#include <QtGui>
#include <QtConcurrent>
#endif

#ifdef HAS_MPI
#include <boost/mpi.hpp>
using namespace boost;
#endif

using namespace std;

/*
#if __GNUC__ > 4
    #define shared_ptr boost::shared_ptr
#endif

using namespace boost;
*/

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

#ifdef HAS_GUI
	QtWriter* QWriter;
	TEMP_SCHEME(PARSER* Input, QTextEdit* Editor);
	void evaluation(PARSER* Input, QProgressBar* progressbar);
	void sa_run(PARSER* Input);
    void run_dock_gui(PARSER *Input, QProgressBar *progressbar);
    void dock_parallel_gui(PARSER *Input, QProgressBar *progressbar);
#endif

/*!
 *
 * Input is defined as a PARSER class to parse the parameters
 * required to McLibela to run.
 *
 */

	PARSER* Input;


/*! REC is defined as a PRMTOP class to include parameters and
 * coordinates of the receptor atoms.
 */

	Mol2 *REC;


/*!
 *  LIG is defined as a PRMTOP class to include parameters and
 * coordinates of the ligand atoms.
 *
 */


	Mol2 *LIG;
	Mol2 *RefLig;

/*!
 * Ene is defined as a ENERGY class to include methods for binding
 * energy computation
 *
 */

    //ENERGY Ene;
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

	/*!
	 *
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
	 *
	 */
	void dock_run();
	int dock_parallel_function(Mol2* Lig2);
	int dock_concurrent_function(string mylig);
    int dock_serial(vector<string> ligand_list, int count, int chunck_size);
	void dock_parallel(void);
    void dock_mpi();

    void eq_run();
    void mcr_run();
#ifdef HAS_GUI
	void dock_concurrent(void);
#endif
};

#endif

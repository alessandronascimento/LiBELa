/*
 * QtWriter.h
 *
 *  Created on: 22/07/2011
 *      Author: Nascimento
 */

#ifndef QTWRITER_H_
#define QTWRITER_H_

#include "../PARSER.h"
#include "../Mol2.h"
#include "../Mol2.h"
#include <QTextEdit>
#include <QString>
#include<QApplication>
#include"zlib.h"

using namespace std;

class QtWriter {
public:
	//! logfile is defined as a C file to generate a written report (log) of the program run.
	FILE *logfile;
	//! Qt TextEdit widget;
	QTextEdit* Editor;
	//! Mol2 Object.
	Mol2 *Cmol;

    gzFile outmol2;



	//! Class constructor. Opens the logfile with the "output" name defined in the input file.
	//! @param Input Pointer to the Parser class.
	QtWriter(PARSER *Input, QTextEdit* Ed);

	//! This function writes down the parameters given by the user as input.
	//! @param Pointer to the Parser class
	void write_params(PARSER *Input);
	void print_param(string p1, double p2, string p3);
	void print_param(string p1, int p2, string p3);
	void print_param(string p1, string p2, string p3);
    void print_param(void);
	//! This function writes a "welcome" message
	void write_welcome(void);
	//! This message shows a char array in the screen and writes it to the logfile.
	//! @param  info char array (98 char long) with a message to be written/shown.
	void write_to_log(char info[98]);
	//! This message shows a char array in the screen and writes it to the logfile.
	//! @param  info char array (98 char long) with a message to be written/shown.
	void print_info(char info[98]);
	//! This function prints a line of stars '*'.
	void print_line(void);
	//! This overloaded function writes a blank line in the logfile and in the screen.
	void write_to_log(void);
	//! This classe writes a gzipped pdb file with dock/MC results.
	//! @param Cmol Mol2 Object.
	//! @param xyz atomic coordinates.
	//! @param energy Energy (kcal/mol) of the accepted pose.
	//! @oaram rmsd RMSD of the accepted pose compared to the original (crystallographic) pose.
	//! @param outname std string used to define the pdb filename.
	void write_pdb(Mol2 *Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname);
	void writeMol2(Mol2 *Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname);
    void writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double energy, double rmsd);
	//! Class destructor. Closes the logfile.
	~QtWriter(void);
    //!
    void write_box(vector<double>center, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z);
};

#endif /* QTWRITER_H_ */

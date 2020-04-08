/*
 * QtWriter.cpp
 *
 *  Created on: 22/07/2011
 *      Author: Nascimento
 */

#include "QtWriter.h"

QtWriter::QtWriter(PARSER *_Input, QTextEdit* Ed) {
	Editor = Ed;
    Input = _Input;
	logfile= fopen((Input->output+".log").c_str(), "w");
	Editor->ensureCursorVisible();

    if (Input->write_mol2 and Input->dock_mode){
        outmol2 = gzopen((Input->output+"_dock.mol2.gz").c_str(), "w");
    }

    this->write_welcome();
    this->write_params();
}

void QtWriter::write_params(){
    fprintf(logfile,"****************************************************************************************************\n");
    fprintf(logfile, "* %-98s *\n", "Parsed parameters:");

	Editor->append(QString("%1\n").arg("Parsed parameters:"));
    fprintf(logfile,"* %98s *\n", "");

    if (Input->dock_mode){
        this->print_param("Mode", "dock", "");
    }
    this->print_param("dock parallel", Input->dock_parallel, "");
    if (Input->dock_parallel){
        this->print_param("parallel jobs", int(Input->parallel_jobs), "");
    }

    this->print_param();

    this->print_param("rec_mol2", Input->rec_mol2, "");
    this->print_param("lig_mol2", Input->lig_mol2, "");
    this->print_param("reflig_mol2", Input->reflig_mol2, "");
    this->print_param("multifile", Input->multifile, "");
    this->print_param("mol2_aa", int(Input->mol2_aa), "");

    this->print_param();

    this->print_param("scoring function", Input->scoring_function, "");
    this->print_param("dielectric_model", Input->dielectric_model, "");
    this->print_param("diel", Input->diel, "");
    this->print_param("solvation sigma", Input->sigma, "");
    this->print_param("solvation_alpha", Input->solvation_alpha, "");
    this->print_param("solvation_beta", Input->solvation_beta, "");
    this->print_param("use grids?", Input->use_grids, "");

    this->print_param();

    this->print_param("only_score?", Input->only_score, "");
    this->print_param("search_box", Input->search_box_x, "");
    this->print_param("minimization_tolerance", Input->min_tol, "");
    this->print_param("minimization delta", Input->min_delta, "");
    this->print_param("dock_min_tol", Input->dock_min_tol, "");
    this->print_param("Minimization timeout", Input->min_timeout, "");
    this->print_param("sort_by_energy", Input->sort_by_energy, "");
    this->print_param("elec_scale", Input->elec_scale, "");
    this->print_param("vdw_scale", Input->vdw_scale, "");
    this->print_param("overlay_optimizer", Input->overlay_optimizer, "");
    this->print_param("energy_optimizer", Input->energy_optimizer, "");
    this->print_param("ignore_h", int(Input->dock_no_h), "");
    this->print_param("deal", Input->deal, "");

    this->print_param();

    this->print_param("generate_conformers", Input->generate_conformers, "");
    if (Input->generate_conformers){
        this->print_param("number_of_conformers", Input->lig_conformers, "");
        this->print_param("conformers_to_rank", Input->conformers_to_evaluate, "");
    }

    this->print_param();

    this->print_param("output_prefix", Input->output, "");
    this->print_param("write_mol2", Input->write_mol2, "");

    this->print_param();

    this->print_line();
	QApplication::processEvents();
}

void QtWriter::print_param(string p1, double p2, string p3){
    Editor->append(QString("%1: %2 %3").arg(QString::fromStdString(p1)).arg(QString::number(p2)).arg(QString::fromStdString(p3)));
    fprintf(logfile, "* %-30s %-46.2f %20.20s*\n", p1.c_str(), p2, p3.c_str());
	QApplication::processEvents();
}

void QtWriter::print_param(string p1, int p2, string p3){
    Editor->append(QString("%1: %2 %3").arg(QString::fromStdString(p1)).arg(QString::number(p2)).arg(QString::fromStdString(p3)));
    fprintf(logfile, "* %-30s %-46d %20s*\n", p1.c_str(), p2, p3.c_str());
	QApplication::processEvents();
}

void QtWriter::print_param(string p1, string p2, string p3){
    Editor->append(QString("%1: %2 %3").arg(QString::fromStdString(p1)).arg(QString::fromStdString(p2)).arg(QString::fromStdString(p3)));
    fprintf(logfile, "* %-30s %-46s %20s*\n", p1.c_str(), p2.c_str(), p3.c_str());
	QApplication::processEvents();
}

void QtWriter::print_param(void){
    Editor->append("");
    fprintf(logfile, "* %-30s %-46s %20s*\n", "", "", "");
    QApplication::processEvents();
}



void QtWriter::write_welcome(void){
	Editor->clear();
    Editor->append("****************************************************************************************************");
    Editor->append("****************************************************************************************************");
    Editor->append("*                  MCLiBELa - Monte Carlo-based Ligand Binding Energy Landscape");
    Editor->append(QString("*                                      Version 0.1  - Build %1").arg(BUILD));
    Editor->append("*                                                                                                  ");
    Editor->append("* University of Sao Paulo                                                                          ");
    Editor->append("* More Info:                                                                                       ");
    Editor->append("*      http://www.biotechmol.ifsc.usp.br/                                                          ");
    Editor->append("*                                                                                                  ");
    Editor->append("****************************************************************************************************");
    Editor->append("****************************************************************************************************");
	Editor->append("");

	fprintf(logfile, "**************************************************************************************\n");
	fprintf(logfile, "                                                                                      \n");
	fprintf(logfile, "                                  iMcLiBELa v. 1.0                                    \n");
	fprintf(logfile, "                    Monte Carlo-Based Ligand Binding Energy Landscape                 \n");
	fprintf(logfile, "                                                                                      \n");
    fprintf(logfile, "                Written by Alessandro S. Nascimento - asnascimento@ifsc.usp.br        \n");
	fprintf(logfile, "                                                                                      \n");
    fprintf(logfile, "                              University of Sao Paulo, Brazil                         \n");
	fprintf(logfile, "                                                                                      \n");
    fprintf(logfile, "                           http://www.biotechmol.ifsc.usp.br/                         \n");
	fprintf(logfile, "                                                                                      \n");
	fprintf(logfile, "**************************************************************************************\n");

	Editor->update();
	QApplication::processEvents();
}

void QtWriter::write_to_log(char info[98]){
	Editor->append(QString("* %1 *").arg(info, -98));
	Editor->update();
    fprintf(logfile, "* %-98.98s *\n", info);
	QApplication::processEvents();
}

void QtWriter::write_to_log(void){
	Editor->append(QString("* %1 *").arg("", -98));
    Editor->update();
    fprintf(logfile, "* %-98.s *\n", "");
	QApplication::processEvents();
}

void QtWriter::print_info(char info[98]){
	Editor->append(QString("* %1").arg(info, -98));
	Editor->update();
    fprintf(logfile, "* %-98.98s *\n", info);
	QApplication::processEvents();
}

void QtWriter::print_line(void){
	Editor->append("****************************************************************************************************************");
	Editor->update();
	fprintf(logfile, "* %-98s *\n", "**************************************************************************************");
	QApplication::processEvents();
}


void QtWriter::write_pdb(Mol2 *Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname){
	gzFile outpdb;
	outpdb = gzopen((outname+".pdb.gz").c_str(), "w");
	gzprintf(outpdb, "MDL\n");
	gzprintf(outpdb, "REMARK\n");
	gzprintf(outpdb, "REMARK %-9s energy = %9.2f rmsd   = %12.3f\n", outname.c_str(), energy, rmsd);
	int i=0;
	unsigned resn=0;

	while (resn < Cmol->residue_pointer.size()-1){
		while(i < Cmol->residue_pointer[resn+1]-1){
			gzprintf(outpdb, "%6s%5d %-4s%1s%3.3s%1s%4d%5s%8.3f%8.3f%8.3f%6.2f%6.2f\n", "ATOM  ", i+1, Cmol->atomnames[i].c_str(), " ",Cmol->resnames[resn].c_str()," ",resn+1," ",xyz[i][0], xyz[i][1], xyz[i][2], 1.0, Cmol->charges[i]);
			i++;
		}
		resn++;
	}

	while (i < Cmol->N){
		gzprintf(outpdb, "%6s%5d %-4s%1s%3.3s%1s%4d%5s%8.3f%8.3f%8.3f%6.2f%6.2f\n", "ATOM  ", i+1, Cmol->atomnames[i].c_str(), " ",Cmol->resnames[resn].c_str()," ",resn+1," ",xyz[i][0], xyz[i][1], xyz[i][2], 1.0, Cmol->charges[i]);
		i++;
	}
	gzprintf(outpdb, "TER\n");
	gzprintf(outpdb, "ENDMDL\n");
	gzclose(outpdb);
	QApplication::processEvents();
}

void QtWriter::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double energy, double rmsd, string outname){
	gzFile outmol2;
	outmol2 = gzopen((outname+".mol2.gz").c_str(), "a");
	gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
	gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
	gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
	gzprintf(outmol2, "SMALL\n");
	gzprintf(outmol2, "USER_CHARGES\n");
	gzprintf(outmol2, "Energy: %7.3f\n", energy);
	gzprintf(outmol2, "RMSD: %7.3f\n", rmsd);
	gzprintf(outmol2, "@<TRIPOS>ATOM\n");
	int i=0;
	unsigned resn=0;

	while(resn < Cmol->residue_pointer.size()-1){
		while(i < Cmol->residue_pointer[resn+1]-1){
			if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
				gzprintf(outmol2, "%7d %-3.3s       %9.4f %9.4f %9.4f %3.3s %5d %5s %13.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			else {
				gzprintf(outmol2, "%7d %-3.3s       %9.4f %9.4f %9.4f %3.3s %5d %5s %13.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			i++;
		}
		resn++;
	}
	while(i < Cmol->N){
		if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
			gzprintf(outmol2, "%7d %-3.3s       %9.4f %9.4f %9.4f %3.3s %5d %5s %13.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		else{
			gzprintf(outmol2, "%7d %-3.3s       %9.4f %9.4f %9.4f %3.3s %5d %5s %13.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		i++;
	}

	gzprintf(outmol2, "@<TRIPOS>BOND\n");

	for (unsigned j=0; j<Cmol->bonds.size(); j++){
		gzprintf(outmol2, "%6d%6s%6s%5s\n", j+1, Cmol->bonds[j][0].c_str(),Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
	}
	gzclose(outmol2);
	QApplication::processEvents();
}

void QtWriter::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double energy, double rmsd){
    gzprintf(outmol2, "\n");
    gzprintf(outmol2, "########## %15.15s: %19.19s\n", "Name", Cmol->molname.c_str());
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "Energy Score", energy);
    gzprintf(outmol2, "########## %15.15s: % 19.6f\n", "RMSD", rmsd);
    gzprintf(outmol2, "\n");
    gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
    gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
    gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
    gzprintf(outmol2, "SMALL\n");
    gzprintf(outmol2, "USER_CHARGES\n");
    gzprintf(outmol2, "Energy: %7.3f\n", energy);
    gzprintf(outmol2, "RMSD/OVERLAY: %7.3f\n", rmsd);
    gzprintf(outmol2, "@<TRIPOS>ATOM\n");
    int i=0;
    unsigned resn=0;

    if (int(xyz.size()) != Cmol->N){
        printf("Mismatch in atom number while writting mol2 file. Please check!\n");
        printf("Cmol->N = %d    xyz.size() = %d\n", Cmol->N, int(xyz.size()));
        exit(1);
    }

    while(resn < Cmol->residue_pointer.size()-1){
        while(i < Cmol->residue_pointer[resn+1]-1){
            if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
            }
            else {
                gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
            }
            i++;
        }
        resn++;
    }
    while(i < Cmol->N){
        if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
        }
        else{
            gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s%4d %18.18s %9.4f \n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
        }
        i++;
    }

    gzprintf(outmol2, "@<TRIPOS>BOND\n");

    for (unsigned j=0; j<Cmol->bonds.size(); j++){
        gzprintf(outmol2, "%6d%6.6s%6.6s%5.5s\n", j+1, Cmol->bonds[j][0].c_str(), Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
    }
}

QtWriter::~QtWriter(void){
	fprintf(logfile,"**************************************************************************************\n");
    if (Input->write_mol2 and Input->dock_mode){
        gzclose(outmol2);
    }
    fclose(logfile);
	QApplication::processEvents();
}

void QtWriter::write_box(vector<double>center, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z){
    FILE *box;
    box = fopen("box.pdb", "w");
    fprintf (box, "REMARK    CENTER OF THE BOX  %2.3f  %2.3f  %2.3f\n", center[0], center[1], center[2]);
    fprintf (box, "ATOM      1  DUA BOX     1     % 7.3f % 7.3f % 7.3f\n", min_x, min_y, min_z);
    fprintf (box, "ATOM      2  DUB BOX     1     % 7.3f % 7.3f % 7.3f\n", max_x, min_y, min_z);
    fprintf (box, "ATOM      3  DUC BOX     1     % 7.3f % 7.3f % 7.3f\n", max_x, min_y, max_z);
    fprintf (box, "ATOM      4  DUD BOX     1     % 7.3f % 7.3f % 7.3f\n", min_x, min_y, max_z);
    fprintf (box, "ATOM      5  DUE BOX     1     % 7.3f % 7.3f % 7.3f\n", min_x, max_y, min_z);
    fprintf (box, "ATOM      6  DUF BOX     1     % 7.3f % 7.3f % 7.3f\n", max_x, max_y, min_z);
    fprintf (box, "ATOM      7  DUG BOX     1     % 7.3f % 7.3f % 7.3f\n", max_x, max_y, max_z);
    fprintf (box, "ATOM      8  DUH BOX     1     % 7.3f % 7.3f % 7.3f\n", min_x, max_y, max_z);
    fprintf (box, "CONECT    1    2    4    5\n");
    fprintf (box, "CONECT    2    1    3    6\n");
    fprintf (box, "CONECT    3    2    4    7\n");
    fprintf (box, "CONECT    4    1    3    8\n");
    fprintf (box, "CONECT    5    1    6    8\n");
    fprintf (box, "CONECT    6    2    5    7\n");
    fprintf (box, "CONECT    7    3    6    8\n");
    fprintf (box, "CONECT    8    4    5    7\n");
    fprintf (box, "END\n");

    fclose(box);
}

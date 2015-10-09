/*
 * ENERGY.cpp
 *
 *  Created on: 23/03/2010
 *      Author: UFABC
 */

#include "ENERGY.h"

using namespace std;


ENERGY::ENERGY(){
}


double ENERGY::distance(double x1, double x2, double y1, double y2, double z1, double z2) {
	return ( sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) ); }

double ENERGY::distance_squared (double x1, double x2, double y1, double y2, double z1, double z2) {
	return (((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))); }


double ENERGY::compute_elec(vector<double> rec_charges, vector<double> lig_charges, vector<double> rec_crd, vector<double> lig_crd ){
	elec=0.00;
    for (unsigned i=0; i<(rec_charges.size()); i++){
        for (unsigned j=0; j<(lig_charges.size()); j++) {
			elec+= (rec_charges[i]*lig_charges[j])/(distance(rec_crd[3*i], lig_crd[3*j],rec_crd[(3*i)+1], lig_crd[(3*j)+1], rec_crd[(3*i)+2], lig_crd[(3*j)+2])); } }
	return (elec);
}

double ENERGY::compute_vdw(vector<double> rec_epsilons, vector<double> lig_epsilons, vector<double> rec_radii, vector<double> lig_radii, vector<double> rec_crd, vector<double> lig_crd){
	vdw=0.00;
    for (unsigned i=0; i<(rec_epsilons.size()); i++){
        for (unsigned j=0; j<(lig_epsilons.size()); j++){
			rij = rec_radii[i]+lig_radii[j];
			eij = sqrt(rec_epsilons[i]*lig_epsilons[j]);
			dij2 = distance_squared(rec_crd[3*i], lig_crd[3*j],rec_crd[(3*i)+1], lig_crd[(3*j)+1], rec_crd[(3*i)+2], lig_crd[(3*j)+2]);
			bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
			acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
			vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));
		}
	}
	return (vdw);
}

double ENERGY::compute_ene(vector<double> rec_epsilons, vector<double> lig_epsilons, vector<double> rec_radii, vector<double> lig_radii, vector<double> rec_charges, vector<double> lig_charges, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd){
	vdw=0.00;
	elec=0.00;
    for (unsigned i=0; i<(rec_charges.size()); i++){
        for (unsigned j=0; j<(lig_charges.size()); j++){
			rij = rec_radii[i]+lig_radii[j];
			eij = sqrt(rec_epsilons[i]*lig_epsilons[j]);
			dij = distance(rec_crd[i][0], lig_crd[j][0],rec_crd[i][1], lig_crd[j][1], rec_crd[i][2], lig_crd[j][2]);
            elec+= (332.0*rec_charges[i]*lig_charges[j])/dij;
			dij2 = dij*dij;
			bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
			acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
			vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));
		}
	}
//	printf("ENERGY DECOMP:: Elec: %.4f VDW: %.4f Total: %.4f\n", elec, vdw, elec+vdw);
	return (vdw+elec);
}

double ENERGY::compute_ene_solvation(vector<double> rec_epsilons, vector<double> lig_epsilons, vector<double> rec_radii, vector<double> lig_radii, vector<double> rec_charges, vector<double> lig_charges, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double deltaij_es6, double sigma){
	vdw=0.00;
	elec=0.00;
	rec_solv=0.0;
	lig_solv=0.0;
    for (unsigned i=0; i<(rec_charges.size()); i++){
		rec_solv_affinity = (0.25*((rec_charges[i])*(rec_charges[i])))-0.005;
        for (unsigned j=0; j<(lig_charges.size()); j++){
			rij = rec_radii[i]+lig_radii[j];
			eij = sqrt(rec_epsilons[i]*lig_epsilons[j]);
			dij = distance(rec_crd[i][0], lig_crd[j][0],rec_crd[i][1], lig_crd[j][1], rec_crd[i][2], lig_crd[j][2]);
            elec+= (332.0*rec_charges[i]*lig_charges[j])/(dij*dij);
			dij2 = dij*dij;
			bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
			acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
			vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));

			deff = pow(((dij2*dij2*dij2)+deltaij_es6), (1.0/3.0));
			lig_solv_affinity = (0.25*((lig_charges[j])*(lig_charges[j])))-0.005;
			lig_solv_distf = ((4.0/3.0)*PI*(rec_radii[i]*rec_radii[i]*rec_radii[i]));
			lig_solv_distf = (lig_solv_distf * exp(-deff/(2*(sigma*sigma))))/(sigma*sigma*sigma);
			rec_solv_distf = ((4.0/3.0)*PI*(lig_radii[j]*lig_radii[j]*lig_radii[j]));
			rec_solv_distf = (rec_solv_distf * exp(-deff/(2*(sigma*sigma))))/(sigma*sigma*sigma);
			rec_solv+= rec_solv_affinity*rec_solv_distf;
			lig_solv+= lig_solv_affinity*lig_solv_distf;
		}
	}
//	printf("Elec: %7.3f   VDW: %7.3f    Solv: %7.3f\n", elec, vdw, rec_solv+lig_solv);
	return (vdw+elec+rec_solv+lig_solv);
}

// Soft-core part

double ENERGY::compute_vdw_softcore(vector<double> rec_epsilons, vector<double> lig_epsilons, vector<double> rec_radii, vector<double> lig_radii, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double deltaij6){
	vdw=0.00;
    for (unsigned i=0; i<(rec_epsilons.size()); i++){
        for (unsigned j=0; j<(lig_epsilons.size()); j++){
			rij = rec_radii[i]*lig_radii[j];
			eij = sqrt(rec_epsilons[i]*lig_epsilons[j]);
			dij2 = distance_squared(rec_crd[i][0], lig_crd[j][0],rec_crd[i][1], lig_crd[j][1], rec_crd[i][2], lig_crd[j][2]);
			acoef = (4096.00)*eij*(rij*rij*rij*rij*rij*rij);
			bcoef = (128.00)*eij*(rij*rij*rij);
			vdw+= ((acoef/(((dij2*dij2*dij2)+deltaij6)*((dij2*dij2*dij2)+deltaij6))) - (bcoef/((dij2*dij2*dij2)+deltaij6)));
		}
	}
	return (vdw);
}

double ENERGY::compute_elec_softcore(vector<double> rec_charges, vector<double> lig_charges, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double diel, double deltaij_es6 ){
	elec=0.00;
    for (unsigned i=0; i<(rec_charges.size()); i++){
        for (unsigned j=0; j<(lig_charges.size()); j++) {
			dij2 = distance_squared(rec_crd[i][0], lig_crd[j][0],rec_crd[i][1], lig_crd[j][1], rec_crd[i][2], lig_crd[j][2]);
			elec+= (332.0*rec_charges[i]*lig_charges[j])/(diel*(pow(((dij2*dij2*dij2)+deltaij_es6),(1.0/3.0))));
		}
	}
	return (elec);
}

double ENERGY::compute_solvation(vector<double>rec_charges, vector<double> lig_charges, vector<double>rec_radii, vector<double>lig_radii, vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double deltaij_es6, double sigma){
	rec_solv = 0.00;
	lig_solv = 0.00;
	for (unsigned i=0; i< rec_charges.size();i++){
		rec_solv_affinity = (0.25*((rec_charges[i])*(rec_charges[i])))-0.005;
		for (unsigned j=0; j< lig_charges.size(); j++){
			dij2 = distance_squared(rec_crd[i][0], lig_crd[j][0],rec_crd[i][1], lig_crd[j][1], rec_crd[i][2], lig_crd[j][2]);
			deff = pow(((dij2*dij2*dij2)+deltaij_es6), (1.0/3.0));
			lig_solv_affinity = (0.25*((lig_charges[j])*(lig_charges[j])))-0.005;
			lig_solv_distf = ((4.0/3.0)*PI*(rec_radii[i]*rec_radii[i]*rec_radii[i]));
			lig_solv_distf = (lig_solv_distf * exp(-deff/(2*(sigma*sigma))))/(sigma*sigma*sigma);
			rec_solv_distf = ((4.0/3.0)*PI*(lig_radii[j]*lig_radii[j]*lig_radii[j]));
			rec_solv_distf = (rec_solv_distf * exp(-deff/(2*(sigma*sigma))))/(sigma*sigma*sigma);
			rec_solv+= rec_solv_affinity*rec_solv_distf;
			lig_solv+= lig_solv_affinity*lig_solv_distf;
		}
	}
//	printf("Receptor Solvation Energy: %.4f  Ligand Solvation Energy: %.4f \n", rec_solv, lig_solv);
	return(rec_solv+lig_solv);
}

double ENERGY::compute_ene_softcore_solvation(vector<double>rec_charges, vector<double> lig_charges, vector<double>rec_radii, vector<double>lig_radii,vector<double> rec_epsilons, vector<double> lig_epsilons,vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double diel, double deltaij_es6, double deltaij6, double sigma){
	rec_solv = 0.00;
	lig_solv = 0.00;
	vdw=0.0;
	elec=0.0;
    for (unsigned i=0; i<(rec_charges.size()); i++){
		rec_solv_affinity = (0.25*((rec_charges[i])*(rec_charges[i])))-0.005;
		for (unsigned j=0; j< lig_charges.size(); j++){

			//! electrostactic energy (softcore)

			dij2 = this->distance_squared(rec_crd[i][0], lig_crd[j][0],rec_crd[i][1], lig_crd[j][1], rec_crd[i][2], lig_crd[j][2]);
			deff = pow(((dij2*dij2*dij2)+deltaij_es6), (1.0/3.0));
			elec+= (332.0*rec_charges[i]*lig_charges[j])/(diel*deff);

			//! VDW energy (softcore)

			rij = rec_radii[i]*lig_radii[j];
			eij = sqrt(rec_epsilons[i]*lig_epsilons[j]);
			acoef = (4096.00)*eij*(rij*rij*rij*rij*rij*rij);
			bcoef = (128.00)*eij*(rij*rij*rij);
			vdw+= ((acoef/(((dij2*dij2*dij2)+deltaij6)*((dij2*dij2*dij2)+deltaij6))) - (bcoef/((dij2*dij2*dij2)+deltaij6)));

			//! Solvation energy

			lig_solv_affinity = (0.25*((lig_charges[j])*(lig_charges[j])))-0.005;
			lig_solv_distf = ((4.0/3.0)*PI*(rec_radii[i]*rec_radii[i]*rec_radii[i]));
			lig_solv_distf = (lig_solv_distf * exp(-deff/(2*(sigma*sigma))))/(sigma*sigma*sigma);
			rec_solv_distf = ((4.0/3.0)*PI*(lig_radii[j]*lig_radii[j]*lig_radii[j]));
			rec_solv_distf = (rec_solv_distf * exp(-deff/(2*(sigma*sigma))))/(sigma*sigma*sigma);

			rec_solv+= rec_solv_affinity*rec_solv_distf;
			lig_solv+= lig_solv_affinity*lig_solv_distf;
		}
	}
//	printf("ELEC:%7.3f   VDW:%7.3f   Rec_Solv:%7.3f   Lig_Solv:%7.3f\n", elec, vdw, rec_solv, lig_solv);
	return(elec+vdw+rec_solv+lig_solv);
}

double ENERGY::compute_ene_softcore(vector<double>rec_charges, vector<double> lig_charges, vector<double>rec_radii, vector<double>lig_radii,vector<double> rec_epsilons, vector<double> lig_epsilons,vector<vector<double> >rec_crd, vector<vector<double> >lig_crd, double diel, double deltaij_es6, double deltaij6){
	vdw=0.0;
	elec=0.0;
    for (unsigned i=0; i<(rec_charges.size()); i++){
		rec_solv_affinity = (0.25*((rec_charges[i])*(rec_charges[i])))-0.005;
		for (unsigned j=0; j< lig_charges.size(); j++){

			//! electrostactic energy (softcore)

			dij2 = distance_squared(rec_crd[i][0], lig_crd[j][0],rec_crd[i][1], lig_crd[j][1], rec_crd[i][2], lig_crd[j][2]);
			deff = pow(((dij2*dij2*dij2)+deltaij_es6), (1.0/3.0));
			elec+= (332.0*rec_charges[i]*lig_charges[j])/(diel*deff);

			//! VDW energy (softcore)

			rij = rec_radii[i]*lig_radii[j];
			eij = sqrt(rec_epsilons[i]*lig_epsilons[j]);
			acoef = (4096.00)*eij*(rij*rij*rij*rij*rij*rij);
			bcoef = (128.00)*eij*(rij*rij*rij);
			vdw+= ((acoef/(((dij2*dij2*dij2)+deltaij6)*((dij2*dij2*dij2)+deltaij6))) - (bcoef/((dij2*dij2*dij2)+deltaij6)));
		}
	}
//	printf("Ele: %7.3f    VDW: %7.3f\n", elec, vdw);
	return(elec+vdw);
}




double ENERGY::compute_ene_softcore_solvation(Mol2* REC, Mol2* LIG,PARSER* Input){
	return(this->compute_ene_softcore_solvation(REC->charges,LIG->charges,REC->radii,LIG->radii,REC->epsilons,LIG->epsilons,REC->xyz,LIG->xyz, Input->diel, Input->deltaij_es6, Input->deltaij6, Input->sigma));
}

double ENERGY::compute_ene_softcore_solvation(Mol2* REC, Mol2* LIG,PARSER* Input, vector<vector<double> >coord){
	return(this->compute_ene_softcore_solvation(REC->charges,LIG->charges,REC->radii,LIG->radii,REC->epsilons,LIG->epsilons,REC->xyz,coord, Input->diel, Input->deltaij_es6, Input->deltaij6, Input->sigma));
}

double ENERGY::compute_ene_softcore_solvation(Mol2* REC, Mol2* LIG,PARSER* Input, RAND *Rand){
	return(this->compute_ene_softcore_solvation(REC->charges,LIG->charges,REC->radii,LIG->radii,REC->epsilons,LIG->epsilons,REC->xyz,LIG->new_mcoords[Rand->lign], Input->diel, Input->deltaij_es6, Input->deltaij6, Input->sigma));
}

double ENERGY::compute_ene_softcore_solvation(Mol2* REC, vector<vector<double> >coord, Mol2* LIG,PARSER* Input, RAND *Rand){
	return(this->compute_ene_softcore_solvation(REC->charges,LIG->charges,REC->radii,LIG->radii,REC->epsilons,LIG->epsilons, coord,LIG->new_mcoords[Rand->lign], Input->diel, Input->deltaij_es6, Input->deltaij6, Input->sigma));
}

double ENERGY::compute_ene_softcore_solvation(Mol2* REC, vector<vector<double> >coord, vector<vector<double> >lig_coords, Mol2* LIG,PARSER* Input){
	return(this->compute_ene_softcore_solvation(REC->charges,LIG->charges,REC->radii,LIG->radii,REC->epsilons,LIG->epsilons, coord, lig_coords, Input->diel, Input->deltaij_es6, Input->deltaij6, Input->sigma));
}


double ENERGY::compute_ene(Mol2* REC, vector<vector<double> >coords, Mol2* LIG, RAND* Rand){
	return(this->compute_ene(REC->epsilons, LIG->epsilons, REC->radii,  LIG->radii, REC->charges, LIG->charges, coords, LIG->new_mcoords[Rand->lign]));
}

double ENERGY::compute_ene(Mol2* REC, Mol2* LIG,RAND *Rand){
	return(this->compute_ene(REC->epsilons, LIG->epsilons,REC->radii, LIG->radii, REC->charges, LIG->charges, REC->xyz, LIG->new_mcoords[Rand->lign]));
}

double ENERGY::compute_ene(Mol2* REC, Mol2* LIG,vector<vector<double> >coord){
	return(this->compute_ene(REC->epsilons, LIG->epsilons, REC->radii, LIG->radii, REC->charges, LIG->charges, REC->xyz, coord));
}

double ENERGY::compute_ene(Mol2* REC, vector<vector<double> >coord, vector<vector<double> >lig_coords, Mol2* LIG){
	return(this->compute_ene(REC->epsilons, LIG->epsilons, REC->radii, LIG->radii, REC->charges, LIG->charges, coord, lig_coords));
}

double ENERGY::compute_ene_from_grids(Grid* Grids, Mol2* Lig){
	double elec =0.0, vdwA=0.0, vdwB = 0.00, rec_solv=0.00, lig_solv=0.00, lig_affinity;
	int a, b, c;
	double af, bf, cf;
	for (int i=0; i< Lig->N; i++){
		af = (Lig->xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
		bf = (Lig->xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
		cf = (Lig->xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
		a = round(af);
		b = round(bf);
		c = round(cf);

		if (a > 0 and b > 0 and c > 0 and a < int(Grids->elec_grid.size()) and b < int(Grids->elec_grid[0].size()) and c < int(Grids->elec_grid[0][0].size())){

			elec += Lig->charges[i]* Grids->elec_grid[a][b][c];

			vdwA += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 6) * Grids->vdwA_grid[a][b][c];
			vdwB += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 3) * Grids->vdwB_grid[a][b][c];

			lig_affinity = (0.25*Lig->charges[i]*Lig->charges[i]) - 0.005;
            rec_solv += Grids->solv_gauss[a][b][c]* (4.0/3.0) * PI * pow(Lig->radii[i], 3);
			lig_solv += lig_affinity*Grids->rec_solv_gauss[a][b][c];
		}
		else {
#ifdef DEBUG
			printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
		}
	}
	return(elec+vdwA-vdwB+rec_solv+lig_solv);
}

double ENERGY::compute_ene_from_grids(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
	double elec =0.0, vdwA=0.0, vdwB = 0.00, rec_solv=0.00, lig_solv=0.00, lig_affinity;
	int a, b, c;
	double af, bf, cf;
	for (int i=0; i< Lig->N; i++){
		af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
		bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
		cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
		a = round(af);
		b = round(bf);
		c = round(cf);

		if (a > 0 and b > 0 and c > 0 and a < int(Grids->elec_grid.size()) and b < int(Grids->elec_grid[0].size()) and c< int(Grids->elec_grid[0][0].size())){

			elec += Lig->charges[i]* Grids->elec_grid[a][b][c];

			vdwA += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 6) * Grids->vdwA_grid[a][b][c];
			vdwB += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 3) * Grids->vdwB_grid[a][b][c];

			lig_affinity = (0.25*Lig->charges[i]*Lig->charges[i]) - 0.005;
			rec_solv += Grids->solv_gauss[a][c][c]* (4.0/3.0) * PI * pow(Lig->radii[i], 3);
			lig_solv += lig_affinity*Grids->rec_solv_gauss[a][b][c];
		}
		else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
#ifdef DEBUG
			printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
		}
	}
//	printf("Elec: %7.3f VDW: %7.3f Rec_Solv: %7.3f Lig_solv:%7.3f\n", elec, vdwA-vdwB, rec_solv, lig_solv);
	return(elec+vdwA-vdwB+rec_solv+lig_solv);
}

double ENERGY::compute_ene_from_grids_nosolv(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
	double elec =0.0, vdwA=0.0, vdwB = 0.00;
	int a, b, c;
	double af, bf, cf;
	for (int i=0; i< Lig->N; i++){
		af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
		bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
		cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
		a = round(af);
		b = round(bf);
		c = round(cf);

		if (a > 0 and b > 0 and c > 0 and a < int(Grids->elec_grid.size()) and b < int(Grids->elec_grid[0].size()) and c< int(Grids->elec_grid[0][0].size())){

			elec += Lig->charges[i]* Grids->elec_grid[a][b][c];
			vdwA += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 6) * Grids->vdwA_grid[a][b][c];
			vdwB += Lig->epsilons_sqrt[i] * pow(Lig->radii[i], 3) * Grids->vdwB_grid[a][b][c];
		}
		else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
#ifdef DEBUG
			printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
		}
	}
//	printf("Elec: %7.3f VDW: %7.3f Rec_Solv: %7.3f Lig_solv:%7.3f\n", elec, vdwA-vdwB, rec_solv, lig_solv);
	return(elec+vdwA-vdwB);
}

double ENERGY::compute_ene_from_grids_hardcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
	double elec =0.0, vdwA=0.0, vdwB = 0.00, rec_solv=0.00, lig_solv=0.00, lig_affinity;
	int a, b, c;
	double af, bf, cf;
	for (int i=0; i< Lig->N; i++){
		af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
		bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
		cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
		a = round(af);
		b = round(bf);
		c = round(cf);

		if (a > 0 and b > 0 and c > 0 and a < int(Grids->elec_grid.size()) and b < int(Grids->elec_grid[0].size()) and c< int(Grids->elec_grid[0][0].size())){

			elec += Lig->charges[i]* Grids->elec_grid[a][b][c];

			vdwA += sqrt(Lig->epsilons[i]*pow((2*Lig->radii[i]), 12))* Grids->vdwA_grid[a][b][c];
			vdwB += sqrt(2 * Lig->epsilons[i]*pow((2*Lig->radii[i]), 6)) * Grids->vdwB_grid[a][b][c];

			lig_affinity = (0.25*Lig->charges[i]*Lig->charges[i]) - 0.005;
            rec_solv += Grids->solv_gauss[a][b][c]* (4.0/3.0) * PI * pow(Lig->radii[i], 3);
			lig_solv += lig_affinity*Grids->rec_solv_gauss[a][b][c];
		}
		else {
            elec += 999999.9;
            vdwA += 999999.9;
            vdwB += 999999.9;
#ifdef DEBUG
			printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
		}
	}
//	printf("Elec: %7.3f VDW: %7.3f Rec_Solv: %7.3f Lig_solv:%7.3f\n", elec, vdwA-vdwB, rec_solv, lig_solv);
	return(elec+(vdwA-vdwB)+rec_solv+lig_solv);
}

double ENERGY::compute_ene_from_grids_hardcore_nosolv(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz){
	double elec =0.0, vdwA=0.0, vdwB = 0.00;
	int a, b, c;
	double af, bf, cf;
	for (int i=0; i< Lig->N; i++){
		af = (xyz[i][0] - Grids->xbegin)/Grids->grid_spacing;
		bf = (xyz[i][1] - Grids->ybegin)/Grids->grid_spacing;
		cf = (xyz[i][2] - Grids->zbegin)/Grids->grid_spacing;
		a = round(af);
		b = round(bf);
		c = round(cf);

		if (a > 0 and b > 0 and c > 0 and a < int(Grids->elec_grid.size()) and b < int(Grids->elec_grid[0].size()) and c< int(Grids->elec_grid[0][0].size())){

			elec += Lig->charges[i]* Grids->elec_grid[a][b][c];

			vdwA += sqrt(Lig->epsilons[i]*pow((2*Lig->radii[i]), 12))* Grids->vdwA_grid[a][b][c];
			vdwB += sqrt(2 * Lig->epsilons[i]*pow((2*Lig->radii[i]), 6)) * Grids->vdwB_grid[a][b][c];
		}
		else {
			elec += 999999.9;
			vdwA += 999999.9;
			vdwB += 999999.9;
#ifdef DEBUG
			printf("Ligand %s is slipping from computation box. Inaccurate energies will be estimated.\n", Lig->molname.c_str());
#endif
		}
	}
//	printf("Elec: %7.3f VDW: %7.3f\n", elec, vdwA-vdwB);
	return(elec+vdwA-vdwB);
}

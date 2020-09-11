/*
 * Grid.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: asn
 */

#include "Grid.h"
#include "iMcLiBELa.h"

#define MAX_ALJ 100.0
#define MAX_BLJ 50.0

using namespace std;

double Grid::distance(double x1, double x2, double y1, double y2, double z1, double z2) {
	return ( sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}

double Grid::distance_squared(double x1, double x2, double y1, double y2, double z1, double z2) {
	return ((((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}


Grid::Grid(PARSER* _Input, WRITER* _Writer) {
    this->Input = _Input;
    this->Writer = _Writer;
	grid_spacing = Input->grid_spacing;
    this->pbsa_loaded = false;
    this->delphi_loaded = false;
}

Grid::Grid(PARSER* _Input, WRITER* _Writer, Mol2* Rec, vector<double> com){
	this->Input = _Input;
    this->Writer = _Writer;
    this->pbsa_loaded = false;
    this->delphi_loaded = false;
    if (Input->pbsa_grid != "" and Input->use_pbsa){
        this->load_Ambergrids_from_file();
    }
    else if (Input->delphi_grid != "" and Input->use_delphi){
        this->load_phimap_from_file(Input->delphi_gsize);
    }
    else if (Input->delphi_cube_grid != "" and Input->use_delphi){
        if (Input->delphi_cube_grid.substr(Input->delphi_cube_grid.size()-3, 3) == ".gz"){
            this->load_delphi_gzcube();
        }
        else{
            this->load_delphi_cube();
        }
    }
    else {
        this->grid_spacing = Input->grid_spacing;
        this->generate_points(com);
    }

    switch (Input->scoring_function) {
    case 0:
        if (Input->dock_parallel and Input->parallel_jobs > 1){
            this->compute_grid_softcore_omp(Rec);
        }
        else {
            this->compute_grid_softcore(Rec);
        }
        break;
    case 1:
        if (Input->dock_parallel and Input->parallel_jobs > 1){
            this->compute_grid_softcore_omp(Rec);
        }
        else {
            this->compute_grid_softcore(Rec);
        }
        break;
    case 2:
        if (Input->dock_parallel and Input->parallel_jobs > 1){
            this->compute_grid_hardcore_omp(Rec);
        }
        else{
            this->compute_grid_hardcore(Rec);
        }
        break;
    case 3:
        if (Input->dock_parallel and Input->parallel_jobs > 1){
            this->compute_grid_hardcore_omp(Rec);
        }
        else{
            this->compute_grid_hardcore(Rec);
        }
        break;
    case 4:
        if (Input->dock_parallel and Input->parallel_jobs > 1){
            this->compute_grid_hardcore_omp_Gaussian(Rec);
        }
        else{
            this->compute_grid_hardcore_Gaussian(Rec);
        }
        break;
    case 5:
        if (Input->dock_parallel and Input->parallel_jobs > 1){
            this->compute_grid_hardcore_omp_Gaussian(Rec);
        }
        else{
            this->compute_grid_hardcore_Gaussian(Rec);
        }
        break;
    }

	if (Input->write_grids){
		this->write_grids_to_file();
	}
}

void Grid::generate_points(vector<double> ref_com){
    this->npointsx = round(Input->x_dim/this->grid_spacing);
    this->npointsy = round(Input->y_dim/this->grid_spacing);
    this->npointsz = round(Input->z_dim/this->grid_spacing);
	this->xbegin = ref_com[0] - Input->x_dim/2;
	this->ybegin = ref_com[1] - Input->y_dim/2;
	this->zbegin = ref_com[2] - Input->z_dim/2;
	this->xend = (npointsx*grid_spacing)+xbegin;
	this->yend = (npointsy*grid_spacing)+ybegin;
	this->zend = (npointsz*grid_spacing)+zbegin;

#ifdef DEBUG
	printf("Npoints: %d %d %d\n", npointsx, npointsy, npointsz);
	printf("Beginning of grid coordinates: %7.3f %7.3f %7.3f\n", xbegin, ybegin, zbegin);
	printf("Grid box goes until %7.3f %7.3f %7.3f\n", xend, yend, zend);
#endif

}

void Grid::compute_grid_softcore(Mol2* Rec){
	if (Input->scoring_function < 2){
		vector<double> elec_t1(npointsz), vdwA_t1(npointsz), vdwB_t1(npointsz), solv_t1(npointsz),rec_solv_t1(npointsz);
		vector<vector<double> > elec_t2, vdwA_t2, vdwB_t2, solv_t2, rec_solv_t2;

        double elec, d, d2, d3, x, y, z, vdwA, vdwB, denom, d6, solv, rec_solv;

		for(int a=0; a< this->npointsx; a++){
			x = (a*grid_spacing) + this->xbegin;

			for (int b=0; b< this->npointsy; b++){
				y = (b*this->grid_spacing) + this->ybegin;
				for (int c=0; c<this->npointsz; c++){
					z = (c*this->grid_spacing) + this->zbegin;
					elec = 0.0;
					vdwA = 0.0;
					vdwB = 0.0;
					solv=0.0;
					rec_solv=0.0;

					for (int i=0; i< Rec->N; i++){
                        d = this->distance(x, Rec->xyz[i][0], y,Rec->xyz[i][1], z, Rec->xyz[i][2]);
                        d2 = d*d;
                        d3 = d*d*d;
						d6 = d2*d2*d2;

                        denom = pow((d3 + Input->deltaij_es3), (1.0/3.0));

                        if (Input->dielectric_model == "constant"){
                            elec += (332.0*Rec->charges[i])/(Input->diel*denom);
                            denom = pow((d6 + Input->deltaij_es6), (1.0/3.0));
                            solv += ((Input->solvation_alpha * Rec->charges[i] * Rec->charges[i])+ Input->solvation_beta) *  exp((-denom/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                            rec_solv += (4.0/3.0) * PI * pow(Rec->radii[i], 3) * exp((-denom/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                        }

                        else {                      // dielectric model == r
                            denom = pow((d6 + Input->deltaij_es6), (1.0/3.0));
                            elec += (332.0* Rec->charges[i]/(Input->diel*denom));
                            solv += ((Input->solvation_alpha * Rec->charges[i] * Rec->charges[i])+ Input->solvation_beta) *  exp((-denom/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                            rec_solv += (4.0/3.0) * PI * pow(Rec->radii[i], 3) * exp((-denom/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                        }

						denom = (d6 + Input->deltaij6);
                        vdwA += (4096.0 * Rec->epsilons_sqrt[i] * pow(Rec->radii[i], 6)) / (denom*denom);
                        vdwB += ( 128.0 * Rec->epsilons_sqrt[i] * pow(Rec->radii[i], 3)) / denom;
					}
					elec_t1[c] = (elec);
					vdwA_t1[c] = (vdwA);
					vdwB_t1[c] = (vdwB);
					solv_t1[c] = (solv);
					rec_solv_t1[c] = (rec_solv);
				}
				elec_t2.push_back(elec_t1);
				vdwA_t2.push_back(vdwA_t1);
				vdwB_t2.push_back(vdwB_t1);
				solv_t2.push_back(solv_t1);
				rec_solv_t2.push_back(rec_solv_t1);

			}
			this->elec_grid.push_back(elec_t2);
			this->vdwA_grid.push_back(vdwA_t2);
			this->vdwB_grid.push_back(vdwB_t2);
			this->solv_gauss.push_back(solv_t2);
			this->rec_solv_gauss.push_back(rec_solv_t2);
			elec_t2.clear();
			vdwA_t2.clear();
			vdwB_t2.clear();
			solv_t2.clear();
			rec_solv_t2.clear();
		}

		this->rec_si = 0.00;
		for(int i=0; i<Rec->N; i++){
            this->rec_si += (Input->solvation_alpha*Rec->charges[i]*Rec->charges[i]) + Input->solvation_beta;
		}
	}
	else {
		this->compute_grid_hardcore(Rec);
	}
}

void Grid::compute_grid_softcore_omp(Mol2* Rec){
    if (Input->scoring_function < 2){
        vector<double> elec_t1(npointsz), vdwA_t1(npointsz), vdwB_t1(npointsz), solv_t1(npointsz),rec_solv_t1(npointsz);
        vector<vector<double> > elec_t2, vdwA_t2, vdwB_t2, solv_t2, rec_solv_t2;

        // initializing the vectors;

        for (int i=0; i<this->npointsy; i++){
            elec_t2.push_back(elec_t1);
            vdwA_t2.push_back(vdwA_t1);
            vdwB_t2.push_back(vdwB_t1);
            solv_t2.push_back(solv_t1);
            rec_solv_t2.push_back(rec_solv_t1);
        }

        for (int i=0; i<this->npointsx; i++){
            this->elec_grid.push_back(elec_t2);
            this->vdwA_grid.push_back(vdwA_t2);
            this->vdwB_grid.push_back(vdwB_t2);
            this->solv_gauss.push_back(solv_t2);
            this->rec_solv_gauss.push_back(rec_solv_t2);
        }

        // Now, starting OMP...
#pragma omp parallel num_threads(Input->parallel_jobs)
        {
#pragma omp for schedule(static, 1)
            for(int a=0; a< this->npointsx; a++){
                double x = (a*grid_spacing) + this->xbegin;

                for (int b=0; b< this->npointsy; b++){
                    double y = (b*this->grid_spacing) + this->ybegin;
                    for (int c=0; c<this->npointsz; c++){
                        double z = (c*this->grid_spacing) + this->zbegin;
                        double elec = 0.0;
                        double vdwA = 0.0;
                        double vdwB = 0.0;
                        double solv=0.0;
                        double rec_solv=0.0;

                        for (int i=0; i< Rec->N; i++){
                            double d = this->distance(x, Rec->xyz[i][0], y,Rec->xyz[i][1], z, Rec->xyz[i][2]);
                            double d2 = d*d;
                            double d3 = d*d*d;
                            double d6 = d2*d2*d2;

                            double denom = pow((d3 + Input->deltaij_es3), (1.0/3.0));

                            if (Input->dielectric_model == "constant"){
                                elec += (332.0*Rec->charges[i])/(Input->diel*denom);
                                denom = pow((d6 + Input->deltaij_es6), (1.0/3.0));
                                solv += ((Input->solvation_alpha * Rec->charges[i] * Rec->charges[i])+ Input->solvation_beta) *  exp((-denom/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                                rec_solv += (4.0/3.0) * PI * pow(Rec->radii[i], 3) * exp((-denom/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                            }

                            else {                      // dielectric model == r
                                denom = pow((d6 + Input->deltaij_es6), (1.0/3.0));
                                elec += (332.0* Rec->charges[i]/(Input->diel*denom));
                                solv += ((Input->solvation_alpha * Rec->charges[i] * Rec->charges[i])+ Input->solvation_beta) *  exp((-denom/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                                rec_solv += (4.0/3.0) * PI * pow(Rec->radii[i], 3) * exp((-denom/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                            }

                            denom = (d6 + Input->deltaij6);
                            vdwA += (4096.0 * Rec->epsilons_sqrt[i] * pow(Rec->radii[i], 6)) / (denom*denom);
                            vdwB += ( 128.0 * Rec->epsilons_sqrt[i] * pow(Rec->radii[i], 3)) / denom;
                        }
                        this->elec_grid[a][b][c] = elec;
                        this->vdwA_grid[a][b][c] = vdwA;
                        this->vdwB_grid[a][b][c] = vdwB;
                        this->solv_gauss[a][b][c] = solv;
                        this->rec_solv_gauss[a][b][c] = rec_solv;
                    }
                }
            }

            this->rec_si = 0.00;
            for(int i=0; i<Rec->N; i++){
                this->rec_si += (Input->solvation_alpha*Rec->charges[i]*Rec->charges[i]) + Input->solvation_beta;
            }
        }
    }
    else {
        this->compute_grid_hardcore(Rec);
    }
}

void Grid::write_grids_to_file(){
	FILE* outgrid;
    double elec, vdwA, vdwB, solv, rec_solv, pbsa, delphi;
    int pbsa_flag = 0;
    if (this->pbsa_loaded){
        pbsa_flag = 1;
    }
    else if (this->delphi_loaded){
        pbsa_flag = 2;
    }
	outgrid = fopen((Input->grid_prefix + ".grid").c_str(), "wb");
	if (outgrid == NULL){
		printf("Could not open McGrid file. Please check");
		exit(1);
	}
	fwrite(&this->npointsx, sizeof(int), 1, outgrid);
	fwrite(&this->npointsy, sizeof(int), 1, outgrid);
	fwrite(&this->npointsz, sizeof(int), 1, outgrid);

	fwrite(&this->grid_spacing, sizeof(double), 1, outgrid);

	fwrite(&this->xbegin, sizeof(double), 1, outgrid);
	fwrite(&this->ybegin, sizeof(double), 1, outgrid);
	fwrite(&this->zbegin, sizeof(double), 1, outgrid);

    fwrite(&pbsa_flag, sizeof(int), 1, outgrid);

	for (int a=0; a<npointsx; a++){
		for (int b=0; b<npointsy; b++){
			for (int c=0; c<npointsz; c++){
				elec = this->elec_grid[a][b][c];
				vdwA = this->vdwA_grid[a][b][c];
				vdwB = this->vdwB_grid[a][b][c];
				solv = this->solv_gauss[a][b][c];
				rec_solv = this->rec_solv_gauss[a][b][c];
                if (this->pbsa_loaded){
                    pbsa = this->pbsa_grid[a][b][c];
                }
                else if (this->delphi_loaded){
                    delphi = this->delphi_grid[a][b][c];
                }
                else {
                    pbsa=0.0;
                    delphi=0.0;
                }
				fwrite(&elec, sizeof(double), 1, outgrid);
				fwrite(&vdwA, sizeof(double), 1, outgrid);
				fwrite(&vdwB, sizeof(double), 1, outgrid);
				fwrite(&rec_solv, sizeof(double), 1, outgrid);
				fwrite(&solv, sizeof(double), 1, outgrid);
                fwrite(&pbsa, sizeof(double), 1, outgrid);
                fwrite(&delphi, sizeof(double), 1, outgrid);
			}
		}
	}
	fclose(outgrid);
}

void Grid::load_grids_from_file(){
	FILE* ingrid;
    double elec, vdwA, vdwB, solv, rec_solv, pbsa, delphi;
    size_t garbage;
    int pbsa_flag;
	ingrid = fopen((Input->grid_prefix + ".grid").c_str(), "rb");
	if (ingrid == NULL){
		printf("Could not open McGrid file. Please check");
		exit(1);
	}

	garbage = fread(&this->npointsx, sizeof(int), 1, ingrid);
	garbage = fread(&this->npointsy, sizeof(int), 1, ingrid);
	garbage = fread(&this->npointsz, sizeof(int), 1, ingrid);

	garbage = fread(&this->grid_spacing, sizeof(double), 1, ingrid);

	garbage = fread(&this->xbegin, sizeof(double), 1, ingrid);
	garbage = fread(&this->ybegin, sizeof(double), 1, ingrid);
	garbage = fread(&this->zbegin, sizeof(double), 1, ingrid);

    garbage = fread(&pbsa_flag, sizeof(int), 1, ingrid);

    switch (pbsa_flag) {
    case 1:
        this->pbsa_loaded = true;
        break;
    case 2:
        this->delphi_loaded = true;
        break;
    }

    vector<double> elec_t1(npointsz), vdwA_t1(npointsz), vdwB_t1(npointsz), solv_t1(npointsz),rec_solv_t1(npointsz), pbsa_t1(npointsz), delphi_t1(npointsz);
    vector<vector<double> > elec_t2, vdwA_t2, vdwB_t2, solv_t2, rec_solv_t2, pbsa_t2, delphi_t2;

		for(int a=0; a< this->npointsx; a++){
			for (int b=0; b< this->npointsy; b++){
				for (int c=0; c<this->npointsz; c++){
					garbage = fread(&elec, sizeof(double), 1, ingrid);
					garbage = fread(&vdwA, sizeof(double), 1, ingrid);
					garbage = fread(&vdwB, sizeof(double), 1, ingrid);
					garbage = fread(&rec_solv, sizeof(double), 1, ingrid);
					garbage = fread(&solv, sizeof(double), 1, ingrid);
                    garbage = fread(&pbsa, sizeof(double), 1, ingrid);
                    garbage = fread(&delphi, sizeof(double), 1, ingrid);

					elec_t1[c] = (elec);
					vdwA_t1[c] = (vdwA);
					vdwB_t1[c] = (vdwB);
					solv_t1[c] = (solv);
					rec_solv_t1[c] = (rec_solv);
                    pbsa_t1[c] = (pbsa);
                    delphi_t1[c] = (delphi);
				}
				elec_t2.push_back(elec_t1);
				vdwA_t2.push_back(vdwA_t1);
				vdwB_t2.push_back(vdwB_t1);
				solv_t2.push_back(solv_t1);
				rec_solv_t2.push_back(rec_solv_t1);
                pbsa_t2.push_back(pbsa_t1);
                delphi_t2.push_back(delphi_t1);
			}
			this->elec_grid.push_back(elec_t2);
			this->vdwA_grid.push_back(vdwA_t2);
			this->vdwB_grid.push_back(vdwB_t2);
			this->solv_gauss.push_back(solv_t2);
			this->rec_solv_gauss.push_back(rec_solv_t2);
            this->pbsa_grid.push_back(pbsa_t2);
            this->delphi_grid.push_back(delphi_t2);
			elec_t2.clear();
			vdwA_t2.clear();
			vdwB_t2.clear();
			solv_t2.clear();
			rec_solv_t2.clear();
            pbsa_t2.clear();
            delphi_t2.clear();
		}
	fclose(ingrid);
}

Grid::~Grid() {
	this->elec_grid.clear();
	this->vdwA_grid.clear();
	this->vdwB_grid.clear();
	this->rec_solv_gauss.clear();
	this->solv_gauss.clear();
}

void Grid::compute_grid_hardcore(Mol2* Rec){
	vector<double> elec_t1(npointsz), vdwA_t1(npointsz), vdwB_t1(npointsz), solv_t1(npointsz),rec_solv_t1(npointsz);
	vector<vector<double> > elec_t2, vdwA_t2, vdwB_t2, solv_t2, rec_solv_t2;

    double elec, d, d2, d6, x, y, z, vdwA, vdwB, solv, rec_solv, deff;
    double sqrt2 = sqrt(2.0);

	for(int a=0; a< this->npointsx; a++){
        x = (a*this->grid_spacing) + this->xbegin;

		for (int b=0; b< this->npointsy; b++){
			y = (b*this->grid_spacing) + this->ybegin;

			for (int c=0; c<this->npointsz; c++){
				z = (c*this->grid_spacing) + this->zbegin;
				elec = 0.0;
				vdwA = 0.0;
				vdwB = 0.0;
				solv=0.0;
				rec_solv=0.0;

				for (int i=0; i< Rec->N; i++){
					d2 = this->distance_squared(x, Rec->xyz[i][0], y,Rec->xyz[i][1], z, Rec->xyz[i][2]);
					d6 = d2 * d2 * d2;
                    if (Input->dielectric_model == "constant"){
                        d = sqrt(d2);
                        elec += 332.0*((Rec->charges[i])/(d*Input->diel));
                    }
                    else if (Input->dielectric_model == "4r") {     // epsilon = 4r
                        elec += 332.0 * (Rec->charges[i]/(4*d2));
                    }
                    else {                                          // Input->dielectric_model = "r"
                        elec += 332.0 * (Rec->charges[i]/d2);
                    }

                    vdwA += Rec->epsilons_sqrt[i]*64.0*pow(Rec->radii[i], 6) / (d6*d6);
                    vdwB += sqrt2*Rec->epsilons_sqrt[i]*8.0*pow(Rec->radii[i], 3) / d6;

                    deff = (d2);

                    solv += ((Input->solvation_alpha * Rec->charges[i] * Rec->charges[i])+ Input->solvation_beta) *  exp((-(deff)/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                    rec_solv += (4.0/3.0) * PI * pow(Rec->radii[i], 3) * exp((-(deff)/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
				}
                elec_t1[c] = (elec);
				vdwA_t1[c] = (vdwA);
				vdwB_t1[c] = (vdwB);
				solv_t1[c] = (solv);
				rec_solv_t1[c] = (rec_solv);
			}
			elec_t2.push_back(elec_t1);
			vdwA_t2.push_back(vdwA_t1);
			vdwB_t2.push_back(vdwB_t1);
			solv_t2.push_back(solv_t1);
			rec_solv_t2.push_back(rec_solv_t1);

		}
		this->elec_grid.push_back(elec_t2);
		this->vdwA_grid.push_back(vdwA_t2);
		this->vdwB_grid.push_back(vdwB_t2);
		this->solv_gauss.push_back(solv_t2);
		this->rec_solv_gauss.push_back(rec_solv_t2);
		elec_t2.clear();
		vdwA_t2.clear();
		vdwB_t2.clear();
		solv_t2.clear();
		rec_solv_t2.clear();
	}
	this->rec_si = 0.00;
	for(int i=0; i<Rec->N; i++){
        this->rec_si += (Input->solvation_alpha*Rec->charges[i]*Rec->charges[i]) + Input->solvation_beta;
	}
}

void Grid::load_Ambergrids_from_file(){
    FILE *pbsa_map;
    char str[80];

    if (Input->pbsa_grid.substr(Input->pbsa_grid.size()-3, 3) == ".gz"){
        this->load_gzAmbergrids_from_file();
    }
    else {

        pbsa_map=fopen(Input->pbsa_grid.c_str(),"r");

        if (pbsa_map == NULL){
            printf("Could not open PBSA Grid file. Please check");
            exit(1);
        }

        str[0] = '#';
        while(str[0] =='#'){
            fgets(str, 80, pbsa_map);
        }

        float space, x0, y0, z0;
        sscanf(str, "%f %f %f %f", &space, &x0, &y0, &z0);
        fscanf(pbsa_map, "%d %d %d", &this->npointsx, &this->npointsy, &this->npointsz);

        this->grid_spacing = double(space);
        this->xbegin = double(x0);
        this->ybegin = double(y0);
        this->zbegin = double(z0);
        this->xend = (npointsx*grid_spacing)+xbegin;
        this->yend = (npointsy*grid_spacing)+ybegin;
        this->zend = (npointsz*grid_spacing)+zbegin;

        //    vector initialization

        vector<double> vz(npointsz);
        vector<vector<double> > vtmp;
        for (int i=0; i< npointsy; i++){
            vtmp.push_back(vz);
        }

        for (int i=0; i< npointsx; i++){
            this->pbsa_grid.push_back(vtmp);
        }

        float phi;

        for (int x=0; x<npointsx; x++){
            for (int y=0; y<npointsy; y++){
                for (int z=0; z< npointsz; z++){
                    fscanf(pbsa_map, "%f", &phi);
                    this->pbsa_grid[x][y][z] = double(phi);
                }
            }
        }

        fclose(pbsa_map);

        sprintf(info, "PBSA Grid file %s read!", Input->pbsa_grid.c_str());
        Writer->print_info(info);

        this->pbsa_loaded = true;
    }
}

void Grid::load_gzAmbergrids_from_file(){
    char str[80];

    gzFile pbsa_map=gzopen(Input->pbsa_grid.c_str(),"r");

    if (pbsa_map == NULL){
        printf("Could not open PBSA Grid file. Please check");
        exit(1);
    }

    str[0] = '#';
    while(str[0] =='#'){
        gzgets(pbsa_map, str, 80);
    }

    float space, x0, y0, z0;
    sscanf(str, "%f %f %f %f", &space, &x0, &y0, &z0);
    gzgets(pbsa_map, str, 80);
    sscanf(str, "%d %d %d", &this->npointsx, &this->npointsy, &this->npointsz);

    this->grid_spacing = double(space);
    this->xbegin = double(x0);
    this->ybegin = double(y0);
    this->zbegin = double(z0);
    this->xend = (npointsx*grid_spacing)+xbegin;
    this->yend = (npointsy*grid_spacing)+ybegin;
    this->zend = (npointsz*grid_spacing)+zbegin;

//    vector initialization

    vector<double> vz(npointsz);
    vector<vector<double> > vtmp;
    for (int i=0; i< npointsy; i++){
        vtmp.push_back(vz);
    }

    for (int i=0; i< npointsx; i++){
        this->pbsa_grid.push_back(vtmp);
    }

    float phi;
    int count_line;

    for (int x=0; x<npointsx; x++){
        for (int y=0; y<npointsy; y++){
            count_line = 0;
            for (int z=0; z< npointsz; z++){
                count_line++;
               gzgets(pbsa_map, str, 19);
               phi = atof(str);
               this->pbsa_grid[x][y][z] = double(phi);
               if (count_line == 6){
                   count_line = 0;
                   gzgets(pbsa_map, str, 80);               // go to the end of line
               }
            }
        }
    }
    gzclose(pbsa_map);

    sprintf(info, "PBSA Grid file %s read!", Input->pbsa_grid.c_str());
    Writer->print_info(info);

    this->pbsa_loaded = true;
}

void Grid::compute_grid_hardcore_omp(Mol2* Rec){
    vector<double> elec_t1(npointsz), vdwA_t1(npointsz), vdwB_t1(npointsz), solv_t1(npointsz),rec_solv_t1(npointsz);
    vector<vector<double> > elec_t2, vdwA_t2, vdwB_t2, solv_t2, rec_solv_t2;

    double sqrt2 = sqrt(2.0);

// initializing the vectors;

    for (int i=0; i<this->npointsy; i++){
        elec_t2.push_back(elec_t1);
        vdwA_t2.push_back(vdwA_t1);
        vdwB_t2.push_back(vdwB_t1);
        solv_t2.push_back(solv_t1);
        rec_solv_t2.push_back(rec_solv_t1);
    }

    for (int i=0; i<this->npointsx; i++){
        this->elec_grid.push_back(elec_t2);
        this->vdwA_grid.push_back(vdwA_t2);
        this->vdwB_grid.push_back(vdwB_t2);
        this->solv_gauss.push_back(solv_t2);
        this->rec_solv_gauss.push_back(rec_solv_t2);
    }

// Now, starting OMP...

#pragma omp parallel num_threads(Input->parallel_jobs)
{
#pragma omp for schedule(static, 1)
    for(int a=0; a< this->npointsx; a++){

        double x = (a*this->grid_spacing) + this->xbegin;

        for (int b=0; b< this->npointsy; b++){
            double y = (b*this->grid_spacing) + this->ybegin;

            for (int c=0; c<this->npointsz; c++){
                double z = (c*this->grid_spacing) + this->zbegin;
                double elec = 0.0;
                double vdwA = 0.0;
                double vdwB = 0.0;
                double solv=0.0;
                double rec_solv=0.0;

                for (int i=0; i< Rec->N; i++){
                    double d2 = this->distance_squared(x, Rec->xyz[i][0], y,Rec->xyz[i][1], z, Rec->xyz[i][2]);
                    double d6 = d2 * d2 * d2;
                    if (Input->dielectric_model == "constant"){
                        double d = sqrt(d2);
                        elec += 332.0*((Rec->charges[i])/(d*Input->diel));
                    }
                    else if (Input->dielectric_model == "4r") {     // epsilon = 4r
                        elec += 332.0 * (Rec->charges[i]/(4*d2));
                    }
                    else {                                          // Input->dielectric_model = "r"
                        elec += 332.0 * (Rec->charges[i]/d2);
                    }

                    vdwA += Rec->epsilons_sqrt[i]*64.0*pow(Rec->radii[i], 6) / (d6*d6);
                    vdwB += sqrt2*Rec->epsilons_sqrt[i]*8.0*pow(Rec->radii[i], 3) / d6;

                    double deff = (d2);

                    solv += ((Input->solvation_alpha * Rec->charges[i] * Rec->charges[i])+ Input->solvation_beta) *  exp((-(deff)/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                    rec_solv += (4.0/3.0) * PI * pow(Rec->radii[i], 3) * exp((-(deff)/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                }
                this->elec_grid[a][b][c] = elec;
                this->vdwA_grid[a][b][c] = vdwA;
                this->vdwB_grid[a][b][c] = vdwB;
                this->solv_gauss[a][b][c] = solv;
                this->rec_solv_gauss[a][b][c] = rec_solv;
            }
        }
    }
}                       // end of pragma
    this->rec_si = 0.00;
    for(int i=0; i<Rec->N; i++){
        this->rec_si += (Input->solvation_alpha*Rec->charges[i]*Rec->charges[i]) + Input->solvation_beta;
    }
}

void Grid::load_Delphi_Grid_from_file(){
    FILE *phimap;
    char *title;
    int ivary, nbyte, inddat, intx, inty, intz, tmpi;
    double xang, yang, zang, xstart, xend, ystart, yend, zstart, zend;
    double extent;

    phimap = fopen(Input->delphi_grid.c_str(), "rb");

    fread(&tmpi, sizeof(int), 1, phimap);

    title = (char *) malloc(sizeof(char) * 62);
    for (int i=0; i<60; i++) {
        title[i] = fgetc(phimap);
    }
    title[60] = '\n';
    title[61] = (char) 0;

    fread(&tmpi, sizeof(int), 1, phimap);
    fread(&tmpi, sizeof(int), 1, phimap);

    fread(&ivary, sizeof(int), 1, phimap);
    fread(&nbyte, sizeof(int), 1, phimap);
    fread(&inddat, sizeof(int), 1, phimap);

    fread(&extent, sizeof(double), 1, phimap);
    fread(&extent, sizeof(double), 1, phimap);
    fread(&extent, sizeof(double), 1, phimap);

    fread(&xang, sizeof(double), 1, phimap);
    fread(&yang, sizeof(double), 1, phimap);
    fread(&zang, sizeof(double), 1, phimap);

    fread(&xstart, sizeof(double), 1, phimap);
    fread(&xend, sizeof(double), 1, phimap);

    fread(&ystart, sizeof(double), 1, phimap);
    fread(&yend, sizeof(double), 1, phimap);

    fread(&zstart, sizeof(double), 1, phimap);
    fread(&zend, sizeof(double), 1, phimap);

    fread(&intx, sizeof(int), 1, phimap);
    fread(&inty, sizeof(int), 1, phimap);
    fread(&intz, sizeof(int), 1, phimap);

    this->npointsx = intx+1;
    this->npointsy = inty+1;
    this->npointsz = intz+1;

    this->xbegin = extent*xstart;
    this->xend = extent*xend;

    this->ybegin = extent*ystart;
    this->yend = extent*yend;

    this->zbegin = extent*ystart;
    this->zend = extent*yend;

    this->grid_spacing = 0.25;

#ifdef DEBUG
    printf("%s\n", title);
    printf("ivary: %d\n", ivary);
    printf("nbyte: %d\n", nbyte);
    printf("inddat: %d\n", inddat);
    printf("Box angles: %f %f %f\n", xang, yang, zang);
    printf("Box limits: %f -  %f   %f - %f   %f - %f\n", xstart, xend, ystart, yend, zstart, zend);
    printf("%d %d %d\n", intx, inty, intz);
#endif

    double phi;

    vector<double> vz(intz+1);
    vector<vector<double> > vtmp;
    for (int i=0; i< inty+1; i++){
        vtmp.push_back(vz);
    }

    for (int i=0; i< intx+1; i++){
        this->delphi_grid.push_back(vtmp);
    }

    for (int z=0; z<intz+1; z++){
        for(int y=0; y<inty+1; y++){
            for (int x=0; x< intx+1; x++){
                fread(&phi, sizeof(double), 1, phimap);
                this->delphi_grid[x][y][z] = (0.593f*phi); // potential is given in KT/e
            }
        }
    }

    fclose(phimap);

    sprintf(info, "DelPhi Grid file %s read!", Input->delphi_grid.c_str());
    Writer->print_info(info);

    this->delphi_loaded = true;
}

void Grid::load_phimap_from_file(int gsize){
    FILE *phimap;
    char *uplbl, *nxtlbl, *toplbl, *botlbl;
    double scale, oldmid_x, oldmid_y, oldmid_z;

    phimap = fopen(Input->delphi_grid.c_str(), "rb");

    if (phimap == NULL){
        printf("Could not open file %s. Please check!\n", Input->delphi_grid.c_str());
        exit(1);
    }

    uplbl = (char *) malloc(sizeof(char) * 23);
    for (int i=0; i<21; i++) {
        uplbl[i] = fgetc(phimap);
    }
    uplbl[21] = '\n';
    uplbl[22] = (char) 0;

    nxtlbl = (char *) malloc(sizeof(char) * 12);
    for (int i=0; i<10; i++) {
        nxtlbl[i] = fgetc(phimap);
    }
    nxtlbl[10] = '\n';
    nxtlbl[11] = (char) 0;

    toplbl = (char *) malloc(sizeof(char) * 56);
    for (int i=0; i<54; i++) {
        toplbl[i] = fgetc(phimap);
    }
    toplbl[54] = '\n';
    toplbl[55] = (char) 0;


#ifdef DEBUG
    printf("%s\n", uplbl);
    printf("%s\n", nxtlbl);
    printf("%s\n", toplbl);
#endif

    vector<double> vz(gsize);
    vector<vector<double> > vtmp;
    for (int i=0; i<gsize; i++){
        vtmp.push_back(vz);
    }

    for (int i=0; i<gsize; i++){
        this->delphi_grid.push_back(vtmp);
    }

    double kt_phi;
    for (int nz=0; nz < gsize; nz++){
        for (int ny=0; ny < gsize; ny++){
            for (int nx=0; nx < gsize; nx++){
                fread(&kt_phi, sizeof(double), 1, phimap);
                this->delphi_grid[nx][ny][nz] = 0.593f*kt_phi;  //converting kt units to kcal/mol
            }
        }
    }

    botlbl = (char *) malloc(sizeof(char) * 18);
    for (int i=0; i<16; i++) {
        botlbl[i] = fgetc(phimap);
    }
    botlbl[16] = '\n';
    botlbl[17] = (char) 0;

#ifdef DEBUG
    printf("%s\n", botlbl);
#endif

    int igrid;
    fread(&scale, sizeof(double), 1, phimap);
    fread(&oldmid_x, sizeof(double), 1, phimap);
    fread(&oldmid_y, sizeof(double), 1, phimap);
    fread(&oldmid_z, sizeof(double), 1, phimap);
    fread(&igrid, sizeof(int), 1, phimap);

    if (igrid =! gsize){
        printf("Mismatch between GSIZE given and found in phimap file. Please check.\n");
        exit(1);
    }

    this->grid_spacing = 1/scale;

    this->npointsx = gsize;
    this->npointsy = gsize;
    this->npointsz = gsize;

    this->xbegin = (((1)-((gsize+1)/2))/scale)+oldmid_x;
    this->xend = (((gsize+1)-((gsize+1)/2))/scale)+oldmid_x;

    this->ybegin = (((1)-((gsize+1)/2))/scale)+oldmid_y;
    this->yend = (((gsize+1)-((gsize+1)/2))/scale)+oldmid_y;

    this->zbegin = (((1)-((gsize+1)/2))/scale)+oldmid_z;
    this->zend = (((gsize+1)-((gsize+1)/2))/scale)+oldmid_z;

    fclose(phimap);
    sprintf(info, "DelPhi Grid file %s read!", Input->delphi_grid.c_str());
    Writer->print_info(info);
    sprintf(info, "DelPhi Grid Scale / Spacing: %7.4f %7.4f", scale, this->grid_spacing);
    Writer->print_info(info);
    sprintf(info, "DelPhi Grid GSIZE: %5d. Mid XYZ: %7.4f %7.4f %7.4f", gsize, oldmid_x, oldmid_y, oldmid_z);
    Writer->print_info(info);
    sprintf(info, "Grid Dimensions XYZ: %7.4f - %7.4f %7.4f - %7.4f %7.4f - %7.4f", this->xbegin, this->xend, this->ybegin, this->yend,
            this->zbegin, this->zend);
    Writer->print_info(info);

    this->delphi_loaded = true;
}

void Grid::load_delphi_cube(){
    FILE* phimap;
    int igrid, itmp;
    float scale, centerx, centery, centerz, dtmp;
    char str[200];

    phimap = fopen(Input->delphi_cube_grid.c_str(), "r");

    fscanf(phimap, "%f %d %f %f %f", &scale, &igrid, &centerx, &centery, &centerz);
    this->grid_spacing = double(1.0/scale);
    this->xbegin = double(centerx - ((igrid-1)/2)*this->grid_spacing);
    this->xend = double(centerx + (((igrid+1)/2)*this->grid_spacing));
    this->npointsx = igrid;

    this->ybegin = double(centery - ((igrid-1)/2)*this->grid_spacing);
    this->yend = double(centery + (((igrid+1)/2)*this->grid_spacing));
    this->npointsy = igrid;

    this->zbegin = double(centerz - ((igrid-1)/2)*this->grid_spacing);
    this->zend = double(centerz + (((igrid+1)/2)*this->grid_spacing));
    this->npointsz = igrid;

    fgets(str, 80, phimap);
    fgets(str, 80, phimap);
    fscanf(phimap, "%d %f %f %f", &itmp, &dtmp, &dtmp, &dtmp);
    fscanf(phimap, "%d %f %f %f", &itmp, &dtmp, &dtmp, &dtmp);
    fscanf(phimap, "%d %f %f %f", &itmp, &dtmp, &dtmp, &dtmp);
    fscanf(phimap, "%d %f %f %f", &itmp, &dtmp, &dtmp, &dtmp);
    fscanf(phimap, "%d %f %f %f %f", &itmp, &dtmp, &dtmp, &dtmp, &dtmp);

    vector<double> vz(igrid);
    vector<vector<double> > vtmp;
    for (int i=0; i<igrid; i++){
        vtmp.push_back(vz);
    }

    for (int i=0; i<igrid; i++){
        this->delphi_grid.push_back(vtmp);
    }

    float phi;

    for (int nx=0; nx < igrid; nx++){
        for (int ny=0; ny < igrid; ny++){
            for (int nz=0; nz < igrid; nz++){
                fscanf(phimap, "%f", &phi);
                this->delphi_grid[nx][ny][nz] = double(0.593f*phi);  //converting kt units to kcal/mol
            }
        }
    }

    fclose(phimap);

    sprintf(info, "DelPhi CUBE Grid file %s read!", Input->delphi_cube_grid.c_str());
    Writer->print_info(info);
    sprintf(info, "DelPhi CUBE Grid Spacing: %7.4f", this->grid_spacing);
    Writer->print_info(info);
    sprintf(info, "DelPhi CUBE Grid GSIZE: %5d. Center: %7.4f %7.4f %7.4f", igrid, centerx, centery, centerz);
    Writer->print_info(info);
    sprintf(info, "Grid Dimensions XYZ: [%7.4f - %7.4f] [%7.4f - %7.4f] [%7.4f - %7.4f]", this->xbegin, this->xend, this->ybegin, this->yend,
            this->zbegin, this->zend);
    Writer->print_info(info);
    this->delphi_loaded = true;
}

void Grid::load_delphi_gzcube(){
    gzFile phimap = gzopen(Input->delphi_cube_grid.c_str(),"r");

    if (phimap == NULL){
        printf("Could not open DelPhi Grid file %s. Please check", Input->delphi_cube_grid.c_str());
        exit(1);
    }

    int igrid;
    float scale, centerx, centery, centerz;
    char str[200];

    gzgets(phimap, str, 80);
    sscanf(str, "%10f%6d%9f%9f%9f", &scale, &igrid, &centerx, &centery, &centerz);
    this->grid_spacing = double(1.0/scale);
    this->xbegin = double(centerx - ((igrid-1)/2)*this->grid_spacing);
    this->xend = double(centerx + (((igrid+1)/2)*this->grid_spacing));
    this->npointsx = igrid;

    this->ybegin = double(centery - ((igrid-1)/2)*this->grid_spacing);
    this->yend = double(centery + (((igrid+1)/2)*this->grid_spacing));
    this->npointsy = igrid;

    this->zbegin = double(centerz - ((igrid-1)/2)*this->grid_spacing);
    this->zend = double(centerz + (((igrid+1)/2)*this->grid_spacing));
    this->npointsz = igrid;

    gzgets(phimap, str, 80);
    gzgets(phimap, str, 80);
    gzgets(phimap, str, 80);
    gzgets(phimap, str, 80);
    gzgets(phimap, str, 80);
    gzgets(phimap, str, 80);

    vector<double> vz(igrid);
    vector<vector<double> > vtmp;
    for (int i=0; i<igrid; i++){
        vtmp.push_back(vz);
    }

    for (int i=0; i<igrid; i++){
        this->delphi_grid.push_back(vtmp);
    }

    float phi;
    int count=0;

    for (int nx=0; nx < igrid; nx++){
        for (int ny=0; ny < igrid; ny++){
            for (int nz=0; nz < igrid; nz++){
                count++;
                gzgets(phimap, str, 14);
                phi = atof(str);
                this->delphi_grid[nx][ny][nz] = double(0.593f*phi);  //converting kt units to kcal/mol
                if (count == 6){
                    gzgets(phimap, str, 14);
                    count = 0;
                }
            }
            gzgets(phimap, str, 14);
            count = 0;
        }
    }

    gzclose(phimap);

    sprintf(info, "DelPhi CUBE Grid file %s read!", Input->delphi_cube_grid.c_str());
    Writer->print_info(info);
    sprintf(info, "DelPhi CUBE Grid Spacing: %7.4f", this->grid_spacing);
    Writer->print_info(info);
    sprintf(info, "DelPhi CUBE Grid GSIZE: %5d. Center: %7.4f %7.4f %7.4f", igrid, centerx, centery, centerz);
    Writer->print_info(info);
    sprintf(info, "Grid Dimensions XYZ: [%7.4f - %7.4f] [%7.4f - %7.4f] [%7.4f - %7.4f]", this->xbegin, this->xend, this->ybegin, this->yend,
            this->zbegin, this->zend);
    Writer->print_info(info);
    this->delphi_loaded = true;

}

void Grid::compute_grid_hardcore_Gaussian(Mol2* Rec){
    vector<double> elec_t1(npointsz), vdwA_t1(npointsz), vdwB_t1(npointsz), solv_t1(npointsz),rec_solv_t1(npointsz);
    vector<vector<double> > elec_t2, vdwA_t2, vdwB_t2, solv_t2, rec_solv_t2;

    double elec, d, d2, d6, x, y, z, vdwA, vdwB, solv, rec_solv, deff, LJ_gauss_weight, Coulomb_gauss_weight;
    double sqrt2 = sqrt(2.0);

    for(int a=0; a< this->npointsx; a++){
        x = (a*this->grid_spacing) + this->xbegin;

        for (int b=0; b< this->npointsy; b++){
            y = (b*this->grid_spacing) + this->ybegin;

            for (int c=0; c<this->npointsz; c++){
                z = (c*this->grid_spacing) + this->zbegin;
                elec = 0.0;
                vdwA = 0.0;
                vdwB = 0.0;
                solv=0.0;
                rec_solv=0.0;

                for (int i=0; i< Rec->N; i++){
                    d2 = this->distance_squared(x, Rec->xyz[i][0], y,Rec->xyz[i][1], z, Rec->xyz[i][2]);
                    d = sqrt(d2);
                    d6 = d2 * d2 * d2;

                    double dw = d-(1.2*Rec->radii[i]);
                    dw > 0. ? LJ_gauss_weight = exp(-( (d-(Rec->radii[i]+Rec->radii[i])) * (d-(Rec->radii[i]+Rec->radii[i])) ) /(2.0*Input->LJ_sigma*Input->LJ_sigma)) : -1.0;

                    Coulomb_gauss_weight = exp(-( (d-(Rec->radii[i]+Rec->radii[i])) * (d-(Rec->radii[i]+Rec->radii[i])) ) /(2.0*Input->coulomb_sigma*Input->coulomb_sigma));

                    // Electrostatic Potential..

                    if (Input->dielectric_model == "constant"){
                        if (Input->use_GW_Coulomb){
                            elec += (332.0*((Rec->charges[i])/(d*Input->diel)))*Coulomb_gauss_weight;
                        }
                        else {
                            elec += (332.0*((Rec->charges[i])/(d*Input->diel)));
                        }
                    }
                    else if (Input->dielectric_model == "4r") {     // epsilon = 4r
                        if (Input->use_GW_Coulomb){
                            elec += (332.0 * (Rec->charges[i]/(4*d2)))*Coulomb_gauss_weight;
                        }
                        else {
                            elec += (332.0 * (Rec->charges[i]/(4*d2)));
                        }
                    }
                    else {                                          // Input->dielectric_model = "r"
                        if (Input->use_GW_Coulomb){
                            elec += (332.0 * (Rec->charges[i]/d2))*Coulomb_gauss_weight;
                        }
                        else{
                            elec += (332.0 * (Rec->charges[i]/d2));
                        }
                    }

                    // VDW Repulsive Potential

                    if (Input->use_GW_LJ12){
                        LJ_gauss_weight > 0 ? vdwA += (Rec->epsilons_sqrt[i]*64.0*pow(Rec->radii[i], 6) / (d6*d6))*LJ_gauss_weight: vdwA += MAX_ALJ ;
                    }
                    else {
                        vdwA += (Rec->epsilons_sqrt[i]*64.0*pow(Rec->radii[i], 6) / (d6*d6));
                    }


                    // VDW Attractive Potential

                    if (Input->use_GW_LJ6){
                        LJ_gauss_weight > 0 ? vdwB += (sqrt2*Rec->epsilons_sqrt[i]*8.0*pow(Rec->radii[i], 3) / d6)*LJ_gauss_weight: vdwB += MAX_BLJ;
                    }
                    else {
                        vdwB += (sqrt2*Rec->epsilons_sqrt[i]*8.0*pow(Rec->radii[i], 3) / d6);
                    }

                    deff = (d2);

                    solv += ((Input->solvation_alpha * Rec->charges[i] * Rec->charges[i])+ Input->solvation_beta) *  exp((-(deff)/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                    rec_solv += (4.0/3.0) * PI * pow(Rec->radii[i], 3) * exp((-(deff)/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                }
                elec_t1[c] = (elec);
                vdwA_t1[c] = (vdwA);
                vdwB_t1[c] = (vdwB);
                solv_t1[c] = (solv);
                rec_solv_t1[c] = (rec_solv);
            }
            elec_t2.push_back(elec_t1);
            vdwA_t2.push_back(vdwA_t1);
            vdwB_t2.push_back(vdwB_t1);
            solv_t2.push_back(solv_t1);
            rec_solv_t2.push_back(rec_solv_t1);

        }
        this->elec_grid.push_back(elec_t2);
        this->vdwA_grid.push_back(vdwA_t2);
        this->vdwB_grid.push_back(vdwB_t2);
        this->solv_gauss.push_back(solv_t2);
        this->rec_solv_gauss.push_back(rec_solv_t2);
        elec_t2.clear();
        vdwA_t2.clear();
        vdwB_t2.clear();
        solv_t2.clear();
        rec_solv_t2.clear();
    }
    this->rec_si = 0.00;
    for(int i=0; i<Rec->N; i++){
        this->rec_si += (Input->solvation_alpha*Rec->charges[i]*Rec->charges[i]) + Input->solvation_beta;
    }
}

void Grid::compute_grid_hardcore_omp_Gaussian(Mol2* Rec){
    vector<double> elec_t1(npointsz), vdwA_t1(npointsz), vdwB_t1(npointsz), solv_t1(npointsz),rec_solv_t1(npointsz);
    vector<vector<double> > elec_t2, vdwA_t2, vdwB_t2, solv_t2, rec_solv_t2;

    double sqrt2 = sqrt(2.0);

// initializing the vectors;

    for (int i=0; i<this->npointsy; i++){
        elec_t2.push_back(elec_t1);
        vdwA_t2.push_back(vdwA_t1);
        vdwB_t2.push_back(vdwB_t1);
        solv_t2.push_back(solv_t1);
        rec_solv_t2.push_back(rec_solv_t1);
    }

    for (int i=0; i<this->npointsx; i++){
        this->elec_grid.push_back(elec_t2);
        this->vdwA_grid.push_back(vdwA_t2);
        this->vdwB_grid.push_back(vdwB_t2);
        this->solv_gauss.push_back(solv_t2);
        this->rec_solv_gauss.push_back(rec_solv_t2);
    }

// Now, starting OMP...

#pragma omp parallel num_threads(Input->parallel_jobs)
{
#pragma omp for schedule(static, 1)
    for(int a=0; a< this->npointsx; a++){

        double x = (a*this->grid_spacing) + this->xbegin;

        for (int b=0; b< this->npointsy; b++){
            double y = (b*this->grid_spacing) + this->ybegin;

            for (int c=0; c<this->npointsz; c++){
                double z = (c*this->grid_spacing) + this->zbegin;
                double elec = 0.0;
                double vdwA = 0.0;
                double vdwB = 0.0;
                double solv=0.0;
                double rec_solv=0.0;

                for (int i=0; i< Rec->N; i++){
                    double d2 = this->distance_squared(x, Rec->xyz[i][0], y,Rec->xyz[i][1], z, Rec->xyz[i][2]);
                    double d6 = d2 * d2 * d2;
                    double d = sqrt(d2);

                    double LJ_gauss_weight;

                    double dw = d-(1.2*Rec->radii[i]);
                    dw > 0. ? LJ_gauss_weight = exp(-( (d-(Rec->radii[i]+Rec->radii[i])) * (d-(Rec->radii[i]+Rec->radii[i])) ) /(2.0*Input->LJ_sigma*Input->LJ_sigma)) : -1.0;
                    double Coulomb_gauss_weight = exp(-( (d-(Rec->radii[i]+Rec->radii[i])) * (d-(Rec->radii[i]+Rec->radii[i])) ) /(2.0*Input->coulomb_sigma*Input->coulomb_sigma));

                    if (Input->dielectric_model == "constant"){
                        double d = sqrt(d2);
                        if (Input->use_GW_Coulomb){
                            elec += (332.0*((Rec->charges[i])/(d*Input->diel)))*Coulomb_gauss_weight;
                        }
                        else{
                            elec += (332.0*((Rec->charges[i])/(d*Input->diel)));
                        }
                    }
                    else if (Input->dielectric_model == "4r") {     // epsilon = 4r
                        if (Input->use_GW_Coulomb){
                            elec += (332.0 * (Rec->charges[i]/(4*d2)))*Coulomb_gauss_weight;
                        }
                        else{
                            elec += (332.0 * (Rec->charges[i]/(4*d2)));
                        }
                    }
                    else {                                          // Input->dielectric_model = "r"
                        if (Input->use_GW_Coulomb){
                            elec += (332.0 * (Rec->charges[i]/d2))*Coulomb_gauss_weight;
                        }
                        else{
                            elec += (332.0 * (Rec->charges[i]/d2));
                        }
                    }

                    if (Input->use_GW_LJ12){
                        LJ_gauss_weight > 0 ? vdwA += (Rec->epsilons_sqrt[i]*64.0*pow(Rec->radii[i], 6) / (d6*d6))*LJ_gauss_weight: vdwA += MAX_ALJ ;
                    }
                    else {
                        vdwA += (Rec->epsilons_sqrt[i]*64.0*pow(Rec->radii[i], 6) / (d6*d6));
                    }

                    if (Input->use_GW_LJ6){
                        LJ_gauss_weight > 0 ? vdwB += (sqrt2*Rec->epsilons_sqrt[i]*8.0*pow(Rec->radii[i], 3) / d6)*LJ_gauss_weight: vdwB += MAX_BLJ;
                    }
                    else {
                        vdwB += (sqrt2*Rec->epsilons_sqrt[i]*8.0*pow(Rec->radii[i], 3) / d6);
                    }


                    double deff = (d2);

                    solv += ((Input->solvation_alpha * Rec->charges[i] * Rec->charges[i])+ Input->solvation_beta) *  exp((-(deff)/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                    rec_solv += (4.0/3.0) * PI * pow(Rec->radii[i], 3) * exp((-(deff)/(2*Input->sigma*Input->sigma))) / (Input->sigma*Input->sigma*Input->sigma);
                }
                this->elec_grid[a][b][c] = elec;
                this->vdwA_grid[a][b][c] = vdwA;
                this->vdwB_grid[a][b][c] = vdwB;
                this->solv_gauss[a][b][c] = solv;
                this->rec_solv_gauss[a][b][c] = rec_solv;
            }
        }
    }
}                       // end of pragma
    this->rec_si = 0.00;
    for(int i=0; i<Rec->N; i++){
        this->rec_si += (Input->solvation_alpha*Rec->charges[i]*Rec->charges[i]) + Input->solvation_beta;
    }
}

/*
 * Grid.cpp
 *
 *  Created on: Jan 7, 2013
 *      Author: asn
 */

#include "Grid.h"
#include "iMcLiBELa.h"

using namespace std;

double Grid::distance(double x1, double x2, double y1, double y2, double z1, double z2) {
	return ( sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}

double Grid::distance_squared(double x1, double x2, double y1, double y2, double z1, double z2) {
	return ((((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}


Grid::Grid(PARSER* _Input) {
	this->Input = _Input;
	grid_spacing = Input->grid_spacing;
}

Grid::Grid(PARSER* _Input, Mol2* Rec, vector<double> com){
	this->Input = _Input;
    if (Input->pbsa_grid != ""){
        this->load_Ambergrids_from_file();
    }
    else {
        this->grid_spacing = Input->grid_spacing;
        this->generate_points(com);
    }

    if (Input->scoring_function < 2){
        this->compute_grid_softcore(Rec);
    }
    else {
            this->compute_grid_hardcore(Rec);
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

void Grid::write_grids_to_file(){
	FILE* outgrid;
	double elec, vdwA, vdwB, solv, rec_solv;
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

	for (int a=0; a<npointsx; a++){
		for (int b=0; b<npointsy; b++){
			for (int c=0; c<npointsz; c++){
				elec = this->elec_grid[a][b][c];
				vdwA = this->vdwA_grid[a][b][c];
				vdwB = this->vdwB_grid[a][b][c];
				solv = this->solv_gauss[a][b][c];
				rec_solv = this->rec_solv_gauss[a][b][c];
				fwrite(&elec, sizeof(double), 1, outgrid);
				fwrite(&vdwA, sizeof(double), 1, outgrid);
				fwrite(&vdwB, sizeof(double), 1, outgrid);
				fwrite(&rec_solv, sizeof(double), 1, outgrid);
				fwrite(&solv, sizeof(double), 1, outgrid);
			}
		}
	}
	fclose(outgrid);
}

void Grid::load_grids_from_file(){
	FILE* ingrid;
	double elec, vdwA, vdwB, solv, rec_solv;
    size_t garbage;
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

	vector<double> elec_t1(npointsz), vdwA_t1(npointsz), vdwB_t1(npointsz), solv_t1(npointsz),rec_solv_t1(npointsz);
	vector<vector<double> > elec_t2, vdwA_t2, vdwB_t2, solv_t2, rec_solv_t2;

		for(int a=0; a< this->npointsx; a++){
			for (int b=0; b< this->npointsy; b++){
				for (int c=0; c<this->npointsz; c++){
					garbage = fread(&elec, sizeof(double), 1, ingrid);
					garbage = fread(&vdwA, sizeof(double), 1, ingrid);
					garbage = fread(&vdwB, sizeof(double), 1, ingrid);
					garbage = fread(&rec_solv, sizeof(double), 1, ingrid);
					garbage = fread(&solv, sizeof(double), 1, ingrid);

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
    double sqrt2 = sqrt(2);

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
//					vdwA += sqrt(Rec->epsilons[i]*(pow((2.0*Rec->radii[i]), 12.0))) / (d6*d6) ;
//                  vdwB += sqrt(2.0 * Rec->epsilons[i]*(pow((2.0*Rec->radii[i]), 6.0))) / (d6) ;

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
    FILE* ingrid;
    char line [256];

//    double grid_spacing, xbegin, ybegin, zbegin;
//    int npointsx, npointsy, npointsz;
    vector <double > values;
    double tempor;


    ingrid=fopen(Input->pbsa_grid.c_str(),"r");

    if (ingrid == NULL){
        printf("Could not open PBSA Grid file. Please check");
        exit(1);
    }

    int row = 0;
    vector<double> vtmp;
    vector<vector<double> > vvtmp;

    while(fgets(line,sizeof(line),ingrid)){
        ++row;

        if (row==9){
            sscanf(line,"%lf %lf %lf %lf",&this->grid_spacing, &this->xbegin, &this->ybegin, &this->zbegin);
        }
        else if (row==10){
            sscanf(line,"%i %i %i",&this->npointsx,&this->npointsy,&this->npointsz);
            this->xend = (npointsx*grid_spacing)+xbegin;
            this->yend = (npointsy*grid_spacing)+ybegin;
            this->zend = (npointsz*grid_spacing)+zbegin;
            for (int a=0; a<this->npointsz; a++){
                vtmp.push_back(0.0);
            }
            for (int a=0; a<this->npointsy; a++){
                vvtmp.push_back(vtmp);
            }
            for (int a=0; a<this->npointsz; a++){
                this->pbsa_grid.push_back(vvtmp);
            }
        }
        else if (row>10){
            istringstream iss(line);
            while (iss>>tempor){
                values.push_back(tempor);
            }
        }
    }
    int counter=0;

    for(int c=0; c<npointsz; c++){
        for (int b=0; b<npointsy; b++){
            for (int a=0; a<npointsx; a++){
                this->pbsa_grid[a][b][c] = values[counter];
                counter++;
            }
        }
    }
    fclose(ingrid);
}


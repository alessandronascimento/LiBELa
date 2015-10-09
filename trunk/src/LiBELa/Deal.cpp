/*
 * Deal.cpp
 *
 *  Created on: 01/06/2012
 *      Author: Nascimento
 */

#include "Deal.h"

Deal::Deal(Mol2* _Rec, Mol2* _Lig, PARSER* _Input) {
	population_size = _Input->deal_pop_size;
	slices = this->define_chunck_sizes();
	chunck_size=int(population_size/QThread::idealThreadCount())+1;
	Rec = _Rec;
	Lig = _Lig;
	Input = _Input;

	int seed = rand() % 100 + 1;
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, seed);

	mu = new double[6];
	sigma = new double[6];

	mu[0]=0.0;
	mu[1]=0.0;
	mu[2]=0.0;
	mu[3]=0.0;
	mu[4]=0.0;
	mu[5]=0.0;
	sigma[0]=60.0;
	sigma[1]=30.0;
	sigma[2]=60.0;
	sigma[3]=1.0;
	sigma[4]=1.0;
	sigma[5]=1.0;

	time_i=clock();

	this->iterate(Input);

	time_f=clock();
	printf("Time elapsed for evaluation: %7.2f msecs.\n", double(((time_f-time_i)*1000)/CLOCKS_PER_SEC));

	//preparing final coordinates
	final_coords = this->rototranslate(Lig, rotrans[0][0], rotrans[0][1], rotrans[0][2], rotrans[0][3], rotrans[0][4], rotrans[0][5]);
}

void Deal::iterate(PARSER* Input){
	for (int i=1; i <= Input->deal_generations; i++){
		this->sort_and_evaluate(this->mu, i);
		this->update_mu(Input->deal_top_ranked);
	}
}


void Deal::sort_array(void){
	double temp;
	vector<double> temp2;
	for (unsigned i=0; i< energies.size()-1 ; i++){
		for (unsigned j=i; j<energies.size(); j++){
			if (energies[j] < energies[i]){
				temp = energies[i];
				temp2 = rotrans[i];

				energies[i]=energies[j];
				rotrans[i] = rotrans[j];

				energies[j]=temp;
				rotrans[j]=temp2;
			}
		}
	}
}


void Deal::sort_and_evaluate(double *mu, int iteration){

	energies.clear();
	rotrans.clear();

	QVector<double> temp;
	QVector<double> Qenergies;
	QVector<QVector<double> > Qrotrans;
	double sqrt_iter = sqrt(iteration);

	for (int i=0; i< population_size; i++){
		temp.append((gsl_ran_gaussian_ziggurat(r, sigma[0]/sqrt_iter)) + mu[0]);
		temp.append((gsl_ran_gaussian_ziggurat(r, sigma[1]/sqrt_iter)) + mu[1]);
		temp.append((gsl_ran_gaussian_ziggurat(r, sigma[2]/sqrt_iter)) + mu[2]);
		temp.append((gsl_ran_gaussian_ziggurat(r, sigma[3]/sqrt_iter)) + mu[3]);
		temp.append((gsl_ran_gaussian_ziggurat(r, sigma[4]/sqrt_iter)) + mu[4]);
		temp.append((gsl_ran_gaussian_ziggurat(r, sigma[5]/sqrt_iter)) + mu[5]);
		Qrotrans.append(temp);
		temp.clear();
	}

	Qenergies = QtConcurrent::blockingMapped(Qrotrans, std::bind1st(std::mem_fun(&Deal::evaluate_energy), this));
	energies = Qenergies.toStdVector();
	for (int i=0; i<population_size; i++){
		rotrans.push_back(Qrotrans.at(i).toStdVector());
	}

	this->sort_array();

/*
	Energy = new ENERGY;
	energies.clear();
	rotrans.clear();
	vector<double> temp;
	double sqrt_iter = sqrt(iteration);
	time_i=clock();
	for (int i=0; i< population_size; i++){
		temp.push_back((gsl_ran_gaussian_ziggurat(r, sigma[0]/sqrt_iter)) + mu[0]);
		temp.push_back((gsl_ran_gaussian_ziggurat(r, sigma[1]/sqrt_iter)) + mu[1]);
		temp.push_back((gsl_ran_gaussian_ziggurat(r, sigma[2]/sqrt_iter)) + mu[2]);
		temp.push_back((gsl_ran_gaussian_ziggurat(r, sigma[3]/sqrt_iter)) + mu[3]);
		temp.push_back((gsl_ran_gaussian_ziggurat(r, sigma[4]/sqrt_iter)) + mu[4]);
		temp.push_back((gsl_ran_gaussian_ziggurat(r, sigma[5]/sqrt_iter)) + mu[5]);
		rotrans.push_back(temp);
		new_coord = Coord.rototranslate(Lig->xyz, Lig, temp[0], temp[1], temp[2], temp[3], temp[4], temp[5]);
		temp.clear();
		energies.push_back(Energy->compute_ene_softcore_solvation(Rec, Lig, Input, new_coord));
	}
	delete Energy;
	this->sort_array();
	time_f=clock();
	printf("Time elapsed for evaluation: %7.2f msecs.\n", double(((time_f-time_i)*1000)/CLOCKS_PER_SEC));
*/
}


double Deal::evaluate_energy(QVector<double> temp){
	double energy;
	vector<vector<double> > new_coord;
	new_coord = this->rototranslate(Lig, temp.at(0), temp.at(1), temp.at(2), temp.at(3), temp.at(4), temp.at(5));
	energy = this->compute_ene_softcore_solvation(Rec, Lig, Input, new_coord);
	return (energy);
}


void Deal::update_mu(int top_ranked){
	for (int i=0; i< top_ranked; i++){
			mu[0] += rotrans[i][0];
			mu[1] += rotrans[i][1];
			mu[2] += rotrans[i][2];
			mu[3] += rotrans[i][3];
			mu[4] += rotrans[i][4];
			mu[5] += rotrans[i][5];
		}
		for (int i=0; i<6; i++){
			mu[i] = mu[i]/((top_ranked+1)*1.0);
		}
}

Deal::~Deal() {
	gsl_rng_free(r);
}

vector<int> Deal::define_chunck_sizes(void){
	int slice = int(population_size/QThread::idealThreadCount());
	int count=0;
	vector<int> chuncks;
	chuncks.push_back(0);
	while(count < population_size){
		count += (slice+1);
		if (count > population_size){
			count=population_size;
		}
		chuncks.push_back(count);
	}
	return(chuncks);
}

vector<vector<double> >Deal::rototranslate(Mol2* Lig, double alpha, double beta, double gamma, double transx, double transy, double transz){
	vector<vector<double> >new_coordinates;
	vector<double> txyz;
	vector<double> COM = this->compute_com(Lig);
	double x, y, z;
	for(int i=0; i < Lig->N ; i++){
		x=Lig->xyz[i][0]-COM[0];
		y=Lig->xyz[i][1]-COM[1];
		z=Lig->xyz[i][2]-COM[2];
		txyz.push_back((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx+COM[0]);
		txyz.push_back((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy + COM[1]);
		txyz.push_back((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz + COM[2]);
		new_coordinates.push_back(txyz);
		txyz.clear();
	}
	return(new_coordinates);
}

vector<double> Deal::compute_com(Mol2 *Cmol){
	double centerx=0.0;
	double centery=0.0;
	double centerz=0.0;
	double totalmass=0.0;
	vector<double>com;
	for(int i=0; i<Cmol->N; i++){
		centerx+= (Cmol->masses[i]*Cmol->xyz[i][0]);
		centery+= (Cmol->masses[i]*Cmol->xyz[i][1]);
		centerz+= (Cmol->masses[i]*Cmol->xyz[i][2]);
		totalmass+= Cmol->masses[i];
	}
	com.push_back(centerx/totalmass);
	com.push_back(centery/totalmass);
	com.push_back(centerz/totalmass);
	return(com);
}

double Deal::compute_ene_softcore_solvation(Mol2* Rec, Mol2* Lig, PARSER* Input, vector<vector<double> >coords){
	double rec_solv = 0.00;
	double lig_solv = 0.00;
	double vdw=0.0;
	double elec=0.0;
	double energy=0.0;
	double dij, dij2, deff, rij, eij, acoef, bcoef, lig_solv_affinity, lig_solv_distf, rec_solv_affinity, rec_solv_distf;
	switch(Input->scoring_function){
	case 0:
		for (int i=0; i<Rec->N; i++){
				rec_solv_affinity = (0.25*((Rec->charges[i])*(Rec->charges[i])))-0.005;
				for (int j=0; j< Lig->N; j++){

					//! electrostactic energy (softcore)

					dij2 = distance_squared(Rec->xyz[i][0], coords[j][0],Rec->xyz[i][1], coords[j][1], Rec->xyz[i][2], coords[j][2]);
					deff = pow(((dij2*dij2*dij2)+Input->deltaij_es6), (1.0/3.0));
					elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*deff);

					//! VDW energy (softcore)

					rij = Rec->radii[i]*Lig->radii[j];
					eij = sqrt(Rec->epsilons[i]*Lig->epsilons[j]);
					acoef = (4096.00)*eij*(rij*rij*rij*rij*rij*rij);
					bcoef = (128.00)*eij*(rij*rij*rij);
					vdw+= ((acoef/(((dij2*dij2*dij2)+Input->deltaij6)*((dij2*dij2*dij2)+Input->deltaij6))) - (bcoef/((dij2*dij2*dij2)+Input->deltaij6)));

					//! Solvation energy

					lig_solv_affinity = (0.25*((Lig->charges[j])*(Lig->charges[j])))-0.005;
					lig_solv_distf = ((4.0/3.0)*PI*(Rec->radii[i]*Rec->radii[i]*Rec->radii[i]));
					lig_solv_distf = (lig_solv_distf * exp(-deff/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
					rec_solv_distf = ((4.0/3.0)*PI*(Lig->radii[j]*Lig->radii[j]*Lig->radii[j]));
					rec_solv_distf = (rec_solv_distf * exp(-deff/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);

					rec_solv+= rec_solv_affinity*rec_solv_distf;
					lig_solv+= lig_solv_affinity*lig_solv_distf;
				}
			}
		energy=elec+vdw+rec_solv+lig_solv;
		break;

	case 1:
		for (int i=0; i<(Rec->N); i++){
			rec_solv_affinity = (0.25*((Rec->charges[i])*(Rec->charges[i])))-0.005;
			for (int j=0; j< Lig->N; j++){
				//! electrostactic energy (softcore)

				dij2 = distance_squared(Rec->xyz[i][0], coords[j][0],Rec->xyz[i][1], coords[j][1], Rec->xyz[i][2], coords[j][2]);
				deff = pow(((dij2*dij2*dij2)+Input->deltaij_es6), (1.0/3.0));
				elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(Input->diel*deff);

				//! VDW energy (softcore)

				rij = Rec->radii[i]*Lig->radii[j];
				eij = sqrt(Rec->epsilons[i]*Lig->epsilons[j]);
				acoef = (4096.00)*eij*(rij*rij*rij*rij*rij*rij);
				bcoef = (128.00)*eij*(rij*rij*rij);
				vdw+= ((acoef/(((dij2*dij2*dij2)+Input->deltaij6)*((dij2*dij2*dij2)+Input->deltaij6))) - (bcoef/((dij2*dij2*dij2)+Input->deltaij6)));
				}
			}
		energy=elec+vdw;
		break;

	case 2:
		for (int i=0; i<Rec->N; i++){
			rec_solv_affinity = (0.25*((Rec->charges[i])*(Rec->charges[i])))-0.005;
			for (int j=0; j<Lig->N; j++){
				rij = Rec->radii[i]+Lig->radii[j];
				eij = sqrt(Rec->epsilons[i]*Lig->epsilons[j]);
				dij = distance(Rec->xyz[i][0], coords[j][0],Rec->xyz[i][1], coords[j][1], Rec->xyz[i][2], coords[j][2]);
				elec+= (332.0*Rec->charges[i]*Lig->charges[j])/(dij*Input->diel);
				dij2 = dij*dij;
				bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
				acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
				vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));

				deff = pow(((dij2*dij2*dij2)+Input->deltaij_es6), (1.0/3.0));
				lig_solv_affinity = (0.25*((Lig->charges[j])*(Lig->charges[j])))-0.005;
				lig_solv_distf = ((4.0/3.0)*PI*(Rec->radii[i]*Rec->radii[i]*Rec->radii[i]));
				lig_solv_distf = (lig_solv_distf * exp(-deff/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
				rec_solv_distf = ((4.0/3.0)*PI*(Lig->radii[j]*Lig->radii[j]*Lig->radii[j]));
				rec_solv_distf = (rec_solv_distf * exp(-deff/(2*(Input->sigma*Input->sigma))))/(Input->sigma*Input->sigma*Input->sigma);
				rec_solv+= rec_solv_affinity*rec_solv_distf;
				lig_solv+= lig_solv_affinity*lig_solv_distf;
			}
		}
		energy=elec+vdw+rec_solv+lig_solv;
		break;

	case 3:
		for (int i=0; i<(Rec->N); i++){
				for (int j=0; j<(Lig->N); j++){
					rij = Rec->radii[i]+Lig->radii[j];
					eij = sqrt(Rec->epsilons[i]*Lig->epsilons[j]);
					dij = distance(Rec->xyz[i][0], coords[j][0],Rec->xyz[i][1], coords[j][1], Rec->xyz[i][2], coords[j][2]);
					elec+= (332.0*Rec->charges[i]*Lig->charges[j])/dij;
					dij2 = dij*dij;
					bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
					acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
					vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));
				}
			}
		energy = elec+vdw;
		break;
	}
	return(energy);
}

double Deal::distance_squared (double x1, double x2, double y1, double y2, double z1, double z2) {
	return (((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))); }

double Deal::distance (double x1, double x2, double y1, double y2, double z1, double z2) {
	return ( sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) ); }


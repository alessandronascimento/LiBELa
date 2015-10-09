#include<iostream>
#include<zlib.h>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<string>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

  string infile;
  int c;
  int nbins=1000;
  double ene_lower_limit=-30.0;
  double ene_upper_limit=30.0;
  double rmsd_limit = 6.0;
  double Temp = 100.0;
  int N=0;
  

  if (argc < 2){
    printf("Usage %s -i <inputfile> -n <number_of_bins> -t <temperature> -l <energy_lower_limit> -u <energy_upper_limit> -r rmsd_limit [-h]\n", argv[0]);
	exit(1);
  }

  while ((c = getopt(argc, argv, "n:i:t:l:u:r:h")) != -1)
	switch (c){
	case 'n':
		nbins = atoi(optarg);
		break;
	case 'i':
		infile = string(optarg);
		break;
	case 'h':
		printf("Usage %s -i <inputfile> -n <number_of_bins> [-h]\n", argv[0]);
		break;
		exit(1);
	case '?':
		printf("Usage %s -i <inputfile> -n <number_of_bins> [-h]\n", argv[0]);
		break;
		exit(1);
	case 't':
		Temp = double(atof(optarg));
		break;
	case 'l':
		ene_lower_limit = double(atof(optarg));
		break;
	case 'u':
		ene_upper_limit = double(atof(optarg));
		break;
	case 'r':
		rmsd_limit = double(atof(optarg));
		break;	
	}

  gzFile inpfile = gzopen(infile.c_str(), "r");

  
  char str[150];

  int matrix[nbins][nbins];
  float energies[nbins][nbins];
  float binding_energies[nbins][nbins];

  for (int i=0; i<nbins; i++){
	for (int j=0; j<nbins; j++){
		matrix[i][j] = 0;
		energies[i][j] = 9999.;
	}
  }

 double rmsd_step= rmsd_limit/nbins;
 double ene_step = (ene_upper_limit-ene_lower_limit)/nbins;
  
  for (int i=0; i<4; i++){ 
	  gzgets(inpfile, str, 150);
  }

  float ene, rmsd, iene, dtmp, energy, crmsd;
  int tint;
  vector<float> conformer_energies;
  bool found_conf;
  float bene;

  while (! gzeof(inpfile)){
	gzgets(inpfile, str, 150);
	sscanf(str, "%d %f %f %f %f %f %f %f %f %d %f", &tint, &ene, &rmsd, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &tint, &iene);
	bene = ene-iene; 								// binding energy = total_energy - internal energy;

	energy=((ene)-ene_lower_limit)/ene_step;
	crmsd=(rmsd/rmsd_step);
	matrix[int(crmsd)][int(energy)]++;
	if ((ene) < energies[int(crmsd)][int(energy)]){
		energies[int(crmsd)][int(energy)] = ene;
		binding_energies[int(crmsd)][int(energy)] = bene;
	}
    N++;

	found_conf = false;
	for (unsigned i=0; i<conformer_energies.size(); i++){
		if (conformer_energies[i] == iene){
			found_conf = true;
		}
	}

	if (!found_conf){
		conformer_energies.push_back(iene);
	}
  }
  
  printf("Found %d sampled conformations in input file.\n", N);

  gzclose(inpfile);

  int count=0;
  long double Sum_pi_ln_pi=0.0;
  long double Zcomp = 0.0;
  long double dH = 0.0;
  for (int i=0; i<nbins; i++){
        for (int j=0; j<nbins; j++){
		Zcomp += exp(-energies[i][j]/(0.001987*Temp)) * (1.0*matrix[i][j]/N);
		dH += binding_energies[i][j]*(1.0*matrix[i][j]/N);
        if (matrix[i][j] > 0){
            count++;
            Sum_pi_ln_pi += (1.0*matrix[i][j]/N) * log(1.0*matrix[i][j]/N);
        }
	}
  }
	
 long double Zlig = 0.0;
 for (unsigned i=0; i<conformer_energies.size(); i++){
	Zlig += exp(-conformer_energies[i]/(0.001987*Temp));
 }
 long double dF = -0.001987*Temp*log(Zcomp/Zlig);

 printf("Z_comp = %10.6LG\n", Zcomp);
 printf("Z_lig  = %10.6LG for %d ligand conformers found.\n", Zlig, conformer_energies.size());
 printf("F_lig  = %10.6LG kcal/mol\n", dF);
 printf("dH     = %10.6LG kcal/mol\n", dH);
 printf("TdS    = %10.6LG kcal/mol\n", -(dF-dH));
 printf("S      = %10.6LG kcal/mol\n", -0.001987*Temp*Sum_pi_ln_pi);
 printf("Number of states: %6d out of %6d possible states (%5.3f%%).\n", count, (nbins*nbins), (count*1.0/(nbins*nbins)));
 
 return 0;
}

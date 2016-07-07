#include<cstdlib>
#include<cstdio>
#include<zlib.h>
#include<string>
#include<cmath>

using namespace std;

int main(int argc, char* argv[]){

    string infile;
    int c;

    if (argc < 2){
      printf("Usage %s -i <inputfile> [-h]\n", argv[0]);
      exit(1);
    }

    while ((c = getopt(argc, argv, "i:h")) != -1)
      switch (c){
      case 'i':
          infile = string(optarg);
          break;
      case 'h':
          printf("Usage %s -i <inputfile> [-h]\n", argv[0]);
          break;
          exit(1);
      case '?':
          printf("Usage %s -i <inputfile> [-h]\n", argv[0]);
          break;
          exit(1);
      }

    gzFile inpfile = gzopen(infile.c_str(), "r");

    char str[200];
    int Nrot=0;
    char tstr[10];
    float ene, rmsd, dtmp, iene, T;
    int step, tint, nsteps=0;
    double average = 0.0, fluct = 0.0;
    double ftmp, entropy;
    string parsing_string = "%d %f %f %f %f %f %f %f %f %d %f";

// Parsing two first line (comments)
    gzgets(inpfile, str, 200);
    gzgets(inpfile, str, 200);
    sscanf(str, "%s %s %s %s %s %s %s %f", tstr,tstr, tstr, tstr, tstr, tstr, tstr, &T);
    printf("Temperature: %10.3f K\n", T);
    gzgets(inpfile, str, 200);
    sscanf(str, "%s %s %d", tstr, tstr, &Nrot);
    printf("Rotatable bonds: %10d\n", Nrot);

    for (int i=0; i< Nrot; i++){
        parsing_string += " %f";
    }

    gzgets(inpfile, str, 200);
    gzgets(inpfile, str, 200);

//    printf("Computing average energies. This may take a while...\n");

    while (! gzeof(inpfile)){
        gzgets(inpfile, str, 200);
        sscanf(str, parsing_string.c_str(), &step, &ene, &rmsd, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &tint, &iene, &dtmp, &dtmp, &dtmp);
        average += (ene-iene);
        nsteps++;
    }
    gzclose(inpfile);

    average = average/nsteps;
    printf("%20s %10.3f kcal/mol\n", "Average energy:", average);
//    printf("Computing flucutation. This may take a while...\n");

    inpfile = gzopen(infile.c_str(), "r");

    for (int i=0; i<5; i++){
        gzgets(inpfile, str, 200);
    }

    while (! gzeof(inpfile)){
        gzgets(inpfile, str, 200);
        sscanf(str, parsing_string.c_str(), &step, &ene, &rmsd, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &tint, &iene, &dtmp, &dtmp, &dtmp);
        ftmp = (ene-iene) - average;
        fluct += exp(ftmp/(0.0019858775203792202*T));
    }

    fluct = fluct/nsteps;
//    printf("<e^DE/kT> = %10.3f\n", fluct);
    entropy = 0.0019858775203792202*T*log(fluct);
    printf("%20s %10.3f kcal/mol.\n", "-TDS:", entropy);
    printf("%20s %10.3f kcal/mol.\n", "DG:", average+entropy);
    printf("%20s %10.3f nM.\n", "Ligand Affinity:", (1./exp(-(average+entropy)/(0.0019858775203792202*T)))*1E9);

    return 0;

}

//###################################################################################################################################################################

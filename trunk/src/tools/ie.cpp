#include<cstdlib>
#include<cstdio>
#include<zlib.h>
#include<string>
#include<cmath>

using namespace std;

int main(int argc, char* argv[]){

    string compfile, ligfile;
    int c;

    if (argc <= 3){
      printf("Usage %s -c <complex_file> -l <ligand_file> [-h]\n", argv[0]);
      exit(1);
    }

    while ((c = getopt(argc, argv, "c:l:h")) != -1)
      switch (c){
      case 'c':
          compfile = string(optarg);
          break;
      case 'l':
          ligfile = string(optarg);
          break;
      case 'h':
          printf("Usage %s -c <complex_file> -l <ligand_file> [-h]\n", argv[0]);
          break;
          exit(1);
      case '?':
          printf("Usage %s -c <complex_file> -l <ligand_file> [-h]\n", argv[0]);
          break;
          exit(1);
      }

    gzFile inpfile = gzopen(compfile.c_str(), "r");

    char str[200];
    int Nrot=0;
    char tstr[10];
    float ene, rmsd, dtmp, iene, T;
    int step, tint, nsteps=0;
    double average_comp = 0.0, average_lig=0.0,  average_bind=0.0, fluct = 0.0, k=0.0019858775203792202;
    double entropy;
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
        average_comp += ene;
        nsteps++;
    }
    gzclose(inpfile);

    average_comp = average_comp/nsteps;
    printf("%20s %10.3f kcal/mol\n", "Average complex energy:", average_comp);

//    printf("Computing flucutation. This may take a while...\n");

    inpfile = gzopen(ligfile.c_str(), "r");
    nsteps=0;
    for (int i=0; i<5; i++){
        gzgets(inpfile, str, 200);
    }

    while (! gzeof(inpfile)){
        gzgets(inpfile, str, 200);
        sscanf(str, parsing_string.c_str(), &step, &ene, &rmsd, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &tint, &iene, &dtmp, &dtmp, &dtmp);
        average_lig += ene ;
        nsteps++;
    }

    gzclose(inpfile);

    average_lig = average_lig/nsteps;
    printf("%20s %10.3f kcal/mol\n", "Average ligand energy:", average_lig);

    average_bind = average_comp-average_lig;
    printf("%20s %10.3f kcal/mol\n", "Average binding energy:", average_bind);

    inpfile = gzopen(compfile.c_str(), "r");
    nsteps=0;
    for (int i=0; i<5; i++){
        gzgets(inpfile, str, 200);
    }

    while (! gzeof(inpfile)){
        gzgets(inpfile, str, 200);
        sscanf(str, parsing_string.c_str(), &step, &ene, &rmsd, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &dtmp, &tint, &iene, &dtmp, &dtmp, &dtmp);
        fluct += exp((ene-iene-average_bind)/(k*T));
        nsteps++;
    }
    gzclose(inpfile);


    fluct = fluct/nsteps;
    entropy = k*T*log(fluct);
    printf("%20s %10.3f kcal/mol.\n", "-TDS:", entropy);
    printf("%20s %10.3f kcal/mol.\n", "DG:", average_bind+entropy);
    printf("%20s %10.3f nM.\n", "Ligand Affinity:", (1./exp(-(average_bind+entropy)/(k*T)))*1E9);

    return 0;

}

//###################################################################################################################################################################

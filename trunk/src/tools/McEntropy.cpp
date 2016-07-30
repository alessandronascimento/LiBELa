#include<iostream>
#include<zlib.h>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<string>
#include <vector>
#include "../LiBELa/PARSER.h"
#include "../LiBELa/Mol2.h"
#include "../LiBELa/COORD_MC.h"

using namespace std;

int main(int argc, char* argv[]){

    int MAX_STEPS = 10000000; // 10^6
    string infile, ligfile;
    int c;
    int rot_bins=360, trans_bins=60, translation_window=30;
    double Temp = 100.0;
    PARSER* Input = new PARSER;

    if (argc < 2){
        printf("Usage %s -i <inputfile> -t <number of bins for translation> -r <number of bins for rotation> -l <ligand mol2 file> [-h]\n", argv[0]);
        exit(1);
    }

    while ((c = getopt(argc, argv, "i:t:r:l:h")) != -1)
        switch (c){
        case 'i':
            infile = string(optarg);
            break;
        case 'h':
            printf("Usage %s -i <inputfile> -t <number of bins for translation> -r <number of bins for rotation> [-h]\n", argv[0]);
            break;
            exit(1);
        case '?':
            printf("Usage %s -i <inputfile> -t <number of bins for translation> -r <number of bins for rotation> [-h]\n", argv[0]);
            break;
            exit(1);
        case 't':
            trans_bins = atoi(optarg);
            break;
        case 'r':
            rot_bins = atoi(optarg);
            break;
        case 'l':
            ligfile = string(optarg);
            break;
        }

    Mol2* Lig = new Mol2(Input, ligfile);
    COORD_MC* Coord = new COORD_MC;
    vector<double> com = Coord->compute_com(Lig);

    delete Coord;
    delete Lig;
    delete Input;

    float translation_step = translation_window*1.0/translation_bins; //typically 0.5 Ang
    float rotation_step = 360.0/rot_bins;                             //typically 1.0 degree

    float hist_x[trans_bins];
    float hist_y[trans_bins];
    float hist_z[trans_bins];

    for (unsigned i=0; i< trans_bins; i++){
        hist_x[i]= 0.0;
        hist_y[i]= 0.0;
        hist_z[i]= 0.0;
    }

    float hist_alpha[rot_bins];
    float hist_beta[rot_bins];
    float hist_gamma[rot_bins];

    for (unsigned i=0; i< rot_bins; i++){
        hist_alpha[i] = 0.0;
        hist_beta[i] = 0.0;
        hist_gamma[i] = 0.0;
    }

    gzFile inpfile = gzopen(infile.c_str(), "r");
    char str[150];
    char tstr[20];
    float tfloat;
    int tint, n_rot;
    float x, y, z, alpha, beta, gamma;

    gzgets(inpfile, str, 150);
    gzgets(inpfile, str, 150);
    sscanf(str, "%s %f %s %f %s %f %s %f", tstr, &tfloat, tstr, &tfloat, tstr, &tfloat, tstr, &Temp);
    gzgets(inpfile, str, 150);
    sscanf(str, "%s %s %d", tstr, tstr, &n_rot);
    gzgets(inpfile, str, 150);
    gzgets(inpfile, str, 150);

    float** torsion = new float*[MAX_STEPS];

    int count=0;

    while (! gzeof(inpfile)){
        gzgets(inpfile, str, 150);
        count++;
        if (count > MAX_STEPS){
            printf("Maximum number of steps exceeded. Please change this parameter in the source file and recompile.\n");
            exit(1);
        }

        torsion[count-1] = new float[n_rot];

        sscanf(str, "%d %f %f %f %f %f %f %f %f %f %d %f %{%f%}", &tint, &tfloat, &tfloat, &x, &y, &z, &alpha,
               &beta, &gamma, &tint, &tfloat, torsion[count-1]);

        hist_x[int((x-com[0]-(translation_window*1.0/2))/(translation_step))] += 1.0;
        hist_y[int((y-com[1]-(translation_window*1.0/2))/(translation_step))] += 1.0;
        hist_z[int((z-com[2]-(translation_window*1.0/2))/(translation_step))] += 1.0;

        hist_alpha[int(alpha/rotation_step)] += 1.0;
        hist_gamma[int(alpha/rotation_step)] += 1.0;
        hist_beta[int(alpha/rotation_step)] += 1.0;

        /*
         * process torsions here
         */
    }

    double Strans=0.0;
    double Srot=0.0;

    for (unsigned i=0; i< trans_bins; i++){
        hist_x[i] = hist_x/count;
        Strans += hist_x[i] * log(hist_x[i]);
        hist_y[i] = hist_y/count;
        Strans += hist_y[i] * log(hist_y[i]);
        hist_z[i] = hist_z/count;
        Strans += hist_z[i] * log(hist_z[i]);
    }

    for (unsigned i=0; i< rot_bins; i++){
        hist_alpha[i] = hist_alpha[i]/count;
        Srot += hist_alpha[i] * log(hist_alpha[i]);
        hist_beta[i] = hist_beta[i]/count;
        Srot += hist_beta[i] * log(hist_beta[i]);
        hist_gamma[i] = hist_gamma[i]/count;
        Srot += hist_gamma[i] * log(hist_gamma[i]);
    }

    printf("Strans: *10.3f\n", Strans);
    printf("Srot: %10.3f\n", Srot);

    return 0;
}

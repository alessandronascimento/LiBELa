#include<iostream>
#include<zlib.h>
#include<cstdlib>
#include<cstdio>
#include<stdio.h>
#include<cmath>
#include<string>
#include <vector>
#include<sstream>
#include "../LiBELa/PARSER.cpp"
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/COORD_MC.cpp"

using namespace std;

int main(int argc, char* argv[]){

    string infile, ligfile;
    int c;
    int rot_bins=360, trans_bins=60, translation_window=30;
    double Temp = 100.0;
    double k=0.001987;
    PARSER* Input = new PARSER;
    FILE* hist_output_trans;
    FILE* hist_output_rot;

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
    printf("Ligand COM: %7.3f %7.3f %7.3f\n", com[0], com[1], com[2]);

    delete Coord;
    delete Lig;
    delete Input;

    float translation_step = translation_window*1.0/trans_bins;       //typically 0.5 Ang
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


    gzFile inpfile = gzopen(infile.c_str(), "r");
    char str[250];
    char tstr[20];
    float tfloat;
    int tint, n_rot;
    float x, y, z, alpha, beta, gamma;

    gzgets(inpfile, str, 250);
    gzgets(inpfile, str, 250);
    sscanf(str, "%s %f %s %f %s %f %s %f", tstr, &tfloat, tstr, &tfloat, tstr, &tfloat, tstr, &tfloat);
    Temp = double(tfloat);
    gzgets(inpfile, str, 250);
    sscanf(str, "%s %s %d", tstr, tstr, &n_rot);
    gzgets(inpfile, str, 250);
    gzgets(inpfile, str, 250);

    float** hist_torsions = new float*[n_rot];

    for (int i=0; i < n_rot; i++){
        hist_torsions[i] = new float[rot_bins];
        for (unsigned j=0; j< rot_bins; j++){
            hist_torsions[i][j] = 0.0;
            hist_alpha[j] = 0.0;
            hist_beta[j] = 0.0;
            hist_gamma[j] = 0.0;
        }
    }


    printf("Parsing file %s. Temp = %7.3f K. N_rot = %3d\n", infile.c_str(), Temp, n_rot);


    int count=0;
    string line;


    while (!gzeof(inpfile)){
        gzgets(inpfile, str, 250);
        count++;

        float* torsion = new float[n_rot];

        line = string(str);
        istringstream ss(line);
        ss >> tint >> tfloat >> tfloat >> x >> y >> z >> alpha >> beta >> gamma >> tint >> tfloat;
        for (int i=0; i< n_rot; i++){
            ss >> torsion[i];
        }

        hist_x[int(round((x-com[0]+(translation_window*1.0/2.))/(translation_step)))] += 1.0;
        hist_y[int(round((y-com[1]+(translation_window*1.0/2.))/(translation_step)))] += 1.0;
        hist_z[int(round((z-com[2]+(translation_window*1.0/2.))/(translation_step)))] += 1.0;

        hist_alpha[int(round(alpha/rotation_step))] += 1.0;
        hist_gamma[int(round(alpha/rotation_step))] += 1.0;
        hist_beta[int(round(alpha/rotation_step))] += 1.0;

        for (int i=0; i< n_rot; i++){
            hist_torsions[i][int(round(torsion[i]/rotation_step))] += 1.0;
        }
        delete [] torsion;
    }

    gzclose(inpfile);

    hist_output_trans = fopen("histogram_translation.dat","w");
    hist_output_rot = fopen("histogram_rotation.dat","w");

    fprintf(hist_output_trans, "%7.7s %7.7s %7.7s %7.7s %7.7s %7.7s\n", "x", "h_x", "y", "h_y", "z", "h_z");
    fprintf(hist_output_rot,   "%7.7s %7.7s %7.7s %7.7s %7.7s %7.7s ","alpha", "h_alpha", "beta", "h_beta", "gamma", "h_gamma");
    for (int i=0; i<n_rot; i++){
        sprintf(str, "tor[%2d]", i+1);
        fprintf(hist_output_rot, " %7.7s ", str);
    }
    fprintf(hist_output_rot, "\n");


    double Strans=0.0;
    double Srot=0.0;
    double Stor=0.0;

    for (unsigned i=0; i< trans_bins; i++){
        hist_x[i] = hist_x[i]/count;
        if (hist_x[i] > 0.0){
            Strans += hist_x[i] * log(hist_x[i]);
        }
         fprintf(hist_output_trans, "%7.4f %7.4f ", i+1, (i*1.0*translation_step)+com[0]-(translation_window*1.0/2.), hist_x[i]);


        hist_y[i] = hist_y[i]/count;
        if (hist_y[i] > 0.0){
            Strans += hist_y[i] * log(hist_y[i]);
        }
        fprintf(hist_output_trans, "%7.4f %7.4f ", (i*1.0*translation_step)+com[1]-(translation_window*1.0/2.), hist_y[i]);


        hist_z[i] = hist_z[i]/count;
        if (hist_z[i] > 0.0){
            Strans += hist_z[i] * log(hist_z[i]);
        }
        fprintf(hist_output_trans, "%7.4f %7.4f\n", (i*1.0*translation_step)+com[2]-(translation_window*1.0/2.), hist_z[i]);
    }

    for (unsigned i=0; i< rot_bins; i++){
        hist_alpha[i] = hist_alpha[i]/count;
        if (hist_alpha[i]> 0.0){
            Srot += hist_alpha[i] * log(hist_alpha[i]);
        }
        fprintf(hist_output_rot, "%7.1f %7.4f " , i*rotation_step, hist_alpha[i]);



        hist_beta[i] = hist_beta[i]/count;
        if (hist_beta[i]> 0.0){
            Srot += hist_beta[i] * log(hist_beta[i]);
        }
        fprintf(hist_output_rot, "%7.1f %7.4f " , i*rotation_step, hist_beta[i]);

        hist_gamma[i] = hist_gamma[i]/count;
        if (hist_gamma[i]> 0.0){
            Srot += hist_gamma[i] * log(hist_gamma[i]);
        }
        fprintf(hist_output_rot, "%7.1f %7.4f " , i*rotation_step, hist_gamma[i]);

        for (int j=0; j< n_rot; j++){
            hist_torsions[j][i] = hist_torsions[j][i]/count;
            if (hist_torsions[j][i] > 0.0){
                Stor += hist_torsions[j][i] * log(hist_torsions[j][i]);
            }
            fprintf(hist_output_rot, "%7.1f %7.4f" , i*rotation_step, hist_torsions[j][i]);
        }
        fprintf(hist_output_rot, "\n");
    }



    /*
     * Cleaning up....
     */

    fclose(hist_output_trans);
    fclose(hist_output_rot);


    for (int i=0; i < n_rot; i++){
         delete [] hist_torsions[i];
    }
    delete [] hist_torsions;

    Strans = -k*Strans;
    Srot = -k*Srot;
    Stor = -k*Stor;
    printf("Strans: %10.3f\n", Strans);
    printf("Srot: %10.3f\n", Srot);
    printf("Stor: %10.3f\n", Stor);
    double S = Strans + Srot + Stor;
    printf("-TS = %10.3f\n", -Temp*S);

    return 0;
}

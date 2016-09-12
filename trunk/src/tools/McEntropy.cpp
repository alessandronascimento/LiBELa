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

    string infile, ligfile, prefix="McEntropy";
    int c;
    int rot_bins=360, trans_bins=60, translation_window=30;
    double Temp = 100.0;
    double k=0.001987;
    PARSER* Input = new PARSER;
    FILE* hist_output_trans;
    FILE* hist_output_rot;

    printf("******************************************************************************************\n");
    printf("*                                                                                        *\n");
    printf("*              McEntropy - Entropy Estimation From 1st Order Approximation               *\n");
    printf("*        Written by Alessandro S. Nascimento and Joao Victor S. Cunha   /  2016          *\n");
    printf("*                      University of Sao Paulo - USP - Brazil                            *\n");
    printf("*                                                                                        *\n");
    printf("*                             asnascimento@ifsc.usp.br                                   *\n");
    printf("*                                                                                        *\n");
    printf("******************************************************************************************\n");


    if (argc < 2){
        printf("Usage %s -i <inputfile> -t <number of bins for translation> -r <number of bins for rotation> -l <ligand mol2 file> -o <output_prefix>[-h]\n", argv[0]);
        exit(1);
    }

    while ((c = getopt(argc, argv, "i:t:r:l:o:h")) != -1)
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
        case 'o':
            prefix = string(optarg);
            break;
        }

    Mol2* Lig = new Mol2(Input, ligfile);
    COORD_MC* Coord = new COORD_MC;
    vector<double> com = Coord->compute_com(Lig);
    printf("* Ligand COM: %7.3f %7.3f %7.3f                                                    *\n", com[0], com[1], com[2]);
    delete Coord;
    delete Lig;
    delete Input;

    double translation_step = translation_window*1.0/trans_bins;       //typically 0.5 Ang
    double rotation_step = 360.0/rot_bins;                             //typically 1.0 degree

    vector<double> hist_x(trans_bins);
    vector<double> hist_y(trans_bins);
    vector<double> hist_z(trans_bins);

    for (unsigned i=0; i< trans_bins; i++){
        hist_x.push_back(0.0);
        hist_y.push_back(0.0);
        hist_z.push_back(0.0);
    }

    vector<double> hist_alpha(rot_bins);
    vector<double> hist_beta(rot_bins);
    vector<double> hist_gamma(rot_bins);


/*
 *
 * Parsing data from trajectory header...
 *
 */

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
    for (int i=0; i< n_rot; i++){
        hist_torsions[i] = new float[rot_bins];
    }

    for (int i=0; i<n_rot; i++){
        for (int j=0; j< rot_bins; j++){
            hist_torsions[i][j] = 0.0;
        }
    }

    for (int i=0; i< rot_bins; i++){
        hist_alpha.push_back(0.0);
        hist_beta.push_back(0.0);
        hist_gamma.push_back(0.0);
    }

    printf("* Parsing file %s. Temp = %7.3f K. N_rot = %3d                      *\n", infile.c_str(), Temp, n_rot);


    int count=0;
    string line;

/*
 *
 * Parsing trajectory data
 *
 */

   while (!gzeof(inpfile)){
        gzgets(inpfile, str, 250);
        count++;

        float* torsion = new float[n_rot];

        line = string(str);
        istringstream ss(line);
        ss >> tint >> tfloat >> tfloat >> x >> y >> z >> alpha >> beta >> gamma >> tint >> tfloat;
        for (int i=0; i< n_rot; i++){
            ss >> torsion[i];
            if (torsion[i] < 0){
                torsion[i] = torsion[i]+360.;
            }
        }

        hist_x[int(round((x-com[0]+(translation_window*1.0/2.))/(translation_step)))] += 1.0;
        hist_y[int(round((y-com[1]+(translation_window*1.0/2.))/(translation_step)))] += 1.0;
        hist_z[int(round((z-com[2]+(translation_window*1.0/2.))/(translation_step)))] += 1.0;

        hist_alpha[int(round(alpha/rotation_step))] += 1.0;
        hist_gamma[int(round(beta/rotation_step))] += 1.0;
        hist_beta[int(round(gamma/rotation_step))] += 1.0;

        for (int i=0; i< n_rot; i++){
            int angle = round(torsion[i]/rotation_step);
            if (angle < 0 or angle > rot_bins){
                printf("ANGLE OFFSET: %d!\n", angle);
                printf("Check McEntropy.cpp file and recompile.\n");
                exit(1);
            }
            else{
                hist_torsions[i][angle] += 1.0;
            }
        }
        delete [] torsion;
    }


    gzclose(inpfile);

    hist_output_trans = fopen((prefix + "_histogram_translation.dat").c_str(),"w");
    hist_output_rot = fopen((prefix + "_histogram_rotation.dat").c_str(),"w");

    fprintf(hist_output_trans, "#%7.7s %7.7s %7.7s %7.7s %7.7s %7.7s\n", "x", "h_x", "y", "h_y", "z", "h_z");
    fprintf(hist_output_rot,   "#%7.7s %7.7s %7.7s %7.7s %7.7s %7.7s ","alpha", "h_alpha", "beta", "h_beta", "gamma", "h_gamma");

    for (int i=0; i<n_rot; i++){
        sprintf(str, "tor[%2d]", i+1);
        fprintf(hist_output_rot, " %7.7s ", str);
    }
    fprintf(hist_output_rot, "\n");

    double Strans=0.0;
    double Srot=0.0;
    double Stor=0.0;

/*
 *
 * Histograming....
 *
 */

    for (unsigned i=0; i< trans_bins; i++){
        hist_x[i] = hist_x[i]/count;
        if (hist_x[i] > 0.0){
            Strans += hist_x[i] * log(hist_x[i]);
        }
         fprintf(hist_output_trans, "%7.4f %7.4f ", (i*1.0*translation_step)+com[0]-(translation_window*1.0/2.), hist_x[i]);


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
    hist_x.clear();
    hist_y.clear();
    hist_z.clear();
    hist_alpha.clear();
    hist_beta.clear();
    hist_gamma.clear();

    for (int i=0; i< n_rot; i++){
        delete [] hist_torsions[i];
    }

    delete [] hist_torsions;

/*
 * Computing Shannon Entropies
*/

    Strans = -k*Strans;
    Srot = -k*Srot;
    Stor = -k*Stor;
    printf("*                                                                                        *\n");
    printf("* Translational Entropy: %10.4f kcal/(mol.K)                                         *\n", Strans);
    printf("* Rotational    Entropy: %10.4f kcal/(mol.K)                                         *\n", Srot);
    printf("* Torsional     Entropy: %10.4f kcal/(mol.K)                                         *\n", Stor);
    double S = Strans + Srot + Stor;
    printf("* Total         Entropy: %10.4f kcal/(mol.K)                                         *\n", S);
    printf("*                                                                                        *\n");
    printf("* -TS = %10.4f                                                                       *\n", -Temp*S);
    printf("*                                                                                        *\n");
    printf("******************************************************************************************\n");

    return 0;
}

#include<cstdio>
#include<cstdlib>
#include<getopt.h>
#include<string>
#include<sstream>
#include<zlib.h>
//#include<vector>
//#include<cmath>
//#include<time.h>
#include <gsl/gsl_statistics.h>

using namespace std;

int main (int argc, char *argv[]){

    int c;
    int ncycles=1E6;
    int stride=1;
    string infile;


    if (argc < 2){
        printf("Usage %s -i <inputfile> -n <data_size> -s <stride> [-h]\n", argv[0]);
        exit(1);
    }

    while ((c = getopt(argc, argv, "n:i:s:h")) != -1)
        switch (c){
        case 'n':
            ncycles= atoi(optarg);
            break;
        case 'i':
            infile = string(optarg);
            break;
        case 'h':
            printf("Usage %s -i <inputfile> -n <data_size> -s <stride> [-h]\n", argv[0]);
            break;
            exit(1);
        case '?':
            printf("Usage %s -i <inputfile> -n <data_size> -s <stride> [-h]\n", argv[0]);
            break;
            exit(1);
        case 's':
            stride = atoi(optarg);
            break;
        }


    // Reading data


    gzFile datain = gzopen(infile.c_str(), "r");
    char str[250];

    for (int i=0; i<5; i++){
        gzgets(datain, str, 250);
    }

    int count=0;
    string line;
    int tint;
    double ene, dtmp, x, y, z, alpha, beta, gamma;
    double data_ene[ncycles], data_x[ncycles], data_y[ncycles], data_z[ncycles], data_alpha[ncycles], data_beta[ncycles], data_gamma[ncycles];


    while (count < ncycles && !gzeof(datain)){
        gzgets(datain, str, 250);
        line = string(str);
        istringstream ss(line);
        ss >> tint >> ene >> dtmp >> x >> y >> z >> alpha >> beta >> gamma;
        data_ene[count] = ene;
        data_x[count] = x;
        data_y[count] = y;
        data_z[count] = z;
        data_alpha[count] = alpha;
        data_beta[count] = beta;
        data_gamma[count] = gamma;
        count++;
    }

    gzclose(datain);

    double autocorr_ene = gsl_stats_lag1_autocorrelation(data_ene, stride, count);

    double autocorr_x = gsl_stats_lag1_autocorrelation(data_x, stride, count);
    double autocorr_y = gsl_stats_lag1_autocorrelation(data_y, stride, count);
    double autocorr_z = gsl_stats_lag1_autocorrelation(data_z, stride, count);

    double autocorr_alpha = gsl_stats_lag1_autocorrelation(data_alpha, stride, count);
    double autocorr_beta = gsl_stats_lag1_autocorrelation(data_beta, stride, count);
    double autocorr_gamma = gsl_stats_lag1_autocorrelation(data_gamma, stride, count);

    printf("Autocorrelation for Energy: %7.4f\n", autocorr_ene);
    printf("Autocorrelation for X Translation: %7.4f\n", autocorr_x);
    printf("Autocorrelation for Y Translation: %7.4f\n", autocorr_y);
    printf("Autocorrelation for Z Translation: %7.4f\n", autocorr_z);
    printf("Autocorrelation for Rotation in alpha: %7.4f\n", autocorr_alpha);
    printf("Autocorrelation for Rotation in beta : %7.4f\n", autocorr_beta);
    printf("Autocorrelation for Rotation in gamma: %7.4f\n", autocorr_gamma);

    return 0;

}

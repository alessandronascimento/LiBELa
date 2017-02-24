#include<cstdio>
#include<cstdlib>
#include<getopt.h>
#include<string>
#include<sstream>
#include<zlib.h>
#include <gsl/gsl_statistics.h>

using namespace std;

double compute_autocorr(double* data, int n, int lag){
    double num=0.0, denom=0.0;
    double mean = gsl_stats_mean(data, 1, n);
    for (int i=0; i<(n-lag); i++){
        num += (data[i]-mean)*(data[i+lag]-mean);
        denom += (data[i]-mean)*(data[i]-mean);
    }
    return(num/denom);
}


int main (int argc, char *argv[]){

    int c;
    int ncycles=1E6;
    int lag = 1;
    int stride=1;
    string infile;


    if (argc < 2){
        printf("Usage %s -i <inputfile> -n <data_size> -s <stride> [-h] -k <lag>\n", argv[0]);
        exit(1);
    }

    while ((c = getopt(argc, argv, "n:i:s:k:h")) != -1)
        switch (c){
        case 'n':
            ncycles= atoi(optarg);
            break;
        case 'i':
            infile = string(optarg);
            break;
        case 'h':
            printf("Usage %s -i <inputfile> -n <data_size> -s <stride> [-h] -k <lag>\n", argv[0]);
            break;
            exit(1);
        case '?':
            printf("Usage %s -i <inputfile> -n <data_size> -s <stride> [-h] -k <lag>\n", argv[0]);
            break;
            exit(1);
        case 's':
            stride = atoi(optarg);
            break;
        case 'k':
            lag = atoi(optarg);
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
    double* data_ene = new double[ncycles];
    double* data_x = new double[ncycles];
    double* data_y = new double[ncycles];
    double* data_z = new double[ncycles];
    double* data_alpha = new double[ncycles];
    double* data_beta = new double[ncycles];
    double* data_gamma = new double[ncycles];


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


    double autocorr_ene = compute_autocorr(data_ene, count, lag);

    double autocorr_x = compute_autocorr(data_x, count, lag);
    double autocorr_y = compute_autocorr(data_y, count, lag);
    double autocorr_z = compute_autocorr(data_z, count, lag);

    double autocorr_alpha = compute_autocorr(data_alpha, count, lag);
    double autocorr_beta = compute_autocorr(data_beta, count, lag);
    double autocorr_gamma = compute_autocorr(data_gamma, count, lag);

    printf("Autocorrelation for Energy: %7.4f\n", autocorr_ene);
    printf("Autocorrelation for X Translation: %7.4f\n", autocorr_x);
    printf("Autocorrelation for Y Translation: %7.4f\n", autocorr_y);
    printf("Autocorrelation for Z Translation: %7.4f\n", autocorr_z);
    printf("Autocorrelation for Rotation in alpha: %7.4f\n", autocorr_alpha);
    printf("Autocorrelation for Rotation in beta : %7.4f\n", autocorr_beta);
    printf("Autocorrelation for Rotation in gamma: %7.4f\n", autocorr_gamma);

    return 0;

}

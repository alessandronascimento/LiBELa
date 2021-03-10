#include<cstdio>
#include<cstdlib>
#include<getopt.h>
#include<string>
#include<vector>
#include<cmath>
#include<time.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


using namespace std;


int main (int argc, char *argv[]){

    int c;
    int ncycles=1E6;
    string infile;
    double sigma;
    vector<double> datax, datay;
    float x,y;
    char str[80];


    if (argc < 2){
      printf("Usage %s -i <inputfile> -n <number_cycles> -e <estimated_error> [-h]\n", argv[0]);
      exit(1);
    }

    while ((c = getopt(argc, argv, "n:i:e:h")) != -1)
      switch (c){
      case 'n':
          ncycles= atoi(optarg);
          break;
      case 'i':
          infile = string(optarg);
          break;
      case 'h':
          printf("Usage %s -i <inputfile> -n <number_cycles> -e <estimated_error> [-h]\n", argv[0]);
          break;
          exit(1);
      case '?':
          printf("Usage %s -i <inputfile> -n <number_cycles> -e <estimated_error> [-h]\n", argv[0]);
          break;
          exit(1);
      case 'e':
          sigma = double(atof(optarg));
          break;
      }


    printf("*****************************************************************************\n");
    printf("*                                                                           *\n");
    printf("*                               McBootstrapper                              *\n");
    printf("* A sery simple program to compute averages/standard deviations in datasets *\n");
    printf("* using the Bootstrap method. Written by Alessandro Nascimento @ IFSC/USP   *\n");
    printf("*                                                                           *\n");
    printf("*****************************************************************************\n");


// Reading data



    FILE *datafile = fopen(infile.c_str(), "r");
    while(!feof(datafile)){
        fgets(str, 80, datafile);
        sscanf(str, "%f %f", &x, &y);
        datax.push_back(x);
        datay.push_back(y);
    }

    fclose(datafile);

    printf("The dataset has %d data points.\n", int(datax.size()));
    printf("Read dataset:\n");
    printf("*****************************************************************************\n");
    printf("%7.7s %7.7s\n", "[X]", "[Y]");
    for (unsigned i=0; i<datax.size()-1; i++){
        printf("%7.4f %7.4f\n", datax[i], datay[i]);
    }
    printf("*****************************************************************************\n");

    if (datax.size() != datay.size()){
        printf("Data inconsistency. Different number of values for sets. Please check!\n");
        exit(1);
    }
    int N = int(datax.size());

    double dataX[datax.size()];
    double dataY[datay.size()];

    double new_dataX[datax.size()];

    copy(datax.begin(), datax.end(), dataX);
    copy(datay.begin(), datay.end(), dataY);

    double original_r = gsl_stats_correlation(dataX, 1, dataY, 1, N);

    printf("Correlation Coefficient: %7.4f\n", original_r);

    double new_r;
    double sum_r = 0.0;
    double sum_r_squared = 0.0;


// Preparing for Bootstrap

//    const gsl_rng_type *T;
    gsl_rng *r;
    srand(rand());
    r = gsl_rng_alloc (gsl_rng_ranlxs2);

    FILE* output = fopen("Bootstrap.dat", "w");
    fprintf(output, "%10.10s %7.7s\n", "#cycle", "r");

    for (int i=0; i< ncycles; i++){
        for (int j=0; j< N; j++){
            new_dataX[j] = dataX[j] + gsl_ran_gaussian(r, sigma);
        }
        new_r = gsl_stats_correlation(new_dataX, 1, dataY, 1, N);
        sum_r += new_r;
        sum_r_squared += (new_r*new_r);
        fprintf(output, "%10d %7.4f\n", i, new_r);
    }

    gsl_rng_free(r);

    double average_r = sum_r / ncycles;
    double stdev = (sum_r_squared - ((sum_r*sum_r))/ncycles)/(ncycles-1);
    stdev = sqrt(stdev);

    fprintf(output, "#Average r = %7.4f +/- %7.4f\n", average_r, stdev);
    fclose(output);

    printf("Average correlation coefficient: %7.4f +/- %7.4f (N=%d).\n", average_r, stdev, ncycles);
    printf("*****************************************************************************\n");

    return 0;



}

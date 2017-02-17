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
    string infile;


    if (argc < 2){
      printf("Usage %s -i <inputfile> -n <data_size> [-h]\n", argv[0]);
      exit(1);
    }

    while ((c = getopt(argc, argv, "n:i:h")) != -1)
      switch (c){
      case 'n':
          ncycles= atoi(optarg);
          break;
      case 'i':
          infile = string(optarg);
          break;
      case 'h':
          printf("Usage %s -i <inputfile> -n <data_size> [-h]\n", argv[0]);
          break;
          exit(1);
      case '?':
          printf("Usage %s -i <inputfile> -n <data_size> [-h]\n", argv[0]);
          break;
          exit(1);
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
    double ene;
    double data[ncycles];


    while (count < ncycles && !gzeof(datain)){
        gzgets(datain, str, 250);
        line = string(str);
        istringstream ss(line);
        ss >> tint >> ene;
        data[count] = ene;
        count++;
    }

    gzclose(datain);

    double autocorr = gsl_stats_lag1_autocorrelation(data, 1, count);

    printf("Autocorrelation for data: %f\n", autocorr);

    return 0;

}

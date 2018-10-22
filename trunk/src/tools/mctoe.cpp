#include<iostream>
#include<zlib.h>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<string>
#include<sstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

    string line;
    float Temp;
    double beta;
    int nrot, tint;
    float size_x, size_y, size_z;
    char str[300];
    char stmp[80];
    int count=0;
    float tenergy, tconf_energy,  tfloat;
    long double sum_ene=0.0;
    double k=0.0019858775203792202;
    vector<double> energies;
    vector<double> conf_energies;
    vector<double> delta_energies;
    double average_delta;

    gzFile inpfile = gzopen(argv[1], "r");
    gzgets(inpfile, str, 250);
    gzgets(inpfile, str, 250);
    sscanf(str, "%s %f %s %f %s %f %s %f ", stmp, &size_x, stmp, &size_y, stmp, &size_z, stmp, &Temp);
    printf("Temperature: %.3f\n", Temp);

    beta=1/(k*Temp);

    gzgets(inpfile, str, 200);
    sscanf(str, "%s %s %d", stmp, stmp, &nrot);
    printf("Number of Rotatable Bonds: %d.\n", nrot);

    gzgets(inpfile, str, 200);
    gzgets(inpfile, str, 200);

    while (!gzeof(inpfile)){
        gzgets(inpfile, str, 250);
        count++;

        float* torsion = new float[nrot];

        line = string(str);
        istringstream ss(line);
        ss >> tint >> tenergy >> tfloat >> tfloat >> tfloat >> tfloat >> tfloat >>
                tfloat >> tfloat >> tint >> tconf_energy;
        for (int i=0; i< nrot; i++){
            ss >> torsion[i];
        }

        sum_ene+= double(tenergy-tconf_energy);
        energies.push_back(double(tenergy));
        conf_energies.push_back(double(tconf_energy));
        delta_energies.push_back(double(tenergy-tconf_energy));
    }
        gzclose(inpfile);

        average_delta = sum_ene/count;
        printf("Average Delta_Energy: %10.4f\n", average_delta);

        double deltaU;
        double second_order=0.0, third_order=0.0;
        for (int i=0; i< delta_energies.size(); i++){
            deltaU=delta_energies[i]-average_delta;
            second_order+= (deltaU*deltaU);
            third_order+= (deltaU*deltaU*deltaU);
        }
        second_order = ((beta/2)*(second_order/delta_energies.size()));
        third_order = ((beta*beta/6)*(third_order/delta_energies.size()));

        printf("Second-Order Term: %10.4f kcal/mol\n", second_order);
        printf("Third-Order Term:  %10.4f kcal/mol\n", third_order);

        double DG=average_delta + second_order + third_order;
        printf("DeltaG = %10.3f kcal/mol\n", DG);



    }

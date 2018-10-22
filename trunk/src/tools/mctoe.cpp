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
    double average_RL;
    double average_L;

    gzFile inpfile = gzopen(argv[1], "r");
    gzgets(inpfile, str, 250);
    gzgets(inpfile, str, 250);
    sscanf(str, "%s %f %s %f %s %f %s %f ", stmp, &size_x, stmp, &size_y, stmp, &size_z, stmp, &Temp);
    printf("Temperature: %.3f\n", Temp);

    beta=1.0/(k*Temp);

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

        average_L+= double(tconf_energy);
        average_RL+= double(tenergy);
        average_delta+= double(tenergy-tconf_energy);
        energies.push_back(double(tenergy));
        conf_energies.push_back(double(tconf_energy));
        delta_energies.push_back(double(tenergy-tconf_energy));
    }
        gzclose(inpfile);


        average_L = average_L / count;
        average_RL = average_RL / count;
        average_delta = average_delta / count;
        printf("Average Delta_Energy: %10.4f\n", average_delta);

        double deltaU_L, deltaU_RL;;
        double RL_2=0.0, RL_3=0.0, L_2=0, L_3=0.0;
        for (int i=0; i< delta_energies.size(); i++){
            deltaU_L=conf_energies[i]-average_L;
            deltaU_RL = energies[i]-average_RL;

            L_2+= (deltaU_L*deltaU_L);
            L_3+= (deltaU_L*deltaU_L*deltaU_L);

            RL_2+= (deltaU_RL*deltaU_RL);
            RL_3+= (deltaU_RL*deltaU_RL*deltaU_RL);
        }

        double second_order = ((beta*beta/2)*(RL_2 - L_2));
        double third_order = ((beta*beta*beta/6)*(RL_3-L_3));

        printf("Second-Order Term: %10.4f kcal/mol\n", second_order);
        printf("Third-Order Term:  %10.4f kcal/mol\n", third_order);

        double DG=average_RL - average_L + second_order + third_order;
        printf("DeltaG = %10.3f kcal/mol\n", DG);



    }

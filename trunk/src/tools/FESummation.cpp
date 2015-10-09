#include <iostream>
#include <iostream>
#include <fstream>
#include <iostream>
#include<stdio.h>
#include<string>
#include<vector>
#include<cmath>
#include <cstring>
#include<cstdlib>
#include "zlib.h"

const double K = 0.001987;
const double PI = 3.14159265359;


using namespace std;

void FreeEnergyLogGZ(string file_in_log, string file_in_gz, double T, int TotalStepstoUse){

    vector< long double> enerconf;
    int nconf;
    double volx,voly,volz;


    ifstream cllog;
    string s = "";
    string k = "";
    char fileinlog[1000];
    char buffer[1000];
    char buffer2[1000];
    strcpy(fileinlog, file_in_log.c_str());
    cllog.open(fileinlog);

    while(s !="* generate_conformers            1                                                                 *"){
        getline(cllog,s);
    }

    cllog >> s >> s >> nconf;

    while(s !="Conformer energies (GAFF):"){
        getline(cllog,s);

    }

    long double energycf;

    for(int i=0; i<=nconf+1; i++){
        getline(cllog,s);

        if (s == "****************************************************************************************************"){
            break;
        }

        strcpy(buffer, s.c_str());
        sscanf(buffer,"%s %Lf %s", &buffer2, &energycf, &buffer2);
        enerconf.push_back(energycf);
        printf("Conformer Energy [%2d]: %5.2Lf\n", i, enerconf[i]);
    }

    getline(cllog,s);

    cllog >> volx >> voly >> volz;
    printf("Binding site dimensions: %5.2f %5.2f %5.2f\n", volx, voly, volz);
    cllog.close();

    char fileingz[1000];
    long double tener = 0.0;
    strcpy(fileingz, file_in_gz.c_str());
    gzFile fillegz  = gzopen(fileingz, "r");
    int nlines =0;
    long double norm = 0.0;
    long double Dl, dX, dY, dZ, dA, dB, dG,energy,confenergy,rmsd,Dt;
    long double Angular_Factor, Translational_Factor, FE, Mean_energy;

    for (int i=0; i<enerconf.size(); i++){
        norm = norm + exp(-((enerconf[i])/(K*T)));
    }

    int count = 0;
    while (! gzeof(fillegz)){

        if (nlines==0)
        {
            gzgets(fillegz, buffer,1000);
            gzgets(fillegz, buffer,1000);
            gzgets(fillegz, buffer,1000);
            gzgets(fillegz, buffer,1000);
            nlines++;
        }

        else{

            gzgets(fillegz, buffer,1000);
            sscanf(buffer, "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf ", &Dl, &energy, &rmsd, &dX, &dY, &dZ, &dA, &dB, &dG, &Dt, &confenergy);
            Angular_Factor = (abs(dA)*abs(dB)*abs(dG)*PI*PI*PI)/(8*PI*PI*180*180*180);
            Translational_Factor = (abs(dX)*abs(dY)*abs(dZ))/(volx*voly*volz);
            tener = tener + exp(-((energy)/(K*T)))*(Angular_Factor*Translational_Factor)/(norm);
            count++;

            //Mean internal energy calculation
            Mean_energy = Mean_energy + (energy-confenergy);
        }
    }


    FE = -K*T*log(tener);
    Mean_energy = Mean_energy/count;

    printf("FE = %6.3Lf kcal/mol\n", FE);
    printf("Mean Energy (dH) = %6.3Lf kcal/mol\n", Mean_energy);
    printf("TdS = %6.3Lf kcal/mol\n", Mean_energy-FE);
}




void FreeEnergyLog(string file_in_log, string file_in, double T, int TotalStepstoUse){

    vector< long double> enerconf;
    int nconf;
    double volx,voly,volz;


    ifstream cllog;
    string s = "";
    string k = "";
    char fileinlog[1000];
    char buffer[1000];
    char buffer2[1000];
    strcpy(fileinlog, file_in_log.c_str());
    cllog.open(fileinlog);

    while(s !="* generate_conformers            1                                                                 *"){
        getline(cllog,s);
    }
    cllog >> s >> s >> nconf;

    while(s !="Conformer energies (GAFF):"){
        getline(cllog,s);

    }
    long double energycf;



    for(int i=0; i<= nconf+1; i++){
        getline(cllog,s);

        if (s == "****************************************************************************************************"){
            break;
        }

        strcpy(buffer, s.c_str());
        sscanf(buffer,"%s %Lf %s", &buffer2, &energycf, &buffer2);
        enerconf.push_back(energycf);
        printf("Conformer Energy [%2d]: %5.2Lf\n", i, enerconf[i]);


    }
    getline(cllog,s);

    cllog >> volx >> voly >> volz;
    printf("Binding site dimensions: %5.2f %5.2f %5.2f\n", volx, voly, volz);
    cllog.close();

    char filein[1000];
    long double tener = 0.0;
    strcpy(filein, file_in.c_str());
    FILE * fille  = fopen(filein, "r");
    int nlines =0;
    long double norm = 0.0;
    long double Dl, dX, dY, dZ, dA, dB, dG,energy,confenergy,rmsd,Dt;
    long double Angular_Factor, Translational_Factor, FE, Mean_energy;
    for (int i=0; i<enerconf.size(); i++){
        norm = norm + exp(-((enerconf[i])/(K*T)));
    }
    
    int count=0;
    while (fgetc(fille) != EOF){

        if (nlines==0)
        {
            fgets(buffer,1000,fille);
            fgets(buffer,1000,fille);
            fgets(buffer,1000,fille);
            fgets(buffer,1000,fille);

            nlines++;
        }

        else{

            fgets(buffer,1000,fille);
            sscanf(buffer, "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf ", &Dl, &energy, &rmsd, &dX, &dY, &dZ, &dA, &dB, &dG, &Dt, &confenergy);
            Angular_Factor = (abs(dA)*abs(dB)*abs(dG)*PI*PI*PI)/(8*PI*PI*180*180*180);
            Translational_Factor = (abs(dX)*abs(dY)*abs(dZ))/(volx*voly*volz);
            tener = tener + exp(-((energy)/(K*T)))*(Angular_Factor*Translational_Factor)/(norm);

            //Mean internal energy calculation
            Mean_energy = Mean_energy + (energy-confenergy) ;


        }
    }


    FE = -K*T*log(tener);
    Mean_energy = Mean_energy/count;

    printf("FE = %6.3Lf kcal/mol\n", FE);
    printf("Mean Energy (dH) = %6.3Lf kcal/mol\n", Mean_energy);
    printf("TdS = %6.3Lf kcal/mol\n", Mean_energy-FE);
}





int main(int argc, char *argv[])
{

    string input;
    double T;
    int c;
    int simsteps;
    if (argc != 7){
        printf("Usage: ./%s -i <file.dat.gz> -t <temperature> -s <simulation_steps> \n" , argv[0]);
        exit(1);
    }

    while ((c = getopt(argc, argv, "i:t:s:")) != -1)
        switch (c){
        case 'i':
            input = string(optarg);
            break;
        case 'h':
            printf("Usage %s -i <inputfile> -T <temperature> -s <Simulation Steps>\n", argv[0]);
            break;
            exit(1);
        case '?':
            printf("Usage %s -i <inputfile> -T <temperature> -s <Simulation Steps>\n", argv[0]);
            break;
            exit(1);
        case 't':
            T = double(atof(optarg));

            break;
        case 's':
            simsteps = double(atof(optarg));
            break;
        }

    printf("#*****************************************************************************************\n");
    printf("#                                                                                        *\n");
    printf("#              Free_Energy LiBELa - Free Energy Calculation for iMcLiBELa                *\n");
    printf("#        Written by Alessandro S. Nascimento and Joao Victor S. Cunha   /  2015          *\n");
    printf("#                      University of Sao Paulo - USP - Brazil                            *\n");
    printf("#                                                                                        *\n");
    printf("#                             asnascimento@ifsc.usp.br                                   *\n");
    printf("#                             joao.victor.cunha@usp.br                                   *\n");
    printf("#                                                                                        *\n");
    printf("#*****************************************************************************************\n");


    printf("Input file = %s\n", input.c_str());
    printf("Temperature = %5.1f K\n", T);

    string file_in = input;
    int n = input.length()-10;
    string file_in_log_gz = input.substr(0,n)+".log";
    string file_in_log = input.substr(0,n+3)+".log";
    if (input.substr(input.size()-3, 3) == ".gz")
    {
        FreeEnergyLogGZ(file_in_log_gz,file_in,T,simsteps);
    }
    else{
        FreeEnergyLog(file_in_log,file_in,T,simsteps);
    }
}



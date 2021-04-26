#ifndef IMCLIBELA_H
#define IMCLIBELA_H

#define ORGANIZATION "Universidade de Sao Paulo"
#define NAME "iMcLiBELa"
#define LONG_NAME "iMcLiBELa - Monte Carlo-Based Ligand Binding Energy Landscape"
#define VERSION "1.0"
#define PI 3.14159265359


struct energy_result_t{
    double vdw= 0.0;
    double elec= 0.0;
    double rec_solv= 0.0;
    double lig_solv= 0.0;
    double hb_donor = 0.0;
    double hb_acceptor = 0.0;
    double restraints = 0.0;
    double total = 0.0;

};

#endif

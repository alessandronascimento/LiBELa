#ifndef ENERGY2_H
#define ENERGY2_H

#include<vector>
#include "PARSER.h"
#include "Mol2.h"
#include "iMcLiBELa.h"
#include "Grid.h"

class Energy2
{
public:
    struct GridInterpol{
        double vdwA;
        double vdwB;
        double elec;
        double pbsa;
        double solv_gauss;
        double rec_solv_gauss;
    };

    PARSER* Input;
    Energy2(PARSER* _Input);
    double distance(double x1, double x2, double y1, double y2, double z1, double z2);
    double distance_squared(double x1, double x2, double y1, double y2, double z1, double z2);

    double compute_energy_softcore_solvation(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz);
    double compute_energy_softcore(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz);
    double compute_energy_hardcore_solvation(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz);
    double compute_energy_hardcore(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz);

    double compute_ene_from_grids_softcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
    double compute_ene_from_grids_softcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
    double compute_ene_from_grids_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
    double compute_ene_from_grids_hardcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);


    double compute_energy_softcore_solvation(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);
    double compute_energy_softcore(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);
    double compute_energy_hardcore_solvation(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);
    double compute_energy_hardcore(Mol2* Rec, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);

    double compute_ene(Mol2* Rec,   Mol2* Lig, vector<vector<double> >lig_xyz);
    double compute_ene(Grid *Grids, Mol2* Lig, vector<vector<double> >lig_xyz);
    double compute_ene(Grid* Grids, Mol2* Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);
    double compute_ene(Mol2 *Rec, Mol2 *Lig, vector<vector<double> > lig_xyz, energy_result_t* energy_result);

    double compute_ene_from_grids_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result);
    double compute_ene_from_grids_hardcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result);
    double compute_ene_from_grids_softcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result);
    double compute_ene_from_grids_softcore(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz, energy_result_t* energy_result);

    double trilinear_interpolation(vector<vector<vector<double> > > grid, double x, double y, double z, int x0, int y0, int z0, int x1, int y1, int z1);
    void trilinear_interpolation(Grid* Grids, double x, double y, double z, int x0, int y0, int z0, int x1, int y1, int z1, GridInterpol* GI);

    double evaluate_forces_hardcore_solvation(Grid* Grids, Mol2* Lig, vector<vector<double> > xyz);
};

#endif // ENERGY2_H

#include<iostream>
#include<zlib.h>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<string>
#include <vector>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/rotor.h>
#include <zlib.h>
#include "../LiBELa/iMcLiBELa.h"
#include "../LiBELa/main.h"
#include "../LiBELa/Mol2.h"
#include "../LiBELa/COORD_MC.h"
#include "../LiBELa/PARSER.h"
#include "../LiBELa/Optimizer.h"
#include "../LiBELa/McEntropy.h"


#define k 0.0019872041

using namespace std;
using namespace OpenBabel;

OBMol* GetMol(const string &molfile){
    OBMol* mol = new OBMol;

    OBConversion conv;
    OBFormat *format = conv.FormatFromExt(molfile.c_str());
    if (!format || !conv.SetInFormat(format)) {
        printf("Could not find input format for file\n");
        return mol;
    }

    ifstream ifs(molfile.c_str());
    if (!ifs) {
        printf("Could not open %s for reading.\n", molfile.c_str());
        return mol;
    }

    if (!conv.Read(mol, &ifs)) {
        printf("Could not read molecule from file\n");
        return mol;
    }

    return mol;
}

vector<vector<double> > copy_from_obmol(OBMol mymol){
    double* myxyz = new double[mymol.NumAtoms()*3];
    vector<vector<double > > vec_xyz;
    vector<double> tmp(3);
    myxyz = mymol.GetCoordinates();
    for (unsigned i=0; i < mymol.NumAtoms(); i++){
        tmp[0] = (myxyz[3*i]);
        tmp[1] = (myxyz[(3*i)+1]);
        tmp[2] = (myxyz[(3*i)+2]);
        vec_xyz.push_back(tmp);
    }

    tmp.clear();
    return vec_xyz;
}

void copy_to_obmol(vector<vector<double> > vec_xyz, OBMol* mol){
    double* dxyz = new double[vec_xyz.size()*3];
    for (unsigned i=0; i<vec_xyz.size(); i++){
        dxyz[3*i] = vec_xyz[i][0];
        dxyz[(3*i)+1] = vec_xyz[i][1];
        dxyz[(3*i)+2] = vec_xyz[i][2];
    }
    mol->SetCoordinates(dxyz);
    delete [] dxyz;
}

double* copy_to_obmol(vector<vector<double> > vec_xyz){
    double* dxyz = new double[vec_xyz.size()*3];
    for (unsigned i=0; i<vec_xyz.size(); i++){
        dxyz[3*i] = vec_xyz[i][0];
        dxyz[(3*i)+1] = vec_xyz[i][1];
        dxyz[(3*i)+2] = vec_xyz[i][2];
    }
    return(dxyz);
}

double get_dihedral(vector<double> a, vector<double> b, vector<double> c, vector<double> d){
    // Computing vectors C,B,C
    double xij = a[0]-b[0];
    double yij = a[1]-b[1];
    double zij = a[2]-b[2];
    double xkj = c[0]-b[0];
    double ykj = c[1]-b[1];
    double zkj = c[2]-b[2];
    double xkl = c[0]-d[0];
    double ykl = c[1]-d[1];
    double zkl = c[2]-d[2];

    //     Calculate the normals to the two planes n1 and n2
    //      this is given as the cross products:
    //       AB x BC
    //      --------- = n1
    //      |AB x BC|
    //
    //       BC x CD
    //      --------- = n2
    //      |BC x CD|
    //
    double dxi = (yij * zkj) - (zij * ykj);     // Normal to plane 1
    double dyi = (zij * xkj) - (xij * zkj);
    double dzi = (xij * ykj) - (yij * xkj);
    double gxi = (zkj * ykl) - (ykj * zkl);    // Normal to plane 2
    double gyi = (xkj * zkl) - (zkj * xkl);
    double gzi = (ykj * xkl) - (xkj * ykl);

    //  Calculate the length of the two normals

    double bi = (dxi * dxi) + (dyi * dyi) + (dzi * dzi);
    double bk = (gxi * gxi) + (gyi * gyi) + (gzi * gzi);
    double ct = (dxi * gxi) + (dyi * gyi) + (dzi * gzi);

    double boi2 = 1./bi;
    double boj2 = 1./bk;
    bi   = sqrt(bi);
    bk   = sqrt(bk);

    double z1   = 1./bi;
    double z2   = 1./bk;
    double bioj = bi * z2;
    double bjoi = bk * z1;


    ct   = ct * z1 * z2;


    if (ct >  1.0){
        ct = 1.0;
    }
    else if (ct < (-1.0)){
        ct = -1.0;
    }
    double ap = acos(ct);

    double s = xkj * (dzi * gyi - dyi * gzi) + ykj * (dxi * gzi - dzi * gxi) + zkj * (dyi * gxi - dxi * gyi);

    if (s < 0.0){
        ap = -ap;
    }

    else if (ap > 0.0){
        ap = PI - ap;
    }
    else{
        ap = -(PI + ap);
    }
    ap = ap*180/PI;       // convert to degrees
    return (ap);
}

double get_OB_dihedral(OBMol* mol, int a, int b, int c, int d){
    return mol->GetTorsion(a, b, c, d);
}

double check_angle(double angle){
    if ( angle > 360.0){
        while (angle > 360.0){
            angle -= 360.0;
        }
    }
    else if (angle < 0.0){
        while (angle < 0.0){
            angle += 360.0;
        }
    }
    return angle;
}

int main(int argc, char* argv[]){

    string trajfile;
    string ref_mol2;
    int c, nrot, nframes=1E7, nthreads=1;
    double dx, dy, dz, dalpha, dbeta, dgamma, angle, rmsdi, rmsdf, T=300.0;
    double energy;
    vector<vector<int> > atoms_in_dihedrals;

    if (argc < 2){
        printf("Usage %s -t <traj_file> -n <number_of_frames> -r <ref_mol2> -T <temperature> -p <nthreads> [-h]\n", argv[0]);
        exit(1);
    }

    while ((c = getopt(argc, argv, "t:r:T:n:p:h")) != -1)
        switch (c){
        case 't':
            trajfile = string(optarg);
            break;
        case 'h':
            printf("Usage %s -t <traj_file> -n <number_of_frames> -r <ref_mol2> -T <temperature> -p <nthreads> [-h]\n", argv[0]);
            exit(1);
            break;
        case '?':
            printf("Usage %s -t <traj_file> -n <number_of_frames> -r <ref_mol2> -T <temperature> -p <nthreads> [-h]\n", argv[0]);
            exit(1);
            break;
        case 'r':
            ref_mol2 = string(optarg);
            break;
        case 'T':
            T = double(atof(optarg));
            break;
        case 'n':
            nframes = atoi(optarg);
            break;
        case 'p':
            nthreads = atoi(optarg);
            break;
        }

    printf("#*****************************************************************************************\n");
    printf("#                                                                                        *\n");
    printf("#         McEnergy - A program for Binding Energy computation within McLiBELa            *\n");
    printf("#                    Written by Alessandro S. Nascimento - 2018                          *\n");
    printf("#                      University of Sao Paulo - USP - Brazil                            *\n");
    printf("#                                                                                        *\n");
    printf("#                             asnascimento@ifsc.usp.br                                   *\n");
    printf("#                                                                                        *\n");
    printf("#*****************************************************************************************\n");
    printf("# Trajectory file:  %-74.74s\n", trajfile.c_str());
    printf("# Reference file: %-74.74s\n", ref_mol2.c_str());
    printf("# Temperature: %-74.3f\n", T);
    printf("#*****************************************************************************************\n");


    PARSER* Input = new PARSER;
    Input->temp = T;
    Input->sample_torsions = true;
    COORD_MC* Coord = new COORD_MC;

    Mol2* RefMol = new Mol2(Input, ref_mol2);                   // read the initial coordinates of the ligand

    OBMol* mol = new OBMol;
    mol = GetMol(ref_mol2);

    OBRotorList RotorList;
    OBRotorIterator RotorIterator;
    OBRotor *Rotor;

    RotorList.Setup(*mol);
    Rotor = RotorList.BeginRotor(RotorIterator);
    mol->ToInertialFrame();
    vector<int> tmp(4);

    nrot = int(RotorList.Size());
    printf("Found %d rotatable bonds in ligand %s.\n", nrot, RefMol->molname.c_str());

    for (unsigned i = 0; i < RotorList.Size(); ++i, Rotor = RotorList.NextRotor(RotorIterator)) {
        tmp = Rotor->GetDihedralAtoms();
        for (unsigned j=0; j<4; j++){
            tmp[j] = tmp[j]-1;      // now indexes start with 0 instead of the OB default 1.
        }
        angle = get_dihedral(RefMol->xyz[tmp[0]], RefMol->xyz[tmp[1]], RefMol->xyz[tmp[2]], RefMol->xyz[tmp[3]]);
        printf("Torsion [%2d]: %2d %2d %2d %2d = %8.4f\n", i+1, tmp[0], tmp[1], tmp[2], tmp[3],angle);
        atoms_in_dihedrals.push_back(tmp);
        tmp.clear();
    }

    printf("#*****************************************************************************************\n");

    delete mol;

    vector<double> torsions(nrot);

    unique_ptr<Mol2> TrajMol2(new Mol2(Input, ref_mol2));
    vector<double> com = Coord->compute_com(RefMol);

    unique_ptr<McEntropy> Entropy(new McEntropy(Input, Coord, com, nrot));

    printf("#%10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s", "Frame", "DX", "DY", "DZ", "DALPHA", "DBETA", "DGAMMA", "RMSDi", "RMSDf", "OptStatus");

    for (int j=0; j< nrot; j++){
        printf(" ROT[%3d] ", j);
    }

    printf("\n");

    energy = 0.0;
    Optimizer::align_t align_data;
    align_data.ref_xyz = RefMol->xyz;

    gzFile trajectory = gzopen(trajfile.c_str(), "r");

    if (trajectory == NULL){
        printf("# Could not open file %s...\n", trajfile.c_str());
        exit(1);
    }

    int count=0;
    int opt_status;

    for (int i=0; i< nframes; i++){
        if (! gzeof(trajectory)){
            unique_ptr<Optimizer> opt(new Optimizer(RefMol, RefMol, Input));
            Optimizer::align_result_t  opt_result;
            for (unsigned a=0; a<3; a++){
                opt_result.rotation.push_back(0.0);
                opt_result.translation.push_back(0.0);
            }
            align_data.current_xyz.clear();
            align_data.current_xyz = TrajMol2->get_next_xyz(Input, trajectory);

            if (align_data.ref_xyz.size() == align_data.current_xyz.size()){
                count++;

                opt->minimize_alignment_nlopt_simplex(&align_data, &opt_result);

                dx = opt_result.translation[0];
                dy = opt_result.translation[1];
                dz = opt_result.translation[2];
                dalpha = opt_result.rotation[0];
                dbeta = opt_result.rotation[1];
                dgamma = opt_result.rotation[2];
                rmsdi = Coord->compute_rmsd(align_data.ref_xyz, align_data.current_xyz, RefMol->N);
                rmsdf = opt_result.rmsd;
                energy+= TrajMol2->ensemble_energies[i];


                for (unsigned j=0; j< unsigned(nrot); j++){
                    angle = get_dihedral(align_data.current_xyz[atoms_in_dihedrals[j][0]], align_data.current_xyz[atoms_in_dihedrals[j][1]], align_data.current_xyz[atoms_in_dihedrals[j][2]], align_data.current_xyz[atoms_in_dihedrals[j][3]]);
                    angle = check_angle(angle);
                    torsions[j] = (angle);
                }

                Entropy->update(dx, dy, dz, dalpha, dbeta, dgamma, torsions);

                printf("%10d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10d", i+1, dx, dy, dz, dalpha, dbeta, dgamma, rmsdi, rmsdf, opt_result.opt_status);

                for (unsigned j=0; j< RotorList.Size(); j++){
                    printf("%10.4f ", torsions[j]);
                }
                printf("\n");

            }
            else{
                printf("Sizes of reference molecule (%3lu) and trajectory (%3lu) don't match. Please, check!",
                align_data.ref_xyz.size(), align_data.current_xyz.size());
                exit(1);
            }
        }
    }

    gzclose(trajectory);

    McEntropy::entropy_t McEnt;
    McEntropy::entropy_t Max_Ent;
    McEnt.Srot = 0.0; McEnt.Storsion = 0.0; McEnt.Strans = 0.0;

    Entropy->get_results(&McEnt, &Max_Ent, count);

    energy = energy / count;

    printf("#*****************************************************************************************\n");

    printf("First-Order Approximation Translation Entropy (TS): %10.4g kcal/mol @ %7.2f K\n", McEnt.Strans*T, T);
    printf("First-Order Approximation Rotation Entropy (TS):    %10.4g kcal/mol @ %7.2f K\n", McEnt.Srot*T, T);
    printf("First-Order Approximation Torsion Entropy (TS):     %10.4g kcal/mol @ %7.2f K\n", McEnt.Storsion*T, T);
    printf("First-Order Approximation Total Entropy (S):        %10.4g kcal/(mol.K)@ %7.2f K\n", McEnt.S, T);
    printf("First-Order Approximation -TS (-TS):                %10.4g kcal/mol @ %7.2f K\n", -McEnt.TS, T);
    printf("First-Order Approximation -TS @ 300K:               %10.4g kcal/mol @ %7.2f K\n", -McEnt.S*300., 300.);

    printf("#*****************************************************************************************\n");

    printf("Maximal Entropies Computed for this System:\n");
    printf("First-Order Approximation Translation Entropy (TS): %10.4g kcal/mol @ %7.2f K\n", Max_Ent.Strans*T, T);
    printf("First-Order Approximation Rotation Entropy (TS):    %10.4g kcal/mol @ %7.2f K\n", Max_Ent.Srot*T, T);
    printf("First-Order Approximation Torsion Entropy (TS):     %10.4g kcal/mol @ %7.2f K\n", Max_Ent.Storsion*T, T);
    printf("First-Order Approximation Total Entropy (S):        %10.4g kcal/(mol.K)@ %7.2f K\n", Max_Ent.S, T);
    printf("First-Order Approximation -TS (-TS):                %10.4g kcal/mol @ %7.2f K\n", -Max_Ent.TS, T);
    printf("First-Order Approximation -TS @ 300K:               %10.4g kcal/mol @ %7.2f K\n", -Max_Ent.S*300., 300.);

    printf("#*****************************************************************************************\n");

    printf("Entropy loss (-TdS): %10.4g kcal/mol (%10.4f %s) @ %7.2f K\n", (-McEnt.TS - (-Max_Ent.TS)), ((-McEnt.TS/-Max_Ent.TS)*100), "%", T);

    printf("#*****************************************************************************************\n");

    printf("Ebind = <E> -TS\n");
    printf("Ebind = %8.4f - %8.4f = %8.4f kcal/mol\n", energy, McEnt.TS, (energy-McEnt.TS));

    printf("#*****************************************************************************************\n");

    printf("# Optimization Status:\n");
    printf("# Failure = -1, Out of memory = -3, Roundoff limited = -4, Forced stop = -5,\n");
    printf("# Reached Stopval = 2, Ftol = 3, Xtol = 4, Maxeval = 5, Maxtime =   6\n");

    printf("#*****************************************************************************************\n");


    delete RefMol;
    delete Input;
    delete Coord;

    return 0;
}


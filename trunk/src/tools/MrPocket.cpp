#include <iostream>
#include <vector>
#include <cstdio>
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/PARSER.cpp"
#include "../LiBELa/COORD_MC.cpp"
#include "../LiBELa/WRITER.cpp"
using namespace std;

double distance_squared(double x1, double x2, double y1, double y2, double z1, double z2){
    return (((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1)));
}

void append_atom(Mol2* Pocket, double x, double y, double z, double charge){
    Pocket->N +=1;
    Pocket->radii.push_back(1.9080);
    Pocket->epsilons.push_back(0.0860);
    vector<double> xyz;
    xyz.push_back(x);
    xyz.push_back(y);
    xyz.push_back(z);
    Pocket->xyz.push_back(xyz);
    xyz.clear();
    Pocket->charges.push_back(charge);
    Pocket->residue_pointer.push_back(1);
    Pocket->resnames.push_back("DUM");
    Pocket->sybyl_atoms.push_back("C.3");
    Pocket->atomnames.push_back("C");
}

int main (int argc, char* argv[]){

    double probe_radii = 1.9080;          //gaff c
    double probe_epsilon = 0.0860;        //gaff c
    double box_size = 20.0;
    int res_site = 0;
    int res_number;
    char *ligfile;
    char *inputfile;
    int c;

    if (argc < 2){
        printf("Usage: McPocket -i iMcLiBELa.inp\n");
        printf("\t -l <ligand.mol> to define the pocket using the reference ligand\n");
        printf("\t -r to define the pochet using the receptor center of mass\n");
        printf("\t -n <resnumber> to define residue <resnumber> as an anchor on the pocket\n");
        printf("\t -b <box_size> to define the box size. Default is 15.0 Angstrom\n");
        exit(1);
    }

    while ((c = getopt (argc, argv, "i:rln:b:")) != -1)
        switch (c)
        {
        case 'i':
            inputfile = optarg;
            break;
        case 'r':
            break;
        case 'n':
            res_site = 1;
            res_number = atoi(optarg);
            break;
        case 'l':
            res_site = -1;
            ligfile = optarg;
            break;
        case 'b':
            box_size = atof(optarg);
            break;
        case '?':
            if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt)){
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                printf("Usage: McPocket -i iMcLiBELa.inp\n");
                printf("\t -l to define the pocket using the reference ligand\n");
                printf("\t -r to define the pochet using the receptor center of mass\n");
                printf("\t -n <resnumber> to define residue <resnumber> as an anchor on the pocket\n");
                printf("\t -b <box_size> to define the box size. Default is 15.0 Angstrom\n");
            }
            else{
                fprintf (stderr, "Unknown option character `\\x%x'.\n",optopt);
                printf("Usage: McPocket -i iMcLiBELa.inp\n");
                printf("\t -l to define the pocket using the reference ligand\n");
                printf("\t -r to define the pochet using the receptor center of mass\n");
                printf("\t -n <resnumber> to define residue <resnumber> as an anchor on the pocket\n");
                printf("\t -b <box_size> to define the box size. Default is 15.0 Angstrom\n");
            }
            return 1;
            break;
        default:
            printf("Usage: McPocket -i iMcLiBELa.inp\n");
            printf("\t -l <ligand.mol2> to define the pocket using the reference ligand\n");
            printf("\t -r to define the pochet using the receptor center of mass\n");
            printf("\t -n <resnumber> to define residue <resnumber> as an anchor on the pocket\n");
            printf("\t -b <box_size> to define the box size. Default is 15.0 Angstrom\n");
        }


    PARSER* Input = new PARSER;
    Input->set_parameters(inputfile);
    Input->write_mol2 = true;

    Mol2* Rec = new Mol2(Input, Input->rec_mol2);
    printf("Molecule %s read!\n", Rec->molname.c_str());

    if (res_site > 0) {
        printf("Residue %d selected: %s\n", res_site, Rec->resnames[res_number-1].c_str());
    }

    double dij, dij2, rij, eij, acoef, bcoef, vdw, elec;


    Mol2* Pocket = new Mol2();
    Pocket->N=0;

    COORD_MC* Coord = new COORD_MC;

    WRITER* Writer = new WRITER("pocket", Input);

    vector<double> rec_com;

    if (res_site == 0){
        rec_com = Coord->compute_com(Rec);
        printf("Receptor center of mass: %.3f %.3f %.3f\n", rec_com[0], rec_com[1], rec_com[2]);
    }
    else if (res_site > 0){
        rec_com.push_back(Rec->xyz[Rec->residue_pointer[res_site-1]-1][0]);
        rec_com.push_back(Rec->xyz[Rec->residue_pointer[res_site-1]-1][1]);
        rec_com.push_back(Rec->xyz[Rec->residue_pointer[res_site-1]-1][2]);
        printf("Selection center of mass: %.3f %.3f %.3f\n", rec_com[0], rec_com[1], rec_com[2]);
    }
    else {
        Mol2* Lig = new Mol2(Input, Input->reflig_mol2);
        rec_com = Coord->compute_com(Lig);
        delete Lig;
        printf("Ligand center of mass: %.3f %.3f %.3f\n", rec_com[0], rec_com[1], rec_com[2]);
    }

//    printf("Receptor center of mass: %.3f %.3f %.3f\n", rec_com[0], rec_com[1], rec_com[2]);
    printf("Defining a search box with %.1f x %.1f x %.1f Angstroms for cavity searching....\n", box_size, box_size, box_size);

    double x_min = rec_com[0]-(box_size/2);
    double x_max = rec_com[0]+(box_size/2);
    double y_min = rec_com[1]-(box_size/2);
    double y_max = rec_com[1]+(box_size/2);
    double z_min = rec_com[2]-(box_size/2);
    double z_max = rec_com[2]+(box_size/2);

    Writer->write_box(rec_com, x_min, y_min, z_min, x_max, y_max, z_max);

    for (double x=x_min; x<=x_max; x+=1.5){
        for (double y=y_min; y<=y_max; y+=1.5){
            for (double z=z_min; z<=z_max; z+=1.5){
                vdw=0.0;
                elec=0.0;
                for (int i=0; i<Rec->N; i++){
                    rij = Rec->radii[i] + probe_radii;
                    eij = sqrt(Rec->epsilons[i]*probe_radii);
                    dij2 = distance_squared(Rec->xyz[i][0], x, Rec->xyz[i][1], y, Rec->xyz[i][2], z);
                    dij = sqrt(dij2);
                    bcoef = 2*eij*(rij*rij*rij*rij*rij*rij);
                    acoef = eij*(rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij*rij);
                    vdw+= ((acoef/(dij2*dij2*dij2*dij2*dij2*dij2)) - (bcoef/(dij2*dij2*dij2)));
                    if (Input->dielectric_model == "constant"){
                        (332.0*Rec->charges[i])/(Input->diel*dij);
                    }
                    else if (Input->dielectric_model == "4r") {
                        elec += (332.0*Rec->charges[i])/(4*dij2);
                    }
                    else {                                          // Input->dielectric_model = "r"
                        elec += (332.0*Rec->charges[i])/(dij2);
                    }
                }
                if (vdw < -5.0){
                    //        printf("%.3f %.3f %.3f: %.4f\n", x, y, z, vdw);
                    append_atom(Pocket, x, y, z, -elec);
                }
            }
        }
    }


    Mol2* Pocket2 = new Mol2();
    Pocket2->N=0;

    int neighbors=0;
    for (unsigned i=0; i<Pocket->xyz.size()-1; i++){
        neighbors=0;
        for (unsigned j=i+1; j<Pocket->xyz.size(); j++){
            if (distance_squared(Pocket->xyz[i][0], Pocket->xyz[j][0], Pocket->xyz[i][1], Pocket->xyz[j][1], Pocket->xyz[i][2], Pocket->xyz[j][2]) <= 12.25){
                neighbors++;
            }
        }
        //      printf("Dummy atom %d has %d neighbors.\n", i, neighbors);
        if (neighbors > 3){
            append_atom(Pocket2, Pocket->xyz[i][0],Pocket->xyz[i][1], Pocket->xyz[i][2], Pocket->charges[i]);
        }
    }


    Writer->writeMol2(Pocket, Pocket->xyz, 0.00, 0.00);
    Writer->writeMol2(Pocket2, Pocket2->xyz, 0.00, 0.00, "pocket2");
    delete Writer;
    delete Coord;
    delete Pocket;
    delete Rec;
    delete Input;
    return 0;
}

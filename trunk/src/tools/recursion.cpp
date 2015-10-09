
#include "MCRecursion.h"
#include "../WRITER.h"
#include "../Mol2.h"
#include "../PARSER.h"
#include "../Grid.h"
#include "../Conformer.h"
#include "../COORD_MC.h"
#include "../Energy2.h"


using namespace std;

int main(int argc, char *argv[]){

    char info[98];
    clock_t start, end;
    COORD_MC* Coord = new COORD_MC;
    PARSER* Input = new PARSER;
    Input->set_parameters(argv[1]);

    WRITER* Writer = new WRITER("MCRecursion", Input);
    MC* MCR = new MC(Writer);

    Energy2* Ene = new Energy2(Input);

    Mol2* Lig = new Mol2(Input, Input->lig_mol2);
    Mol2* Rec = new Mol2(Input, Input->rec_mol2);
    Grid* Grids = new Grid(Input);
    vector<double> center;

    if (Input->use_grids){
        center = Coord->compute_com(Lig);
        if (Input->load_grid_from_file){
            sprintf(info,"Loading grids from file %s.grid...", Input->grid_prefix.c_str());
            Writer->print_info(info);
            Grids->load_grids_from_file();
            sprintf(info, "Loaded energy grids with %d points spaced by %5.3f Angstroms in each directon.", Grids->npointsx*Grids->npointsy*Grids->npointsz, Grids->grid_spacing);
            Writer->print_info(info);
        }
        else{
            sprintf(info,"Generating energy grids. It can take a couple of minutes. Coffee time maybe ?");
            Writer->print_info(info);
            start = clock();
            Grids = new Grid(Input, Rec, center);
            end = clock();
            sprintf(info, "Computed energy grids with %d points spaced by %5.3f Angstroms in each directon.", Grids->npointsx*Grids->npointsy*Grids->npointsz, Grids->grid_spacing);
            Writer->print_info(info);
            sprintf(info, "Grid computation took %d seconds.", int((end-start)/CLOCKS_PER_SEC));
            Writer->print_info(info);
        }
        double grid_energy = Ene->compute_ene(Grids, Lig, Lig->xyz);
        sprintf(info,"Original Grid energy: %.4f kcal/mol.", grid_energy);
        Writer->print_info(info);
    }

    if (Input->generate_conformers){
        if (Input->conformer_generator == "GA"){
            Conformer* Conf = new Conformer;
            Conf->generate_conformers_GA(Input, Lig, Input->lig_mol2);
            delete Conf;
        }
        else {
            Conformer* Conf = new Conformer;
            Conf->generate_conformers_WRS(Input, Lig, Input->lig_mol2);
            delete Conf;
        }
    }

    if (Input->use_grids){
        MCR->run(Grids, Lig , Lig, Lig->xyz, Input, 1.189207115, 150.0);
    }
    else {
        //            EqMC->run(REC, RefLig , LIG, LIG->xyz, Input);
    }

    delete Coord;
    delete Ene;
    delete Input;
    delete Lig;
    delete Rec;
    delete Writer;
    delete MCR;


    return 0;

}

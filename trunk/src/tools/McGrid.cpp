#include <iostream>
#include <cstdlib>
#include <vector>
#include "../LiBELa/Grid.cpp"
#include "../LiBELa/PARSER.cpp"
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/COORD_MC.cpp"

using namespace std;

int main(int argc, char* argv[]){

        if (argc < 2){
                printf("Usage: %s input_file\n", argv[0]);
                exit(1);
        }

        printf("****************************************************************************************************\n");
        printf("****************************************************************************************************\n");
        printf("*                  McGrid - Energy Grid Computation for Molecular Docking with McLiBELa            *\n");
        printf("*                                                                                                  *\n");
        printf("* University of Sao Paulo                                                                          *\n");
        printf("* More Info:                                                                                       *\n");
        printf("*      http://www.biotechmol.ifsc.usp.br/                                                          *\n");
        printf("*                                                                                                  *\n");
        printf("****************************************************************************************************\n");
        printf("****************************************************************************************************\n");


        clock_t start, end;

        start = clock();

        PARSER* Input = new PARSER();
        Input->set_parameters(argv[1]);
        Input->write_grids = true;

        Mol2* Rec = new Mol2(Input, Input->rec_mol2);
        Mol2* Lig = new Mol2(Input, Input->reflig_mol2);

        vector<double> com;
        COORD_MC* Coord = new COORD_MC();
        com = Coord->compute_com(Lig);

        delete Lig;
        delete Coord;

        Grid* Grids = new Grid(Input, Rec, com);
        printf("Computed energy grids with %d points spaced by %5.3f Angstroms in each directon.\n", Grids->npointsx*Grids->npointsy*Grids->npointsz, Grids->grid_spacing);
        printf("Grid X limits: %10.4f %10.4f.\n", Grids->xbegin, Grids->xend);
        printf("Grid Y limits: %10.4f %10.4f.\n", Grids->ybegin, Grids->yend);
        printf("Grid Z limits: %10.4f %10.4f.\n", Grids->zbegin, Grids->zend);
        end = clock();
        printf("Grid computation took %d seconds.\n", int((end-start)/CLOCKS_PER_SEC));

        delete Rec;
        delete Input;

        return 0;
}

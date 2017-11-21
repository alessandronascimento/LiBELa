#include <cstdio>
#include <cstdlib>
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/COORD_MC.cpp"
#include "../LiBELa/PARSER.cpp"
#include <vector>
using namespace std;

int main(int argc, char* argv[]){
    if (argc < 2){
        printf("Usage: %s ligand.mol2.\n", argv[0]);
        exit(1);
    }
    vector<double>com;
    PARSER* Input = new PARSER;
    Mol2* Lig = new Mol2(Input, string(argv[1]));
    COORD_MC* Coord = new COORD_MC;
    com = Coord->compute_com(Lig);
    printf("Ligand COM: %10.5f %10.5f %10.5f\n.", com[0], com[1], com[2]);
    return 0;
}



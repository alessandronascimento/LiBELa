#include <cstdlib>
#include <vector>
#include "../Mol2.cpp"
#include "../COORD_MC.cpp"
#include "../PARSER.cpp"
#include "../WRITER.cpp"


using namespace std;

int main (int argc, char* argv[]) {

    printf("****************************************************************************************************\n");
    printf("****************************************************************************************************\n");
    printf("*                       LigShift - Reorients a Ligand for Benchmark Purposes                       *\n");
    printf("*                                                                                                  *\n");
    printf("* University of Sao Paulo                                                                          *\n");
    printf("* More Info:                                                                                       *\n");
    printf("*      http://www.biotechmol.ifsc.usp.br/                                                          *\n");
    printf("*                                                                                                  *\n");
    printf("****************************************************************************************************\n");
    printf("****************************************************************************************************\n");

    if (argc < 7){
        printf("Usage: %s file.mol2 alpha beta gamma dx dy dz\n", argv[0]);
        exit(1);
    }

    PARSER* Input = new PARSER;
    COORD_MC* Coord = new COORD_MC;
    WRITER* Writer = new WRITER("LigShifted", Input);
    Mol2* Mol = new Mol2(Input, string(argv[1]));

    vector<vector<double> > new_xyz;
    double alpha, beta, gamma, dx, dy, dz;

    alpha=atof(argv[1]);
    beta=atof(argv[2]);
    gamma=atof(argv[3]);
    dx=atof(argv[4]);
    dy=atof(argv[5]);
    dz=atof(argv[6]);

    new_xyz = Coord->rototranslate(Mol->xyz, Mol, alpha, beta, gamma, dx, dy, dz);
    double rmsd = Coord->compute_rmsd(Mol->xyz, new_xyz, Mol->N);
    Writer->writeMol2(Mol, new_xyz, 0.0, rmsd);

    printf("* Molecule was rotated and translated. New pose written to LigShift.mol2.gz.                *\n");
    printf("* RMSD after rotation/translation: %7.3f                                                  *\n", rmsd);

    delete Writer;
    delete Coord;
    delete Mol;
    delete Input;

    printf("* Finishing LigShift...                                                                     *\n");
    return 0;
}

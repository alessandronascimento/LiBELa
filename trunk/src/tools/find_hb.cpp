#include <iostream>
#include <cstdio>
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/PARSER.cpp"
#include "../LiBELa/COORD_MC.cpp"
#include "../LiBELa/FindHB.cpp"

using namespace std;


void usage(){
    printf("Usage: Find_HB -l <lig.mol2> \n");
}

int main(int argc, char* argv[]){

    int c;
    char *ligfile;
    double cutoff=10.0;
    PARSER* Input = new PARSER;
    vector<int> hb_donors;
    vector<int> hb_acceptors;

    if (argc < 2){
        usage();
        exit(1);
    }

    while ((c = getopt (argc, argv, "l:")) != -1)
        switch (c)
        {
        case 'l':
            ligfile = optarg;
            Input->reflig_mol2 = string(ligfile);
            Input->lig_mol2 = string(ligfile);
            break;
        case '?':
            if (optopt == 'c')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt)){
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                usage();
            }
            else{
                fprintf (stderr, "Unknown option character `\\x%x'.\n",optopt);
                usage();
            }
            return 1;
            break;
        default:
            usage();
        }


    Mol2* Lig = new Mol2(Input, Input->lig_mol2);
    FindHB* HBLig = new FindHB;
    HBLig->find_ligandHB(Input->lig_mol2, Lig);
    printf("Number of donors / acceptors found: %5lu / %5lu\n", Lig->HBdonors.size(), Lig->HBacceptors.size());

    delete HBLig;
    delete Lig;


    return  0;

}

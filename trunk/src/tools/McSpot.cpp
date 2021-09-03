#include <cstdio>
#include <string>
#include "../LiBELa/Mol2.h"
#include "../LiBELa/PARSER.h"

using namespace std;

int main(int argc, char* argv[]){

    PARSER* Input = new PARSER;
    Mol2* Rec;
    Mol2* Lig;
    string receptor, ligand;
    int c;

    while ((c = getopt (argc, argv, "r:l:")) != -1)
        switch (c)
        {
        case 'r':
            receptor = (string(optarg));
            Rec = new Mol2(Input, receptor);
            break;
        case 'l':
            ligand = (string(optarg));
            Lig = new Mol2(Input, ligand);
            break;
        case '?':
            if (optopt == 'c') {
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                printf("Usage: %s -r <receptor> -l <ligand>\n", argv[0]);
            }
            else if (isprint (optopt)) {
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                printf("Usage: %s -r <receptor> -l <ligand>\n", argv[0]);
            }
            else {
                fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
                printf("Usage: %s -r <receptor> -l <ligand>\n", argv[0]);
            }
            return 1;
        default:
            abort ();
        }
}

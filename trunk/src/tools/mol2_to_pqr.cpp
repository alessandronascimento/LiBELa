#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<string>
#include <vector>
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/WRITER.cpp"

using namespace std;

int main(int argc, char* argv[]){
    if (argc < 3){
        printf("Usage: %s mymol.mol2 mymol.pqr.\n", argv[0]);
        exit(1);
    }
    PARSER* Input = new PARSER;
    Input->output = "mol2_to_pqr";
    Mol2* Lig = new Mol2(Input, string(argv[1]));
    WRITER* Writer = new WRITER(Input);
    Writer->write_pqr(Lig, string(argv[2]));
    return 0;
}

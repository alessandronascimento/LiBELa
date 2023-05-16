#include<zlib.h>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<string>
#include <vector>
#include "../LiBELa/Mol2.cpp"
#include "../LiBELa/PARSER.cpp"

using namespace std;

double distance(double x1, double x2, double y1, double y2, double z1, double z2) {
    double d = ((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1));
    if (d != 0.0) {
        d = sqrt(d);
    }
    return (d);
}

int main(int argc, char* argv[]){

    if (argc < 2){
        printf("Usage: %s file.mol2", argv[0]);
        exit(1);
    }

    PARSER* Input = new PARSER;
    string filename = string(argv[1]);
    Mol2* mol = new Mol2(Input, filename);
    double max_dist=0.0, dist;

    for (int i=0; i<mol->N-1; i++){
        for (int j=i+1; j<mol->N; j++){
            dist = distance(mol->xyz[i][0], mol->xyz[j][0],mol->xyz[i][1], mol->xyz[j][1],mol->xyz[i][2], mol->xyz[j][2]);
            if (dist > max_dist){
                max_dist = dist;
            }
        }
    }
    delete mol;
    delete Input;
    printf("Maximal distance: %8.3f", max_dist);
    return 0;
}

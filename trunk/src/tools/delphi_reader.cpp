#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

void parse_chunck(int i, FILE* file){
    fread(&i, sizeof(int), 1, file);
}

int main(int argc, char* argv[]){

    FILE *phimap;
    char header[132];
    char *title;
    int ivary, nbyte, inddat, intx, inty, intz, tmpi;
    double xang, yang, zang, xstart, xend, ystart, yend, zstart, zend;
    double extent;


    phimap = fopen(argv[1], "rb");

    parse_chunck(tmpi, phimap);

    title = (char *) malloc(sizeof(char) * 62);
    for (int i=0; i<60; i++) {
        title[i] = fgetc(phimap);
    }
    title[60] = '\n';
    title[61] = (char) 0;
    printf("%s\n", title);

    parse_chunck(tmpi, phimap);
    parse_chunck(tmpi, phimap);

    fread(&ivary, sizeof(int), 1, phimap);
    fread(&nbyte, sizeof(int), 1, phimap);
    fread(&inddat, sizeof(int), 1, phimap);

    fread(&extent, sizeof(double), 1, phimap);
    fread(&extent, sizeof(double), 1, phimap);
    fread(&extent, sizeof(double), 1, phimap);

    fread(&xang, sizeof(double), 1, phimap);
    fread(&yang, sizeof(double), 1, phimap);
    fread(&zang, sizeof(double), 1, phimap);

    fread(&xstart, sizeof(double), 1, phimap);
    fread(&xend, sizeof(double), 1, phimap);

    fread(&ystart, sizeof(double), 1, phimap);
    fread(&yend, sizeof(double), 1, phimap);

    fread(&zstart, sizeof(double), 1, phimap);
    fread(&zend, sizeof(double), 1, phimap);

    fread(&intx, sizeof(int), 1, phimap);
    fread(&inty, sizeof(int), 1, phimap);
    fread(&intz, sizeof(int), 1, phimap);

    printf("ivary: %d\n", ivary);
    printf("nbyte: %d\n", nbyte);
    printf("inddat: %d\n", inddat);

    printf("Box angles: %f %f %f\n", xang, yang, zang);
    printf("Box limits: %f -  %f   %f - %f   %f - %f\n", xstart, xend, ystart, yend, zstart, zend);

    printf("%d %d %d\n", intx, inty, intz);



    vector<vector<vector<double> > > delphi_grid;
    double phi;

    vector<double> vz(intz+1);
    vector<vector<double> > vtmp;
    for (int i=0; i< inty+1; i++){
        vtmp.push_back(vz);
    }

    for (int i=0; i< intx+1; i++){
        delphi_grid.push_back(vtmp);
    }

    for (int z=0; z<intz+1; z++){
        for(int y=0; y<inty+1; y++){
            for (int x=0; x< intx+1; x++){
                fread(&phi, sizeof(double), 1, phimap);
                delphi_grid[x][y][z] = phi;
            }
        }
    }

    for (int a=0; a< intx+1; a++){
        for (int b=0; b< inty+1; b++){
            for (int c=0; c<intz+1; c++){
                printf("%10.6f ", delphi_grid[a][b][c]);
            }
        }
        printf("\n");
    }

    return 0;
}

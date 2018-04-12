#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

void parse_chunck(int i, FILE* file){
    fread(&tmpi, sizeof(int), 1, phimap);
}

int main(int argc, char* argv[]){

    int tmpi;
    FILE *phimap;
    char *uplbl, *nxtlbl, *toplbl, *botlbl;
    double scale, oldmid_x, oldmid_y, oldmid_z;

    phimap = fopen(argv[1], "rb");

    int gsize = atoi(argv[2]);
    printf("Expecting a %d size grid\n", gsize);

    fread(&tmpi, sizeof(int), 1, phimap);

    uplbl = (char *) malloc(sizeof(char) * 22);
    for (int i=0; i<20; i++) {
        uplbl[i] = fgetc(phimap);
    }
    uplbl[20] = '\n';
    uplbl[21] = (char) 0;
    printf("%s\n", uplbl);

    fread(&tmpi, sizeof(int), 1, phimap);
    fread(&tmpi, sizeof(int), 1, phimap);

    nxtlbl = (char *) malloc(sizeof(char) * 12);
    for (int i=0; i<10; i++) {
        nxtlbl[i] = fgetc(phimap);
    }
    nxtlbl[10] = '\n';
    nxtlbl[11] = (char) 0;
    printf("%s\n", nxtlbl);


    toplbl = (char *) malloc(sizeof(char) * 62);
    for (int i=0; i<60; i++) {
        toplbl[i] = fgetc(phimap);
    }
    toplbl[60] = '\n';
    toplbl[61] = (char) 0;
    printf("%s\n", toplbl);

    vector<double> vz(gsize);
    vector<vector<double> > vtmp;
    for (int i=0; i<gsize; i++){
        vtmp.push_back(vz);
    }


    vector<vector<vector<double> > > phi;
    for (int i=0; i<gsize; i++){
        phi.push_back(vtmp);
    }

    fread(&tmpi, sizeof(int), 1, phimap);
    fread(&tmpi, sizeof(int), 1, phimap);

    double kt_phi;
    for (int nx=0; nx < gsize; nx++){
        for (int ny=0; ny < gsize; ny++){
            for (int nz=0; nz < gsize; nz++){
                fread(&kt_phi, sizeof(double), 1, phimap);
                phi[nx][ny][nz] = 0.593*kt_phi;
            }
        }
    }

    fread(&tmpi, sizeof(int), 1, phimap);
    fread(&tmpi, sizeof(int), 1, phimap);

    botlbl = (char *) malloc(sizeof(char) * 18);
    for (int i=0; i<16; i++) {
        botlbl[i] = fgetc(phimap);
    }
    botlbl[16] = '\n';
    botlbl[17] = (char) 0;
    printf("%s\n", botlbl);

    fread(&tmpi, sizeof(int), 1, phimap);
    fread(&tmpi, sizeof(int), 1, phimap);

    fread(&scale, sizeof(double), 1, phimap);
    fread(&oldmid_x, sizeof(double), 1, phimap);
    fread(&oldmid_y, sizeof(double), 1, phimap);
    fread(&oldmid_z, sizeof(double), 1, phimap);

    printf("Scale: %10.4f\n", 1./scale);
    printf("oldmid: %10.4f %10.4f %10.4f\n", oldmid_x, oldmid_y, oldmid_z);

    for (int nx=0; nx < gsize; nx++){
        for (int ny=0; ny < gsize; ny++){
            for (int nz=0; nz < gsize; nz++){
                double x = (((nx+1)-((gsize+1)/2))/scale)+oldmid_x;
                double y = (((ny+1)-((gsize+1)/2))/scale)+oldmid_y;
                double z = (((nz+1)-((gsize+1)/2))/scale)+oldmid_z;
                printf("Potential at %10.4f %10.4f %10.4f: %10.4f\n", x, y, z, phi[nx][ny][nz]);
            }
        }
    }

    return 0;
}

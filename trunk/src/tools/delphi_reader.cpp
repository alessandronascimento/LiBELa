#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

int main(int argc, char* argv[]){

    int tmpi, igrid;
    FILE *phimap;
    char *uplbl, *nxtlbl, *toplbl, *botlbl;
    float scale, oldmid_x, oldmid_y, oldmid_z;

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

    fread(&tmpi, sizeof(int), 1, phimap);
    fread(&tmpi, sizeof(int), 1, phimap);

    vector<float> vz(gsize);
    vector<vector<float> > vtmp;
    for (int i=0; i<gsize; i++){
        vtmp.push_back(vz);
    }

    vector<vector<vector<float> > > phi;
    for (int i=0; i<gsize; i++){
        phi.push_back(vtmp);
    }

    float kt_phi;
    for (int nx=0; nx < gsize; nx++){
        for (int ny=0; ny < gsize; ny++){
            for (int nz=0; nz < gsize; nz++){
                fread(&kt_phi, sizeof(float), 1, phimap);
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

    fread(&scale, sizeof(float), 1, phimap);
    fread(&oldmid_x, sizeof(float), 1, phimap);
    fread(&oldmid_y, sizeof(float), 1, phimap);
    fread(&oldmid_z, sizeof(float), 1, phimap);
    fread(&igrid, sizeof(int), 1, phimap);
    fread(&tmpi, sizeof(int), 1, phimap);

    printf("Scale: %10.4f\n", 1./scale);
    printf("oldmid: %10.4f %10.4f %10.4f\n", oldmid_x, oldmid_y, oldmid_z);
    printf("igrid: %d\n", igrid);

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

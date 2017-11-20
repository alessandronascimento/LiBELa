#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;

int main(int argc, char* argv[]){

    double grid_spacing, xbegin, ybegin, zbegin;
    int npointsx, npointsy, npointsz;

    FILE *pbsa_map;
    char str[80];

    pbsa_map=fopen("pbsa_phi.phi","r");

    if (pbsa_map == NULL){
        printf("Could not open PBSA Grid file. Please check");
        exit(1);
    }

    str[0] = '#';
    while(str[0] =='#'){
        fgets(str, 80, pbsa_map);
    }

    fscanf(pbsa_map, "%f %f %f %f", &grid_spacing, &xbegin, &ybegin, &zbegin);
    fscanf(pbsa_map, "%d %d %d", &npointsx, &npointsy, &npointsz);

    double ***pbsa_grid = new double**[npointsx];
    for (int i=0; i<npointsx; i++){
        pbsa_grid[i] = new double*[npointsy];
    }
    for (int i=0; i<npointsx; i++){
        for (int j=0; j<npointsy; j++){
            pbsa_grid[i][j] = new double[npointsz];
        }
    }

    printf("Reading grids for %d x %d x %d map...\n", npointsx, npointsy, npointsz);

    double phi;

    for (int z=0; z<npointsz; z++){
        for (int y=0; y<npointsy; y++){
            for (int x=0; x< npointsx; x++){
               fscanf(pbsa_map, "%f", &phi);
               pbsa_grid[x][y][z] = phi;
            }
        }
    }

    vector<double> vz(npointsz);
    vector<vector<double> > vtmp;
    for (int i=0; i< npointsy; i++){
        vtmp.push_back(vz);
    }
    vector<vector<vector<double> > > pbsa_map2;
    for (int i=0; i< npointsx; i++){
        pbsa_map2.push_back(vtmp);
    }


    for (int x=0; x<npointsx; x++){
        for (int y=0; y<npointsy; y++){
            for (int z=0; z< npointsz; z++){
               pbsa_map2[x][y][z] = pbsa_grid[x][y][z];
            }
        }
    }

    return 0;
}


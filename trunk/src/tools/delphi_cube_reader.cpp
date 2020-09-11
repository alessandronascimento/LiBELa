#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]){
    FILE* phimap;
    int igrid, itmp;
    float scale, centerx, centery, centerz, dtmp;
    char str[200];

    if (argc < 2 ){
        printf("Usage: %s <phimap.cube> \n", argv[0]);
        exit(1);
    }

    phimap = fopen(argv[1], "r");
    fscanf(phimap, "%10f%5d%10f%8f%10f", &scale, &igrid, &centerx, &centery, &centerz);
    printf("Grid Size = %f\n", 1.0/scale);
    printf("igrid = %10d\n", igrid);
    printf("Grid center: %.6f %.6f %.6f\n", centerx, centery, centerz);
    fgets(str, 80, phimap);
    fgets(str, 80, phimap);
//    fscanf(phimap, "%s", str);
//    printf("%s\n", str);
    fscanf(phimap, "%d %f %f %f", &itmp, &dtmp, &dtmp, &dtmp);
//    printf("%d %f\n", itmp, dtmp);
    fscanf(phimap, "%d %f %f %f", &itmp, &dtmp, &dtmp, &dtmp);
//    printf("%d %f\n", itmp, dtmp);
    fscanf(phimap, "%d %f %f %f", &itmp, &dtmp, &dtmp, &dtmp);
//    printf("%d %f\n", itmp, dtmp);
    fscanf(phimap, "%d %f %f %f", &itmp, &dtmp, &dtmp, &dtmp);
//    printf("%d %f\n", itmp, dtmp);
    fscanf(phimap, "%d %f %f %f %f", &itmp, &dtmp, &dtmp, &dtmp, &dtmp);
//    printf("%d %f\n", itmp, dtmp);
//    printf("Grid Size = %6.2f\n", 1.0/scale);
//    printf("Grid center: %f %f %f\n", centerx, centery, centerz);

    float phi;

    for (int i=0; i<((igrid*igrid*igrid)+1); i++){
        fscanf(phimap, "%f", &phi);
        printf("%6.5f\n", phi);

    }

    return 0;
}

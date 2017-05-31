#include <cstring>
#include<cstdlib>
#include <string>
#include <vector>


using namespace std;

struct atom_param{
    string type;
    double radius;
    double epsilon;
};

atom_param get_gaff_atomic_parameters(string gaff_atom, vector<atom_param> vgaff){
    atom_param atomic_parameters;
    for (unsigned i=0; i<vgaff.size(); i++){
        if (vgaff[i].type == gaff_atom){
            atomic_parameters = vgaff[i];
        }
    }
    return atomic_parameters;
}

int main (int argc, char *argv[]){
    FILE *gaff_file;
    char str[80];
    char at[3];
    float r, e;
    char filename[80];

    char* dir_path = getenv("LIBELA");
    if (dir_path== NULL){
        printf("Environment variable LIBELA is not set.\n");
        exit(1);
    }
    else {
        strcpy(filename, dir_path);
        strcat(filename, "/src/tools/gaff2_vdw.dat");
    }

    gaff_file = fopen(filename, "r");
    vector<atom_param> vgaff;
    if (gaff_file!= NULL){
        while (!feof(gaff_file)){
            fgets(str, 80, gaff_file);
            if (str[0] != '#'){
                sscanf(str, "%s %f %f", at, &r, &e);
                atom_param v;
                v.type = string(at);
                v.radius = double(r);
                v.epsilon = double(e);

                vgaff.push_back(v);
            }
        }
    }

    fclose(gaff_file);

    for (unsigned i=0; i< vgaff.size(); i++){
        printf("%10s %10.4f %10.4f\n", vgaff[i].type.c_str(), vgaff[i].radius, vgaff[i].epsilon);
    }

    atom_param test = get_gaff_atomic_parameters("hc", vgaff);

    printf("Parameters for %s: %10.4f %10.4f\n", test.type.c_str(), test.radius, test.epsilon);

    return 0;

}

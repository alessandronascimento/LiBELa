/*
 * Mol2.cpp
 *
 *  Created on: 10/10/2011
 *      Author: Nascimento
 */

#include "Mol2.h"

using namespace std;

Mol2::Mol2(){

}

Mol2::Mol2(PARSER *Input, ifstream &mol2file) {

	if (mol2file.is_open()){
		line="";
		while(line != "@<TRIPOS>MOLECULE"){
			getline(mol2file, line);
		}

		getline(mol2file, line);		// Molecule name

		mol2file >> this->N;
		mol2file >> this->Nbonds;
		mol2file >> this->Nres;
		str[0]  = '#';
		while (line != "@<TRIPOS>ATOM"){
			getline(mol2file, line);
		}

		double tx, ty, tz;
		vector<double> txyz;
		int tres;
		double tcharge;
		int count=0;

		if (Input->mol2_aa){
			this->get_gaff_parameters();
		}

		for (int i=0; i<N; i++){
			mol2file >> line; // atom number (i+1);

			mol2file >> line;
			this->atomnames.push_back(line);

			mol2file >> tx >> ty >> tz;
			txyz.push_back(tx);
			txyz.push_back(ty);
			txyz.push_back(tz);
			this->xyz.push_back(txyz);
			txyz.clear();

			mol2file >> line; //atom type (SYBYL)
			if (Input->mol2_aa){
				this->amberatoms.push_back(line);
/*
 * If atomtypes are given as gaff types, we *still* need to get VDW parameter
 * in the vdw.param file.
 * For SYBYL atom types there is no need to read the file, since the parameters
 * are parsed in the convert2gaff method.
 */

				this->get_epsilon(this->atomtypes_prm,this->amberatoms[i],this->welldepth);
				this->get_radius(this->atomtypes_prm, this->amberatoms[i], this->radius);
			}
			else {
                this->amberatoms.push_back(this->convert2gaff2(line));
			}

			this->get_masses(amberatoms[i]);

			mol2file >> tres; // residue number;
			if (tres > count){
				this->residue_pointer.push_back(i+1);
				count = tres;
				mol2file >> line;
				this->resnames.push_back(line);
			}
			else {
				mol2file >> line; // residue name
			}

			mol2file >> tcharge;
			this->charges.push_back(tcharge);
		}

		while (line != "@<TRIPOS>BOND"){
			getline(mol2file, line);
		}
/*		vector<int> bond;
		for (int i=0; i<this->Nbonds; i++){
			mol2file >> line >> tx >> ty >> tz;
			bond.push_back(i+1);
			bond.push_back(tx);
			bond.push_back(ty);
			bond.push_back(tz);
			this->bonds.push_back(bond);
			bond.clear();
		}
*/
	}
	else {
		printf("Mol2 file could not be opened! Please check!\n");
		exit(1);
	}
}



Mol2::Mol2(PARSER *Input, string molfile) {
    if ((molfile.substr(molfile.size()-3, 3) == ".gz") or (molfile.substr(molfile.size()-2, 2) == ".z")){
        this->parse_gzipped_file(Input, molfile);
    }
    else {
        FILE *mol2file;
        mol2file = fopen(molfile.c_str(), "r");

        int tint;
        float tx, ty, tz;
        vector<double> txyz;
        int tres;
        float tcharge;
        int count=0;
        char tatomtype[10];
        char resname[10];
        string cpstr;

        if (Input->mol2_aa){
            this->get_gaff_parameters();
        }

        if (mol2file !=NULL){
            str[0]='#';
            while(str[0] !='@'){
                fgets(str, 80, mol2file);

            }
            fgets(str, 80, mol2file);
            this->molname = str;
            this->molname = this->molname.substr(0,this->molname.size()-2);
            fscanf(mol2file, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);

#ifdef DEBUG
            printf("Number of atoms: %d\n", this->N);
#endif

            cpstr = string(str);
            while (cpstr.substr(0,13) != "@<TRIPOS>ATOM"){
                fgets(str, 80, mol2file);
                cpstr = string(str);
            }

            for (int i=0; i<this->N; i++){
                fscanf(mol2file, "%d %s %f %f %f %5s%d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
                txyz.push_back(tx);
                txyz.push_back(ty);
                txyz.push_back(tz);
                this->xyz.push_back(txyz);
                txyz.clear();

                this->charges.push_back(tcharge);
                this->atomnames.push_back(str);

                if (Input->mol2_aa){
                    this->amberatoms.push_back(tatomtype);
                    this->get_epsilon(this->atomtypes_prm,this->amberatoms[i],this->welldepth);
                    this->get_radius(this->atomtypes_prm, this->amberatoms[i], this->radius);
                }
                else{
                    this->amberatoms.push_back(this->convert2gaff2(string(tatomtype)));
                    this->sybyl_atoms.push_back(string(tatomtype));
                }

                this->get_masses(amberatoms[i]);

                if (tres > count){
                    this->residue_pointer.push_back(i+1);
                    count = tres;
                    this->resnames.push_back(string(resname));
                }
            }

            fscanf(mol2file, "%s\n", str);
            if (str[0] != '@'){
                while (str[0] != '@'){
                    fgets(str, 80, mol2file);
                }
            }

            vector<string> bond;
            char s1[6], s2[6], s3[5];
            for (int i=0; i<this->Nbonds; i++){
                fscanf(mol2file, "%d%s%s%s\n", &tint, s1, s2, s3);
                bond.push_back(string(s1));
                bond.push_back(string(s2));
                bond.push_back(string(s3));
                this->bonds.push_back(bond);
                bond.clear();
            }
        }
        else {
            printf("Mol2 file could not be opened! Please check!\n");
            exit(1);
        }
#ifdef DEBUG
        printf("atomnames size: %d\ncharges size:%d\nradii size: %d\n", this->atomnames.size(), this->charges.size(), this->radii.size());
#endif
        fclose(mol2file);
    }
}

bool Mol2::parse_mol2file(PARSER *Input, string molfile) {
	FILE *mol2file;
	int tint;
	float tx, ty, tz;
	vector<double> txyz;
	int tres;
	float tcharge;
	int count=0;
	char tatomtype[10];
	char resname[10];
	string cpstr;
	bool bret = false;

	if (Input->mol2_aa){
		this->get_gaff_parameters();
	}

	mol2file = fopen(molfile.c_str(), "r");

	if (mol2file !=NULL){
		str[0]='#';
		while(str[0] !='@'){
			fgets(str, 80, mol2file);

		}
		fgets(str, 80, mol2file);
		this->molname = str;
        this->molname = this->molname.substr(0,this->molname.size()-1);
		fscanf(mol2file, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);

		cpstr = string(str);
		while (cpstr.substr(0,13) != "@<TRIPOS>ATOM"){
			fgets(str, 80, mol2file);
			cpstr = string(str);
		}

		for (int i=0; i<this->N; i++){
            fscanf(mol2file, "%d %s %f %f %f %5s%d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
			txyz.push_back(tx);
			txyz.push_back(ty);
			txyz.push_back(tz);
			this->xyz.push_back(txyz);
			txyz.clear();

			this->charges.push_back(tcharge);
			this->atomnames.push_back(str);

			if (Input->mol2_aa){
				this->amberatoms.push_back(tatomtype);
				this->get_epsilon(this->atomtypes_prm,this->amberatoms[i],this->welldepth);
				this->get_radius(this->atomtypes_prm, this->amberatoms[i], this->radius);
			}
			else{
                this->amberatoms.push_back(this->convert2gaff2(string(tatomtype)));
				this->sybyl_atoms.push_back(string(tatomtype));
			}

			this->get_masses(amberatoms[i]);

			if (tres > count){
				this->residue_pointer.push_back(i+1);
				count = tres;
				this->resnames.push_back(string(resname));
			}
		}

		fscanf(mol2file, "%s\n", str);
		if (str[0] != '@'){
			while (str[0] != '@'){
				fgets(str, 80, mol2file);
			}
		}

		vector<string> bond;
		char s1[6], s2[6], s3[5];
		for (int i=0; i<this->Nbonds; i++){
			fscanf(mol2file, "%d%s%s%s\n", &tint, s1, s2, s3);
			bond.push_back(string(s1));
			bond.push_back(string(s2));
			bond.push_back(string(s3));
			this->bonds.push_back(bond);
			bond.clear();
		}
		bret = true;
	}
	else {
		bret = false;
		printf("Skipping file %s\n", molfile.c_str());
	}
	fclose(mol2file);
	return (bret);
}

bool Mol2::parse_gzipped_file(PARSER* Input, string molfile){
    bool bret = false;
    int tint;
    float tx, ty, tz;
    vector<double> txyz;
    int tres;
    float tcharge;
    int count=0;
    char tatomtype[10];
    char resname[10];
    string cpstr;

    if (Input->mol2_aa){
//        this->get_gaff_parameters();
        this->initialize_gaff();
    }

    gzFile mol2file = gzopen(molfile.c_str(), "r");
    if (mol2file != NULL){
        str[0]='#';
        while(str[0] !='@'){
            gzgets(mol2file, str, 80);

        }
        gzgets(mol2file, str, 80);
        this->molname = str;
        this->molname = this->molname.substr(0,this->molname.size()-1);
        gzgets(mol2file, str, 80);
        sscanf(str, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);

#ifdef DEBUG
        printf("Number of atoms: %d\n", this->N);
#endif

        cpstr = string(str);
        while (cpstr.substr(0,13) != "@<TRIPOS>ATOM"){
            gzgets(mol2file, str, 80);
            cpstr = string(str);
        }

        for (int i=0; i<this->N; i++){
            gzgets(mol2file, str, 80);
            sscanf(str, "%d %s %f %f %f %5s%d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
            txyz.push_back(tx);
            txyz.push_back(ty);
            txyz.push_back(tz);
            this->xyz.push_back(txyz);
            txyz.clear();

            this->charges.push_back(tcharge);
            this->atomnames.push_back(str);

            if (Input->mol2_aa){
                this->amberatoms.push_back(tatomtype);
                atom_param* at = new atom_param;
                this->get_gaff_atomic_parameters(string(tatomtype), at);
                this->radii.push_back(at->radius);
                this->epsilons.push_back(at->epsilon);
                this->epsilons_sqrt.push_back(sqrt(at->epsilon));
                delete at;
            }
            else{
                this->amberatoms.push_back(this->convert2gaff2(string(tatomtype)));
                this->sybyl_atoms.push_back(string(tatomtype));
            }

            this->get_masses(amberatoms[i]);

            if (tres > count){
                this->residue_pointer.push_back(i+1);
                count = tres;
                this->resnames.push_back(string(resname));
            }
        }

        gzgets(mol2file, str, 80);
        if (str[0] != '@'){
            while (str[0] != '@'){
                gzgets(mol2file, str, 80);
            }
        }

        vector<string> bond;
        char s1[6], s2[6], s3[5];
        for (int i=0; i<this->Nbonds; i++){
            gzgets(mol2file, str, 80);
            sscanf(str, "%d%s%s%s\n", &tint, s1, s2, s3);
            bond.push_back(string(s1));
            bond.push_back(string(s2));
            bond.push_back(string(s3));
            this->bonds.push_back(bond);
            bond.clear();
        }
        bret = true;
    }
    else {
        printf("Skipping file %s...\n", molfile.c_str());
        bret = false;
    }
    gzclose(mol2file);
    return (bret);
}

string Mol2::convert2gaff(string atom){
	string gaff_atom;
	if (atom == "C.3"){
		gaff_atom = "c3";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.1094);
		this->epsilons_sqrt.push_back(sqrt(0.1094));
	}
	else if (atom =="C.2"){
		gaff_atom = "c2";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.0860);
		this->epsilons_sqrt.push_back(sqrt(0.0860));
	}

	else if (atom =="C.1"){
		gaff_atom = "c1";
		this->radii.push_back(1.9080);
        this->epsilons.push_back(0.2100);
        this->epsilons_sqrt.push_back(sqrt(0.2100));

	}

	else if (atom =="C.ar"){
		gaff_atom = "ca";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.0860);
		this->epsilons_sqrt.push_back(sqrt(0.0860));
	}

	else if (atom =="C.cat"){
		gaff_atom = "c";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.0860);
		this->epsilons_sqrt.push_back(sqrt(0.0860));
	}

	else if (atom =="N.3"){
		gaff_atom = "n3";
        this->radii.push_back(1.7685);
        this->epsilons.push_back(0.1405);
        this->epsilons_sqrt.push_back(sqrt(0.1405));
	}

	else if (atom =="N.2"){
		gaff_atom = "n2";
        this->radii.push_back(1.7685);
        this->epsilons.push_back(0.1405);
        this->epsilons_sqrt.push_back(sqrt(0.1405));
	}

	else if (atom =="N.1"){
		gaff_atom = "n1";
        this->radii.push_back(1.7685);
        this->epsilons.push_back(0.1405);
        this->epsilons_sqrt.push_back(sqrt(0.1405));
	}

	else if (atom =="N.ar"){
		gaff_atom = "nh";
        this->radii.push_back(1.7685);
        this->epsilons.push_back(0.1405);
        this->epsilons_sqrt.push_back(sqrt(0.1405));
	}

	else if (atom =="N.am"){
		gaff_atom = "n";
        this->radii.push_back(1.7685);
        this->epsilons.push_back(0.1405);
        this->epsilons_sqrt.push_back(sqrt(0.1405));
	}

	else if (atom =="N.4"){
		gaff_atom = "n4";
        this->radii.push_back(1.7685);
        this->epsilons.push_back(0.1405);
        this->epsilons_sqrt.push_back(sqrt(0.1405));
	}

	else if (atom =="N.pl3"){
		gaff_atom = "na";
        this->radii.push_back(1.7685);
        this->epsilons.push_back(0.1405);
        this->epsilons_sqrt.push_back(sqrt(0.1405));
	}

	else if (atom =="N.p"){
		gaff_atom = "na";
        this->radii.push_back(1.7685);
        this->epsilons.push_back(0.1405);
        this->epsilons_sqrt.push_back(sqrt(0.1405));
	}

	else if (atom =="O.3"){
		gaff_atom = "oh";
		this->radii.push_back(1.7210);
		this->epsilons.push_back(0.2104);
		this->epsilons_sqrt.push_back(sqrt(0.2104));
	}

	else if (atom =="O.2"){
		gaff_atom = "o";
        this->radii.push_back(1.8518);
        this->epsilons.push_back(0.3108);
        this->epsilons_sqrt.push_back(sqrt(0.3108));
	}

	else if (atom =="O.co2"){
		gaff_atom = "o";
        this->radii.push_back(1.8518);
        this->epsilons.push_back(0.3108);
        this->epsilons_sqrt.push_back(sqrt(0.3108));
	}

	else if (atom =="O.spc" or atom == "O.t3p"){
		gaff_atom = "ow";
		this->radii.push_back(1.7683);
		this->epsilons.push_back(0.1520);
		this->epsilons_sqrt.push_back(sqrt(0.1520));
	}

	else if (atom =="S.3"){
		gaff_atom = "sh"; //not sure... sh or ss
		this->radii.push_back(2.0000);
		this->epsilons.push_back(0.2500);
		this->epsilons_sqrt.push_back(sqrt(0.2500));
	}

	else if (atom =="S.2"){
		gaff_atom = "s2";
		this->radii.push_back(2.0000);
		this->epsilons.push_back(0.2500);
		this->epsilons_sqrt.push_back(sqrt(0.2500));
	}

	else if (atom =="S.O" or atom == "S.o"){
		gaff_atom = "s4";
		this->radii.push_back(2.0000);
		this->epsilons.push_back(0.2500);
		this->epsilons_sqrt.push_back(sqrt(0.2500));
	}

	else if (atom =="S.O2" or atom == "S.o2"){
		gaff_atom = "s6";
		this->radii.push_back(2.0000);
		this->epsilons.push_back(0.2500);
		this->epsilons_sqrt.push_back(sqrt(0.2500));
	}

	else if (atom =="P.3"){
		gaff_atom = "p3";
		this->radii.push_back(2.1000);
		this->epsilons.push_back(0.2000);
		this->epsilons_sqrt.push_back(sqrt(0.2000));
	}

	else if (atom =="F"){
		gaff_atom = "f";
		this->radii.push_back(1.75);
		this->epsilons.push_back(0.0610);
		this->epsilons_sqrt.push_back(sqrt(0.0610));
	}

	else if (atom =="H"){
		gaff_atom = "hc";
        this->radii.push_back(0.600);      //        this->radii.push_back(1.4870);
		this->epsilons.push_back(0.0157);
		this->epsilons_sqrt.push_back(sqrt(0.0157));
	}

	else if (atom =="H.spc" or atom=="H.t3p"){
		gaff_atom = "hw";
		this->radii.push_back(0.0000);
        this->epsilons.push_back(0.0000);
		this->epsilons_sqrt.push_back(0.0000);
	}

	else if (atom =="Cl"){
		gaff_atom = "cl";
        this->radii.push_back(1.9452);
        this->epsilons.push_back(0.2638);
        this->epsilons_sqrt.push_back(sqrt(0.2638));
	}

	else if (atom =="Br"){
		gaff_atom = "br";
        this->radii.push_back(2.0275);
        this->epsilons.push_back(0.3932);
        this->epsilons_sqrt.push_back(sqrt(0.3932));
	}

	else if (atom =="I"){
		gaff_atom = "i";
        this->radii.push_back(2.1558);
        this->epsilons.push_back(0.4955);
        this->epsilons_sqrt.push_back(sqrt(0.4955));
	}

	else if (atom =="Mg"){
		gaff_atom = "MG";
		this->radii.push_back(1.5545);
		this->epsilons.push_back(0.00295);
		this->epsilons_sqrt.push_back(sqrt(0.00295));
	}

    else if (atom =="LP" or atom == "Lp"){
		gaff_atom = "DU";
		this->radii.push_back(0.00);
		this->epsilons.push_back(0.00);
		this->epsilons_sqrt.push_back(0.00);
	}

	else if (atom == "Fe"){
		gaff_atom = "FE";
		this->radii.push_back(1.200);
		this->epsilons.push_back(0.05000);
		this->epsilons_sqrt.push_back(0.05000);
	}
	else if (atom == "Zn"){
		gaff_atom = "Zn";
		this->radii.push_back(1.10);
		this->epsilons.push_back(0.0125);
		this->epsilons_sqrt.push_back(sqrt(0.0125));
	}
    else if (atom == "Cu"){
        gaff_atom = "Cu";
        this->radii.push_back(2.20);
        this->epsilons.push_back(0.200);
        this->epsilons_sqrt.push_back(sqrt(0.200));
    }
    else if (atom == "Ca"){
        gaff_atom = "Ca";
        this->radii.push_back(1.790);
        this->epsilons.push_back(0.0140);
        this->epsilons_sqrt.push_back(sqrt(0.0140));
    }

	else{
		printf("Atom type %s not found among GAFF parameters.\nPlease check Mol2.h source file.\n", atom.c_str());
		exit(1);
	}

	return(gaff_atom);
}


string Mol2::convert2gaff2(string atom){
    string gaff_atom;
    if (atom == "C.3"){
        gaff_atom = "c3";
        this->radii.push_back(1.9069);
        this->epsilons.push_back(0.1078);
        this->epsilons_sqrt.push_back(sqrt(0.1078));
    }
    else if (atom =="C.2"){
        gaff_atom = "c2";
        this->radii.push_back(1.8606);
        this->epsilons.push_back(0.0988);
        this->epsilons_sqrt.push_back(sqrt(0.0988));
    }

    else if (atom =="C.1"){
        gaff_atom = "c1";
        this->radii.push_back(1.9525);
        this->epsilons.push_back(0.1596);
        this->epsilons_sqrt.push_back(sqrt(0.1596));

    }

    else if (atom =="C.ar"){
        gaff_atom = "ca";
        this->radii.push_back(1.8606);
        this->epsilons.push_back(0.0988);
        this->epsilons_sqrt.push_back(sqrt(0.0988));
    }

    else if (atom =="C.cat"){
        gaff_atom = "c";
        this->radii.push_back(1.8606);
        this->epsilons.push_back(0.0988);
        this->epsilons_sqrt.push_back(sqrt(0.0988));
    }

    else if (atom =="N.3"){
        gaff_atom = "n3";
        this->radii.push_back(1.8886);
        this->epsilons.push_back(0.0858);
        this->epsilons_sqrt.push_back(sqrt(0.0858));
    }

    else if (atom =="N.2"){
        gaff_atom = "n2";
        this->radii.push_back(1.8993);
        this->epsilons.push_back(0.0941);
        this->epsilons_sqrt.push_back(sqrt(0.0941));
    }

    else if (atom =="N.1"){
        gaff_atom = "n1";
        this->radii.push_back(1.8372);
        this->epsilons.push_back(0.1098);
        this->epsilons_sqrt.push_back(sqrt(0.1098));
    }

    else if (atom =="N.ar"){
        gaff_atom = "nh";
        this->radii.push_back(1.7903);
        this->epsilons.push_back(0.2150);
        this->epsilons_sqrt.push_back(sqrt(0.2150));
    }

    else if (atom =="N.am"){
        gaff_atom = "n";
        this->radii.push_back(1.7852);
        this->epsilons.push_back(0.1636);
        this->epsilons_sqrt.push_back(sqrt(0.1636));
    }

    else if (atom =="N.4"){
        gaff_atom = "n4";
        this->radii.push_back(1.4028);
        this->epsilons.push_back(3.8748);
        this->epsilons_sqrt.push_back(sqrt(3.8748));
    }

    else if (atom =="N.pl3"){
        gaff_atom = "na";
        this->radii.push_back(1.7992);
        this->epsilons.push_back(0.2042);
        this->epsilons_sqrt.push_back(sqrt(0.2042));
    }

    else if (atom =="N.p"){
        gaff_atom = "na";
        this->radii.push_back(7992);
        this->epsilons.push_back(0.2042);
        this->epsilons_sqrt.push_back(sqrt(0.2042));
    }

    else if (atom =="O.3"){
        gaff_atom = "oh";
        this->radii.push_back(1.8200);
        this->epsilons.push_back(0.0930);
        this->epsilons_sqrt.push_back(sqrt(0.0930));
    }

    else if (atom =="O.2"){
        gaff_atom = "o";
        this->radii.push_back(1.7107);
        this->epsilons.push_back(0.1463);
        this->epsilons_sqrt.push_back(sqrt(0.1463));
    }

    else if (atom =="O.co2"){
        gaff_atom = "o";
        this->radii.push_back(1.7107);
        this->epsilons.push_back(0.1463);
        this->epsilons_sqrt.push_back(sqrt(0.1463));
    }

    else if (atom =="O.spc" or atom == "O.t3p"){ //GAFF 1
        gaff_atom = "ow";
        this->radii.push_back(1.7683);
        this->epsilons.push_back(0.1520);
        this->epsilons_sqrt.push_back(sqrt(0.1520));
    }

    else if (atom =="S.3"){
        gaff_atom = "sh"; //not sure... sh or ss
        this->radii.push_back(1.9825);
        this->epsilons.push_back(0.2824);
        this->epsilons_sqrt.push_back(sqrt(0.2824));
    }

    else if (atom =="S.2"){
        gaff_atom = "s2";
        this->radii.push_back(1.9825);
        this->epsilons.push_back(0.2824);
        this->epsilons_sqrt.push_back(sqrt(0.2824));
    }

    else if (atom =="S.O" or atom == "S.o"){
        gaff_atom = "s4";
        this->radii.push_back(1.9825);
        this->epsilons.push_back(0.2824);
        this->epsilons_sqrt.push_back(sqrt(0.2824));
    }

    else if (atom =="S.O2" or atom == "S.o2"){
        gaff_atom = "s6";
        this->radii.push_back(1.9825);
        this->epsilons.push_back(0.2824);
        this->epsilons_sqrt.push_back(sqrt(0.2824));
    }

    else if (atom =="P.3"){
        gaff_atom = "p3";
        this->radii.push_back(2.0732);
        this->epsilons.push_back(0.2295);
        this->epsilons_sqrt.push_back(sqrt(0.2295));
    }

    else if (atom =="F"){
        gaff_atom = "f";
        this->radii.push_back(1.7029);
        this->epsilons.push_back(0.0832);
        this->epsilons_sqrt.push_back(sqrt(0.0832));
    }

    else if (atom =="H"){
        gaff_atom = "hc";
        this->radii.push_back(1.4593);      //        this->radii.push_back(1.4870);
        this->epsilons.push_back(0.0208);
        this->epsilons_sqrt.push_back(sqrt(0.0208));
    }

    else if (atom =="H.spc" or atom=="H.t3p"){
        gaff_atom = "hw";
        this->radii.push_back(0.0000);
        this->epsilons.push_back(0.0000);
        this->epsilons_sqrt.push_back(0.0000);
    }

    else if (atom =="Cl"){
        gaff_atom = "cl";
        this->radii.push_back(1.9452);
        this->epsilons.push_back(0.2638);
        this->epsilons_sqrt.push_back(sqrt(0.2638));
    }

    else if (atom =="Br"){
        gaff_atom = "br";
        this->radii.push_back(2.0275);
        this->epsilons.push_back(0.3932);
        this->epsilons_sqrt.push_back(sqrt(0.3932));
    }

    else if (atom =="I"){
        gaff_atom = "i";
        this->radii.push_back(2.1558);
        this->epsilons.push_back(0.4955);
        this->epsilons_sqrt.push_back(sqrt(0.4955));
    }

/*
 * From now on, parameters are not taken from GAFF.
 */

    else if (atom =="Mg"){
        gaff_atom = "MG";
        this->radii.push_back(1.5545);
        this->epsilons.push_back(0.00295);
        this->epsilons_sqrt.push_back(sqrt(0.00295));
    }

    else if (atom =="LP" or atom == "Lp"){
        gaff_atom = "DU";
        this->radii.push_back(0.00);
        this->epsilons.push_back(0.00);
        this->epsilons_sqrt.push_back(0.00);
    }

    else if (atom == "Fe"){
        gaff_atom = "FE";
        this->radii.push_back(1.200);
        this->epsilons.push_back(0.05000);
        this->epsilons_sqrt.push_back(0.05000);
    }
    else if (atom == "Zn"){
        gaff_atom = "Zn";
        this->radii.push_back(1.10);
        this->epsilons.push_back(0.0125);
        this->epsilons_sqrt.push_back(sqrt(0.0125));
    }
    else if (atom == "Cu"){
        gaff_atom = "Cu";
        this->radii.push_back(2.20);
        this->epsilons.push_back(0.200);
        this->epsilons_sqrt.push_back(sqrt(0.200));
    }
    else if (atom == "Ca"){
        gaff_atom = "Ca";
        this->radii.push_back(1.790);
        this->epsilons.push_back(0.0140);
        this->epsilons_sqrt.push_back(sqrt(0.0140));
    }

    else{
        printf("Atom type %s not found among GAFF parameters.\nPlease check Mol2.h source file.\n", atom.c_str());
        exit(1);
    }

    return(gaff_atom);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Mol2::read_atomtypes_prm() {
	string atom;
	double rad, well;
	char vdw[100];

	elsa_dir_path = getenv("LIBELA_DIR");
	if (elsa_dir_path== NULL){
		printf("Environment variable LIBELA_DIR is not set\n");
		printf("Trying to open vdw.param in local folder...\n");
		strcpy(vdw, "vdw.param");
	}

	else {

		strcpy(vdw, elsa_dir_path);
		strcat(vdw, "/vdw.param");
	}

	ifstream vdwprm(vdw);
		if (vdwprm.is_open()){
			getline(vdwprm, line);
			while (!vdwprm.eof()) {
				vdwprm >> atom >> rad >> well;
				this->atomtypes_prm.push_back(atom);
				this->radius.push_back(rad);
				this->welldepth.push_back(well);
			}
		}
		else {
			printf("Could not open vdw.param file %s. Please check.\n", vdw);
			exit(1);
		}
}

void Mol2::get_epsilon(vector<string>atomtypes_prm, string amberatom, vector<double>welldepth) {
	for (unsigned i=0; i<atomtypes_prm.size(); i++) {
		if (atomtypes_prm[i] == amberatom) {
			this->epsilons.push_back(welldepth[i]);
			this->epsilons_sqrt.push_back(sqrt(welldepth[i]));
		}
	}
}

void Mol2::get_radius(vector<string>atomtypes_prm, string amberatom, vector<double>radius) {
	bool test = false;
	for (unsigned i=0; i<atomtypes_prm.size(); i++) {
		if (atomtypes_prm[i] == amberatom) {

#ifdef DEBUG
			printf("Atom %s == %s. Radii: %.5f\n", amberatom.c_str(), atomtypes_prm[i].c_str(), radius[i]);
#endif

			this->radii.push_back(radius[i]);
			test = true;
		}
	}
	if (! test){
		printf("Radius for atom %s was not found! Please check....\n", amberatom.c_str());
		exit(1);
	}
}

void Mol2::get_masses(string atomname){
	if(atomname.substr(0,1) == "h"){
		this->masses.push_back(1.008);
	}
	else if (atomname.substr(0,1) == "c" and atomname.substr(1,1) != "l"){ // carbon but not clorine
		this->masses.push_back(12.01);
	}
	else if (atomname.substr(0,1) == "n"){
		this->masses.push_back(14.01);
	}
	else if (atomname.substr(0,1) == "o"){
		this->masses.push_back(16.00);
	}
	else if (atomname.substr(0,1) == "p"){
		this->masses.push_back(30.97);
	}
	else if (atomname.substr(0,1) == "s"){
		this->masses.push_back(32.06);
	}
	else if (atomname.substr(0,1) == "f"){
		this->masses.push_back(19.00);
	}
	else if (atomname.substr(0,2) == "cl"){
		this->masses.push_back(32.06);
	}
	else if (atomname.substr(0,2) == "br"){
		this->masses.push_back(79.90);
	}
	else if (atomname.substr(0,1) == "i"){
		this->masses.push_back(126.9);
	}
	else if (atomname.substr(0,1) == "M" and atomname.substr(1,1) == "G"){
		this->masses.push_back(24.305);
	}
	else if (atomname.substr(0,2) == "DU"){
		this->masses.push_back(0.00);
	}
	else if (atomname == "FE"){
		this->masses.push_back(55.0);
	}
	else if (atomname == "C" or atomname == "CA" or atomname == "CB" or atomname == "CC" or atomname == "CD" or atomname == "CK" or atomname == "CM" or atomname == "CN" or atomname == "CQ" or atomname == "CR" or atomname == "CT" or atomname == "CV" or atomname == "CW" or atomname == "CY" or atomname == "CZ" or atomname == "C*"){
		this->masses.push_back(12.01);
	}
	else if (atomname == "N" or atomname == "NA" or atomname == "NB" or atomname == "NC" or atomname == "N2" or atomname == "N3" or atomname == "NT" or atomname == "NY"){
		this->masses.push_back(14.01);
	}
	else if (atomname == "O" or atomname == "O2" or atomname == "OW" or atomname == "OH" or atomname == "OS"){
		this->masses.push_back(16.00);
	}
    else if (atomname == "H" or atomname == "HC" or atomname == "H1" or atomname == "H2" or atomname == "H3" or atomname == "H4" or atomname == "H5" or atomname == "HO" or atomname == "HS" or atomname == "HW" or atomname == "HP" or atomname == "HZ" or atomname == "HA" or atomname == "H0"){
		this->masses.push_back(1.008);
	}
	else if (atomname == "S" or atomname == "SH"){
		this->masses.push_back(32.06);
	}
	else if (atomname == "Zn"){
		this->masses.push_back(53.38);
	}
    else if (atomname == "Cu"){
        this->masses.push_back(63.55);
    }
    else if (atomname == "Ca"){
        this->masses.push_back(40.00);
    }
	else {
		printf("Could not find atomic mass for atom %s\n. Please check Mol2.cpp file.\n", atomname.c_str());
		exit(1);
	}
}

void Mol2::get_trajectory(ifstream &mdcrd){
	if (mdcrd.is_open()){
		getline(mdcrd, atm);
		vector<double> temp;
		vector<vector<double> > xyz;
		double tcoord;
		while(!mdcrd.eof()){
			for (int i=0; i<N; i++){
				for (int j=0; j<3; j++){
					mdcrd >> tcoord;
					temp.push_back(tcoord);
				}
				xyz.push_back(temp);
				temp.clear();
			}
			this->mcoords.push_back(xyz);
			xyz.clear();
		}
	}
	else {
		printf("Could not open trajectory file. Please check!\n");
		exit(1);
	}
}

Mol2::~Mol2(){
	this->xyz.clear();
	this->charges.clear();
	this->radii.clear();
	this->epsilons.clear();
	this->epsilons_sqrt.clear();
	this->resnames.clear();
	this->bonds.clear();
	this->sybyl_atoms.clear();
	this->mcoords.clear();
	this->new_mcoords.clear();
	this->new_xyz.clear();
	this->atomnames.clear();
	this->masses.clear();
	this->amberatoms.clear();

}

void Mol2::get_gaff_parameters(){
	this->atomtypes_prm.clear();
	this->radius.clear();
	this->welldepth.clear();

// Reading GAFF atomtypes;
	this->atomtypes_prm.push_back("h1");		//0
	this->atomtypes_prm.push_back("h2");		//1
	this->atomtypes_prm.push_back("h3");		//2
	this->atomtypes_prm.push_back("h4");		//3
	this->atomtypes_prm.push_back("h5");		//4
	this->atomtypes_prm.push_back("ha");		//5
	this->atomtypes_prm.push_back("hc");		//6
	this->atomtypes_prm.push_back("hn");		//7
	this->atomtypes_prm.push_back("ho");		//8
	this->atomtypes_prm.push_back("hp");		//9
	this->atomtypes_prm.push_back("hs");		//10
	this->atomtypes_prm.push_back("hw");		//11
	this->atomtypes_prm.push_back("hx");		//12
	this->atomtypes_prm.push_back("o");			//13
	this->atomtypes_prm.push_back("oh");		//14
	this->atomtypes_prm.push_back("os");		//15
	this->atomtypes_prm.push_back("ow");		//16
	this->atomtypes_prm.push_back("c");			//17
	this->atomtypes_prm.push_back("cz");		//18
	this->atomtypes_prm.push_back("c1");		//19
	this->atomtypes_prm.push_back("c2");		//20
	this->atomtypes_prm.push_back("c3");		//21
	this->atomtypes_prm.push_back("ca");		//22
	this->atomtypes_prm.push_back("cc");		//23
	this->atomtypes_prm.push_back("cd");		//24
	this->atomtypes_prm.push_back("ce");		//25
	this->atomtypes_prm.push_back("cf");		//26
	this->atomtypes_prm.push_back("cg");		//27
	this->atomtypes_prm.push_back("ch");		//28
	this->atomtypes_prm.push_back("cp");		//29
	this->atomtypes_prm.push_back("cq");		//30
	this->atomtypes_prm.push_back("cu");		//31
	this->atomtypes_prm.push_back("cv");		//32
	this->atomtypes_prm.push_back("cx");		//33
	this->atomtypes_prm.push_back("cy");		//34
	this->atomtypes_prm.push_back("n");			//35
	this->atomtypes_prm.push_back("n1");		//36
	this->atomtypes_prm.push_back("n2");		//37
	this->atomtypes_prm.push_back("n3");		//38
	this->atomtypes_prm.push_back("n4");		//39
	this->atomtypes_prm.push_back("na");		//40
	this->atomtypes_prm.push_back("nb");		//41
	this->atomtypes_prm.push_back("nc");		//42
	this->atomtypes_prm.push_back("nd");		//43
	this->atomtypes_prm.push_back("ne");		//44
	this->atomtypes_prm.push_back("nf");		//45
	this->atomtypes_prm.push_back("nh");		//46
	this->atomtypes_prm.push_back("no");		//47
	this->atomtypes_prm.push_back("s");			//48
	this->atomtypes_prm.push_back("s2");		//49
	this->atomtypes_prm.push_back("s4");		//50
	this->atomtypes_prm.push_back("s6");		//51
	this->atomtypes_prm.push_back("sx");		//52
	this->atomtypes_prm.push_back("sy");		//53
	this->atomtypes_prm.push_back("sh");		//54
	this->atomtypes_prm.push_back("ss");		//55
	this->atomtypes_prm.push_back("p2");		//56
	this->atomtypes_prm.push_back("p3");		//57
	this->atomtypes_prm.push_back("p4");		//58
	this->atomtypes_prm.push_back("p5");		//59
	this->atomtypes_prm.push_back("pb");		//60
	this->atomtypes_prm.push_back("px");		//61
	this->atomtypes_prm.push_back("py");		//62
	this->atomtypes_prm.push_back("f");			//63
	this->atomtypes_prm.push_back("cl");		//64
	this->atomtypes_prm.push_back("br");		//65
	this->atomtypes_prm.push_back("i");			//66
	this->atomtypes_prm.push_back("DU");		//67
												//parm99.dat
	this->atomtypes_prm.push_back("H");
	this->atomtypes_prm.push_back("HO");
	this->atomtypes_prm.push_back("HS");
	this->atomtypes_prm.push_back("HC");
	this->atomtypes_prm.push_back("H1");
	this->atomtypes_prm.push_back("H2");
	this->atomtypes_prm.push_back("H3");
	this->atomtypes_prm.push_back("HP");
	this->atomtypes_prm.push_back("HA");
	this->atomtypes_prm.push_back("H4");
	this->atomtypes_prm.push_back("H5");
	this->atomtypes_prm.push_back("HW");
	this->atomtypes_prm.push_back("HZ");
	this->atomtypes_prm.push_back("O");
	this->atomtypes_prm.push_back("O2");
	this->atomtypes_prm.push_back("OW");
	this->atomtypes_prm.push_back("OH");
	this->atomtypes_prm.push_back("OS");
	this->atomtypes_prm.push_back("C*");
	this->atomtypes_prm.push_back("CT");
	this->atomtypes_prm.push_back("C");
	this->atomtypes_prm.push_back("N");
	this->atomtypes_prm.push_back("N3");
	this->atomtypes_prm.push_back("S");
	this->atomtypes_prm.push_back("SH");
	this->atomtypes_prm.push_back("P");
	this->atomtypes_prm.push_back("IM");
	this->atomtypes_prm.push_back("Li");
	this->atomtypes_prm.push_back("IP");
	this->atomtypes_prm.push_back("Na");
	this->atomtypes_prm.push_back("K");
	this->atomtypes_prm.push_back("Rb");
	this->atomtypes_prm.push_back("Cs");
	this->atomtypes_prm.push_back("MG");
	this->atomtypes_prm.push_back("C0");
	this->atomtypes_prm.push_back("Zn");
	this->atomtypes_prm.push_back("F");
	this->atomtypes_prm.push_back("Cl");
	this->atomtypes_prm.push_back("Br");
	this->atomtypes_prm.push_back("I");
	this->atomtypes_prm.push_back("IB");
	this->atomtypes_prm.push_back("LP");

	this->atomtypes_prm.push_back("CA");
	this->atomtypes_prm.push_back("CB");
	this->atomtypes_prm.push_back("CC");
	this->atomtypes_prm.push_back("CD");
	this->atomtypes_prm.push_back("CK");
	this->atomtypes_prm.push_back("CM");
	this->atomtypes_prm.push_back("CN");
	this->atomtypes_prm.push_back("CQ");
	this->atomtypes_prm.push_back("CR");
	this->atomtypes_prm.push_back("CV");
	this->atomtypes_prm.push_back("CW");
	this->atomtypes_prm.push_back("CY");
	this->atomtypes_prm.push_back("CZ");

	this->atomtypes_prm.push_back("NA");
	this->atomtypes_prm.push_back("NB");
	this->atomtypes_prm.push_back("NC");
	this->atomtypes_prm.push_back("N2");
	this->atomtypes_prm.push_back("NT");
	this->atomtypes_prm.push_back("NY");
    this->atomtypes_prm.push_back("Si"); // from dock6 parameters
    this->atomtypes_prm.push_back("Cu");
    this->atomtypes_prm.push_back("Ca"); // Calcium. Taken from http://personalpages.manchester.ac.uk/staff/Richard.Bryce/amber/ion/metals.dat
    this->atomtypes_prm.push_back("H0"); // Using the same as HA


// Reading GAFF radii. Same order as above
    this->radius.push_back(1.3870);     // 0
    this->radius.push_back(1.2870);     // 1
    this->radius.push_back(1.1870);     // 2
	this->radius.push_back(1.4090);
	this->radius.push_back(1.3590);
	this->radius.push_back(1.4590);
	this->radius.push_back(1.4870);
	this->radius.push_back(0.6000);
	this->radius.push_back(0.6000);
	this->radius.push_back(0.6000);
	this->radius.push_back(0.6000);
	this->radius.push_back(0.6000);
	this->radius.push_back(1.1000);
	this->radius.push_back(1.6612);
	this->radius.push_back(1.7210);
	this->radius.push_back(1.6837);
	this->radius.push_back(1.7683);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(1.75);
	this->radius.push_back(1.948);
	this->radius.push_back(2.22);
	this->radius.push_back(2.35);
	this->radius.push_back(0.0000);
										//parm99.dat
	this->radius.push_back(0.6000);
	this->radius.push_back(0.0000);
	this->radius.push_back(0.6000);
	this->radius.push_back(1.4870);
	this->radius.push_back(1.3870);
	this->radius.push_back(1.2870);
	this->radius.push_back(1.1870);
	this->radius.push_back(1.1000);
	this->radius.push_back(1.4590);
	this->radius.push_back(1.4090);
	this->radius.push_back(1.3590);
	this->radius.push_back(0.0000);
	this->radius.push_back(1.4590);
	this->radius.push_back(1.6612);
	this->radius.push_back(1.6612);
	this->radius.push_back(1.7683);
	this->radius.push_back(1.7210);
	this->radius.push_back(1.6837);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.47);
	this->radius.push_back(1.1370);
	this->radius.push_back(1.8680);
	this->radius.push_back(1.8680);
	this->radius.push_back(2.6580);
	this->radius.push_back(2.9560);
	this->radius.push_back(3.3950);
	this->radius.push_back(0.7926);
	this->radius.push_back(1.7131);
	this->radius.push_back(1.10);
	this->radius.push_back(1.75);
	this->radius.push_back(1.948);
	this->radius.push_back(2.22);
	this->radius.push_back(2.35);
	this->radius.push_back(5.0);
	this->radius.push_back(0.00);

	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);

	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
    this->radius.push_back(2.220); // Si radius from dock6 parameters
    this->radius.push_back(2.20);  // Cu radius
    this->radius.push_back(1.7900); // Calcium
    this->radius.push_back(1.4590); //H0


// Reading GAFF welldepth. Same order as above
    this->welldepth.push_back(0.0157);      // 0
    this->welldepth.push_back(0.0157);      // 1
    this->welldepth.push_back(0.0157);      // 2
    this->welldepth.push_back(0.0150);      // 3
    this->welldepth.push_back(0.0150);      // 4
    this->welldepth.push_back(0.0150);      // 5
    this->welldepth.push_back(0.0157);      // 6
    this->welldepth.push_back(0.0157);      // 7
    this->welldepth.push_back(0.0000);      // 8
    this->welldepth.push_back(0.0157);      // 9
    this->welldepth.push_back(0.0157);      // 10
    this->welldepth.push_back(0.0000);      // 11
    this->welldepth.push_back(0.0157);      // 12
    this->welldepth.push_back(0.2100);      // 13
    this->welldepth.push_back(0.2104);      // 14
    this->welldepth.push_back(0.1700);      // 15
    this->welldepth.push_back(0.1520);      // 16
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.1094);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.0610);
	this->welldepth.push_back(0.2650);
	this->welldepth.push_back(0.3200);
	this->welldepth.push_back(0.4000);
	this->welldepth.push_back(0.0000);
										//parm99.dat
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0000);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0150);
	this->welldepth.push_back(0.0150);
	this->welldepth.push_back(0.0150);
	this->welldepth.push_back(0.0000);
	this->welldepth.push_back(0.0150);
	this->welldepth.push_back(0.2100);
	this->welldepth.push_back(0.2100);
	this->welldepth.push_back(0.1520);
	this->welldepth.push_back(0.2104);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.1094);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.1);
	this->welldepth.push_back(0.0183);
	this->welldepth.push_back(0.00277);
	this->welldepth.push_back(0.00277);
	this->welldepth.push_back(0.000328);
	this->welldepth.push_back(0.00017);
	this->welldepth.push_back(0.0000806);
	this->welldepth.push_back(0.8947);
	this->welldepth.push_back(0.459789);
	this->welldepth.push_back(0.0125);
	this->welldepth.push_back(0.061);
	this->welldepth.push_back(0.265);
	this->welldepth.push_back(0.320);
	this->welldepth.push_back(0.40);
	this->welldepth.push_back(0.1);
	this->welldepth.push_back(0.0000);

	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);

	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
    this->welldepth.push_back(0.320); //Si welldepth from dock6 parameters
    this->welldepth.push_back(0.200); //Cu welldepth
    this->welldepth.push_back(0.0140); // Calcium welldepth
    this->welldepth.push_back(0.0150); // H0

	bool fail = true;
	if ((this->atomtypes_prm.size() == this->radius.size()) and (this->atomtypes_prm.size() == this->welldepth.size())){
		fail = false;
	}
	if (fail){
		printf("Size of vectors with GAFF parameters does not match. Please check Mol2.cpp file!\n");
		exit(1);
	}
}


bool Mol2::parse_gzipped_ensemble(string molfile, int skipper=1){
    char tstr[80];
    bool bret = false;
    int tint;
    float tx, ty, tz;
    vector<double> txyz;
    int tres;
    float tcharge;
    int count=0;
    char tatomtype[10];
    char resname[10];
    string cpstr;
    vector<vector<double> > tcoord;
    int trajsize=0;

    gzFile mol2file = gzopen(molfile.c_str(), "r");
    if (mol2file != NULL){
        str[0]='#';
        while(str[0] !='@'){
            gzgets(mol2file, str, 80);

        }
        gzgets(mol2file, str, 80);
        this->molname = str;
        this->molname = this->molname.substr(0,this->molname.size()-1);
        gzgets(mol2file, str, 80);
        sscanf(str, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);


        cpstr = string(str);
        while (cpstr.substr(0,6) != "Energy"){
            gzgets(mol2file, str, 80);
            cpstr = string(str);
        }
        sscanf(str, "%s %f\n", tstr, &tx);              //parsing ensemble energy
        this->ensemble_energies.push_back(double(tx));
        trajsize++;

        while (cpstr.substr(0,13) != "RMSD/OVERLAY:"){
            gzgets(mol2file, str, 80);
            cpstr = string(str);
        }
        sscanf(str, "%s %f\n", tstr, &tx);              //parsing ensemble rmsd
        this->ensemble_rmsd.push_back(double(tx));

        cpstr = string(str);
        while (cpstr.substr(0,13) != "@<TRIPOS>ATOM"){
            gzgets(mol2file, str, 80);
            cpstr = string(str);
        }

        for (int i=0; i<this->N; i++){
            gzgets(mol2file, str, 80);
            sscanf(str, "%d %s %f %f %f %s %d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
            txyz.push_back(tx);
            txyz.push_back(ty);
            txyz.push_back(tz);
            this->xyz.push_back(txyz);
            txyz.clear();

            this->charges.push_back(tcharge);
            this->atomnames.push_back(str);

            this->amberatoms.push_back(this->convert2gaff2(string(tatomtype)));
            this->sybyl_atoms.push_back(string(tatomtype));
            this->get_masses(amberatoms[i]);

            if (tres > count){
                this->residue_pointer.push_back(i+1);
                count = tres;
                this->resnames.push_back(string(resname));
            }
        }

        gzgets(mol2file, str, 80);
        if (str[0] != '@'){
            while (str[0] != '@'){
                gzgets(mol2file, str, 80);
            }
        }

        vector<string> bond;
        char s1[6], s2[6], s3[5];
        for (int i=0; i<this->Nbonds; i++){
            gzgets(mol2file, str, 80);
            sscanf(str, "%d%s%s%s\n", &tint, s1, s2, s3);
            bond.push_back(string(s1));
            bond.push_back(string(s2));
            bond.push_back(string(s3));
            this->bonds.push_back(bond);
            bond.clear();
        }

        int n=0;

        this->mcoords.push_back(this->xyz);

        while (! gzeof(mol2file)){
            gzgets(mol2file, str, 80);
            while ((str[0] != 'E' or str[6] != ':') and (!gzeof(mol2file))){
                gzgets(mol2file, str, 80);
            }

            if (!gzeof(mol2file)){
                sscanf(str, "%s %f\n", tstr, &tx);
                trajsize++;
                if (trajsize % skipper == 0){
                    this->ensemble_energies.push_back(double(tx));
                }
            }

            while ((str[0] != 'R' or str[12] != ':') and (!gzeof(mol2file))){
                gzgets(mol2file, str, 80);
            }

            if (!gzeof(mol2file)){
                sscanf(str, "%s %f\n", tstr, &tx);
                if (trajsize % skipper == 0){
                    this->ensemble_rmsd.push_back(double(tx));
                }
            }

            while ((str[0] != '@' or str[9] != 'A') and (!gzeof(mol2file))){
                gzgets(mol2file, str, 80);
            }
            if (!gzeof(mol2file) and (trajsize % skipper == 0)){
                txyz.clear();
                for (int i=0; i<this->N; i++){
                    gzgets(mol2file, str, 80);
                    sscanf(str, "%d %s %f %f %f %s %d %s %f\n", &tint, tstr, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);
                    txyz.push_back(tx);
                    txyz.push_back(ty);
                    txyz.push_back(tz);
                    tcoord.push_back(txyz);
                    txyz.clear();
                }
                this->mcoords.push_back(tcoord);
                n++;
                tcoord.clear();
            }
        }
#ifdef DEBUG
        printf("Found %d conformations in file %s\n", n, molfile.c_str());
#endif
        bret = true;
    }

    else {
        printf("Skipping file %s...\n", molfile.c_str());
        bret = false;
    }
    gzclose(mol2file);
    return (bret);
}


void Mol2::initialize_gaff(){
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
    if (gaff_file!= NULL){
        while (!feof(gaff_file)){
            fgets(str, 80, gaff_file);
            if (str[0] != '#'){
                sscanf(str, "%s %f %f", at, &r, &e);
                atom_param v;
                v.type = string(at);
                v.radius = double(r);
                v.epsilon = double(e);

                this->gaff_force_field.push_back(v);
            }
        }
    }

    fclose(gaff_file);
}


void Mol2::get_gaff_atomic_parameters(string gaff_atom, atom_param* ap){
    for (unsigned i=0; i<this->gaff_force_field.size(); i++){
        if (gaff_force_field[i].type == gaff_atom){
            ap = &(gaff_force_field[i]);
        }
    }
}

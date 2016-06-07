#include "MC.h"

MC::MC(WRITER* _Writer)
{
    srand(rand());
    r = gsl_rng_alloc (gsl_rng_ranlxs2);
    Writer = _Writer;
}

MC::MC(Mol2* Lig, PARSER* Input){

/*! Here we will use OpenBabel API to generate a OBMol and use it to get and set torsion angles for all
 * rotatable bonds.
*/

    mol = this->GetMol(Input->lig_mol2);
    OBff = OBForceField::FindForceField("GAFF");

#ifdef DEBUG
    OBff->SetLogFile(&cout);
    OBff->SetLogLevel(OBFF_LOGLVL_LOW);
#endif

    if (!OBff){
        cout << "Could not find FF GAFF!" << endl;
        exit(1);
    }

    OBff->Setup(*mol);
    RotorList.setup(*mol);
    Rotor = RotorList.BeginRotor(RotorIterator);
    mol->toInertialFrame();


    vector<double> tmp(4);
    for (int i = 1; i < RotorList.Size() + 1; ++i, Rotor = RotorList.NextRotor(RotorIterator)) {
        for (unsigned j=0; j < 4; j++){
            tmp[j] = Rotor->GetDihedralAtoms()[j];
        }
        atoms_in_dihedrals.push_back(tmp);
        tmp.clear();
    }

}


MC::~MC(){
    gsl_rng_free (r);
    delete OBff;
    delete mol;
}

void MC::write_conformers(Mol2* Lig){
    for (unsigned i=0; i< Lig->mcoords.size(); i++){
        Writer->writeMol2(Lig, Lig->mcoords[i], 0.0, 0.0, "teste");
    }
}

void MC::run(Grid* Grids, Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T){
    this->xyz = xyz;
    int nReject = 0;

    double sum_x = 0.0;
    double sum_xsquared = 0.0;
    long double sum_Boltzmann_ene = 0.0;
    long double sum_Boltzmann2_ene = 0.0;

    double k=0.0019858775203792202;

    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);

    gsl_rng_set(r, Input->seed);

    sprintf(info, "Starting Monte Carlo equilibration simulation with %5d steps...", (Input->eq_steps));
    Writer->print_info(info);
    gzFile mc_output;
    mc_output = gzopen((Input->output + "_mc.dat.gz").c_str(), "w");

    Energy2* Energy = new Energy2(Input);
    COORD_MC* Coord = new COORD_MC;
    double energy, new_energy, p, rnumber, rmsd;
    step_t* step = new step_t;

    if (Input->generate_conformers){
        this->take_step_flex(Input, Lig, step);
    }
    else {
        this->take_step(Input, Lig, step);
    }

    energy = (Energy->compute_ene(Grids, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);



    sprintf(info, "#%10s %10s %10s", "Step", "Energy", "RMSD");
    gzprintf(mc_output,"##################################################################################################################################\n");
    gzprintf(mc_output,"#BoxsideX= %10.4f BoxsideY= %10.4f BoxsideY= %10.4f Temperature= %10.4f \n", Input->x_dim, Input->y_dim, Input->z_dim, Input->temp);
    gzprintf(mc_output,"##################################################################################################################################\n");
    gzprintf(mc_output, "#      Step    Energy      RMSD    DX         DY        DZ          ALPHA     BETA        GAMMA        NCONF     ConfEnergy\n");

    int count=0;
    int eqcount = 0;
    long double xcom = 0;
    long double ycom = 0;
    long double zcom = 0;

    //Equilibration implementation

    while (eqcount <= Input->eq_steps){

        if (Input->generate_conformers){
            this->take_step_flex(Input, Lig, step);
        }
        else {
            this->take_step(Input, Lig, step);
        }

        new_energy = (Energy->compute_ene(Grids, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);

        if (new_energy <= energy){
            Lig->mcoords = Lig->new_mcoords;
            this->xyz = step->xyz;
            energy = new_energy;
            rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
            xcom = xcom + step->dx;
            ycom = ycom + step->dy;
            zcom = zcom + step->dz;
            this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);
        }
        else{
            p = this->Boltzmman(energy, new_energy, T, Input->bi);
            rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
            if (p > rnumber){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                xcom = xcom + step->dx;
                ycom = ycom + step->dy;
                zcom = zcom + step->dz;
                this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);
            }
        }
        eqcount++;
    }

    Writer->print_line();
    sprintf(info, "Equilibration done with %5d steps. Current system energy: %9.3f kcal/mol.", Input->eq_steps, energy);
    Writer->print_info(info);
    Writer->print_line();

    while (count <= Input->number_steps){

        if (Input->generate_conformers){
            this->take_step_flex(Input, Lig, step);
        }
        else {
            this->take_step(Input, Lig, step);
        }

        new_energy = (Energy->compute_ene(Grids, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);

        if (new_energy <= energy){
            Lig->mcoords = Lig->new_mcoords;
            this->xyz = step->xyz;
            energy = new_energy;
            rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
            xcom = xcom + step->dx;
            ycom = ycom + step->dy;
            zcom = zcom + step->dz;
            this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);
            count++;
            sum_x += energy;
            sum_xsquared += (energy*energy);
            sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
            sum_Boltzmann2_ene += exp((-2.0*energy*(Input->bi-1.0))/(k*T));

            if ((Input->write_mol2) and (count % Input->mc_stride == 0)){
                Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, Input->output + "_MC");
            }
            gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f\n", count, energy, rmsd, step->dx, step->dy, step->dz, step->dalpha, step->dbeta, step->dgamma, step->nconf, Lig->conformer_energies[step->nconf]);
//            sprintf(info, "%10d %10.3f %10.3f", count, energy, rmsd);
            if (count % Input->mc_stride == 0){
                sprintf(info, "Accepted steps: %9d. Current energy for the system: %7.3f kcal/mol.",count, energy);
                Writer->print_info(info);
            }
        }
        else{
            p = this->Boltzmman(energy, new_energy, T, Input->bi);
            rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
            if (p > rnumber){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                xcom = xcom + step->dx;
                ycom = ycom + step->dy;
                zcom = zcom + step->dz;
                this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);
                count++;
                sum_x += energy;
                sum_xsquared += (energy*energy);
                sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
                sum_Boltzmann2_ene += exp((-2.0*energy*(Input->bi-1.0))/(k*T));

                if ((Input->write_mol2) and (count % Input->mc_stride == 0)){
                    Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, Input->output + "_MC");
                }
                gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f\n", count, energy, rmsd, step->dx, step->dy, step->dz, step->dalpha, step->dbeta, step->dgamma, step->nconf, Lig->conformer_energies[step->nconf]);
                if (count % Input->mc_stride == 0){
                    sprintf(info, "Accepted steps: %9d. Current energy for the system: %7.3f kcal/mol.",count, energy);
                    Writer->print_info(info);
                }
            }
            else{
                nReject++;
                sum_x += energy;
                sum_xsquared += (energy*energy);
                sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
                sum_Boltzmann2_ene += exp((-2.0*energy*(Input->bi-1.0))/(k*T));
            }
        }
    }

    gzclose(mc_output);
    delete step;
    delete Energy;
    delete Coord;
    Writer->print_line();
    sprintf(info, "Finished MC simulation after %d accepted steps, %d total steps, acceptance rate %5.3f", Input->number_steps, (Input->number_steps+nReject),double((Input->number_steps*1.0)/(Input->number_steps+nReject)));
    Writer->print_info(info);
    Writer->print_line();
    sprintf(info, "Conformer energies (GAFF):");
    Writer->print_info(info);

    for (unsigned i=0; i< Lig->conformer_energies.size(); i++){
        sprintf(info, "%5d %10.3f kcal/mol", i+1, Lig->conformer_energies[i]);
        Writer->print_info(info);
    }
    Writer->print_line();

    this->XSize = this->MaxMin[1] - this->MaxMin[0];
    this->YSize = this->MaxMin[3] - this->MaxMin[2];
    this->ZSize = this->MaxMin[5] - this->MaxMin[4];
    sprintf(info, "Max Dimensions:");
    Writer->print_info(info);
    sprintf(info, "%10.3f %10.3f %10.3f", this->XSize, this->YSize, this->ZSize);
    Writer->print_info(info);
    Writer->print_line();

    this->average_energy = double(sum_x/(Input->number_steps+nReject));
    this->energy_standard_deviation = ((sum_xsquared/(Input->number_steps+nReject)) - (this->average_energy*this->average_energy));
    this->energy_standard_deviation = sqrt(this->energy_standard_deviation/(Input->number_steps+nReject));
    this->Boltzmann_weighted_average_energy = sum_Boltzmann2_ene/sum_Boltzmann_ene;

    sprintf(info, "Average Monte Carlo energy: %10.3f +- %10.3f @ %7.2f K", this->average_energy, this->energy_standard_deviation, T);
    Writer->print_info(info);

    Writer->print_line();

}


void MC::ligand_run(Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T){

    if (Input->ligsim){

        double sum_x = 0.0;
        double sum_xsquared = 0.0;
        long double sum_Boltzmann_ene = 0.0;
        long double sum_Boltzmann2_ene = 0.0;

        double k=0.0019858775203792202;

        for (int i=0; i< 6; i++){
            this->MaxMin.push_back(0.0);
        }

        this->xyz = xyz;
        int nReject = 0;

        gsl_rng_set(r, Input->seed);

        sprintf(info, "Ligand Monte Carlo simulation with %5d steps", (Input->number_steps));
        Writer->print_info(info);

        gzFile mc_output_lig;
        mc_output_lig = gzopen((Input->output + "_mc.ligsim.dat.gz").c_str(), "w");

        Energy2* Energy = new Energy2(Input);
        COORD_MC* Coord = new COORD_MC;
        double energy, new_energy, p, rnumber, rmsd;
        step_t* step = new step_t;

        if (Input->generate_conformers){
            this->take_step_flex(Input, Lig, step);
        }
        else {
            this->take_step(Input, Lig, step);
        }

        energy = Lig->conformer_energies[step->nconf];

        gzprintf(mc_output_lig,"##################################################################################################################################\n");
        gzprintf(mc_output_lig,"#BoxsideX= %10.4f BoxsideY= %10.4f BoxsideY= %10.4f Temperature= %10.4f \n", Input->x_dim, Input->y_dim, Input->z_dim, Input->temp);
        gzprintf(mc_output_lig,"##################################################################################################################################\n");
        gzprintf(mc_output_lig, "#      Step    Energy      RMSD    NCONF\n");

        int count=0;
        double xcom = 0;
        double ycom = 0;
        double zcom = 0;


        Writer->print_line();

        while (count <= Input->number_steps){

            if (Input->generate_conformers){
                this->take_step_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }

            new_energy = Lig->conformer_energies[step->nconf];

            if (new_energy <= energy){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                xcom = xcom + step->dx;
                ycom = ycom + step->dy;
                zcom = zcom + step->dz;
                this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);

                gzprintf(mc_output_lig, "%10d %10.3f %10.3f  %10d %.3lf %.3lf %.3lf \n", count, energy, rmsd, step->nconf, xcom , ycom, zcom);

                if (count % Input->mc_stride == 0){
                    sprintf(info, "Accepted steps: %9d. Current energy for the system: %7.3f kcal/mol.",count, energy);
                    Writer->print_info(info);
                }


                if ((Input->write_mol2) and (count % Input->mc_stride == 0)){
                    Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, Input->output + "_MC.ligsim");
                }

                count++;
                sum_x += energy;
                sum_xsquared += (energy*energy);
                sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
                sum_Boltzmann2_ene += exp((-2.0*energy*(Input->bi-1.0))/(k*T));
            }
            else{
                p = this->Boltzmman(energy, new_energy, T, Input->bi);
                rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
                if (p > rnumber){
                    Lig->mcoords = Lig->new_mcoords;
                    this->xyz = step->xyz;
                    energy = new_energy;
                    rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                    xcom = xcom + step->dx;
                    ycom = ycom + step->dy;
                    zcom = zcom + step->dz;
                    this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);

                    gzprintf(mc_output_lig, "%10d %10.3f %10.3f  %10d %.3lf %.3lf %.3lf \n", count, energy, rmsd, step->nconf, xcom , ycom, zcom);

                    if (count % Input->mc_stride == 0){
                        sprintf(info, "Accepted steps: %9d. Current energy for the system: %7.3f kcal/mol.",count, energy);
                        Writer->print_info(info);
                    }

                    if ((Input->write_mol2) and (count % Input->mc_stride == 0)){
                        Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, Input->output + "_MC.ligsim");
                    }
                    count++;
                    sum_x += energy;
                    sum_xsquared += (energy*energy);
                    sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
                    sum_Boltzmann2_ene += exp((-2.0*energy*(Input->bi-1.0))/(k*T));
                }
                else{
                    nReject++;
                    sum_x += energy;
                    sum_xsquared += (energy*energy);
                    sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
                    sum_Boltzmann2_ene += exp((-2.0*energy*(Input->bi-1.0))/(k*T));
                }
            }
        }

        gzclose(mc_output_lig);
        delete step;
        delete Energy;
        delete Coord;
        Writer->print_line();
        sprintf(info, "Finished MC simulation for Ligand after %d accepted steps, %d total steps, acceptance rate %5.3f", Input->number_steps, (Input->number_steps+nReject),(Input->number_steps*1.0)/(Input->number_steps+nReject));
        Writer->print_info(info);
        Writer->print_line();

        this->XSize = this->MaxMin[1] - this->MaxMin[0];
        this->YSize = this->MaxMin[3] - this->MaxMin[2];
        this->ZSize = this->MaxMin[5] - this->MaxMin[4];
        sprintf(info, "Max Dimensions:");
        Writer->print_info(info);
        sprintf(info, "%10.3f %10.3f %10.3f", this->XSize, this->YSize, this->ZSize);
        Writer->print_info(info);
        Writer->print_line();

        this->average_energy = double(sum_x/(Input->number_steps+nReject));
        this->energy_standard_deviation = ((sum_xsquared/(Input->number_steps+nReject)) - (this->average_energy*this->average_energy));
        this->energy_standard_deviation = sqrt(this->energy_standard_deviation/(Input->number_steps+nReject));
        this->Boltzmann_weighted_average_energy = sum_Boltzmann2_ene/sum_Boltzmann_ene;
    }
}


void MC::run(Mol2* Rec, Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input, double T){
    this->xyz = xyz;
    int nReject = 0;

    double sum_x = 0.0;
    double sum_xsquared = 0.0;
    long double sum_Boltzmann_ene = 0.0;

    double k=0.0019858775203792202;

    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);
    this->MaxMin.push_back(0.0);

    gsl_rng_set(r, Input->seed);

    if (Input->eq_mode){
        sprintf(info, "Starting Monte Carlo equilibration simulation with %5d steps...", (Input->eq_steps));
        Writer->print_info(info);
        gzFile mc_output;
        mc_output = gzopen((Input->output + "_mc.dat.gz").c_str(), "w");

        Energy2* Energy = new Energy2(Input);
        COORD_MC* Coord = new COORD_MC;
        double energy, new_energy, p, rnumber, rmsd;
        step_t* step = new step_t;

        if (Input->generate_conformers){
            this->take_step_flex(Input, Lig, step);
        }
        else {
            this->take_step(Input, Lig, step);
        }

        energy = (Energy->compute_ene(Rec, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);



        sprintf(info, "#%10s %10s %10s", "Step", "Energy", "RMSD");
        gzprintf(mc_output,"##################################################################################################################################\n");
        gzprintf(mc_output,"#BoxsideX= %10.4f BoxsideY= %10.4f BoxsideY= %10.4f Temperature= %10.4f \n", Input->x_dim, Input->y_dim, Input->z_dim, Input->temp);
        gzprintf(mc_output,"##################################################################################################################################\n");
        gzprintf(mc_output, "#      Step    Energy      RMSD    DX         DY        DZ          ALPHA     BETA        GAMMA        NCONF     ConfEnergy\n");

        int count=0;
        int eqcount = 0;
        long double xcom = 0;
        long double ycom = 0;
        long double zcom = 0;

        //Equilibration implementation
        while (eqcount <= Input->eq_steps){

            if (Input->generate_conformers){
                this->take_step_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }

            new_energy = (Energy->compute_ene(Rec, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);

            if (new_energy <= energy){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                xcom = xcom + step->dx;
                ycom = ycom + step->dy;
                zcom = zcom + step->dz;
                this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);
            }
            else{
                p = this->Boltzmman(energy, new_energy, T, Input->bi);
                rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
                if (p > rnumber){
                    Lig->mcoords = Lig->new_mcoords;
                    this->xyz = step->xyz;
                    energy = new_energy;
                    rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                    xcom = xcom + step->dx;
                    ycom = ycom + step->dy;
                    zcom = zcom + step->dz;
                    this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);
                }
            }
            eqcount++;
        }

        Writer->print_line();
        sprintf(info, "Equilibration done with %5d steps. Current system energy: %9.3f kcal/mol.", Input->eq_steps, energy);
        Writer->print_info(info);
        Writer->print_line();

        while (count <= Input->number_steps){

            if (Input->generate_conformers){
                this->take_step_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }

            new_energy = (Energy->compute_ene(Rec, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);

            if (new_energy <= energy){
                Lig->mcoords = Lig->new_mcoords;
                this->xyz = step->xyz;
                energy = new_energy;
                rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                xcom = xcom + step->dx;
                ycom = ycom + step->dy;
                zcom = zcom + step->dz;
                this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);
                count++;
                sum_x += energy;
                sum_xsquared += (energy*energy);
                sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
            }
            else{
                p = this->Boltzmman(energy, new_energy, T, Input->bi);
                rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
                if (p > rnumber){
                    Lig->mcoords = Lig->new_mcoords;
                    this->xyz = step->xyz;
                    energy = new_energy;
                    rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                    xcom = xcom + step->dx;
                    ycom = ycom + step->dy;
                    zcom = zcom + step->dz;
                    this->MaxMinCM(xcom,ycom,zcom,this->MaxMin);
                    count++;
                    sum_x += energy;
                    sum_xsquared += (energy*energy);
                    sum_Boltzmann_ene += exp(((-Input->bi-1.0)*energy)/(k*T));
                }
                else{
                    nReject++;
                    sum_x += energy;
                    sum_xsquared += (energy*energy);
                }
            }

            if ((Input->write_mol2) and (count % Input->mc_stride == 0)){
                Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, Input->output + "_MC");
            }
            gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f\n", count, energy, rmsd, step->dx, step->dy, step->dz, step->dalpha, step->dbeta, step->dgamma, step->nconf, Lig->conformer_energies[step->nconf]);
            sprintf(info, "%10d %10.3f %10.3f", count, energy, rmsd);
            count++;
            if (count % Input->mc_stride == 0){
                sprintf(info, "Accepted steps: %9d. Current energy for the system: %7.3f kcal/mol.",count, energy);
                Writer->print_info(info);
            }
        }

        gzclose(mc_output);
        delete step;
        delete Energy;
        delete Coord;
        Writer->print_line();
        sprintf(info, "Finished MC simulation after %d accepted steps, %d total steps, acceptance rate %5.3f", Input->number_steps, (Input->number_steps+nReject),double((Input->number_steps*1.0)/(Input->number_steps+nReject)));
        Writer->print_info(info);
        Writer->print_line();
        sprintf(info, "Conformer energies (GAFF):");
        Writer->print_info(info);

        for (unsigned i=0; i< Lig->conformer_energies.size(); i++){
            sprintf(info, "%5d %10.3f kcal/mol", i, Lig->conformer_energies[i]);
            Writer->print_info(info);
        }
        Writer->print_line();

        this->XSize = this->MaxMin[1] - this->MaxMin[0];
        this->YSize = this->MaxMin[3] - this->MaxMin[2];
        this->ZSize = this->MaxMin[5] - this->MaxMin[4];
        sprintf(info, "Max Dimensions:");
        Writer->print_info(info);
        sprintf(info, "%10.3f %10.3f %10.3f", this->XSize, this->YSize, this->ZSize);
        Writer->print_info(info);
        Writer->print_line();

        this->average_energy = double(sum_x/(Input->number_steps+nReject));
        this->energy_standard_deviation = ((sum_xsquared/(Input->number_steps+nReject)) - (this->average_energy*this->average_energy));
        this->energy_standard_deviation = sqrt(this->energy_standard_deviation/(Input->number_steps+nReject));
        this->Boltzmann_weighted_average_energy = sum_Boltzmann_ene/(Input->number_steps+nReject);

        sprintf(info, "Average Monte Carlo energy: %10.3f kcal/mol +- (%10.3f kcal/mol) @ %7.2f K", this->average_energy, this->energy_standard_deviation, T);
        Writer->print_info(info);
        Writer->print_line();
    }
}


void MC::MaxMinCM(double XCOM, double YCOM, double ZCOM, vector <double> Max){

    if (XCOM < Max[0]){
        this->MaxMin[0]=XCOM;
    }
    if (XCOM > Max[1]){
        this->MaxMin[1]=XCOM;
    }
    if (YCOM < Max[2]){
        this->MaxMin[2]=YCOM;
    }
    if (YCOM > Max[3]){
        this->MaxMin[3]=YCOM;
    }
    if (ZCOM < Max[4]){
        this->MaxMin[4]=ZCOM;
    }
    if (ZCOM > Max[5]){
        this->MaxMin[5]=ZCOM;
    }
}


void MC::take_step(PARSER* Input, Mol2* Lig, step_t* step){
    COORD_MC* Coord = new COORD_MC;
    double rnumber; //transx, transy, transz, a, b, g, rnumber;

    rnumber = gsl_rng_uniform(r);
    step->dx = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dy = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dz = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

    rnumber = gsl_rng_uniform(r);
    step->dalpha = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dbeta = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dgamma = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

    step->nconf = 0;

    step->xyz = Coord->rototranslate(this->xyz, Lig, step->dalpha, step->dbeta,step->dgamma, step->dx, step->dy, step->dz);


    delete Coord;
}

void MC::take_step_flex(PARSER* Input, Mol2* Lig, step_t* step){
    COORD_MC* Coord = new COORD_MC;
    double rnumber; //transx, transy, transz, a, b, g, rnumber;
    int ln=0;

    rnumber = gsl_rng_uniform(r);
    step->dx = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dy = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dz = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

    rnumber = gsl_rng_uniform(r);
    step->dalpha = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dbeta = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dgamma = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

    Coord->rototranslate_all(Lig, step->dalpha, step->dbeta, step->dgamma, step->dx, step->dy, step->dz);

    ln = gsl_rng_uniform_int(r, Lig->mcoords.size());
    step->nconf = ln;

    delete Coord;
    step->xyz = Lig->new_mcoords[ln];
}

void MC::take_step_torsion(PARSER* Input, Mol2* Lig, step_t* step){

    COORD_MC* Coord = new COORD_MC;
    double rnumber; //transx, transy, transz, a, b, g, rnumber;
    int ln=0;

// Do rotation and translation

    rnumber = gsl_rng_uniform(r);
    step->dx = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dy = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));
    rnumber = gsl_rng_uniform(r);
    step->dz = -Input->cushion + (1.0 * (rnumber*(2*Input->cushion)));

    rnumber = gsl_rng_uniform(r);
    step->dalpha = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dbeta = -Input->rotation_step + (rnumber*(2*Input->rotation_step));
    rnumber = gsl_rng_uniform(r);
    step->dgamma = -Input->rotation_step + (rnumber*(2*Input->rotation_step));

    Coord->rototranslate_all(Lig, step->dalpha, step->dbeta, step->dgamma, step->dx, step->dy, step->dz);


// Copy coordinates to OBMol

    double* xyz = new double[mol->NumAtoms()*3];
    for (unsigned i=0; i <mol->NumAtoms(); i++){
        xyz[3*i] = Lig->new_mcoords[i][0];
        xyz[(3*i)+1] = Lig->new_mcoords[i][1];
        xyz[(3*i)+2] = Lig->new_mcoords[i][2];
    }

    mol->SetCoordinates(xyz);

// Do torsion search

    double current_angle;
    for (int i=0; i< RotorList.Size(); i++){
    rnumber = gsl_rng_uniform(r);
    current_angle = mol->GetTorsion(mol->GetAtom(atoms_in_dihedrals[i][0]), mol->GetAtom(atoms_in_dihedrals[i][1]), mol->GetAtom(atoms_in_dihedrals[i][2]), mol->GetAtom(atoms_in_dihedrals[i][3]));
    mol->SetTorsion(mol->GetAtom(atoms_in_dihedrals[i][0]), mol->GetAtom(atoms_in_dihedrals[i][1]), mol->GetAtom(atoms_in_dihedrals[i][2]), mol->GetAtom(atoms_in_dihedrals[i][3]),
            ((current_angle + (-Input->torsion_step + (rnumber*(2*Input->torsion_step))))*PI/180));
    }
    step->torsion_angles.push_back(((current_angle + (-Input->torsion_step + (rnumber*(2*Input->torsion_step))))*PI/180));

    xyz = mol->GetCoordinates();
}





double MC::Boltzmman(double ene, double new_ene, double t, double b){
    double de = (new_ene - ene);
    double x = (-(b-1.0)*de) / (0.0019858775203792202*t); // k=0.0019858775203792202 kcal/(mol.K)
    return(exp(x));
}

vector<vector<double> > copy_from_obmol(OBMol* mymol){

}

double copy_to_obmol(vector<vector<double> > xyz){

}

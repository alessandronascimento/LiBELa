#include "mcintegrate.h"

MC::MC(WRITER* _Writer)
{
    srand(rand());
    r = gsl_rng_alloc (gsl_rng_ranlxs2);
    Writer = _Writer;
}

MC::~MC(){
    gsl_rng_free (r);
}

void MC::write_conformers(Mol2* Lig){
    for (unsigned i=0; i< Lig->mcoords.size(); i++){
        Writer->writeMol2(Lig, Lig->mcoords[i], 0.0, 0.0, "teste");
    }
}

void MC::run(Grid* Grids, Mol2* RefLig, Mol2* Lig, vector<vector<double> > xyz, PARSER* Input){
    this->xyz = xyz;
    gsl_rng_set(r, Input->seed);

    if (Input->eq_mode){
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
#ifdef DEBUG
        Writer->print_info(info);
#endif
        gzprintf(mc_output, "#      Step    Energy      RMSD    DX         DY        DZ          ALPHA     BETA        GAMMA        NCONF     ConfEnergy\n");

        int count=0;
        int eqcount = 0;

        //Equilibration implementation
        while (eqcount <= Input->eq_steps){

            if (Input->generate_conformers){
                this->take_step_flex(Input, Lig, step);
            }
            else {
                this->take_step(Input, Lig, step);
            }
            new_energy = (Energy->compute_ene(Grids, Lig, step->xyz)+Lig->conformer_energies[step->nconf]);
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

            if (Input->write_mol2){
                Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, Input->output + "_MC");
            }
            gzprintf(mc_output, "%10d %10.3f %10.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6d %10.3f\n", count, energy, rmsd, step->dx, step->dy, step->dz, step->dalpha, step->dbeta, step->dgamma, step->nconf, Lig->conformer_energies[step->nconf]);
            sprintf(info, "%10d %10.3f %10.3f", count, energy, rmsd);
            count++;
            if (count % 1000 == 0){
                sprintf(info, "Accepted steps: %9d. Current energy for the system: %7.3f kcal/mol.",count, energy);
                Writer->print_info(info);
            }
        }

        double sup = (double) Input->number_steps;
        double low = (double)((Input->number_steps+nReject));
        gzclose(mc_output);
        delete step;
        delete Energy;
        delete Coord;
        Writer->print_line();
        sprintf(info, "Finished MC simulation after %d accepted steps, %d total steps, acceptance rate %5.3f", Input->number_steps, (Input->number_steps+nReject),(sup/low));
        Writer->print_info(info);
        Writer->print_line();
        sprintf(info, "Conformer energies (GAFF):");
        Writer->print_info(info);
        for (unsigned i=0; i< Lig->conformer_energies.size(); i++){
            sprintf(info, "%5d %10.3f kcal/mol", i, Lig->conformer_energies[i]);
            Writer->print_info(info);
        }
        Writer->print_line();


/*
        // Ligand Simulation implementation.

        if (Input->ligsim){

            this->xyz = xyz;
            int nReject = 0;

            gsl_rng_set(r, Input->seed);

            if (Input->ligsim){
                double nsweeps = (Input->number_steps/(Input->sweep_steps));
                double ST =  Input->sweep_steps;
                sprintf(info, "Ligand Monte Carlo simulation with sweep size of %5d and %5d sweeps", (ST , nsweeps));
                Writer->print_info(info);
                gzFile mc_output_lig;
                mc_output_lig = gzopen((Input->output + "_mc.ligsim.dat.gz").c_str(), "w");
                Energy2* Energy = new Energy2(Input);
                COORD_MC* Coord = new COORD_MC;
                double energy, new_energy, p, rnumber, rmsd;
                step_t* step = new step_t;
                vector<vector<vector<double> > > first_mcoords;
                if (Input->generate_conformers){
                    this->take_step_flex(Input, Lig, step);
                }
                else {
                    this->take_step(Input, Lig, step);
                }

                energy = Lig->conformer_energies[0];

#ifdef DEBUG
                Writer->print_info(info);
#endif
                gzprintf(mc_output_lig,"##################################################################################################################################\n");
                gzprintf(mc_output_lig,"#BoxsideX= %10.4f BoxsideY= %10.4f BoxsideY= %10.4f Temperature= %10.4f \n", Input->x_dim, Input->y_dim, Input->z_dim, Input->temp);
                gzprintf(mc_output_lig,"##################################################################################################################################\n");
                gzprintf(mc_output_lig, "#      Step    Energy      RMSD    NCONF\n");

                int count=0;
                int eqcount = 0;
                double xcom = 0;
                double ycom = 0;
                double zcom = 0;


                Writer->print_line();
                Writer->print_info(info);

                while (count <= (Input->number_steps/Input->sweep_steps)){


                    int sweepcount = 1;
                    while(sweepcount <= Input->sweep_steps){

                        if (Input->generate_conformers){
                            this->take_step_flex(Input, Lig, step);
                        }
                        else {
                            this->take_step(Input, Lig, step);
                        }
                        int rand;


                        rand = gsl_rng_uniform_int(r,Lig->conformer_energies.size());

                        new_energy = Lig->conformer_energies[rand];
                        if(abs(xcom)>Input->x_dim/2 or abs(ycom)>Input->y_dim/2 or abs(zcom)>Input->z_dim/2){
                            //if((xcom)>this->MaxMin[1] or (xcom)<this->MaxMin[0] or (ycom)>this->MaxMin[3] or (ycom)<this->MaxMin[2] or (zcom)>this->MaxMin[5] or (zcom)<this->MaxMin[4]){
                            Coord->rototranslate_all(Lig, 0, 0, 0, -xcom, -ycom, -zcom);
                            Lig->mcoords = Lig->new_mcoords;
                            step->nconf = rand;
                            step->xyz = Lig->new_mcoords[rand];
                            xcom =0;
                            ycom =0;
                            zcom =0;

                        }

                        else{
                            if (new_energy <= energy){
                                Lig->mcoords = Lig->new_mcoords;
                                this->xyz = step->xyz;
                                energy = new_energy;
                                rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                                xcom = xcom + step->dx;
                                ycom = ycom + step->dy;
                                zcom = zcom + step->dz;
                                gzprintf(mc_output_lig, "%10d %10.3f %10.3f  %10d %.3lf %.3lf %.3lf \n", count, energy, rmsd, rand, xcom , ycom, zcom);
                                if (Input->write_mol2){
                                    Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, Input->output + "_MC.ligsim");
                                }



#ifdef DEBUG

                                Writer->print_info(info);
#endif
                                sweepcount++;
                            }
                            else{
                                p = this->Boltzmman(energy, new_energy, Input->temp);
                                rnumber = gsl_rng_uniform(r) / (gsl_rng_max(r) + 1.0);
                                if (p > rnumber){

                                    Lig->mcoords = Lig->new_mcoords;
                                    this->xyz = step->xyz;
                                    energy = new_energy;
                                    rmsd = Coord->compute_rmsd(RefLig->xyz, step->xyz, Lig->N);
                                    xcom = xcom + step->dx;
                                    ycom = ycom + step->dy;
                                    zcom = zcom + step->dz;
                                    gzprintf(mc_output_lig, "%10d %10.3f %10.3f  %10d %.3lf %.3lf %.3lf \n", count, energy, rmsd, rand, xcom , ycom, zcom);
                                    if (Input->write_mol2){
                                        Writer->writeMol2(Lig, step->xyz, new_energy, rmsd, Input->output + "_MC.ligsim");
                                    }

                                    sweepcount ++;

                                }
                                else{
                                    nReject++;
                                }


                            }
                        }
                    }

#ifdef DEBUG
                    Writer->print_info(info);
#endif

                    sprintf(info, "%10d %10.3f %10.3f", count, energy, rmsd);
                    count++;
                }


                double sup = (double) Input->number_steps;
                double low = (double)((Input->number_steps+nReject));
                gzclose(mc_output_lig);
                delete step;
                delete Energy;
                delete Coord;
                Writer->print_line();
                sprintf(info, "Finished MC simulation for Ligand after %d accepted steps, %d total steps, acceptance rate %5.3f", Input->number_steps, (Input->number_steps+nReject),(sup/low));
                Writer->print_info(info);
                Writer->print_line();


            }

        }
*/
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

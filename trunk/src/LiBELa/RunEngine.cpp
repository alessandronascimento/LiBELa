/*
 * temp_scheme.cpp
 *
 *  Created on: 06/07/2010
 *      Author: Nascimento
 */


#include "RunEngine.h"

using namespace std;

TEMP_SCHEME::TEMP_SCHEME(int ac, char *av[]){

    this->argc = ac;
    this->argv = av;

    Input = new PARSER();
    Input->set_parameters(av[1]);

    Writer = new WRITER(Input);

    REC = new Mol2(Input, Input->rec_mol2);
    LIG = new Mol2(Input, Input->lig_mol2);
    RefLig = new Mol2(Input, Input->reflig_mol2);

//    RAND Rand;

    Ene = new Energy2(Input);
}

#ifdef HAS_GUI
TEMP_SCHEME::TEMP_SCHEME(PARSER* Input, QTextEdit* Editor){

    QWriter = new QtWriter(Input, Editor);

    // Reading Receptor files

    REC = new Mol2(Input, Input->rec_mol2);
    LIG = new Mol2(Input, Input->lig_mol2);
    RefLig = new Mol2(Input, Input->reflig_mol2);
    Ene = new Energy2(Input);
}
#endif

void TEMP_SCHEME::evaluation(){

    if ((REC->N != int(REC->charges.size())) or (REC->N != int(REC->radii.size())) or (REC->N != int(REC->epsilons.size()))){
        cout <<"The number of atomic parameters and the number of atoms doesn't match!" << endl;
        cout << "Exiting..." << endl;
        exit(1);
    }


    //! Evaluation of the original binding energy
    start_energy = Ene->compute_ene(REC, LIG, LIG->xyz);
    sprintf(info,"Original energy: %.4f kcal/mol", start_energy);
    Writer->print_info(info);
    if (Input->reflig_mol2 == ""){
        center = COORD.compute_com(LIG);
    }
    else {
        center = COORD.compute_com(RefLig);
    }
    Writer->write_box(center, center[0]-Input->x_dim/2, center[1]-Input->y_dim/2, center[2]-Input->z_dim/2, center[0]+Input->x_dim/2, center[1]+Input->y_dim/2, center[2]+Input->z_dim/2);

    sprintf(info, "Center of computation box: %.2f %.2f %.2f", center[0], center[1], center[2]);
    Writer->print_info(info);

    if (Input->use_grids){
        if (Input->load_grid_from_file){
            sprintf(info,"Loading grids from file %s.grid...", Input->grid_prefix.c_str());
            Writer->print_info(info);
            Grids = new Grid(Input);
            Grids->load_grids_from_file();
            sprintf(info, "Loaded energy grids with %d points spaced by %5.3f Angstroms in each directon.", Grids->npointsx*Grids->npointsy*Grids->npointsz, Grids->grid_spacing);
            Writer->print_info(info);
        }
        else{
            sprintf(info,"Generating energy grids. It can take a couple of minutes. Coffee time maybe ?");
            Writer->print_info(info);
            start = clock();
            Grids = new Grid(Input, REC, center);
            end = clock();
            sprintf(info, "Computed energy grids with %d points spaced by %5.3f Angstroms in each directon.", Grids->npointsx*Grids->npointsy*Grids->npointsz, Grids->grid_spacing);
            Writer->print_info(info);
            sprintf(info, "Grid computation took %d seconds.", int((end-start)/CLOCKS_PER_SEC));
            Writer->print_info(info);
        }
        double grid_energy = Ene->compute_ene(Grids, LIG, LIG->xyz);
        sprintf(info,"Original Grid energy: %.4f kcal/mol.", grid_energy);
        Writer->print_info(info);
        sprintf(info, "Energy error = %.4f%s.", fabs((grid_energy-start_energy)/start_energy)*100., "%");
        Writer->print_info(info);
    }

    Writer->print_line();

    this->dock_run();
    this->sa_run();
    this->eq_run();
    this->mcr_run();

    sprintf(info, "Finishing McLibela...");
    Writer->print_info(info);
    Writer->print_line();
    delete Writer;
}

#ifdef HAS_GUI

void TEMP_SCHEME::evaluation(PARSER* Input, QProgressBar* progressbar){

    QWriter->write_welcome();
    QWriter->write_params(Input);

    if ((REC->N != int(REC->charges.size())) or (REC->N != int(REC->radii.size())) or (REC->N != int(REC->epsilons.size()))){
        sprintf(info, "The number of atomic parameters and the number of atoms don't match! Exiting...");
        QWriter->print_info(info);
        exit(1);
    }

    //! Evaluation of the original binding energy
    start_energy = Ene->compute_ene(REC, LIG, LIG->xyz);
    sprintf(info,"Original energy: %.4f kcal/mol", start_energy);
    QWriter->print_info(info);

    center = COORD.compute_com(LIG);

    sprintf(info, "Center of computation box: %.2f %.2f %.2f", center[0], center[1], center[2]);
    QWriter->print_info(info);

    QWriter->print_line();

    if (Input->use_grids){
        if (Input->load_grid_from_file){
            sprintf(info,"Loading grids from file %s.grid...", Input->grid_prefix.c_str());
            QWriter->print_info(info);
            Grids = new Grid(Input);
            Grids->load_grids_from_file();
            sprintf(info, "Loaded energy grids with %d points spaced by %5.3f Angstroms in each directon.", Grids->npointsx*Grids->npointsy*Grids->npointsz, Grids->grid_spacing);
            QWriter->print_info(info);
        }
        else{
            sprintf(info,"Generating energy grids. It can take a couple of minutes. Coffee time maybe ?");
            QWriter->print_info(info);
            start = clock();
            Grids = new Grid(Input, REC, center);
            end = clock();
            sprintf(info, "Computed energy grids with %d points spaced by %5.3f Angstroms in each directon.", Grids->npointsx*Grids->npointsy*Grids->npointsz, Grids->grid_spacing);
            QWriter->print_info(info);
            sprintf(info, "Grid computation took %d seconds.", int((end-start)/CLOCKS_PER_SEC));
            QWriter->print_info(info);
        }
        double grid_energy = Ene->compute_ene(Grids, LIG, LIG->xyz);
        sprintf(info,"Original Grid energy: %.4f kcal/mol.", grid_energy);
        QWriter->print_info(info);
        sprintf(info, "Energy error = %.4f%s.", fabs((grid_energy-start_energy)/start_energy)*100., "%");
        QWriter->print_info(info);
    }

    QWriter->print_line();

    this->run_dock_gui(Input, progressbar);
    //    this->sa_run();
    //    this->eq_run();

    sprintf(info, "Finishing McLibela...");
    QWriter->print_info(info);
    QWriter->print_line();

}
#endif


void TEMP_SCHEME::sa_run(){

    if (Input->sa_mode == true){
        sprintf(info,"%s", "Entering SA Scheme....");
        Writer->print_info(info);
        srand(rand());
        gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs2);
        int seed = (rand() % 101);
        sprintf(info, "Seed: %d", seed);
        Writer->print_info(info);
        gsl_rng_set(r, seed);

        //! rmsd_energy is defined as an output file to keep the trial number, the rmsd and the evaluated energy.
        sprintf(rmsd_ene, "%s_sa.dat", Input->output.c_str());
        rmsd_energy = fopen(rmsd_ene, "w" );




        if ((Input->flex_lig == true)){
            Rand.random(Input->cushion, Input->rotation_step, REC, LIG);
        }
        else {
            Rand.random(Input->cushion, Input->rotation_step);
        }

        Rand.print();

        old_coord = COORD.rototranslate(LIG->xyz, LIG, &Rand);

        if (Input->use_grids){
            start_energy = Ene->compute_ene(Grids, LIG, old_coord);
        }
        else {
            start_energy = Ene->compute_ene(REC, LIG, old_coord);
        }

        sprintf(info, "Start energy: %.4f kcal/mol", start_energy);
        Writer->print_info(info);

        rmsd = COORD.compute_rmsd(LIG->xyz, old_coord, LIG->N);
        fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", -1, rmsd, start_energy, 0.00);

        if( (Input->flex_lig == true)){

            for(double t=Input->sa_start_temp; t>= Input->temp; t-=10*log(t)){

                sprintf(info, "SA Temperature: %.2f", t);
                Writer->print_info(info);

                for (int i=1; i<=Input->sa_steps; i++){
                    Rand.random(Input->cushion, Input->rotation_step, REC, LIG);
                    COORD.rototranslate_all(LIG, &Rand);

                    com = COORD.compute_com(LIG->new_mcoords[Rand.lign], LIG);

                    if ( (com[0] >= center[0]-Input->x_dim/2) and (com[0] <= center[0]+Input->x_dim/2) and (com[1] >= center[1]-Input->y_dim/2) and
                         (com[1] <= center[1]+Input->y_dim/2) and (com[2] >= center[2]-Input->z_dim/2) and (com[2] <= center[2]+Input->z_dim/2) ) {

                        if (Input->use_grids){
                            new_energy = Ene->compute_ene(Grids, LIG,LIG->new_mcoords[Rand.lign]);
                        }
                        else {
                            new_energy = Ene->compute_ene(REC, LIG,LIG->new_mcoords[Rand.lign]);
                        }

                        rmsd = COORD.compute_rmsd(LIG->xyz, LIG->new_mcoords[Rand.lign],LIG->N);

                        if (new_energy <= start_energy){
                            sprintf(info, "Trial: %5d New energy: % 7.2f RMSD: %6.2f Temp: %6.2f Conf: %4d", i, new_energy, rmsd, t, Rand.lign);
                            Writer->print_info(info);
                            fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", i, rmsd, new_energy, t);
                            start_energy=new_energy;
                            LIG->mcoords=LIG->new_mcoords;

                            Writer->writeMol2(LIG, LIG->mcoords[Rand.lign], new_energy, rmsd, Input->output);
                        }

                        else { //energy isn't <= start_energy

                            prob = COORD.compute_prob(start_energy, new_energy, Input->temp);
                            rnumber = gsl_rng_uniform(r);

                            if (prob > rnumber ) {
                                fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", i, rmsd, new_energy, t);
                                sprintf(info, "Trial: %5d New energy: % 7.2f RMSD: %6.2f Temp: %6.2f Conf: %4d", i, new_energy, rmsd, t, Rand.lign);
                                Writer->print_info(info);
                                start_energy=new_energy;
                                LIG->mcoords = LIG->new_mcoords;
                                Writer->writeMol2(LIG, LIG->mcoords[Rand.lign], new_energy, rmsd, Input->output);
                            }
                        }
                    }
                }
            }
        }

        else {

            for(double t=Input->sa_start_temp; t>= Input->temp; t-=10*log(t)){
                sprintf(info, "SA Temperature: %.2f", t);
                Writer->print_info(info);
                for (int i=1; i<= Input->sa_steps; i++){
                    Rand.random(Input->cushion, Input->rotation_step, REC, LIG);
                    new_coord = COORD.rototranslate(old_coord, LIG, &Rand);

                    com = COORD.compute_com(new_coord, LIG);

                    //                    printf("COM: %.2f %.2f %.2f\n", com[0], com[1], com[2]);

                    if ( (com[0] >= center[0]-Input->x_dim/2) and (com[0] <= center[0]+Input->x_dim/2) and (com[1] >= center[1]-Input->y_dim/2) and
                         (com[1] <= center[1]+Input->y_dim/2) and (com[2] >= center[2]-Input->z_dim/2) and (com[2] <= center[2]+Input->z_dim/2) ) {

                        if (Input->use_grids){
                            new_energy = Ene->compute_ene(Grids, LIG, new_coord);
                        }
                        else {
                            new_energy = Ene->compute_ene(REC, LIG, new_coord);
                        }

                        //                        cout << "Energy: " << new_energy << endl;

                        rmsd = COORD.compute_rmsd(LIG->xyz, new_coord, LIG->N);

                        if (new_energy <= start_energy){
                            fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", i, rmsd, new_energy, t);
                            sprintf(info, "Trial: %5d New energy: % 9.4f RMSD: %7.3f Temp: %7.3f", i, new_energy, rmsd, t);
                            Writer->print_info(info);
                            start_energy=new_energy;
                            old_coord=new_coord;
                            Writer->writeMol2(LIG, new_coord, new_energy, rmsd, Input->output);
                        }

                        else {
                            prob = COORD.compute_prob(start_energy, new_energy, t);
                            rnumber = gsl_rng_uniform(r);
                            if (prob > rnumber ) {
                                fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", i, rmsd, new_energy, t);
                                sprintf(info, "Trial: %5d New energy: % 9.4f RMSD: %7.3f Temp: %7.3f", i, new_energy, rmsd, t);
                                Writer->print_info(info);
                                start_energy=new_energy;
                                old_coord=new_coord;
                                Writer->writeMol2(LIG, new_coord, new_energy, rmsd, Input->output);
                            }
                        }
                    }
                }
            }
        }
        fclose(rmsd_energy);
        Writer->print_line();
    }
}

#ifdef HAS_GUI

void TEMP_SCHEME::sa_run(PARSER* Input){
    if (Input->sa_mode == true){
        sprintf(info,"%s", "Entering SA Scheme....");
        QWriter->print_info(info);
        sprintf(rmsd_ene, "%s_sa.dat", Input->output.c_str());
        rmsd_energy = fopen(rmsd_ene, "w" );

        srand(rand());
        gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs2);
        int seed = (rand() % 101);
        sprintf(info, "Seed: %d", seed);
        QWriter->print_info(info);
        gsl_rng_set(r, seed);


        if ((Input->flex_lig == true)){
            Rand.random(Input->cushion, Input->rotation_step, REC, LIG);
        }
        else {
            Rand.random(Input->cushion, Input->rotation_step);
        }

        Rand.print();
        old_coord = COORD.rototranslate(LIG->xyz, LIG, &Rand);

        if (Input->use_grids){
            start_energy = Ene->compute_ene(Grids, LIG, old_coord);
        }
        else {
            start_energy = Ene->compute_ene(REC, LIG, old_coord);
        }

        sprintf(info, "Start energy: %.4f kcal/mol", start_energy);
        QWriter->print_info(info);

        rmsd = COORD.compute_rmsd(LIG->xyz, old_coord, LIG->N);
        fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", -1, rmsd, start_energy, 0.00);

        if( (Input->flex_lig == true)){

            for(double t=Input->sa_start_temp; t>= Input->temp; t-=10*log(t)){

                sprintf(info, "SA Temperature: %.2f", t);
                QWriter->print_info(info);

                for (int i=1; i<=Input->sa_steps; i++){
                    Rand.random(Input->cushion, Input->rotation_step, REC, LIG);
                    COORD.rototranslate_all(LIG, &Rand);

                    com = COORD.compute_com(LIG->new_mcoords[Rand.lign], LIG);

                    if ( (com[0] >= center[0]-Input->x_dim/2) and (com[0] <= center[0]+Input->x_dim/2) and (com[1] >= center[1]-Input->y_dim/2) and
                         (com[1] <= center[1]+Input->y_dim/2) and (com[2] >= center[2]-Input->z_dim/2) and (com[2] <= center[2]+Input->z_dim/2) ) {

                        if (Input->use_grids){
                            new_energy = Ene->compute_ene(Grids, LIG,LIG->new_mcoords[Rand.lign]);
                        }
                        else {
                            new_energy = Ene->compute_ene(REC, LIG,LIG->new_mcoords[Rand.lign]);
                        }

                        rmsd = COORD.compute_rmsd(LIG->xyz, LIG->new_mcoords[Rand.lign],LIG->N);

                        if (new_energy <= start_energy){
                            sprintf(info, "Trial: %5d New energy: % 7.2f RMSD: %6.2f Temp: %6.2f Conf: %4d", i, new_energy, rmsd, t, Rand.lign);
                            QWriter->print_info(info);
                            fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", i, rmsd, new_energy, t);
                            start_energy=new_energy;
                            LIG->mcoords=LIG->new_mcoords;

                            QWriter->writeMol2(LIG, LIG->mcoords[Rand.lign], new_energy, rmsd, Input->output);
                            if (Input->flex_rec ==true){
                                QWriter->writeMol2(REC, REC->mcoords[Rand.recn], 0.00, 0.00, Input->output+"_rec");
                            }
                        }

                        else { //energy isn't <= start_energy

                            prob = COORD.compute_prob(start_energy, new_energy, Input->temp);
                            rnumber = gsl_rng_uniform(r);

                            if (prob > rnumber ) {
                                fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", i, rmsd, new_energy, t);
                                sprintf(info, "Trial: %5d New energy: % 7.2f RMSD: %6.2f Temp: %6.2f Conf: %4d", i, new_energy, rmsd, t, Rand.lign);
                                QWriter->print_info(info);
                                start_energy=new_energy;
                                LIG->mcoords = LIG->new_mcoords;
                                QWriter->writeMol2(LIG, LIG->mcoords[Rand.lign], new_energy, rmsd, Input->output);
                                if (Input->flex_rec ==true){
                                    QWriter->writeMol2(REC, REC->mcoords[Rand.recn], 0.00, 0.00, Input->output+"_rec");
                                }
                            }
                        }
                    }
                }
            }
        }

        else {

            for(double t=Input->sa_start_temp; t>= Input->temp; t-=10*log(t)){
                sprintf(info, "SA Temperature: %.2f", t);
                QWriter->print_info(info);
                for (int i=1; i<= Input->sa_steps; i++){
                    Rand.random(Input->cushion, Input->rotation_step, REC, LIG);
                    new_coord = COORD.rototranslate(old_coord, LIG, &Rand);

                    com = COORD.compute_com(new_coord, LIG);

                    if ( (com[0] >= center[0]-Input->x_dim/2) and (com[0] <= center[0]+Input->x_dim/2) and (com[1] >= center[1]-Input->y_dim/2) and
                         (com[1] <= center[1]+Input->y_dim/2) and (com[2] >= center[2]-Input->z_dim/2) and (com[2] <= center[2]+Input->z_dim/2) ) {

                        if (Input->use_grids){
                            new_energy = Ene->compute_ene(Grids, LIG, new_coord);
                        }
                        else{
                            new_energy = Ene->compute_ene(REC, LIG, new_coord);
                        }

                        rmsd = COORD.compute_rmsd(LIG->xyz, new_coord, LIG->N);

                        if (new_energy <= start_energy){
                            fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", i, rmsd, new_energy, t);
                            sprintf(info, "Trial: %5d New energy: % 9.4f RMSD: %7.3f Temp: %7.3f", i, new_energy, rmsd, t);
                            QWriter->print_info(info);
                            start_energy=new_energy;
                            old_coord=new_coord;
                            QWriter->writeMol2(LIG, new_coord, new_energy, rmsd, Input->output);
                            if (Input->flex_rec ==true){
                                QWriter->writeMol2(REC, REC->mcoords[Rand.recn], 0.00, 0.00, Input->output + "_rec");
                            }
                        }

                        else {
                            prob = COORD.compute_prob(start_energy, new_energy, t);
                            rnumber = gsl_rng_uniform(r);
                            if (prob > rnumber ) {
                                fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", i, rmsd, new_energy, t);
                                sprintf(info, "Trial: %5d New energy: % 9.4f RMSD: %7.3f Temp: %7.3f", i, new_energy, rmsd, t);
                                QWriter->print_info(info);
                                start_energy=new_energy;
                                old_coord=new_coord;
                                QWriter->writeMol2(LIG, new_coord, new_energy, rmsd, Input->output);
                                if (Input->flex_rec ==true){
                                    QWriter->writeMol2(REC, REC->mcoords[Rand.recn], 0.00, 0.00, Input->output + "_rec");
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(rmsd_energy);
    QWriter->print_line();
}

#endif

void TEMP_SCHEME::dock_run(){
    clock_t ti, tf;
    if(Input->dock_mode == true){

        ti = clock();

        sprintf(info, "%5s %-12.12s %-4.4s %-10.10s %-8.8s %-8.8s %-8.8s %-8.8s %-8.8s %3.3s %-5.5s %2.2s", "#", "Ligand", "Resi", "Overlay", "Elec", "VDW", "Rec_Solv", "Lig_Solv", "Energy", "Conf", "MS", "SI");
        Writer->print_info(info);
        Writer->print_line();


        center = COORD.compute_com(RefLig);

        unsigned counter=0;
        if (Input->multifile != ""){
            if (!Input->dock_parallel){
                string ligand = "";
                ifstream multifile(Input->multifile.c_str());
                if (! multifile.is_open()){
                    sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
                    Writer->print_info(info);
                    exit(1);
                }
                multifile >> ligand;
                while ((!multifile.eof()) and (ligand != "EOF")){
                    Mol2* Lig2 = new Mol2(Input, ligand);
                    counter++;
                    if (Input->generate_conformers){
                        Conformer* Conf = new Conformer;
                        if (Input->conformer_generator == "GA"){
                            if (Conf->generate_conformers_GA(Input, Lig2, ligand)){
                                if (Input->use_grids){
                                    Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, Grids, counter);
                                    delete Dock;
                                }
                                else {
                                    Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, counter);
                                    delete Dock;
                                }
                            }
                            else{ 																//test: in case generate_conformers does not succeed.
                                Lig2->mcoords.push_back(Lig2->xyz);
                                if (Input->use_grids){
                                    Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, Grids, counter);
                                    delete Dock;
                                }
                                else {
                                    Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, counter);
                                    delete Dock;
                                }

                            }
                        }
                        else {
                            if(Conf->generate_conformers_WRS(Input, Lig2, ligand)){
                                if (Input->use_grids){
                                    Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, Grids, counter);
                                    delete Dock;
                                }
                                else {
                                    Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, counter);
                                    delete Dock;
                                }
                            }
                        }
                        Conf->~Conformer();
                        delete Conf;
                    }
                    else {
                        if (Input->use_grids){
                            Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, Grids, counter);
                            delete Dock;
                        }
                        else {
                            Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, counter);
                            delete Dock;
                        }

                    }
                    delete Lig2;
                    multifile >> ligand;
                }
                multifile.close();
            }
            else {
                this->dock_parallel();
            }
        }
        else {
            counter=1; 														// just one ligand
            if (Input->generate_conformers){
                Conformer* Conf = new Conformer;
                if (Input->conformer_generator == "GA"){
                    Conf->generate_conformers_GA(Input, LIG, Input->lig_mol2);
                    if (Input->use_grids){
                        Docker* Dock = new Docker(REC, LIG, RefLig, center, Input, Writer, Grids, counter);
                        delete Dock;
                    }
                    else {
                        Docker* Dock = new Docker(REC, LIG, RefLig, center, Input, Writer, counter);
                        delete Dock;
                    }
                }
                else {
                    Conf->generate_conformers_WRS(Input, LIG, Input->lig_mol2);
                    if (Input->use_grids){
                        Docker* Dock = new Docker(REC, LIG, RefLig, center, Input, Writer, Grids, counter);
                        delete Dock;
                    }
                    else {
                        Docker* Dock = new Docker(REC, LIG, RefLig, center, Input, Writer, counter);
                        delete Dock;
                    }
                }
                delete Conf;
            }
            else {
                if (Input->use_grids){
                    Docker* Dock = new Docker(REC, LIG, RefLig, center, Input, Writer, Grids, counter);
                    delete Dock;
                }
                else {
                    Docker* Dock = new Docker(REC, LIG, RefLig, center, Input, Writer, counter);
                    delete Dock;
                }
            }
        }

        Writer->print_line();
        sprintf(info, "Min Status:");
        Writer->print_info(info);
        sprintf(info, "    Failure = -1, Out of memory = -3, Roundoff limited = -4, Forced stop = -5,");
        Writer->print_info(info);
        sprintf(info, "    Stopval reached = 2, Ftol reached = 3, Xtol reached = 4, Maxeval reached=5, Maxtime reached=6");
        Writer->print_info(info);
        Writer->print_line();
        tf = clock()-ti;
        sprintf(info, " Docking computations took %d minute(s)", int((tf/CLOCKS_PER_SEC)/60));
        Writer->print_info(info);
        Writer->print_line();
    }
}

int TEMP_SCHEME::dock_parallel_function(Mol2* Lig2){
    if (Input->use_grids){
        Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, Grids, 0);
        delete Dock;
    }
    else {
        Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, 0);
        delete Dock;
    }

    return 0;
}

int TEMP_SCHEME::dock_serial(vector<string> ligand_list, int count, int chunck_size){
    stringstream buffer;
    buffer << (count+1);
    WRITER* Write_lig = new WRITER(Input->output + "_" + buffer.str(), Input);

    for (unsigned i=0; i< ligand_list.size(); i++){
        Mol2* Lig2 = new Mol2;
        bool lig_is_opened = false;

        if (ligand_list[i].substr(ligand_list[i].size()-3, 3) == ".gz"){
            lig_is_opened = Lig2->parse_gzipped_file(Input, ligand_list[i]);
        }
        else {
            lig_is_opened = Lig2->parse_mol2file(Input, ligand_list[i]);
        }

        if (lig_is_opened){
            if (Input->generate_conformers){
                Conformer* Conf = new Conformer;
                if (Input->conformer_generator == "GA"){
                    Conf->generate_conformers_GA(Input, Lig2, ligand_list[i]);
                }
                else {
                    Conf->generate_conformers_WRS(Input, Lig2, ligand_list[i]);
                }
                delete Conf;
            }
            if (Input->use_grids){
                Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Write_lig, Grids, ((count*chunck_size))+i+1);
                delete Dock;
            }
            else {
                Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Write_lig, ((count*chunck_size))+i+1);
                delete Dock;
            }
        }
        delete Lig2;
    }
    delete Write_lig;
    return 0;
}

int TEMP_SCHEME::dock_concurrent_function(string mylig){
    Mol2* Lig2 = new Mol2(Input, mylig);

    if (Input->generate_conformers and Input->conformer_generator == "GA"){
        Conformer* Conf = new Conformer;
#pragma omp critical
        {
            Conf->generate_conformers_GA(Input, Lig2, mylig);
        }
        delete Conf;
    }

    if (Input->use_grids){
        Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, Grids, 0);
        delete Dock;
    }
    else {
        Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, 0);
        delete Dock;
    }

    return 0;
}

void TEMP_SCHEME::dock_parallel(){

#ifdef HAS_MPI

    this->dock_mpi();

#else

    sprintf(info, "Running Dock mode in parallel with up to %d threads...", int(Input->parallel_jobs));
    Writer->print_info(info);
    vector<string> ligand_list;


/* Here, we read the multifile file and parse the files of the molecules
 * to be docked into a std::vector named ligand_list.
 */

    if (Input->multifile != ""){
        ifstream multifile(Input->multifile.c_str());
        if (! multifile.is_open()){
            sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
            Writer->print_info(info);
            exit(1);
        }
        string ligand;
        multifile >> ligand;

        while ((!multifile.eof()) and (ligand != "EOF")){
            ligand_list.push_back(ligand);
            multifile >> ligand;
        }

        multifile.close();

/*
 * Here the parallel processing begins. Firstly, a parallel section is created with a pragma
 * indicating the number of parallel threads as defined in the input file.
 * After, the molecule objects are created and conformers are calculated using OpenBabel.
 * Due to OpenBabel thread unsafety issues, the conformer generation is done in serial, using
 * the sync directive #pragma critical. Finally, the docking objects are created and the molecules
 * are actually docked.
 */

#pragma omp parallel num_threads(Input->parallel_jobs)
        {
#pragma omp for schedule(static,1)
            for (int i=0; i< int(ligand_list.size()); i++){
                Mol2* Lig2 = new Mol2;
                bool lig_is_opened = false;
#pragma omp critical
                {
                    if (ligand_list[i].substr(ligand_list[i].size()-3, 3) == ".gz"){
                        lig_is_opened = Lig2->parse_gzipped_file(Input, ligand_list[i]);
                    }
                    else {
                        lig_is_opened = Lig2->parse_mol2file(Input, ligand_list[i]);
                    }
                }
                if (lig_is_opened){
                    if (Input->generate_conformers){
                        Conformer* Conf = new Conformer;
#pragma omp critical
                        {
                            if (Input->conformer_generator == "GA"){
                                Conf->generate_conformers_GA(Input, Lig2, ligand_list[i]);
                            }
                            else {
                                Conf->generate_conformers_WRS(Input, Lig2, ligand_list[i]);
                            }
                        }
                        delete Conf;
                        if (Input->use_grids){
                            Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, Grids, i+1);
                            delete Dock;
                        }
                        else {
                            Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, Writer, i+1);
                            delete Dock;
                        }
                    }
                }
                delete Lig2;
            }
        }
        ligand_list.clear();
    }
#endif
}
#ifdef HAS_GUI
void TEMP_SCHEME::dock_concurrent(){
    sprintf(info, "Running Dock mode in parallel with up to %d threads...", int(Input->parallel_jobs));
    Writer->print_info(info);
    QVector<string> ligand_list;

    /* Here we read the multifile file and parse the files of the molecules
 * to be docked into a std::vector named ligand_list.
 */

    if (Input->multifile != ""){
        ifstream multifile(Input->multifile.c_str());
        if (! multifile.is_open()){
            sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
            Writer->print_info(info);
            exit(1);
        }
        string ligand;
        multifile >> ligand;
        while ((!multifile.eof()) and (ligand != "EOF")){
            ligand_list.push_back(ligand);
            multifile >> ligand;
        }
    }
    QtConcurrent::blockingMapped(ligand_list, std::bind1st(std::mem_fun(&TEMP_SCHEME::dock_concurrent_function), this));
}
#endif

#ifdef HAS_MPI
void TEMP_SCHEME::dock_mpi(){

    int chunck_size;
    vector<string> ligand_list;

    mpi::environment env(this->argc, this->argv);
    mpi::communicator world;



    vector<string> tmp;

    if (world.rank() == 0 ){                       // Running only on the *master* job;
        sprintf(info, "Running Dock mode over MPI with up to %d parallel jobs...", world.size());
        Writer->print_info(info);

        if (Input->multifile == ""){
            sprintf(info, "Cannot run single molecule docking in parallel");
            Writer->print_info(info);
            exit(1);
        }

        ifstream multifile(Input->multifile.c_str());
        if (! multifile.is_open()){
            sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
            Writer->print_info(info);
            exit(1);
        }
        string ligand;
        multifile >> ligand;
        while ((!multifile.eof()) and (ligand != "EOF")){
            ligand_list.push_back(ligand);
            multifile >> ligand;
        }
        multifile.close();

        // Multifile read!

        /*
         * Here we ought to split the files in the multifile into N chuncks with approximately the same size. N
         * is determined by MPI (mpirun -np N) during program execution. Each chunk will give rise to a different
         * processReflig
        */

        chunck_size = int(ligand_list.size()/world.size());

        vector<vector<string> > chuncks;

        for (int i=0; i<world.size(); i++){
            tmp.clear();
            for (int j=0; j<chunck_size; j++){
                ligand = ligand_list.front();
                tmp.push_back(ligand);
                ligand_list.erase(ligand_list.begin());
            }
            chuncks.push_back(tmp);
        }
        if (ligand_list.size() != 0){
            for (unsigned i=0; i<ligand_list.size(); i++){
                chuncks[i].push_back(ligand_list[i]);
            }
            ligand_list.clear();
        }

        tmp.clear();

        scatter(world,chuncks, tmp, 0);
        this->dock_serial(tmp, world.rank(), chunck_size);
    }                                                           // End of rank 0;
    else {
        scatter(world, tmp, 0);
        this->dock_serial(tmp, world.rank(), tmp.size());
    }
}

#endif

void TEMP_SCHEME::eq_run(){
    if (Input->eq_mode){
        MC* EqMC = new MC(LIG, Input, Writer);
        if (Input->generate_conformers){
            if (Input->conformer_generator == "GA"){
                Conformer* Conf = new Conformer;
                Conf->generate_conformers_GA(Input, LIG, Input->lig_mol2);
                delete Conf;
            }
            else {
                Conformer* Conf = new Conformer;
                Conf->generate_conformers_WRS(Input, LIG, Input->lig_mol2);
                delete Conf;
            }
        }

        if (Input->use_grids){
            EqMC->run(Grids, RefLig , LIG, LIG->xyz, Input, Input->temp);
        }
        else {
            EqMC->run(REC, RefLig , LIG, LIG->xyz, Input, Input->temp);
        }

        if (Input->ligsim){
            EqMC->ligand_run(RefLig, LIG, LIG->xyz, Input, Input->temp);
        }

        delete EqMC;
    }
}

void TEMP_SCHEME::mcr_run(){
    if (Input->mcr_mode){

        // checking consistency
        if (int(Input->mcr_coefficients.size()) != Input->mcr_size){
            sprintf(info, "Inconsistency found in MCR coefficients. Please check!\n");
            Writer->print_info(info);
            exit(1);
        }

        MC* EqMC = new MC(LIG, Input, Writer);
        if (Input->generate_conformers){
            if (Input->conformer_generator == "GA"){
                Conformer* Conf = new Conformer;
                Conf->generate_conformers_GA(Input, LIG, Input->lig_mol2);
                delete Conf;
            }
            else {
                Conformer* Conf = new Conformer;
                Conf->generate_conformers_WRS(Input, LIG, Input->lig_mol2);
                delete Conf;
            }
        }

        sprintf(info, "MCR %7.7s %7.7s %7.7s %7.7s %7.7s %7.7s %7.7s %7.7s %7.7s %7.7s",  "#i", "bi", "bT", "<ene>" , "SD(ene)", "e-(b-1)U/kT)" , "SD(exp)" ,"ln(W)", "E(ln_W)" ,"Vol(A3)");
        Writer->print_info(info);

        double bt;                      // MC Recursion "effective" temperature (bt) fot ith evaluation;
        double k = 0.0019858775203792202;

        double cum_W = 0.0;
        double cum_W_err = 0.0;
        double max_vol=0.0;
        double volume;

        for (int i=0; i<Input->mcr_size; i++){
            Input->bi = Input->mcr_coefficients[i];
            bt = Input->temp;
            for (int j=0; j < i; j++){
                bt = bt*Input->mcr_coefficients[j];
            }

            if (Input->use_grids){
                EqMC->run(Grids, RefLig , LIG, LIG->xyz, Input, bt);
            }
            else {
                EqMC->run(REC, RefLig , LIG, LIG->xyz, Input, bt);
            }

            volume = (EqMC->XSize*EqMC->YSize*EqMC->ZSize);

            sprintf(info, "MCR %7d %7.4f %7.4g %7.4g %7.4g %7.4Lg %7.4Lg %7.4g %7.4Lg %7.4g", i+1, Input->mcr_coefficients[i], bt, EqMC->average_energy, EqMC->energy_standard_deviation, EqMC->MCR_Boltzmann_weighted_average,
                    EqMC->MCR_Boltzmann_weighted_stdev, log(double(EqMC->MCR_Boltzmann_weighted_average)), ((1.0/double(EqMC->MCR_Boltzmann_weighted_average))*EqMC->MCR_Boltzmann_weighted_stdev), volume);
            Writer->print_info(info);
            Writer->print_line();

            cum_W += (log(double(EqMC->MCR_Boltzmann_weighted_average)));
            cum_W_err += double((1.0/double(EqMC->MCR_Boltzmann_weighted_average))*EqMC->MCR_Boltzmann_weighted_stdev);
            if (volume > max_vol){
                max_vol=volume;
            }
        }

        Writer->print_line();
        sprintf(info, "MCR: SUM of { ln [ <e^([-(b-1)]*U/kT)> ] } = %10.4g +/- %10.4g",  cum_W, cum_W_err);
        Writer->print_info(info);

        sprintf(info, "MCR: Volume: %10.4g.  ln(volume) = %10.4g",  volume, log(volume));
        Writer->print_info(info);

        Writer->print_line();


        double cum_W_lig = 0.0;
        double cum_W_lig_err = 0.0;
        double lig_max_vol = 0.0;
        double lig_volume = 0.0;


        bt = Input->temp;

        if (Input->ligsim){
            for (int i=0; i<Input->mcr_size; i++){
                Input->bi = Input->mcr_coefficients[i];
                bt = Input->temp;
                for (int j=0; j < i; j++){
                    bt = bt*Input->mcr_coefficients[j];
                }
                EqMC->ligand_run(RefLig, LIG, LIG->xyz, Input, bt);
                lig_volume = (EqMC->XSize*EqMC->YSize*EqMC->ZSize);
                sprintf(info, "MCR %7d %7.4f %7.4g %7.4g %7.4g %7.4Lg %7.4Lg %7.4g %7.4Lg %7.4g", i+1, Input->mcr_coefficients[i], bt, EqMC->average_energy, EqMC->energy_standard_deviation, EqMC->MCR_Boltzmann_weighted_average,
                        EqMC->MCR_Boltzmann_weighted_stdev, log(double(EqMC->MCR_Boltzmann_weighted_average)), ((1.0/double(EqMC->MCR_Boltzmann_weighted_average))*EqMC->MCR_Boltzmann_weighted_stdev), lig_volume);
                Writer->print_info(info);

                cum_W_lig += (log(double(EqMC->MCR_Boltzmann_weighted_average)));
                cum_W_lig_err += double((1.0/double(EqMC->MCR_Boltzmann_weighted_average))*EqMC->MCR_Boltzmann_weighted_stdev);
                if (volume > max_vol){
                    lig_max_vol=lig_volume;
                }
            }

            Writer->print_line();
            sprintf(info, "MCR: SUM of { ln [ <e^([-(b-1)]*U/kT)> ] } for the ligand = %10.4g +/- %10.4g",  cum_W_lig, cum_W_lig_err);
            Writer->print_info(info);

            sprintf(info, "MCR: Ligand Volume: %10.4g.  ln(volume) = %10.4g",  lig_volume, log(lig_volume));
            Writer->print_info(info);
            Writer->print_line();

            double complex_A = (-k*Input->temp*log(volume))-(k*Input->temp*cum_W);
            double ligand_A = (-k*Input->temp*log(lig_volume))-(k*Input->temp*cum_W_lig);
            double Delta_A = complex_A - ligand_A;

            sprintf(info, "MCR: Complex Free Energy = %10.4g +/- %10.4g", complex_A, (k*Input->temp*cum_W_err));
            Writer->print_info(info);
            sprintf(info, "MCR: Ligand Free Energy = %10.4g +/- %10.4g", ligand_A, (k*Input->temp*cum_W_lig_err));
            Writer->print_info(info);
            sprintf(info, "MCR: Binding Free Energy = %10.4g +/- %10.4g", Delta_A, (k*Input->temp*(cum_W_err+cum_W_lig_err)));
            Writer->print_info(info);
            Writer->print_line();
        }

        delete EqMC;
    }
}

#ifdef HAS_GUI
void TEMP_SCHEME::run_dock_gui(PARSER* Input, QProgressBar* progressbar){
    clock_t ti, tf;

    if(Input->dock_mode == true){

        ti = clock();

        sprintf(info, "%5s %-12.12s %-4.4s %-10.10s %-8.8s %-8.8s %-8.8s %-8.8s %-8.8s %3.3s %-5.5s %2.2s", "#", "Ligand", "Resi", "Overlay", "Elec", "VDW", "Rec_Solv", "Lig_Solv", "Energy", "Conf", "MS", "SI");
        QWriter->print_info(info);
        QWriter->print_line();

        center = COORD.compute_com(RefLig);

        unsigned counter=0;
        if (!Input->dock_parallel){
            string ligand;
            for (unsigned i=0; i< Input->docking_molecules.size(); i++){
                ligand = Input->docking_molecules[i].toStdString();
                Mol2* Lig2 = new Mol2(Input, ligand);
                counter++;
                if (Input->generate_conformers){                             // has conformers
                    Conformer* Conf = new Conformer;
                    if (Input->conformer_generator == "GA"){
                        if (Conf->generate_conformers_GA(Input, Lig2, ligand)){
                            if (Input->use_grids){
                                Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, Grids, i);
                                delete Dock;
                            }
                            else {
                                Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, i);
                                delete Dock;
                            }
                        }
                        else{ 																//test: in case generate_conformers does not succeed.
                            Lig2->mcoords.push_back(Lig2->xyz);
                            if (Input->use_grids){
                                Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, Grids, i);
                                delete Dock;
                            }
                            else {
                                Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, i);
                                delete Dock;
                            }

                        }
                    }
                    else {
                        if(Conf->generate_conformers_WRS(Input, Lig2, ligand)){
                            if (Input->use_grids){
                                Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, Grids, i);
                                delete Dock;
                            }
                            else {
                                Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, i);
                                delete Dock;
                            }
                        }
                    }
                    Conf->~Conformer();
                    delete Conf;
                }
                else {                                            // single conformation
                    if (Input->use_grids){
                        Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, Grids, i);
                        delete Dock;
                    }
                    else {
                        Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, i);
                        delete Dock;
                    }

                }
                delete Lig2;
                progressbar->setValue(int(((i+1)*100)/Input->docking_molecules.size()));
            }
        }
        else {
            this->dock_parallel_gui(Input, progressbar);
        }


        QWriter->print_line();
        sprintf(info, "Min Status:");
        QWriter->print_info(info);
        sprintf(info, "    Failure = -1, Out of memory = -3, Roundoff limited = -4, Forced stop = -5,");
        QWriter->print_info(info);
        sprintf(info, "    Stopval reached = 2, Ftol reached = 3, Xtol reached = 4, Maxeval reached=5, Maxtime reached=6");
        QWriter->print_info(info);
        QWriter->print_line();
        tf = clock()-ti;
        sprintf(info, " Docking computations took %d minute(s)", int((tf/CLOCKS_PER_SEC)/60));
        QWriter->print_info(info);
        QWriter->print_line();
    }
}

void TEMP_SCHEME::dock_parallel_gui(PARSER* Input, QProgressBar *progressbar){
    sprintf(info, "Running Dock mode in parallel with up to %d threads...", int(Input->parallel_jobs));
    QWriter->print_info(info);

    int count=0;


#pragma omp parallel num_threads(Input->parallel_jobs)
    {
#pragma omp for schedule(static,1)
        for (unsigned i=0; i< Input->docking_molecules.size(); i++){
            Mol2* Lig2 = new Mol2;
            bool lig_is_opened = false;
#pragma omp critical
            {
                if (Input->docking_molecules[i].toStdString().substr(Input->docking_molecules[i].toStdString().size()-3, 3) == ".gz"){
                    lig_is_opened = Lig2->parse_gzipped_file(Input, Input->docking_molecules[i].toStdString());
                }
                else {
                    lig_is_opened = Lig2->parse_mol2file(Input, Input->docking_molecules[i].toStdString());
                }
            }
            if (lig_is_opened){
                if (Input->generate_conformers){
                    Conformer* Conf = new Conformer;
#pragma omp critical
                    {
                        if (Input->conformer_generator == "GA"){
                            Conf->generate_conformers_GA(Input, Lig2, Input->docking_molecules[i].toStdString());
                        }
                        else {
                            Conf->generate_conformers_WRS(Input, Lig2, Input->docking_molecules[i].toStdString());
                        }
                    }
                    delete Conf;
                    if (Input->use_grids){
                        Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, Grids, i+1);
                        delete Dock;
                    }
                    else {
                        Docker* Dock = new Docker(REC, Lig2, RefLig, center, Input, QWriter, i+1);
                        delete Dock;
                    }
                    progressbar->setValue(int(((count+1)*100*Input->parallel_jobs)/Input->docking_molecules.size()));
                }
            }
            delete Lig2;
        }
    }
}

#endif

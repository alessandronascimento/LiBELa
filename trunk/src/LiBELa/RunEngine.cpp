/*
 * temp_scheme.cpp
 *
 *  Created on: 06/07/2010
 *      Author: Nascimento
 */


#include "RunEngine.h"

using namespace std;

/*
 * Constructor for McLIBELa compilation without GUI
*/

TEMP_SCHEME::TEMP_SCHEME(char* inputfile){

    Input = new PARSER();
    Input->set_parameters(inputfile);

    Writer = new WRITER(Input);

    REC = new Mol2(Input, Input->rec_mol2);
    LIG = new Mol2(Input, Input->lig_mol2);
    RefLig = new Mol2(Input, Input->reflig_mol2);
    Ene = new Energy2(Input);

    // Finding HB donors/acceptors for receptor. Using a cutoff of 10. Ang distance between the each CA atom and the ligand COM.
    FindHB* HB = new FindHB;

    for (int i=0; i< REC->residue_pointer.size()-1; i++){
        HB->parse_residue(REC->residue_pointer[i]-1, REC->residue_pointer[i+1]-2, REC->resnames[i], REC, RefLig, 9.0);
    }

    // Finding HB donors/acceptors for RefLig;
    HB->find_ligandHB(Input->reflig_mol2, RefLig);

    delete HB;

    /* For debuggind purposes
 * Here we write the HB bond donors
 *
    for (unsigned i=0; i< REC->HBdonors.size(); i++){
        printf("[%5d]: %4s / %4s --> %4d / %4d \n", i, REC->atomnames[REC->HBdonors[i][0]].c_str(), REC->atomnames[REC->HBdonors[i][1]].c_str(), REC->HBdonors[i][0],REC->HBdonors[i][1]);
    }
*/

}

/*
 * Costructor for compilation with GUI. In this case
 * the compile definition HAS_GUI has to to be set.
 */

#ifdef HAS_GUI
TEMP_SCHEME::TEMP_SCHEME(PARSER* _Input, QPlainTextEdit* Editor, QProgressBar* _progressbar){

    this->Input = _Input;
    this->progressbar = _progressbar;
    QWriter = new QtWriter(Input, Editor);

    // Reading Receptor files

    REC = new Mol2(Input, Input->rec_mol2);
    LIG = new Mol2(Input, Input->lig_mol2);
    RefLig = new Mol2(Input, Input->reflig_mol2);
    Ene = new Energy2(Input);

    // Finding HB donors/acceptors for receptor. Using a cutoff of 10. Ang distance between the each CA atom and the ligand COM.
    FindHB* HB = new FindHB;

    for (int i=0; i< REC->residue_pointer.size()-1; i++){
        HB->parse_residue(REC->residue_pointer[i]-1, REC->residue_pointer[i+1]-2, REC->resnames[i], REC, RefLig, 9.0);
    }

    // Finding HB donors/acceptors for RefLig;
    HB->find_ligandHB(Input->reflig_mol2, RefLig);

    delete HB;
}
#endif

void TEMP_SCHEME::evaluation(){
    if ((REC->N != int(REC->charges.size())) or (REC->N != int(REC->radii.size())) or (REC->N != int(REC->epsilons.size()))){
        cout <<"The number of atomic parameters and the number of atoms doesn't match!" << endl;
        cout << "Exiting..." << endl;
        sprintf(info, "Number of atomic parameters and the number of atoms doesn't match!");
        this->print_info(info);
        sprintf(info, "Exiting...");
        this->print_info(info);
        usleep(2000);
        exit(1);
    }

    //! Evaluation of the original binding energy
    sprintf(info, "The receptor has %5lu / %5lu HB donors/acceptors around the active site.", REC->HBdonors.size(), REC->HBacceptors.size());
    this->print_info(info);
    start_energy = Ene->compute_ene(REC, RefLig, RefLig->xyz);
    sprintf(info,"Original Binding Energy: %.4f kcal/mol", start_energy);
    this->print_info(info);

    if (Input->reflig_mol2 == ""){
        center = COORD.compute_com(LIG);
    }
    else {
        center = COORD.compute_com(RefLig);
    }

    sprintf(info, "Center of computation box: %.2f %.2f %.2f", center[0], center[1], center[2]);
    this->print_info(info);

    if (Input->use_grids){
        if (Input->load_grid_from_file){
            sprintf(info,"Loading grids from file %s.grid...", Input->grid_prefix.c_str());
            this->print_info(info);
            Grids = new Grid(Input, Writer);
            Grids->load_grids_from_file();
            sprintf(info, "Loaded energy grids with %d x %d x %d points spaced by %5.3f Angstroms in each directon.", Grids->npointsx, Grids->npointsy, Grids->npointsz, Grids->grid_spacing);
            this->print_info(info);
            sprintf(info, "Grid Origin: %10.5f %10.5f %10.5f.", Grids->xbegin, Grids->ybegin, Grids->zbegin);
            this->print_info(info);
            if (Grids->pbsa_loaded){
                sprintf(info, "Electrostatic Potencial computed with PBSA.");
                this->print_info(info);
            }
            else if (Grids->delphi_loaded){
                sprintf(info, "Electrostatic Potencial computed with DelPhi.");
                this->print_info(info);
            }
            else {
                sprintf(info, "Electrostatic Potencial computed with Coulomb model.");
                this->print_info(info);
            }
        }
        else{
            sprintf(info,"Generating energy grids. It can take a couple of minutes. Coffee time maybe ?");
            this->print_info(info);
            start = clock();
#ifdef HAS_GUI
            Grids = new Grid(Input, QWriter, REC, center);
#else
            Grids = new Grid(Input, Writer, REC, center);
#endif
            end = clock();
            sprintf(info, "Computed energy grids with %d x %d x %d points spaced by %5.3f Angstroms in each directon.", Grids->npointsx, Grids->npointsy, Grids->npointsz, Grids->grid_spacing);
            this->print_info(info);
            sprintf(info, "Grid Origin: %10.5f %10.5f %10.5f.", Grids->xbegin, Grids->ybegin, Grids->zbegin);
            this->print_info(info);
            sprintf(info, "Grid computation took %d seconds.", int((end-start)/CLOCKS_PER_SEC));
            this->print_info(info);

            this->write_box(center, Grids->xbegin, Grids->ybegin, Grids->zbegin, Grids->xend, Grids->yend, Grids->zend);

        }
        double grid_energy = Ene->compute_ene(Grids, RefLig, RefLig->xyz);
        sprintf(info,"Original Grid energy: %.4f kcal/mol.", grid_energy);
        this->print_info(info);
        sprintf(info, "Energy error = %.4f%s.", fabs((grid_energy-start_energy)/start_energy)*100., "%");
        this->print_info(info);
    }

    this->print_line();

    if (Input->dock_mode){
        this->dock_run();
    }
    else if (Input->eq_mode){
        this->eq_run();
    }
    else if (Input->mcr_mode){
        this->mcr_run();
    }
    else if (Input->full_search_mode){
        this->full_search_run();
    }
    else if (Input->sa_mode){
        this->sa_run();
    }

    sprintf(info, "Finishing McLibela...");
    this->print_info(info);
    this->print_line();
}

void TEMP_SCHEME::sa_run(){

    if (Input->sa_mode == true){
        sprintf(info,"%s", "Entering SA Scheme....");
        this->print_info(info);
        srand(rand());
        gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs2);
        int seed = (rand() % 101);
        sprintf(info, "Seed: %d", seed);
        this->print_info(info);
        gsl_rng_set(r, seed);

        //! rmsd_energy is defined as an output file to keep the trial number, the rmsd and the evaluated energy.
        sprintf(rmsd_ene, "%s_sa.dat", Input->output.c_str());
        rmsd_energy = fopen(rmsd_ene, "w" );




        if ((Input->flex_lig = true)){
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
        this->print_info(info);

        rmsd = COORD.compute_rmsd(LIG->xyz, old_coord, LIG->N);
        fprintf(rmsd_energy, "%6d % 7.3f % 7.3f % 7.3f\n", -1, rmsd, start_energy, 0.00);

        if( (Input->flex_lig = true)){

            for(double t=Input->sa_start_temp; t>= Input->temp; t-=10*log(t)){

                sprintf(info, "SA Temperature: %.2f", t);
                this->print_info(info);

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
                            this->print_info(info);
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
                                this->print_info(info);
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
                this->print_info(info);
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
                            this->print_info(info);
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
                                this->print_info(info);
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
        this->print_line();
    }
}


void TEMP_SCHEME::dock_run(){
    clock_t ti, tf;
    int time_elapsed;
    if(Input->dock_mode == true){
        ti = clock();
        sprintf(info, "%5s %-12.12s %-4.4s %-10.10s %-8.8s %-8.8s %-8.8s %-8.8s %-8.8s %3.3s %-5.5s %2.2s", "#", "Ligand", "Resi", "Overlay", "Elec", "VDW", "Solv", "HBond", "Energy", "Conf", "MS", "SI");
        this->print_info(info);
        this->print_line();


        if (Input->dock_parallel){
#ifdef HAS_GUI
            this->dock_parallel();
#else
            this->dock_parallel();
#endif
        }
        else{
#ifdef HAS_GUI
            Docker* Dock = new Docker(QWriter);
#else
            Docker* Dock = new Docker(Writer);
#endif

            center = COORD.compute_com(RefLig);

            unsigned counter=0;
            if (Input->multifile != ""){
                string ligand = "";
                ifstream multifile(Input->multifile.c_str());
                if (! multifile.is_open()){
                    sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
                    this->print_info(info);
                    exit(1);
                }
                multifile >> ligand;
                while ((!multifile.eof()) and (ligand != "EOF")){
                    Mol2* Lig2 = new Mol2(Input, ligand);
                    counter++;
                    FindHB* HB = new FindHB;
                    HB->find_ligandHB(ligand, Lig2);
                    delete HB;
                    if (Input->generate_conformers){
                        Conformer* Conf = new Conformer;
                        if (Conf->generate_conformers_confab(Input, Lig2, ligand)){
                            if (Input->use_grids){
                                Dock->run(REC, Lig2, RefLig, center, Input, Grids, counter);
                            }
                            else {
                                Dock->run(REC, Lig2, RefLig, center, Input, counter);
                            }

                            if (Input->eq_mode){
                                Writer->writeMol2(Lig2, Lig2->xyz, 0.0, 0.0, string(Input->output + "_" + Lig2->molname));
                                Conformer* Conf = new Conformer;
                                Lig2->mcoords.clear();
                                Conf->generate_conformers_confab(Input, Lig2, string(Input->output + "_" + Lig2->molname+".mol2.gz"));
                                MC* EqMC = new MC(Lig2, Input, Writer);
                                if (Input->use_grids){
                                    EqMC->run(Grids, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                                }
                                else {
                                    EqMC->run(REC, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                                }
                                if (Input->ligsim){
                                    EqMC->ligand_run(RefLig, Lig2, Lig2->xyz, Input, Input->temp);
                                }
                                delete EqMC;
                            }
                        }
                        else{ 																//test: in case generate_conformers does not succeed.
                            Lig2->mcoords.push_back(Lig2->xyz);
                            if (Input->use_grids){
                                Dock->run(REC, Lig2, RefLig, center, Input, Grids, counter);
                            }
                            else {
                                Dock->run(REC, Lig2, RefLig, center, Input, counter);
                            }
                            if (Input->eq_mode){
                                MC* EqMC = new MC(Lig2, Input, Writer);
                                if (Input->use_grids){
                                    EqMC->run(Grids, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                                }
                                else {
                                    EqMC->run(REC, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                                }
                                if (Input->ligsim){
                                    EqMC->ligand_run(RefLig, Lig2, Lig2->xyz, Input, Input->temp);
                                }
                                delete EqMC;
                            }
                        }
                        Conf->~Conformer();
                        delete Conf;
                    }
                    else {
                        if (Input->use_grids){
                            Dock->run(REC, Lig2, RefLig, center, Input, Grids, counter);
                        }
                        else {
                            Dock->run(REC, Lig2, RefLig, center, Input, counter);
                        }
                        if (Input->eq_mode){
                            MC* EqMC = new MC(Lig2, Input, Writer);
                            if (Input->use_grids){
                                EqMC->run(Grids, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                            }
                            else {
                                EqMC->run(REC, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                            }
                            if (Input->ligsim){
                                EqMC->ligand_run(RefLig, Lig2, Lig2->xyz, Input, Input->temp);
                            }
                            delete EqMC;
                        }
                    }
                    delete Lig2;
                    multifile >> ligand;
#ifdef HAS_GUI
                    progressbar->setValue((counter*100)/Input->docking_molecules.size());
#endif
                }
                multifile.close();
            }
            else {
                counter=1; 														// just one ligand
                FindHB* HB = new FindHB;
                HB->find_ligandHB(Input->lig_mol2, LIG);
                delete HB;
                if (Input->generate_conformers){
                    Conformer* Conf = new Conformer;
                    Conf->generate_conformers_confab(Input, LIG, Input->lig_mol2);
                    if (Input->use_grids){
                        Dock->run(REC, LIG, RefLig, center, Input, Grids, counter);
                    }
                    else {
                        Dock->run(REC, LIG, RefLig, center, Input, counter);
                    }
                    if (Input->eq_mode){
                        Writer->writeMol2(LIG, LIG->xyz, 0.0, 0.0, string(Input->output + "_" + LIG->molname));
                        Conformer* Conf = new Conformer;
                        Conf->generate_conformers_confab(Input, LIG, string(Input->output + "_" + LIG->molname+".mol2.gz"));
                        MC* EqMC = new MC(LIG, Input, Writer);
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
                        delete Conf;
                    }
                }
                else {
                    if (Input->use_grids){
                        Dock->run(REC, LIG, RefLig, center, Input, Grids, counter);
                    }
                    else {
                        Dock->run(REC, LIG, RefLig, center, Input, counter);
                    }
                    MC* EqMC = new MC(LIG, Input, Writer);
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
            delete Dock;
        }

        this->print_line();
        sprintf(info, "Min Status:");
        this->print_info(info);
        sprintf(info, "  Failure = -1, Out of memory = -3, Roundoff limited = -4, Forced stop = -5,");
        this->print_info(info);
        sprintf(info, "  Stopval reached = 2, Ftol reached = 3, Xtol reached = 4, Maxeval reached=5, Maxtime reached=6");
        this->print_info(info);
        this->print_line();
        tf = clock()-ti;
        sprintf(info, "Docking computations took %d minute(s)", int((tf/CLOCKS_PER_SEC)/60.));
        this->print_info(info);
        this->print_line();
    }
}

int TEMP_SCHEME::dock_parallel_function(Mol2* Lig2){
#ifdef HAS_GUI
    Docker* Dock = new Docker(QWriter);
#else
    Docker* Dock = new Docker(Writer);
#endif
    if (Input->use_grids){
        Dock->run(REC, Lig2, RefLig, center, Input, Grids, 0);
    }
    else {
        Dock->run(REC, Lig2, RefLig, center, Input, 0);
    }
    delete Dock;
    return 0;
}

int TEMP_SCHEME::dock_serial(vector<string> ligand_list, int count, int chunck_size){
    stringstream buffer;
    buffer << (count+1);
    WRITER* Write_lig = new WRITER(Input->output + "_" + buffer.str(), Input);
#ifdef HAS_GUI
    Docker* Dock = new Docker(QWriter);
#else
    Docker* Dock = new Docker(Write_lig);
#endif

    count = atoi(ligand_list[ligand_list.size()-1].c_str());
    for (unsigned i=0; i< ligand_list.size()-1; i++){
        Mol2* Lig2 = new Mol2;

        bool lig_is_opened = false;

        if (ligand_list[i].substr(ligand_list[i].size()-3, 3) == ".gz"){
            lig_is_opened = Lig2->parse_gzipped_file(Input, ligand_list[i]);
        }
        else {
            lig_is_opened = Lig2->parse_mol2file(Input, ligand_list[i]);
        }

        if (lig_is_opened){
            FindHB* HB = new FindHB;
            HB->find_ligandHB(ligand_list[i], Lig2);
            delete HB;
            if (Input->generate_conformers){
                Conformer* Conf = new Conformer;
                Conf->generate_conformers_confab(Input, Lig2, ligand_list[i]);
                delete Conf;
            }
            if (Input->use_grids){
                Dock->run(REC, Lig2, RefLig, center, Input, Grids, ((count*chunck_size))+i+1);
            }
            else {
                Dock->run(REC, Lig2, RefLig, center, Input, ((count*chunck_size))+i+1);
            }
            if (Input->eq_mode){
                if (Input->generate_conformers){
                    Conformer* Conf = new Conformer;
                    Write_lig->writeMol2(Lig2, Lig2->xyz, 0.0, 0.0, string(Input->output + "_" + Lig2->molname));
                    Lig2->mcoords.clear();
                    Conf->generate_conformers_confab(Input, Lig2, string(Input->output +"_" + Lig2->molname+".mol2.gz"));
                    delete Conf;
                }
                MC* EqMC = new MC(Lig2, Input, Write_lig);
                if (Input->use_grids){
                    EqMC->run(Grids, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                }
                else {
                    EqMC->run(REC, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                }
                if (Input->ligsim){
                    EqMC->ligand_run(RefLig, Lig2, Lig2->xyz, Input, Input->temp);
                }
                delete EqMC;
            }
        }
        delete Lig2;
    }
    delete Dock;
    delete Write_lig;
    return 0;
}

int TEMP_SCHEME::dock_concurrent_function(string mylig){
    Mol2* Lig2 = new Mol2(Input, mylig);

    if (Input->generate_conformers){
        Conformer* Conf = new Conformer;
#pragma omp critical
        {
            Conf->generate_conformers_confab(Input, Lig2, mylig);
        }
        delete Conf;
    }
#ifdef HAS_GUI
    Docker* Dock = new Docker(QWriter);
#else
    Docker* Dock = new Docker(Writer);
#endif
    if (Input->use_grids){
        Dock->run(REC, Lig2, RefLig, center, Input, Grids, 0);
    }
    else {
        Dock->run(REC, Lig2, RefLig, center, Input, 0);
    }
    delete Dock;
    return 0;
}

void TEMP_SCHEME::dock_parallel(){

#ifdef HAS_MPI

    this->dock_mpi();

#else

    sprintf(info, "Running Dock mode in parallel with up to %d threads...", int(Input->parallel_jobs));
    this->print_info(info);
    vector<string> ligand_list;


    /* Here, we read the multifile file and parse the files of the molecules
 * to be docked into a std::vector named ligand_list.
 */

    if (Input->multifile != ""){
        ifstream multifile(Input->multifile.c_str());
        if (! multifile.is_open()){
            sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
            this->print_info(info);
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


#pragma omp critical
                    {
                        FindHB* HB = new FindHB;
                        HB->find_ligandHB(ligand_list[i], Lig2);
                        delete HB;
                    }


                    if (Input->generate_conformers){

#pragma omp critical
                        {
                            Conformer* Conf = new Conformer;
                            Conf->generate_conformers_confab(Input, Lig2, ligand_list[i]);
                            delete Conf;
                        }

                    }

#ifdef HAS_GUI
                    Docker* Dock = new Docker(QWriter);
#else
                    Docker* Dock = new Docker(Writer);
#endif
                    if (Input->use_grids){
                        Dock->run(REC, Lig2, RefLig, center, Input, Grids, i+1);
                    }
                    else {
                        Dock->run(REC, Lig2, RefLig, center, Input, i+1);
                    }
                    delete Dock;
                    if (Input->eq_mode){
                        if (Input->generate_conformers){
                            Writer->writeMol2(Lig2, Lig2->xyz, 0.0, 0.0, string(Input->output + "_" + Lig2->molname));
                            Conformer* Conf = new Conformer;
                            Conf->generate_conformers_confab(Input, Lig2, string(Input->output + "_" + Lig2->molname+".mol2.gz"));
                            delete Conf;
                        }
                        MC* EqMC = new MC(Lig2, Input, Writer);
                        if (Input->use_grids){
                            EqMC->run(Grids, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                        }
                        else {
                            EqMC->run(REC, RefLig , Lig2, Lig2->xyz, Input, Input->temp);
                        }
                        if (Input->ligsim){
                            EqMC->ligand_run(RefLig, Lig2, Lig2->xyz, Input, Input->temp);
                        }
                        delete EqMC;
                    }
                }
                delete Lig2;

#ifdef HAS_GUI
                if (omp_get_thread_num() ==0) {
                    progressbar->setValue(round((i+1)*100/int(ligand_list.size())));
                    QApplication::processEvents();
                }
#endif
            }
        }
        ligand_list.clear();

#ifdef HAS_GUI

        progressbar->setValue(100);

#endif

    }

#endif

}

#ifdef HAS_GUI
void TEMP_SCHEME::dock_concurrent(){
    sprintf(info, "Running Dock mode in parallel with up to %d threads...", int(Input->parallel_jobs));
    this->print_info(info);
    QVector<string> ligand_list;

    /* Here we read the multifile file and parse the files of the molecules
 * to be docked into a std::vector named ligand_list.
 */

    if (Input->multifile != ""){
        ifstream multifile(Input->multifile.c_str());
        if (! multifile.is_open()){
            sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
            this->print_info(info);
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
        this->print_info(info);

        if (Input->multifile == ""){
            sprintf(info, "Cannot run single molecule docking in parallel");
            this->print_info(info);
            exit(1);
        }

        ifstream multifile(Input->multifile.c_str());
        if (! multifile.is_open()){
            sprintf(info, "Could not open file %s. Exiting...\n", Input->multifile.c_str());
            this->print_info(info);
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
                chuncks[chuncks.size()-1-i].push_back(ligand_list[i]);
            }
            ligand_list.clear();
        }

        int counter=0;
        for (unsigned i=0; i< chuncks.size(); i++){
            char buffer[10];
            sprintf(buffer, "%d", counter);
            chuncks[i].push_back(string(buffer));
            counter += chuncks[i].size()-1;
        }

        scatter(world,chuncks, tmp,  0);
        this->dock_serial(tmp, world.rank(), 1);
    }                                                           // End of rank 0;
    else {
        scatter(world, tmp, 0);
        this->dock_serial(tmp, world.rank(), 1);
    }
}

#endif

void TEMP_SCHEME::eq_run(){
    if (Input->eq_mode){
        MC* EqMC = new MC(LIG, Input, Writer);
        if (Input->generate_conformers){
            Conformer* Conf = new Conformer;
            if (Input->dock_mode){
                Writer->writeMol2(LIG, LIG->xyz, 0.0, 0.0, "Lig_docked");
                Conf->generate_conformers_confab(Input, LIG, "Lig_docked.mol2.gz");
            }
            else {
                Conf->generate_conformers_confab(Input, LIG, Input->lig_mol2);
            }
            delete Conf;
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

        if (Input->use_only_binding_energy){
            EqMC->average_deltaE = EqMC->average_bound_energy;
        }
        else {
            EqMC->average_deltaE = EqMC->average_bound_energy - EqMC->average_freeligand_energy;
        }
        double dE= EqMC->average_deltaE - (EqMC->boundTS-EqMC->freeTS);
        Writer->print_line();
        Writer->print_line();
        sprintf(info, " Final Binding Energies:");
        Writer->print_info(info);
        sprintf(info, "     DeltaU = <Ebound> - <Efree> = %7.3f kcal/mol", EqMC->average_deltaE);
        Writer->print_info(info);
        sprintf(info, "     TdS = - (<TSbound> - <TSfree>) = %7.3f kcal/mol", (EqMC->boundTS-EqMC->freeTS));
        Writer->print_info(info);
        sprintf(info, "     Binding Energy = <dU> - TdS = %7.3f kcal/mol", dE);
        Writer->print_info(info);
        Writer->print_line();
        Writer->print_line();
        delete EqMC;
    }
}

void TEMP_SCHEME::mcr_run(){
    if (Input->mcr_mode){

        // checking consistency
        if (int(Input->mcr_coefficients.size()) != Input->mcr_size){
            sprintf(info, "Inconsistency found in MCR coefficients. Please check!\n");
            this->print_info(info);
            exit(1);
        }

        MC* EqMC = new MC(LIG, Input, Writer);
        if (Input->generate_conformers){
            Conformer* Conf = new Conformer;
            Conf->generate_conformers_confab(Input, LIG, Input->lig_mol2);
            delete Conf;
        }

        sprintf(info, "MCR %7.7s %7.7s %10.10s %10.10s %10.10s %10.10s %10.10s %7.7s %7.7s",  "#i", "bi", "bT", "<ene>" , "SD(ene)", "W(b,bt)" , "-kTln(W)", "A_ex" , "Vol(A3)");
        this->print_info(info);
        this->print_line();

        double bt;                          // MC Recursion "effective" temperature (bt) fot ith evaluation;
        double k = 0.0019858775203792202;   // Boltzmann constant in kcal/(mol.K)

        long double cum_W = 0.0L;
        long double cum_W_err = 0.0L;
        double max_vol=0.0;
        double volume;
        long double cum_Ln_W = 0.0L;

        // Doing a equilibrium simulation at the default temperature
        // before starting the recursion.

        Input->bi = 2.0;
        if (Input->use_grids){
            EqMC->run(Grids, RefLig , LIG, LIG->xyz, Input, Input->temp);
        }
        else {
            EqMC->run(REC, RefLig , LIG, LIG->xyz, Input, Input->temp);
        }

        volume = (EqMC->XSize*EqMC->YSize*EqMC->ZSize);

        this->print_line();
        sprintf(info, "%s", "Starting equilibrium simulation before recursion");
        this->print_info(info);
        sprintf(info, "MCR %7d %7.4f %10.4g %10.4g %10.4g %10.4g %10.4g %7.3Lf %7.4g", 0, 1.0, Input->temp, EqMC->average_energy, EqMC->energy_standard_deviation, EqMC->average_energy,
                -k*Input->temp*log(double(EqMC->average_energy)), cum_W, volume);
        this->print_info(info);
        this->print_line();

        //
        // Now starting the MC recursion...
        //

        string mcr_output_prefix = Input->output;

        for (int i=0; i<Input->mcr_size; i++){
            Input->bi = Input->mcr_coefficients[i];
            bt = Input->temp;
            for (int j=0; j <= i; j++){
                bt = bt*Input->mcr_coefficients[j];
            }

            char buffer_output [mcr_output_prefix.size()+10];
            sprintf(buffer_output, "%s_MCR_%d", mcr_output_prefix.c_str(), i);
            Input->output = string(buffer_output);

            if (Input->use_grids){
                EqMC->run(Grids, RefLig , LIG, LIG->xyz, Input, bt);
            }
            else {
                EqMC->run(REC, RefLig , LIG, LIG->xyz, Input, bt);
            }

            volume = (EqMC->XSize*EqMC->YSize*EqMC->ZSize);

            cum_Ln_W += -k*bt*log(double(EqMC->average_energy));
//            cum_Ln_W += -k*Input->temp*log(double(EqMC->MCR_Boltzmann_weighted_average));

            cum_W += (log(double(EqMC->average_energy)));
            cum_W_err += double((1.0/double(EqMC->average_energy))*EqMC->energy_standard_deviation);

            this->print_line();
            sprintf(info, "MCR %7d %7.4f %10.4g %10.4g %10.4g %10.4g %10.4g %7.3Lf %7.4g", i+1, Input->bi, bt, EqMC->average_energy, EqMC->energy_standard_deviation, EqMC->average_energy,
                    -k*bt*log(double(EqMC->average_energy)), cum_W, volume);
            this->print_info(info);
            this->print_line();


            if (volume > max_vol){
                max_vol=volume;
            }
        }

        this->print_line();
        sprintf(info, "MCR: A_excess = %10.4Lg +/- %10.4Lg",  -k*Input->temp*cum_W, k*Input->temp*cum_W_err);
        this->print_info(info);

        sprintf(info, "MCR: Volume: %10.4g.  ln(volume) = %10.4g",  volume, log(volume));
        this->print_info(info);

        sprintf(info, "MCR: A_complex = %10.4Lg",  (-k*Input->temp*log(volume))+(-k*Input->temp*cum_W));
        this->print_info(info);

        this->print_line();


        long double cum_W_lig = 0.0L;
        long double cum_W_lig_err = 0.0L;
        double lig_max_vol = 0.0;
        double lig_volume = 0.0;
        long double cum_Ln_W_lig = 0.0L;

        // Equilibration before recursion

        Input->bi = 1.0;
        this->print_line();
        sprintf(info, "%s", "Starting equilibrium simulation before recursion");
        this->print_info(info);

        EqMC->ligand_run(RefLig, LIG, LIG->xyz, Input, Input->temp);

        lig_volume = (EqMC->XSize*EqMC->YSize*EqMC->ZSize);
        sprintf(info, "MCR %7d %7.4f %10.4g %10.4g %10.4g %10.4g %10.4g %7.7Lf %7.4g", 0, 1.0, Input->temp, EqMC->average_energy, EqMC->energy_standard_deviation, EqMC->average_energy,
                -k*Input->temp*log(double(EqMC->average_energy)), cum_W_lig, lig_volume);
        this->print_info(info);
        this->print_line();

        //
        // Now, starting the MC recursion
        //

        if (Input->ligsim){
            for (int i=0; i<Input->mcr_size; i++){
                Input->bi = Input->mcr_coefficients[i];
                bt = Input->temp;
                for (int j=0; j <= i; j++){
                    bt = bt*Input->mcr_coefficients[j];
                }

                EqMC->ligand_run(RefLig, LIG, LIG->xyz, Input, bt);

                lig_volume = (EqMC->XSize*EqMC->YSize*EqMC->ZSize);

                cum_Ln_W_lig += -k*bt*log(double(EqMC->average_energy));
//                cum_Ln_W_lig += -k*Input->temp*log(double(EqMC->MCR_Boltzmann_weighted_average));

                cum_W_lig += (log(double(EqMC->average_energy)));
                cum_W_lig_err += double((1.0/double(EqMC->average_energy))*EqMC->energy_standard_deviation);

                sprintf(info, "MCR %7d %7.4f %10.4g %10.4g %10.4g %10.4Lg %10.4g %7.7Lf %7.4g", i+1, Input->bi, bt, EqMC->average_energy, EqMC->energy_standard_deviation, EqMC->MCR_Boltzmann_weighted_average,
                        -k*bt*log(double(EqMC->average_energy)), cum_W_lig, lig_volume);
                this->print_info(info);
                this->print_line();


                if (volume > max_vol){
                    lig_max_vol=lig_volume;
                }
            }

            this->print_line();
            sprintf(info, "MCR: A_excess for the ligand = %10.4Lg +/- %10.4Lg",  -k*Input->temp*cum_W_lig, k*Input->temp*cum_W_lig_err);
            this->print_info(info);

            sprintf(info, "MCR: Ligand Volume: %10.4g.  ln(volume) = %10.4g",  lig_volume, log(lig_volume));
            this->print_info(info);

            sprintf(info, "MCR: A_lig = %10.4Lg",  (-k*Input->temp*log(lig_volume))+(-k*Input->temp*cum_W_lig));
            this->print_info(info);

            this->print_line();

            long double complex_A = (-k*Input->temp)*(log(volume)+(cum_W));
            long double ligand_A = (-k*Input->temp)*(log(lig_volume)+(cum_W_lig));
            long double Delta_A = complex_A - ligand_A;

            sprintf(info, "MCR: Complex Free Energy = %10.4Lg +/- %10.4Lg", complex_A, (k*Input->temp*cum_W_err));
            this->print_info(info);
            sprintf(info, "MCR: Ligand Free Energy = %10.4Lg +/- %10.4Lg", ligand_A, (k*Input->temp*cum_W_lig_err));
            this->print_info(info);
            sprintf(info, "MCR: Binding Free Energy = %10.4Lg +/- %10.4Lg", Delta_A, (k*Input->temp*(cum_W_err+cum_W_lig_err)));
            this->print_info(info);
            this->print_line();
        }
        delete EqMC;
    }
}

void TEMP_SCHEME::full_search_run(){
    if (Input->full_search_mode){
        if (Input->use_grids){
            FullSearch* Search = new FullSearch(Input, LIG, Writer, Grids);
            Search->do_search();
            delete Search;
        }
        else {
            FullSearch* Search = new FullSearch(Input, LIG, Writer);
            Search->do_search();
            delete Search;
        }    }
}

void TEMP_SCHEME::print_info(char info[98]){
#ifdef HAS_GUI
    QWriter->print_info(info);
#else
    Writer->print_info(info);
#endif
}

void TEMP_SCHEME::print_line(){
#ifdef HAS_GUI
    QWriter->print_line();
#else
    Writer->print_line();
#endif
}

void TEMP_SCHEME::write_box(vector<double>center, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z){
#ifdef HAS_GUI
    QWriter->write_box(center, min_x, min_y, min_z, max_x, max_y, max_z);
#else
    Writer->write_box(center, min_x, min_y, min_z, max_x, max_y, max_z);
#endif
}


TEMP_SCHEME::~TEMP_SCHEME(){
    delete Ene;
    delete RefLig;
    delete LIG;
    delete REC;
#ifdef HAS_GUI
    delete QWriter;
#else
    delete Writer;
    delete Input;
#endif
}


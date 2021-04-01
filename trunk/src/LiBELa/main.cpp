#include "main.h"

using namespace std;


/*********************************************************************************************
 *           iMcLiBELa - Monte Carlo - Based Ligand Binding Energy Landscape                 *
 *                                                                                           *
 *               This program was written by Alessandro Nascimento - 2012-2020               *
 * 					University of Sao Paulo - Sao Carlos, SP. Brazil                         *
 *                      Sao Carlos Institute of Physics - IFSC/USP                           *
 *                                asnascimento@ifsc.usp.br                                   *
 *                                                                                           *
 *           Please visit http://nascimento.ifsc.usp.br/  for detailed instructions          *
 *                  Distributed under the terms of GNU LGPLv3 license                        *
 *********************************************************************************************/


void write_dock_input(){
    printf("# Mode\n");
    printf("mode                    dock\n");
    printf("dock_parallel           no\n");
    printf("parallel_jobs           1\n");
    printf("\n");
    printf("#\n");
    printf("# Input files\n");
    printf("#\n");
    printf("rec_mol2              	../mol2/rec.mol2.gz\n");
    printf("lig_mol2                ../mol2/lig.mol2.gz\n");
    printf("reflig_mol2             ../mol2/lig.mol2.gz\n");
    printf("mol2_aa                 no\n");
    printf("multifile               multimol.dat\n");
    printf("\n");
    printf("#\n");
    printf("# force field parameters\n");
    printf("#\n");
    printf("\n");
    printf("search_box              8.0 8.0 8.0\n");
    printf("scoring_function        0\n");
    printf("dielectric_model        r\n");
    printf("diel                    2.0\n");
    printf("deltaij                 1.75\n");
    printf("deltaij_es              1.50\n");
    printf("use_grids               yes\n");
    printf("grid_spacing            0.4\n");
    printf("grid_box                30.0    30.0    30.0\n");
    printf("write_grids             McGrid\n");
    printf("use_delphi              no\n");
    printf("solvation_alpha         0.1\n");
    printf("solvation_beta          -0.005\n");
    printf("LJ_sigma                3.0\n");
    printf("scale_vdw_energy        1.0\n");
    printf("scale_elec_energy       1.0\n");
    printf("\n");
    printf("#\n");
    printf("# PBSA\n");
    printf("#\n");
    printf("\n");
    printf("pbsa_grid               ../pbsa_grids/pbsa_phi.phi.gz\n");
    printf("use_pbsa                no\n");
    printf("\n");
    printf("#\n");
    printf("# DelPhi Parameters\n");
    printf("#\n");
    printf("\n");
    printf("delphi_grid             ../delphi_grids/phimap.phi\n");
    printf("use_delphi              no\n");
    printf("delphi_gsize            75\n");
    printf("\n");
    printf("#\n");
    printf("# Optimization\n");
    printf("#\n");
    printf("\n");
    printf("search_box              8.0    8.0    8.0\n");
    printf("minimization_tolerance  1.0e-4\n");
    printf("minimization_delta      1.0e-4\n");
    printf("dock_min_tol            1.0e-4\n");
    printf("minimization_timeout    10\n");
    printf("overlay_optimizer       mma\n");
    printf("energy_optimizer      	direct\n");
    printf("ignore_h                no\n");
    printf("deal                    no\n");
    printf("elec_scale              1.0\n");
    printf("vdw_scale               1.0\n");
    printf("sort_by_energy          no\n");
    printf("use_docking_restraints	no\n");
    printf("restraints_weight       0.0\n");
    printf("use_overlay_cutoff      no\n");
    printf("overlay_cutoff          0.6\n");
    printf("\n");
    printf("#\n");
    printf("# SA Options\n");
    printf("#\n");
    printf("\n");
    printf("sa_start_temp           100\n");
    printf("temperature             0.001\n");
    printf("sa_steps                1000\n");
    printf("cushion                 1.0\n");
    printf("rotation_step           10.0\n");
    printf("sa_mu_t                 1.1\n");
    printf("\n");
    printf("#\n");
    printf("# Output\n");
    printf("#\n");
    printf("\n");
    printf("output_prefix           LiBELa\n");
    printf("write_mol2              yes\n");
    printf("\n");
    printf("#\n");
    printf("# Flexible Ligands\n");
    printf("#\n");
    printf("\n");
    printf("generate_conformers     yes\n");
    printf("number_of_conformers	10\n");
    printf("conformers_to_rank      3\n");
    printf("conf_search_trials      10000\n");
}

void write_mc_input(){
    printf("# Mode\n");
    printf("mode                    eq\n");
    printf("dock_parallel           yes\n");
    printf("parallel_jobs           2\n");
    printf("\n");
    printf("#\n");
    printf("# Input files\n");
    printf("#\n");
    printf("rec_mol2              	../mol2/rec.mol2.gz\n");
    printf("lig_mol2                ../mol2/lig.mol2.gz\n");
    printf("reflig_mol2             ../mol2/lig.mol2.gz\n");
    printf("mol2_aa                 no\n");
    printf("multifile               multimol.dat\n");
    printf("\n");
    printf("#\n");
    printf("# force field parameters\n");
    printf("#\n");
    printf("\n");
    printf("scoring_function        0\n");
    printf("dielectric_model        r\n");
    printf("diel                    2.0\n");
    printf("deltaij                 1.75\n");
    printf("deltaij_es              1.50\n");
    printf("use_grids               yes\n");
    printf("grid_spacing            0.4\n");
    printf("grid_box                30.0    30.0    30.0\n");
    printf("write_grids             McGrid\n");
    printf("use_delphi              no\n");
    printf("solvation_alpha         0.1\n");
    printf("solvation_beta          -0.005\n");
    printf("LJ_sigma                3.0\n");
    printf("scale_vdw_energy        1.0\n");
    printf("scale_elec_energy       1.0\n");
    printf("\n");
    printf("#\n");
    printf("# PBSA\n");
    printf("#\n");
    printf("\n");
    printf("pbsa_grid               ../pbsa_grids/pbsa_phi.phi.gz\n");
    printf("use_pbsa                no\n");
    printf("\n");
    printf("#\n");
    printf("# DelPhi Parameters\n");
    printf("#\n");
    printf("\n");
    printf("delphi_grid             ../delphi_grids/phimap.phi\n");
    printf("use_delphi              no\n");
    printf("delphi_gsize            75\n");
    printf("\n");
    printf("#\n");
    printf("# Optimization\n");
    printf("#\n");
    printf("\n");
    printf("search_box              8.0    8.0    8.0\n");
    printf("minimization_tolerance  1.0e-4\n");
    printf("minimization_delta      1.0e-4\n");
    printf("dock_min_tol            1.0e-4\n");
    printf("minimization_timeout    10\n");
    printf("overlay_optimizer       mma\n");
    printf("energy_optimizer      	direct\n");
    printf("ignore_h                no\n");
    printf("deal                    no\n");
    printf("elec_scale              1.0\n");
    printf("vdw_scale               1.0\n");
    printf("sort_by_energy          no\n");
    printf("use_docking_restraints	no\n");
    printf("restraints_weight       0.0\n");
    printf("use_overlay_cutoff      no\n");
    printf("overlay_cutoff          0.6\n");
    printf("\n");
    printf("#\n");
    printf("# MC Options\n");
    printf("#\n");
    printf("\n");
    printf("nsteps                  10000000\n");
    printf("temperature             150.0\n");
    printf("cushion                 0.5\n");
    printf("max_atom_displacement   0.01\n");
    printf("rotation_step			1.25\n");
    printf("torsion_step			1.25\n");
    printf("sample_torsions			yes\n");
    printf("mc_full_flex			no\n");
    printf("compute_rotation_entropy	no\n");
    printf("equilibration_steps		1000000\n");
    printf("ligand_simulation		yes\n");
    printf("mc_stride               5000\n");
    printf("seed                    56784\n");
    printf("#\n");
    printf("# MCR Options\n");
    printf("#\n");
    printf("\n");
    printf("mcr_size                20\n");
    printf("mcr_coefficients	1.090507733 1.090507733 1.090507733 1.189207115 1.189207115 1.189207115 1.189207115 1.189207115 1.189207115 1.189207115 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0\n");
    printf("\n");
    printf("#\n");
    printf("# SA Options\n");
    printf("#\n");
    printf("\n");
    printf("sa_start_temp           100\n");
    printf("temperature             0.001\n");
    printf("sa_steps                1000\n");
    printf("cushion                 1.0\n");
    printf("rotation_step           10.0\n");
    printf("sa_mu_t                 1.1\n");
    printf("\n");
    printf("#\n");
    printf("# Output\n");
    printf("#\n");
    printf("\n");
    printf("output_prefix           LiBELa\n");
    printf("write_mol2              yes\n");
    printf("\n");
    printf("#\n");
    printf("# Flexible Ligands\n");
    printf("#\n");
    printf("\n");
    printf("generate_conformers     no\n");
    printf("number_of_conformers	10\n");
    printf("conformers_to_rank      1\n");
    printf("conf_search_trials      10000\n");
}

void usage(char *executable){
    printf("Usage:\n");
    printf("\t %s -i input_file       OR \n", executable);
    printf("\t %s -d                  to generate a docking input file\n", executable);
    printf("\t %s -m                  to generate a Monte Carlo input file\n", executable);
    printf("\t %s                     to use a GUI version (if compiled with this option)\n", executable);
}


int main(int argc, char *argv[]){

#ifdef HAS_GUI

    if (argc != 1){

        int c;
        char *inputfile;

        while ((c = getopt (argc, argv, "i:dm")) != -1)
            switch (c)
            {
            case 'i':
                inputfile = optarg;
                break;
            case 'd':
                write_dock_input();
                break;
            case 'm':
                write_mc_input();
                break;
            case '?':
                if (optopt == 'c')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt)){
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    usage(argv[0]);
                }
                else{
                    fprintf (stderr, "Unknown option character `\\x%x'.\n",optopt);
                    usage(argv[0]);
                }
                return 1;
                break;

            default:
                usage(argv[0]);
            }

        TEMP_SCHEME* RunEngine = new TEMP_SCHEME(inputfile);
        RunEngine->evaluation();
        delete RunEngine;
        return(0);
    }

    else {			// using GUI
        QApplication app(argc, argv);
        app.setOrganizationName(ORGANIZATION);
        app.setApplicationName(NAME);
        app.setApplicationVersion(VERSION);

        /*
        QFile stylesheetfile(":/css/stylesheet.css");
        if (stylesheetfile.open(QFile::ReadOnly)){
            app.setStyleSheet(stylesheetfile.readAll());
            stylesheetfile.close();
        }
*/
        QSplashScreen *splash = new QSplashScreen;
        splash->setPixmap(QPixmap(":/libel.png"));
        splash->show();
        splash->showMessage(app.organizationName());
        QIcon windowIcon(":/libel.png");

        GUI gui;
        gui.setWindowIcon(windowIcon);

        QTimer::singleShot(1000, splash, SLOT(close()));
        QTimer::singleShot(1000, &gui, SLOT(showMaximized()));

        return app.exec();
    }
#else
    if (argc != 1){

        int c;
        char *inputfile;

        while ((c = getopt (argc, argv, "i:dm")) != -1)
            switch (c)
            {
            case 'i':
                inputfile = optarg;
                break;
            case 'd':
                write_dock_input();
                break;
            case 'm':
                write_mc_input();
                break;
            case '?':
                if (optopt == 'c')
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint (optopt)){
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                    usage(argv[0]);
                }
                else{
                    fprintf (stderr, "Unknown option character `\\x%x'.\n",optopt);
                    usage(argv[0]);
                }
                return 1;
                break;

            default:
                usage(argv[0]);
            }

        TEMP_SCHEME* RunEngine = new TEMP_SCHEME(inputfile);
        RunEngine->evaluation();
        delete RunEngine;
        return(0);
    }
#endif
}

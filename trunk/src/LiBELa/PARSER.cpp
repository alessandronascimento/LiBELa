
/*
 * PARSER.cpp
 *
 *  Created on: 28/05/2010
 *      Author: Nascimento
 */

#include "PARSER.h"

using namespace std;

PARSER::PARSER(){
	this->flex_lig = false;
	this->flex_rec = false;
	this->dynamic_rate = false;
    this->timeout = 120;
    this->scoring_function = 3;
    this->min_delta = 1.0e-6;
    this->min_tol = 1.0e-6;
    this->min_timeout = 30;
	this->elec_scale = 1.0;
	this->vdw_scale = 1.0;
	this->multifile = "";
    this->overlay_optimizer = "ln_auglag";
    this->energy_optimizer = "direct";
	this->deltaij6 = pow(2.75, 6);
	this->deltaij_es6 = pow(1.75, 6);
	this->cushion = 1.5;
	this->sigma = 3.5;
	this->sa_scheme = false;
	this->sa_start_temp = 5000.0;
	this->sa_steps = 500;
	this->sa_mu_t = 1.01;
	this->diel = 1.0;
	this->output = "McLiBELa";
	this->rotation_step = 15.0;
	this->x_dim = 15.0;
	this->y_dim = 15.0;
	this->z_dim = 15.0;
	this->temp = 300.0;
	this->number_steps = 1000;
    this->eq_steps = 0;
	this->dock_no_h = false;
	this->deal_pop_size = 100;
	this->deal_generations = 30;
	this->deal_top_ranked = 20;
	this->deal = false;
    this->generate_conformers = true;
	this->lig_conformers = 10;
	this->conformers_to_evaluate = 1;
	this->mol2_aa = false;
	this->conformer_generator = "GA";
    this->conformer_min_steps = 1000;
	this->dock_min_tol = 0.001;
	this->sa_mode = false;
	this->dock_mode = false;
	this->eq_mode = false;
    this->mcr_mode = false;
	this->dock_parallel = false;
    this->parallel_jobs = 1 ;
	this->write_mol2 = true;
    this->grid_spacing = 0.3;
	this->use_grids = false;
	this->search_box_x = 12.0;
	this->search_box_y = 12.0;
	this->search_box_z = 12.0;
	this->load_grid_from_file = false;
	this->write_grids = false;
	this->show_rmsd = false;
	this->sort_by_energy = false;
    this->dielectric_model = "r";
    this->only_score = false;
    this->solvation_alpha = 0.2;
    this->solvation_beta = -0.005;
    this->seed = (rand()/(RAND_MAX + 1.0));
    this->ligsim = false;
    this->mc_stride = 1;
    this->mcr_size = 1;
    this->bi = 2.0;
    this->torsion_step = 10.0;
    this->sample_torsions = false;
    this->verbose = false;
    this->entropy_rotation_bins = 360;
    this->entropy_translation_bins = this->x_dim*2;
    this->ligand_energy_model = "GAFF";
    this->atomic_model_ff = "GAFF";
    this->elec_weight = 1.0;
    this->vdw_weight = 1.0;
    this->solvation_weight = 1.0;
    this->pbsa_grid = "";
}

void PARSER::comparing (string param, ifstream &input) {
	if (param == "nsteps") {
		input >> PARSER::number_steps;
	}
	else if (param == "temperature") {
		input >> PARSER::temp;
	}
	else if (param=="sa_start_temp"){
			input >> PARSER::sa_start_temp;
	}
	else if (param == "mode"){
		input >> tmp;
		string mode;
		size_t plus = 0;
		while (plus != string::npos){
			plus = tmp.find("+");
			mode = tmp.substr(0,plus);
			if (mode == "dock" or mode == "DOCK" or mode == "Dock"){
				this->dock_mode = true;
			}
			else if (mode == "sa" or mode == "SA" or mode == "Sa"){
				this->sa_mode = true;
			}
			else if (mode == "equilibrium" or mode == "eq" or mode == "EQ"){
				this->eq_mode = true;
			}
            else if (mode == "mcr" or mode == "MCR"){
                this->mcr_mode = true;
            }
			tmp=tmp.substr(plus+1, tmp.size());
		}
	}
	else if (param=="sa_steps"){
		input >> PARSER::sa_steps;
	}
	else if (param == "sa_mu_t"){
		input >> this->sa_mu_t;
	}
	else if (param == "cushion") {
		input >> PARSER::cushion;
	}
	else if (param == "deltaij6") {
		input >> PARSER::deltaij6;
    }
    else if (param == "deltaij"){
        double dij;
        input >> dij;
        this->deltaij6 = pow(dij, 6);
    }
	else if (param == "deltaij_es6"){
		input >> PARSER::deltaij_es6;
        this->deltaij_es3 = sqrt(this->deltaij_es6);
	}
    else if (param == "deltaij_es"){
        double dijes;
        input >> dijes;
        this->deltaij_es6 = pow(dijes, 6);
        this->deltaij_es3 = pow(dijes, 3);
    }
    else if (param == "diel"){
		input >> PARSER::diel;
	}
	else if (param == "sigma"){
		input >> PARSER::sigma;
	}
	else if (param == "rec_mol2"){
		input >> this->rec_mol2;
	}
	else if (param == "lig_mol2"){
		input >> this->lig_mol2;
	}
	else if(param == "mol2_aa"){
		input >> this->sa;
		if (sa == "yes" or sa == "YES" or sa == "Yes"){
			this->mol2_aa = true;
		}
		else {
			this->mol2_aa = false;
		}
	}
	else if (param == "output_prefix"){
			input >> PARSER::output;
	}
	else if (param == "rotation_step"){
				input >> PARSER::rotation_step;
	}
    else if (param == "sample_torsions"){
        input >> tmp;
        if (tmp == "yes" or tmp == "YES" or tmp == "Yes"){
            this->sample_torsions = true;
        }
    }
    else if (param == "torsion_step"){
        input >> this->torsion_step;
    }
	else if (param == "timeout"){
					input >> PARSER::timeout;
	}
	else if (param == "grid_box"){
		input >> PARSER::x_dim >> PARSER::y_dim >> PARSER::z_dim;
	}
    else if (param == "lig_traj"){
        input >> PARSER::lig_traj;
        PARSER::flex_lig = true;
    }
    else if (param == "rec_traj"){
        input >> PARSER::rec_traj;
        PARSER::flex_rec = true;
	}
	else if( param == "scoring_function"){
		input >> PARSER::scoring_function;
	}
	else if ( param == "dyn_rate_inf_limit"){
		this->dynamic_rate = true;
		input >> this->dyn_rate_inf_limit;
	}
	else if (param == "dyn_rate_steps"){
		input >> this->dyn_rate_steps;
	}
	else if (param == "dyn_rate_sup_limit"){
		this->dynamic_rate = true;
		input >> this->dyn_rate_sup_limit;
	}
	else if (param == "elec_scale"){
		input >> this->elec_scale;
	}
	else if (param == "vdw_scale"){
		input >> this->vdw_scale;
	}
	else if (param == "reflig_mol2"){
		input >> this->reflig_mol2;
	}
	else if (param == "minimization_tolerance"){
		input >> this->min_tol;
	}
	else if (param == "minimization_delta"){
		input >> this->min_delta;
	}
	else if (param == "minimization_timeout"){
		input >> this->min_timeout;
	}
	else if (param == "multifile"){
		input >> this->multifile;
	}
	else if (param == "overlay_optimizer"){
		input >> this->overlay_optimizer;
	}
	else if (param == "energy_optimizer"){
		input >> this->energy_optimizer;
	}
	else if (param == "ignore_h"){
		input >> this->tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->dock_no_h = true;
		}
		else {
			this->dock_no_h = false;
		}
	}
	else if (param == "deal"){
		input >> this->tmp;
		if (tmp == "yes" or tmp == "YES" or tmp == "Yes"){
			this->deal = true;
		}
		else {
			this->deal = false;
		}
	}
	else if (param == "deal_population_size"){
		input >> this->deal_pop_size;
	}
	else if (param == "deal_generations"){
		input >> this->deal_generations;
	}
	else if (param == "deal_top_ranked"){
		input >> this->deal_top_ranked;
	}

	else if (param == "generate_conformers"){
		input >> tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->generate_conformers = true;
		}
		else {
			this->generate_conformers = false;
		}
        this->flex_lig = true;
	}

	else if (param == "number_of_conformers"){
		input >> this->lig_conformers;
	}

	else if (param == "conformers_to_rank"){
		input >> this->conformers_to_evaluate;
	}

	else if (param == "conformer_generator"){
		input >> tmp;
		if (tmp == "GA" or tmp == "ga" or tmp == "Ga"){
			this->conformer_generator = "GA";
		}
		else {
			this->conformer_generator = "WRS";
		}
	}

    else if (param == "conformer_min_steps"){
        input >> this->conformer_min_steps;
	}

	else if (param == "WRS_geom_steps"){
		input >> this->WRS_geom_steps;
	}

	else if (param == "dock_min_tol"){
		input >> this->dock_min_tol;
	}
	else if (param == "dock_parallel"){
		input >> tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->dock_parallel = true;
		}
	}
	else if (param == "parallel_jobs"){
		input >> (this->parallel_jobs);
	}
	else if (param == "write_mol2"){
		input >> tmp;
		if ( tmp == "no" or tmp == "No" or tmp == "NO"){
			this->write_mol2 = false;
		}
	}
	else if (param == "use_grids"){
		input >> tmp;
		if (tmp == "Yes" or tmp == "Yes" or tmp == "yes"){
			this->use_grids = true;
		}
	}
	else if (param == "grid_spacing"){
		input >> this->grid_spacing;
	}
	else if (param == "search_box"){
		input >> this->search_box_x >> this->search_box_y >> this->search_box_z;
	}
	else if (param == "load_grids"){
		this->load_grid_from_file = true;
		input >> this->grid_prefix;
	}
	else if (param == "write_grids"){
		this->write_grids = true;
		input >> this->grid_prefix;
	}
	else if (param == "compute_rmsd"){
		input >> tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->show_rmsd = true;
		}
	}
	else if (param == "sort_by_energy"){
		input >> tmp;
		if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
			this->sort_by_energy = true;
		}
	}
    else if (param == "dielectric_model"){
        input >> this->dielectric_model ;
    }
    else if (param == "only_score"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "Yes"){
            this->only_score = true;
        }
    }
    else if (param == "solvation_alpha"){
        input >> this->solvation_alpha;
    }
    else if (param == "solvation_beta"){
        input >> this->solvation_beta;
    }
    else if (param == "seed"){
        input >> this->seed;
    }
    else if (param == "ligand_simulation"){
        input >> tmp;
        if(tmp == "yes" or tmp == "1" or tmp == "y")
            this->ligsim = true;
    }
    else if (param == "equilibration_steps"){
        input >> this->eq_steps;
    }
    else if (param == "mc_stride"){
        input >> this->mc_stride;
    }
    else if (param == "mcr_size"){
        input >> this->mcr_size;
    }
    else if (param == "verbose"){
        input >> tmp;
        if (tmp == "yes" or tmp == "Yes" or tmp == "YES"){
            this->verbose = true;
        }
    }
    else if (param == "mcr_coefficients"){
        double bi;
        for (int i=0; i<mcr_size; i++){
            input >> bi;
            mcr_coefficients.push_back(bi);
        }
    }
    else if (param == "entropy_rotation_bins"){
        input >> this->entropy_rotation_bins;
    }

    else if (param == "entropy_translation_bins"){
        input >> this->entropy_translation_bins;
    }
    else if (param == "ligand_energy_model"){
        input >> this->ligand_energy_model;
    }
    else if (param == "atomic_model_ff"){
        input >> this->atomic_model_ff;
    }
    else if (param == "elec_weight"){
        input >> this->elec_weight;
    }
    else if (param == "vdw_weight"){
        input >> this->vdw_weight;
    }
    else if (param == "solvation_weight"){
        input >> this->solvation_weight;
    }
    else if (param == "pbsa_grid"){
        input >> this->pbsa_grid;
    }

	else {
		cout << "Unknown parameter: " << param << endl;
		exit(1);
	}
}


void PARSER::set_parameters(char* arg){

	ifstream input(arg);
	char line[256];
	while (!input.eof()){
		input >> param;
		if ((param.substr(0,1) == "#") or (param.substr(0,2) == "//")){
			input.getline(line, 256);
		}
		else {
			PARSER::comparing (param, input);
		}
	}
	this->check_parameters_sanity();
}

void PARSER::check_parameters_sanity(void){
	if (this->conformers_to_evaluate > this->lig_conformers){
		this->conformers_to_evaluate = this->lig_conformers;
	}
}


/*
 * Optimizer.cpp
 *
 *  Created on: 21/03/2012
 *      Author: Nascimento
 */

#include "Optimizer.h"

Optimizer::Optimizer(Mol2* _Rec, Mol2* _RefLig, PARSER* _Parser) {
	this->RefLig = _RefLig;
	this->Rec = _Rec;
	this->Parser = _Parser;
}
Optimizer::Optimizer(Mol2* _Rec, Mol2* _RefLig, PARSER* _Parser, Grid* _Grids) {
	this->RefLig = _RefLig;
	this->Rec = _Rec;
	this->Parser = _Parser;
	this->Grids = _Grids;
}

Optimizer::~Optimizer() {
}

vector<vector<double> > Optimizer::update_coords(const std::vector<double> &x, Mol2* Lig2){
	COORD_MC* Coord = new COORD_MC;
	return(Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]));
}

double Optimizer::evaluate_rmsd(Mol2* Lig1, Mol2* Lig2){
	double rmsd=0.0;
	for (int i=0; i<Lig1->N; i++){
		for (int j=0; j< 3; j++){
			rmsd += (Lig1->xyz[i][j] - Lig2->xyz[i][j]) * (Lig1->xyz[i][j] - Lig2->xyz[i][j]);
		}
	}
	rmsd = sqrt(rmsd);
	return(rmsd);
}

double Optimizer::evaluate_energy(Mol2* Lig2, vector<vector<double> > new_xyz){
    Energy2* Ene = new Energy2(Parser);
    energy_result_t* energy_results = new energy_result_t;
    double energy=0.0;
    if (Parser->use_grids){
        energy = Ene->compute_ene(Grids, Lig2, new_xyz, energy_results);
    }
    else {
        energy = Ene->compute_ene(Rec, Lig2, new_xyz, energy_results);
    }
	delete Ene;
    delete energy_results;
	return(energy);
}

void Optimizer::evaluate_energy(Mol2* Lig2, vector<vector<double> > new_xyz, energy_result_t* energy_result){
    Energy2* Ene = new Energy2(Parser);
    if (Parser->use_grids){
        Ene->compute_ene(Grids, Lig2, new_xyz, energy_result);
    }
    else{
        Ene->compute_ene(Rec, Lig2, new_xyz, energy_result);
    }
}

double Optimizer::objective_energy_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
	vector<vector<double> > new_xyz;
	COORD_MC* Coord = new COORD_MC;
	double f, f2;

	Mol2* Lig2 = (Mol2*) data;
	new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]);
	f = evaluate_energy(Lig2, new_xyz);

    Gaussian* Gauss = new Gaussian;
    double t1, t2, t3, si;
    t1 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, RefLig->xyz, new_xyz);
    t2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, RefLig, RefLig->xyz, RefLig->xyz);
    t3 = Gauss->compute_shape_and_charge_density(Parser, Lig2, Lig2, new_xyz,new_xyz);
    si = (2*t1) / (t2+t3);
    delete Gauss;

    if(!grad.empty()){
		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0]+Parser->min_delta, x[1], x[2], x[3], x[4], x[5]);
		f2 = evaluate_energy(Lig2, new_xyz);
		grad[0] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1]+Parser->min_delta, x[2], x[3], x[4], x[5]);
		f2 = evaluate_energy(Lig2, new_xyz);
		grad[1] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2]+Parser->min_delta, x[3], x[4], x[5]);
		f2 = evaluate_energy(Lig2, new_xyz);
		grad[2] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3]+Parser->min_delta, x[4], x[5]);
		f2 = evaluate_energy(Lig2, new_xyz);
		grad[3] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4]+Parser->min_delta, x[5]);
		f2 = evaluate_energy(Lig2, new_xyz);
		grad[4] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]+Parser->min_delta);
		f2 = evaluate_energy(Lig2, new_xyz);
		grad[5] = (f2-f)/Parser->min_delta;
	}
	delete Coord;
//    return (si*f);
    return (f);
}

double Optimizer::objective_overlay_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
	Gaussian* Gauss = new Gaussian;
	COORD_MC* Coord = new COORD_MC;
	vector<vector<double> > new_xyz;
	double f, f2;
	Mol2* Lig2 = (Mol2*) data;

	new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]);

	f = (Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz));

	if (!grad.empty()){
		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0] + Parser->min_delta, x[1], x[2], x[3], x[4], x[5]);
		f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
		grad[0] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2,x[0], x[1] + Parser->min_delta, x[2], x[3], x[4], x[5]);
		f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
		grad[1] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2]+ Parser->min_delta, x[3], x[4], x[5]);
		f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
		grad[2] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3]+ Parser->min_delta, x[4], x[5]);
		f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
		grad[3] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4]+Parser->min_delta, x[5]);
		f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
		grad[4] = (f2-f)/Parser->min_delta;

		new_xyz = Coord->rototranslate(Lig2->xyz, Lig2, x[0], x[1], x[2], x[3], x[4], x[5]+Parser->min_delta);
		f2 = Gauss->compute_shape_and_charge_density(Parser, RefLig, Lig2, new_xyz);
		grad[5] = (f2-f)/Parser->min_delta;
	}
	delete Gauss;
	delete Coord;
	return (f);
}

void  Optimizer::minimize_energy_nlopt_cobyla(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LN_COBYLA,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);
	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void  Optimizer::minimize_energy_nlopt_lbfgs(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LD_LBFGS,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);
	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;

}

void  Optimizer::minimize_energy_nlopt_ln_auglag(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LN_AUGLAG,6);
    nlopt::opt local_opt(nlopt::LD_MMA,6);

    local_opt.set_xtol_rel(Parser->dock_min_tol);
    local_opt.set_maxtime(Parser->min_timeout);
    opt->set_local_optimizer(local_opt);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);
	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_energy_nlopt_ld_auglag(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LD_AUGLAG,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);
	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;

}

void  Optimizer::minimize_energy_nlopt_mma(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LD_MMA,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);

	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_subplex(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LN_SBPLX,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);

	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_simplex(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LN_NELDERMEAD,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);

	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_crs(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::GN_CRS2_LM,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);

	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_direct(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_DIRECT_L_RAND,6);

	vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
	lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
    nlopt::result nres;

    try {
        nres = opt->optimize(x,f_minimum);
    }
    catch(...){
        if (Parser->verbose){
            printf("DIRECT optimization of molecule %s stopped with an exception.\n", Lig2->molname.c_str());
        }
    }

    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);

    delete opt;

    Lig2->xyz = xyz;
    nlopt::opt *opt2 = new nlopt::opt(nlopt::LN_AUGLAG,6);
    opt2->set_lower_bounds(lb);
    opt2->set_upper_bounds(ub);

    opt2->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt2->set_xtol_rel(Parser->dock_min_tol);
    opt2->set_maxtime(Parser->min_timeout);
    vector<double> x2(6);
    x2[0] = 0.0;
    x2[1] = 0.0;
    x2[2] = 0.0;
    x2[3] = 0.0;
    x2[4] = 0.0;
    x2[5] = 0.0;

    f_minimum=0.00;
    nres = opt2->optimize(x2,f_minimum);
    xyz = Optimizer::update_coords(x2, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);

    if (f_minimum > 9999.9) {
        opt_result->f_min = 9999.9;
    }
    else {
        opt_result->f_min = f_minimum;
    }
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

    delete opt2;

}

void Optimizer::minimize_energy_nlopt_direct_only(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_DIRECT_L_RAND,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);

    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);
    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_energy_nlopt_stogo(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::GD_STOGO,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
	opt->set_xtol_rel(Parser->dock_min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double f_minimum=0.00;
	nlopt::result nres = opt->optimize(x,f_minimum);
	vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);

	delete opt;

	opt_result->f_min = f_minimum;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_isres(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_ISRES,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}

void Optimizer::minimize_energy_nlopt_esch(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_ESCH,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
    opt->set_xtol_rel(Parser->dock_min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double f_minimum=0.00;
    nlopt::result nres = opt->optimize(x,f_minimum);
    vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
    Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);

    delete opt;

    opt_result->f_min = f_minimum;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;

}



void Optimizer::minimize_overlay_nlopt_cobyla(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LN_COBYLA,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
	opt->set_xtol_rel(Parser->min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double fo;
	nlopt::result nres = opt->optimize(x,fo);
	vector<vector<double> > xyz = update_coords(x, Lig2);
	double f_minimum = evaluate_energy(Lig2, xyz);
	delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_lbfgs(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LD_LBFGS,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
	opt->set_xtol_rel(Parser->min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double fo;
	nlopt::result nres = opt->optimize(x,fo);
	vector<vector<double> > xyz = update_coords(x, Lig2);
	double f_minimum = evaluate_energy(Lig2, xyz);
	delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_ln_auglag(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LN_AUGLAG,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
	opt->set_xtol_rel(Parser->min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double fo;
	nlopt::result nres = opt->optimize(x,fo);
	vector<vector<double> > xyz = update_coords(x, Lig2);
	double f_minimum = evaluate_energy(Lig2, xyz);
	delete opt;

    opt_result->energy_result->total = f_minimum;
	opt_result->f_min = fo;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_ld_auglag(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LD_AUGLAG,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] =  -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
	opt->set_xtol_rel(Parser->min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double fo;
	nlopt::result nres = opt->optimize(x,fo);
	vector<vector<double> > xyz = update_coords(x, Lig2);
	double f_minimum = evaluate_energy(Lig2, xyz);
	delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_mma(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LD_MMA,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
	opt->set_xtol_rel(Parser->min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double fo;
	nlopt::result nres = opt->optimize(x,fo);
	vector<vector<double> > xyz = update_coords(x, Lig2);
	double f_minimum = evaluate_energy(Lig2, xyz);
	delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_subplex(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::LN_SBPLX,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
	opt->set_xtol_rel(Parser->min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double fo;
	nlopt::result nres = opt->optimize(x,fo);
	vector<vector<double> > xyz = update_coords(x, Lig2);
	double f_minimum = evaluate_energy(Lig2, xyz);
	delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_crs(Mol2* Lig2, opt_result_t* opt_result){
	nlopt::opt *opt = new nlopt::opt(nlopt::GN_CRS2_LM,6);

	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = -(Parser->search_box_x/2.0);
	lb[4] = -(Parser->search_box_y/2.0);
	lb[5] = -(Parser->search_box_z/2.0);
	vector<double> ub(6);
	ub[0] = 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Parser->search_box_x/2.0;
	ub[4] = Parser->search_box_y/2.0;
	ub[5] = Parser->search_box_z/2.0;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
	opt->set_xtol_rel(Parser->min_tol);
	opt->set_maxtime(Parser->min_timeout);

	vector<double> x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;

	double fo;
	nlopt::result nres = opt->optimize(x,fo);
	vector<vector<double> > xyz = update_coords(x, Lig2);
	double f_minimum = evaluate_energy(Lig2, xyz);
	delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
	opt_result->optimization_status = nres;
	opt_result->optimized_xyz = xyz;
}

void Optimizer::minimize_overlay_nlopt_direct(Mol2* Lig2, opt_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::GN_DIRECT_L_RAND,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::objective_overlay_function, Lig2);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    vector<vector<double> > xyz = update_coords(x, Lig2);
    double f_minimum = evaluate_energy(Lig2, xyz);
    delete opt;

    opt_result->energy_result->total = f_minimum;
    opt_result->f_min = fo;
    opt_result->optimization_status = nres;
    opt_result->optimized_xyz = xyz;
}

void Optimizer::Simulated_Annealing(Mol2* Lig, opt_result_t* opt_result){
	gsl_rng* r;
	r = gsl_rng_alloc(gsl_rng_ranlxs2);
	srand(time(NULL));
	gsl_rng_set(r, (rand() % 100));
	vector<vector<double> > new_xyz;

	SA* Anneal = new SA();
	new_xyz = Anneal->optimize(Lig, Parser, Grids, r);
	delete Anneal;
    evaluate_energy(Lig, new_xyz, opt_result->energy_result);

    opt_result->f_min = opt_result->energy_result->total;
    opt_result->optimized_xyz = new_xyz;
	opt_result->optimization_status = 0;
	gsl_rng_free(r);
}

void Optimizer::minimize_energy_adaptative(Mol2* Lig2, opt_result_t* opt_result){
    double tol = 1.0e-4;
    double deltaij = 2.75;
    double deltaij_es = 1.65;
    Parser->deltaij6 = pow(deltaij,6);
    Parser->deltaij_es6 = pow(deltaij_es,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] =  -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    int count = 1;
    while (deltaij >= 0.0){
        cout << "Adaptative cycle " << count << " ..." << endl;
        cout << "deltaij = " << deltaij << endl;
        cout << "deltaij_es = " << deltaij_es << endl;

        nlopt::opt *opt = new nlopt::opt(nlopt::LD_AUGLAG,6);
        opt->set_lower_bounds(lb);
        opt->set_upper_bounds(ub);

        opt->set_min_objective(Optimizer::objective_energy_function, Lig2);
        opt->set_xtol_rel(tol);
        opt->set_maxtime(Parser->min_timeout);

        vector<double> x(6);
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = 0.0;
        x[4] = 0.0;
        x[5] = 0.0;

        double f_minimum=0.00;
        nlopt::result nres = opt->optimize(x,f_minimum);
        vector<vector<double> > xyz = Optimizer::update_coords(x, Lig2);
        Optimizer::evaluate_energy(Lig2, xyz, opt_result->energy_result);
        delete opt;

        opt_result->f_min = f_minimum;
        opt_result->optimization_status = nres;
        opt_result->optimized_xyz = xyz;


        deltaij -= 0.25;
        deltaij_es -=0.15;

        Parser->deltaij6 = pow(deltaij, 6);
        Parser->deltaij_es6 = pow(deltaij_es, 6);
        count++;
    }
}

double Optimizer::superpose_function(const std::vector<double> &x, std::vector<double> &grad, void *data){
    COORD_MC* Coord = new COORD_MC;
    vector<vector<double> > new_xyz;
    double f, f2;
    align_t* align_data= (align_t*) data;

    new_xyz = Coord->rototranslate(align_data->current_xyz, int(align_data->current_xyz.size()), x[0], x[1], x[2], x[3], x[4], x[5]);

    f = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));

    if (!grad.empty()){
        new_xyz = Coord->rototranslate(align_data->current_xyz, int(align_data->current_xyz.size()), x[0] + Parser->min_delta, x[1], x[2], x[3], x[4], x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[0] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, int(align_data->current_xyz.size()), x[0], x[1] + Parser->min_delta, x[2], x[3], x[4], x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[1] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, int(align_data->current_xyz.size()), x[0], x[1], x[2]+ Parser->min_delta, x[3], x[4], x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[2] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, int(align_data->current_xyz.size()), x[0], x[1], x[2], x[3]+ Parser->min_delta, x[4], x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[3] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, int(align_data->current_xyz.size()), x[0], x[1], x[2], x[3], x[4]+Parser->min_delta, x[5]);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[4] = (f2-f)/Parser->min_delta;

        new_xyz = Coord->rototranslate(align_data->current_xyz, int(align_data->current_xyz.size()), x[0], x[1], x[2], x[3], x[4], x[5]+Parser->min_delta);
        f2 = Coord->compute_rmsd(align_data->ref_xyz, new_xyz, int(align_data->ref_xyz.size()));
        grad[5] = (f2-f)/Parser->min_delta;
    }
    delete Coord;
    return (f);
}

void Optimizer::minimize_alignment_nlopt_ln_auglag(align_t* align_data, align_result_t* opt_result){
    nlopt::opt *opt = new nlopt::opt(nlopt::LN_AUGLAG,6);

    vector<double> lb(6);
    lb[0] = -180.0;
    lb[1] = -90.0;
    lb[2] = -180.0;
    lb[3] = -(Parser->search_box_x/2.0);
    lb[4] = -(Parser->search_box_y/2.0);
    lb[5] = -(Parser->search_box_z/2.0);
    vector<double> ub(6);
    ub[0] = 180.0;
    ub[1] = 90.0;
    ub[2] = 180.0;
    ub[3] = Parser->search_box_x/2.0;
    ub[4] = Parser->search_box_y/2.0;
    ub[5] = Parser->search_box_z/2.0;

    opt->set_lower_bounds(lb);
    opt->set_upper_bounds(ub);

    opt->set_max_objective(Optimizer::superpose_function, align_data);
    opt->set_xtol_rel(Parser->min_tol);
    opt->set_maxtime(Parser->min_timeout);

    vector<double> x(6);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;

    double fo;
    nlopt::result nres = opt->optimize(x,fo);
    delete opt;

    opt_result->rmsd = fo;
    opt_result->translation.push_back(x[3]);
    opt_result->translation.push_back(x[4]);
    opt_result->translation.push_back(x[5]);
    opt_result->rotation.push_back(x[0]);
    opt_result->rotation.push_back(x[1]);
    opt_result->rotation.push_back(x[2]);
}

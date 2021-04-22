#include "DockWidget.h"

DockWidget::DockWidget(PARSER *Input, QPlainTextEdit* Ed)
{

	this->Inp = Input;
	this->Editor = Ed;

    mainLayout = new QVBoxLayout(this);
	editLayout = new QGridLayout;
	buttonLayout = new QHBoxLayout;
	messageLayout = new QHBoxLayout;
	progress_layout = new QVBoxLayout;

	mainLayout->addLayout(editLayout);
	mainLayout->addStretch();
	mainLayout->addLayout(buttonLayout);
    mainLayout->setSizeConstraint(QLayout::SetFixedSize);

	rec_mol2 = new QPushButton(tr("Receptor MOL2 file:"));
	choose_rec_mol2 = new QLineEdit(tr("Choose file..."));

	lig_mol2 = new QPushButton(tr("Ligand MOL2 file:"));
	choose_lig_mol2 = new QLineEdit(tr("Choose file..."));

	reflig_mol2 = new QPushButton(tr("Reference MOl2 file:"));
	choose_reflig_mol2 = new QLineEdit(tr("Choose file..."));

	docking_mol2 = new QPushButton(tr("Docking Molecules (MOL2 File): "));
	choose_docking_mol2 = new QLineEdit(tr("Choose files..."));

	mol2_aa = new QComboBox;
	mol2_aa->addItem("False"); // 0
	mol2_aa->addItem("True");  // 1
	mol2_aa->setCurrentIndex(0);
	mol2_aa_lab = new QLabel(tr("MOL2 File has Amber Atom Types? "));

	scoring_function = new QComboBox;
    scoring_function->addItem("Amber Softcore + Desolvation");        //0
    scoring_function->addItem("Amber Softcore");                      //1
    scoring_function->addItem("Amber FF + Desolvation");              //2
    scoring_function->addItem("Amber FF");                            //3
    scoring_function->addItem("GaussSmooth Amber FF + Desolvation");  //4
    scoring_function->addItem("GaussSmooth Amber FF");                //5
    scoring_function->setCurrentIndex(2);
	scoring_function_lab = new QLabel(tr("Scoring function: "));

	diel = new QDoubleSpinBox;
	diel->setMaximum(80.0);
	diel->setMinimum(1.0);
	diel->setDecimals(1);
	diel->setSingleStep(1.0);
    diel->setValue(Inp->diel);
	diel_lab = new QLabel(tr("Dieletric constant: "));

    dielectric_model = new QComboBox;
    dielectric_model->addItem("constant");          //0
    dielectric_model->addItem("r");                 //1
    dielectric_model->setCurrentIndex(1);
    dielectric_model_lab = new QLabel(tr("Dielectric Model: "));

	sigma = new QDoubleSpinBox;
	sigma->setDecimals(1);
	sigma->setSingleStep(0.1);
    sigma->setValue(Inp->sigma);
	sigma->setMinimum(1.0);
    sigma_lab = new QLabel(tr("Solvation sigma: "));

    alpha = new QDoubleSpinBox;
    alpha->setDecimals(2);
    alpha->setSingleStep(0.01);
    alpha->setValue(Inp->solvation_alpha);
    alpha->setMinimum(0.01);
    alpha_lab = new QLabel(tr("Alpha Solvation Parameter: "));

    beta = new QDoubleSpinBox;
    beta->setDecimals(3);
    beta->setMinimum(-10.000);
    beta->setSingleStep(0.001);
    beta->setValue(Inp->solvation_beta);
    beta_lab = new QLabel(tr("Beta Solvation Parameter: "));

	deltaij_vdw = new QDoubleSpinBox;
	deltaij_vdw->setMaximum(10.0);
	deltaij_vdw->setMinimum(0.0);
    deltaij_vdw->setValue(1.75);
    deltaij_vdw->setSingleStep(0.05);
	deltaij_vdw->setDecimals(2);
	deltaij_vdw_lab = new QLabel(tr("Softcore delta for VDW term: "));

	deltaij_elec = new QDoubleSpinBox;
	deltaij_elec->setMaximum(10.0);
	deltaij_elec->setMinimum(0.0);
    deltaij_elec->setValue(1.50);
    deltaij_elec->setSingleStep(0.05);
	deltaij_elec->setDecimals(2);
	deltaij_elec_lab = new QLabel(tr("Softcore delta for electrostatic term: "));

	vdw_scale = new QDoubleSpinBox;
	vdw_scale->setMaximum(100.0);
	vdw_scale->setMinimum(0.0);
	vdw_scale->setDecimals(1);
	vdw_scale->setSingleStep(1.0);
    vdw_scale->setValue(Inp->vdw_scale);
	vdw_scale_lab = new QLabel(tr("Scale for VDW term in overlay optimization: "));

	elec_scale = new QDoubleSpinBox;
	elec_scale->setMaximum(100.0);
	elec_scale->setMinimum(0.0);
	elec_scale->setDecimals(1);
	elec_scale->setSingleStep(1.0);
    elec_scale->setValue(Inp->elec_scale);
	elec_scale_lab = new QLabel(tr("Scale for electrostatic term in overlay optimization: "));

	box_size = new QDoubleSpinBox;
    box_size->setMinimum(0.1);
	box_size->setMaximum(50.0);
    box_size->setSingleStep(0.1);
	box_size->setDecimals(1);
    box_size->setValue(Inp->search_box_x);
    box_size_lab = new QLabel(tr("Search Box : "));

	min_delta = new QDoubleSpinBox;
	min_delta->setMinimum(1E-10);
	min_delta->setMaximum(1.0);
	min_delta->setDecimals(10);
	min_delta->setSingleStep(0.0000000001);
    min_delta->setValue(Inp->min_delta);
	min_delta_lab = new QLabel(tr("Minimization step: "));

	min_tol = new QDoubleSpinBox;
    min_tol->setMinimum(1e-10);
	min_tol->setMaximum(1.0);
    min_tol->setDecimals(10);
    min_tol->setSingleStep(0.0000000001);
    min_tol->setValue(Inp->min_tol);
	min_tol_lab = new QLabel(tr("Minimization tolerance: "));

	min_timeout = new QSpinBox;
	min_timeout->setMinimum(1);
	min_timeout->setMaximum(120);
	min_timeout->setSingleStep(1);
    min_timeout->setValue(Inp->min_timeout);
	min_timeout_lab = new QLabel(tr("Minimization timeout: "));

	overlay_optimizer = new QComboBox;
	overlay_optimizer->addItem("COBYLA");					//0
	overlay_optimizer->addItem("AUGLAG (Gradient-free)");	//1
    overlay_optimizer->addItem("AUGLAG");                   //2
    overlay_optimizer->addItem("LBFGS");					//3
    overlay_optimizer->addItem("MMA");						//4
    overlay_optimizer->addItem("SUBPLEX");                  //5
    overlay_optimizer->addItem("DIRECT");                   //6
    overlay_optimizer->setCurrentIndex(4);
	overlay_optimizer_lab = new QLabel(tr("Overlay optimizer: "));

	energy_optimizer = new QComboBox;
    energy_optimizer->addItem("COBYLA");                    //0
    energy_optimizer->addItem("AUGLAG (Gradient-free)");    //1
    energy_optimizer->addItem("AUGLAG");                    //2
    energy_optimizer->addItem("LBFGS");                     //3
    energy_optimizer->addItem("MMA");                       //4
    energy_optimizer->addItem("SUBPLEX");                   //5
    energy_optimizer->addItem("SIMPLEX");                   //6
    energy_optimizer->addItem("DIRECT");                    //7
    energy_optimizer->addItem("CRS");                       //8
    energy_optimizer->addItem("STOGO");                     //9
    energy_optimizer->addItem("ISRES");                     //10
    energy_optimizer->addItem("ESCH");                      //11
    energy_optimizer->addItem("NONE");                      //12
    energy_optimizer->setCurrentIndex(7);
	energy_optimizer_lab = new QLabel(tr("Binding energy optimizer: "));

	output_prefix = new QLineEdit(tr("iMcLiBELa"));
	output_prefix_lab = new QLabel(tr("Output prefix: "));

    write_mol2 = new QCheckBox (tr("Write Mol2 file? "), this);
    write_mol2->setChecked(Inp->write_mol2);

    use_grids = new QCheckBox(tr("Use grids?"), this);
    use_grids->setChecked(Inp->use_grids);

    grid_spacing = new QDoubleSpinBox;
    grid_spacing->setValue(Inp->grid_spacing);
    grid_spacing->setDecimals(3);
    grid_spacing->setMinimum(0.01);
    grid_spacing->setMaximum(0.5);
    grid_spacing_lab = new QLabel(tr("Grid spacing: "));

    load_write_file = new QComboBox;
    load_write_file->addItem("Choose an option: ");
    load_write_file->addItem("Load File");
    load_write_file->addItem("Write File");

    grid_file = new QLineEdit;
    grid_file->setText("McGrid.grid");

    grid_box = new QDoubleSpinBox;
    grid_box->setValue(30.0);
    grid_box->setDecimals(1);
    grid_box->setMinimum(15.0);
    grid_box->setSingleStep(1.0);
    grid_box_lab = new QLabel(tr("Grid box: "));

	progressbar = new QProgressBar;
	progressbar->setRange(0, 100);
	progressbar->setValue(0);
	progress_label = new QLabel(tr("<h2>Progress</h2>"));
	progress_label->setAlignment(Qt::AlignHCenter);

    advanced_settings = new QCheckBox(tr("&Advanced settings"), this);

    dock_no_h = new QCheckBox(tr("Ignore hydrogens"), this);
	dock_no_h->setChecked(false);

    generate_conformers = new QCheckBox (tr("Generate conformers? "), this);
    generate_conformers->setChecked(Inp->generate_conformers);

    conformers_lab = new QLabel(tr("Number of conformers: "));
    conformers = new QSpinBox;
    conformers->setValue(Inp->lig_conformers);
    conformers->setSingleStep(1);
    conformers->setMinimum(1);

    conformers_to_evaluate_lab = new QLabel(tr("Conformers to optimize: "));
    conformers_to_evaluate = new QSpinBox;
    conformers_to_evaluate->setValue(Inp->conformers_to_evaluate);
    conformers_to_evaluate->setMinimum(1);

    ligand_energy_model_lab = new QLabel(tr("Ligand Energy Model: "));
    ligand_energy_model = new QComboBox;
    ligand_energy_model->addItem("GAFF");
    ligand_energy_model->addItem("MMFF94");

    dock_min_tol_lab = new QLabel(tr("Docking tolerance: "));
    dock_min_tol = new QDoubleSpinBox;
    dock_min_tol->setDecimals(10);
    dock_min_tol->setValue(Inp->dock_min_tol);


    dock_timeout_lab = new QLabel(tr("Docking timeout: "));
    dock_timeout = new QSpinBox;
    dock_timeout->setValue(Inp->min_timeout);

    parallel = new QCheckBox (tr("Parallel execution? "), this);

    parallel_jobs_lab = new QLabel(tr("Parallel jobs: "));
    parallel_jobs = new QSpinBox;
    parallel_jobs->setValue(1);
    parallel_jobs->setMinimum(1);

    use_cutoff_score = new QCheckBox(tr("Use cutoff Score?"));
    use_cutoff_score_lab = new QLabel(tr("Use cutoff for similarity?"));
    cutoff_score = new QDoubleSpinBox;
    cutoff_score->setValue(Inp->writeMol2_score_cutoff);
    cutoff_score->setMaximum(1.0);
    cutoff_score->setMinimum(0.0);
    cutoff_score->setSingleStep(0.05);

    use_cutoff_energy = new QCheckBox(tr("Use cutoff energy"));
    use_cutoff_energy_lab = new QLabel(tr("Use cutoff for binding energy?"));
    cutoff_energy = new QDoubleSpinBox;
    cutoff_energy->setValue(Inp->writeMol2_energy_cutoff);
    cutoff_energy->setSingleStep(0.01);

//    use_docking_restraint_lab = new QLabel(tr("Use harmonic restraints in docking ?"));
    use_docking_restraint = new QCheckBox(tr("Use harmonic restraints in docking ?"));
    docking_restraint_lab = new QLabel(tr("Docking restraint weight"));
    docking_restraint_weight = new QDoubleSpinBox;
    docking_restraint_weight->setValue(Inp->restraints_weight);
    docking_restraint_weight->setSingleStep(0.1);

    scale_ele_energy_lab = new QLabel(tr("Scale for Electrostatic Energy Term:"));
    scale_vdw_energy_lab = new QLabel(tr("Scale for VDW Energy Term:"));
    scale_ele_energy = new QDoubleSpinBox;
    scale_ele_energy->setValue(1.0);
    scale_ele_energy->setSingleStep(0.1);
    scale_vdw_energy = new QDoubleSpinBox;
    scale_vdw_energy->setValue(1.0);
    scale_vdw_energy->setSingleStep(0.1);



    files_label = new QLabel(tr("<h2>Input Files</h2>"));

    int i=0;

    editLayout->addWidget(files_label, i, 0);
    i++;

    editLayout->addWidget(rec_mol2, i, 0);
    editLayout->addWidget(choose_rec_mol2, i, 1);
    i++;

    editLayout->addWidget(reflig_mol2, i, 0);
    editLayout->addWidget(choose_reflig_mol2, i, 1);
    i++;

    editLayout->addWidget(docking_mol2, i, 0);
    editLayout->addWidget(choose_docking_mol2, i, 1);
    i++;

    editLayout->addWidget(mol2_aa_lab, i, 0);
    editLayout->addWidget(mol2_aa, i, 1);
    i++;

    output_label = new QLabel(tr("<h2>Output Files</h2>"));
    editLayout->addWidget(output_label, i, 0);
    i++;

    editLayout->addWidget(output_prefix_lab, i, 0);
    editLayout->addWidget(output_prefix, i, 1);
    i++;

    editLayout->addWidget(write_mol2, i,0);
    i++;

    editLayout->addWidget(use_cutoff_score_lab, i, 0);
    editLayout->addWidget(use_cutoff_score, i, 1);
    editLayout->addWidget(cutoff_score, i, 2);
    i++;

    editLayout->addWidget(use_cutoff_energy_lab, i, 0);
    editLayout->addWidget(use_cutoff_energy, i, 1);
    editLayout->addWidget(cutoff_energy, i, 2);
    i++;

    sf_input = new QLabel(tr("<h2>Scoring options</h2>"));
    editLayout->addWidget(sf_input, i, 0);
    i++;

    editLayout->addWidget(scoring_function_lab, i, 0);
    editLayout->addWidget(scoring_function, i, 1);
    i++;

    editLayout->addWidget(box_size_lab, i, 0);
    editLayout->addWidget(box_size, i, 1);
    i++;

    editLayout->addWidget(diel_lab, i, 0);
    editLayout->addWidget(diel, i, 1);
    i++;

    editLayout->addWidget(dielectric_model_lab, i, 0);
    editLayout->addWidget(dielectric_model, i, 1);
    i++;

    editLayout->addWidget(sigma_lab, i, 0);
    editLayout->addWidget(sigma, i, 1);
    i++;

    editLayout->addWidget(alpha_lab, i, 0);
    editLayout->addWidget(alpha, i, 1);
    i++;

    editLayout->addWidget(beta_lab, i, 0);
    editLayout->addWidget(beta, i, 1);
    i++;

    editLayout->addWidget(deltaij_vdw_lab, i, 0);
    editLayout->addWidget(deltaij_vdw, i, 1);
    i++;

    editLayout->addWidget(deltaij_elec_lab, i, 0);
    editLayout->addWidget(deltaij_elec, i, 1);
    i++;

    editLayout->addWidget(vdw_scale_lab, i, 0);
    editLayout->addWidget(vdw_scale, i, 1);
    i++;

    editLayout->addWidget(elec_scale_lab, i, 0);
    editLayout->addWidget(elec_scale, i, 1);
    i++;

    editLayout->addWidget(scale_ele_energy_lab, i, 0);
    editLayout->addWidget(scale_ele_energy, i, 1);
    i++;

    editLayout->addWidget(scale_vdw_energy_lab, i, 0);
    editLayout->addWidget(scale_vdw_energy, i, 1);
    i++;

    opt_input = new QLabel(tr("<h2>Optimization options</h2>"));
    editLayout->addWidget(opt_input, i, 0);
    i++;

    editLayout->addWidget(min_delta_lab, i, 0);
    editLayout->addWidget(min_delta, i, 1);
    i++;

    editLayout->addWidget(min_tol_lab, i, 0);
    editLayout->addWidget(min_tol, i, 1);
    i++;

    editLayout->addWidget(dock_min_tol_lab, i, 0);
    editLayout->addWidget(dock_min_tol, i, 1);
    i++;

    editLayout->addWidget(min_timeout_lab, i, 0);
    editLayout->addWidget(min_timeout, i, 1);
    i++;

    editLayout->addWidget(dock_timeout_lab, i, 0);
    editLayout->addWidget(dock_timeout, i, 1);
    i++;

    editLayout->addWidget(overlay_optimizer_lab, i, 0);
    editLayout->addWidget(overlay_optimizer, i, 1);
    i++;

    editLayout->addWidget(energy_optimizer_lab, i, 0);
    editLayout->addWidget(energy_optimizer, i, 1);
    i++;

    editLayout->addWidget(box_size_lab, i, 0);
    editLayout->addWidget(box_size, i, 1);
    i++;

    editLayout->addWidget(dock_no_h, i, 0);
    i++;

    editLayout->addWidget(use_docking_restraint, i, 0);
    i++;

    editLayout->addWidget(docking_restraint_lab, i, 0);
    editLayout->addWidget(docking_restraint_weight, i, 1);
    docking_restraint_weight->setDisabled(true);
    i++;

    conf_input = new QLabel(tr("<h2>Conformer options</h2>"));
    editLayout->addWidget(conf_input, i, 0);
    i++;

    editLayout->addWidget(generate_conformers, i, 0);
    i++;

    editLayout->addWidget(conformers_lab, i, 0);
    editLayout->addWidget(conformers, i, 1);
    i++;

    editLayout->addWidget(conformers_to_evaluate_lab, i, 0);
    editLayout->addWidget(conformers_to_evaluate, i, 1);
    i++;

    editLayout->addWidget(ligand_energy_model_lab, i, 0);
    editLayout->addWidget(ligand_energy_model, i, 1);
    i++;

    parallel_label = new QLabel(tr("<h2>Parallel execution</h2>"));
    editLayout->addWidget(parallel_label, i, 0);
    i++;

    editLayout->addWidget(parallel, i, 0);
    i++;

    editLayout->addWidget(parallel_jobs_lab, i, 0);
    editLayout->addWidget(parallel_jobs, i, 1);
    i++;

    grid_label = new QLabel(tr("<h2>Grid Potentials</h2>"));
    editLayout->addWidget(grid_label, i, 0);
    i++;

    editLayout->addWidget(use_grids, i,0);
    i++;

    editLayout->addWidget(grid_spacing_lab, i, 0);
    editLayout->addWidget(grid_spacing, i, 1);
    i++;

    editLayout->addWidget(load_write_file, i, 0);
    editLayout->addWidget(grid_file, i, 1);
    i++;

    editLayout->addWidget(grid_box_lab, i, 0);
    editLayout->addWidget(grid_box, i, 1);
    i++;

	progress_layout->addWidget(advanced_settings);
	progress_layout->addWidget(progress_label, Qt::AlignCenter);
	progress_layout->addWidget(progressbar);

	mainLayout->addStretch();
	mainLayout->addLayout(progress_layout);
	mainLayout->addStretch();

	Run = new QPushButton(tr("Run!"));
	WriteButton = new QPushButton(tr("Write parameter file"));

	buttonLayout->addStretch();
	buttonLayout->addWidget(WriteButton, Qt::AlignHCenter);
	buttonLayout->addWidget(Run, Qt::AlignHCenter);

	Run->setDefault(true);
	setWindowTitle(tr("iMcLiBEla by Alessandro Nascimento"));

	connect(rec_mol2, SIGNAL(clicked()), this, SLOT(choose_rec_file()));
	connect(reflig_mol2, SIGNAL(clicked()), this, SLOT(choose_reflig_file()));
	connect(docking_mol2, SIGNAL(clicked()), this, SLOT(choose_docking_mol2_files()));
	connect(advanced_settings, SIGNAL(toggled(bool)), this, SLOT(show_advanced_setttings(bool)));
	connect(Run, SIGNAL(clicked()), this, SLOT(Start()));
	connect(WriteButton, SIGNAL(clicked()), this, SLOT(write_parameters()));
    connect(scoring_function, SIGNAL(currentIndexChanged(int)), this, SLOT(show_scoring_function_options(int)));

    connect(write_mol2, SIGNAL(stateChanged(int)), this, SLOT(slot_write_mol2(int)));
    connect(dock_no_h, SIGNAL(stateChanged(int)), this, SLOT(slot_ignore_h(int)));
    connect(use_cutoff_score, SIGNAL(stateChanged(int)), this, SLOT(slot_use_cutoff_score(int)));
    connect(use_cutoff_energy, SIGNAL(stateChanged(int)), this, SLOT(slot_use_cutoff_energy(int)));
    connect(generate_conformers, SIGNAL(stateChanged(int)), this, SLOT(slot_generate_conformers(int)));
    connect(parallel, SIGNAL(stateChanged(int)), this, SLOT(slot_parallel(int)));
    connect(parallel_jobs, SIGNAL(valueChanged(int)), this, SLOT(slot_parallel(int)));

    connect(use_docking_restraint, SIGNAL(stateChanged(int)), this, SLOT(slot_use_restraints(int)));

    connect(scale_ele_energy, SIGNAL(valueChanged(double)), this, SLOT(slot_scale_ele(double)));
    connect(scale_vdw_energy, SIGNAL(valueChanged(double)), this, SLOT(slot_scale_vdw(double)));

    connect(load_write_file, SIGNAL(currentIndexChanged(int)), this, SLOT(slot_grid_file(int)));
    connect(use_grids, SIGNAL(stateChanged(int)), this, SLOT(slot_use_grids(int)));

    this->set_initial_parameters();
    this->hide_advanced_options();
}

void DockWidget::choose_rec_file(){
    QString filename = QFileDialog::getOpenFileName(this, tr("Choose a File"), "", tr("MOL2 Files (*.mol2 *.mol2.gz)"));
	choose_rec_mol2->setText(filename.toUtf8());
    Inp->rec_mol2 = string(filename.toStdString());
}

void DockWidget::choose_reflig_file(){
    QString filename = QFileDialog::getOpenFileName(this, tr("Choose a File"), "", tr("MOL2 Files (*.mol2 *.mol2.gz)"));
	choose_reflig_mol2->setText(filename.toUtf8());
    Inp->reflig_mol2 = string(filename.toUtf8().toStdString());
    Inp->lig_mol2 = Inp->reflig_mol2;
}

void DockWidget::choose_docking_mol2_files(){
    QStringList ql = QFileDialog::getOpenFileNames(this, tr("Choose MOL2 Files"), "", tr("MOL2 Files (*.mol2 *.mol2.gz)"));
	Inp->docking_molecules = ql;
	QString tmp = QString("%1 Molecules read!").arg(Inp->docking_molecules.size());
	choose_docking_mol2->setText(tmp);
}

void DockWidget::choose_grid_file(){
    QString filename = QFileDialog::getOpenFileName(this, tr("Choose a File"), "", tr("Grid Files (*.grid)"));
    QStringList list = filename.split(".");
    Inp->grid_prefix = list.at(0).toStdString();
    grid_file->setText(filename);
}

void DockWidget::slot_grid_file(int i){
    switch (i){
    case 1:
        Inp->load_grid_from_file = true;
        Inp->write_grids = false;
        this->choose_grid_file();
        break;
    case 2:
        Inp->write_grids = true;
        Inp->load_grid_from_file = false;
        QStringList prefix = this->grid_file->text().split(".");
        Inp->grid_prefix = prefix[0].toStdString();
    }
}

void DockWidget::show_scoring_function_options(int sf){
    if (sf > 1){
        deltaij_vdw_lab->hide();
        deltaij_vdw->hide();
        deltaij_elec_lab->hide();
        deltaij_elec->hide();
    }
    else {
        deltaij_vdw_lab->show();
        deltaij_vdw->show();
        deltaij_elec_lab->show();
        deltaij_elec->show();
    }
}

void DockWidget::show_advanced_setttings(bool state){
	if (state){
        this->show_advanced_options();
	}
	else {
        this->hide_advanced_options();
	}
}

void DockWidget::slot_write_mol2(int state){
    switch (state) {
    case Qt::Checked:
        Inp->write_mol2 = true;
        break;
    case Qt::Unchecked:
        Inp->write_mol2 = false;
    }
}

void DockWidget::slot_use_cutoff_score(int state){
    switch(state){
    case Qt::Checked:
        Inp->use_writeMol2_score_cutoff = true;
        break;
    case Qt::Unchecked:
        Inp->use_writeMol2_score_cutoff = false;
        break;
    }
}

void DockWidget::slot_use_cutoff_energy(int state){
    switch(state){
    case Qt::Checked:
        Inp->use_writeMol2_energy_cutoff = true;
        break;
    case Qt::Unchecked:
        Inp->use_writeMol2_energy_cutoff = false;
        break;
    }
}

void DockWidget::slot_ignore_h(int state){
    switch(state){
    case Qt::Checked:
        Inp->dock_no_h = true;
        break;
    case Qt::Unchecked:
        Inp->dock_no_h = false;
        break;
    }
}

void DockWidget::slot_generate_conformers(int state){
    switch(state){
    case Qt::Unchecked:
        Inp->generate_conformers = false;
        this->conformers_to_evaluate->setDisabled(true);
        this->conformers->setDisabled(true);
        this->conformers_lab->setDisabled(true);
        break;
    case Qt::Checked:
        Inp->generate_conformers = true;
        this->conformers_to_evaluate->setDisabled(false);
        this->conformers->setDisabled(false);
        this->conformers_lab->setDisabled(false);
        break;
    default:
        Inp->generate_conformers = false;
        break;
    }
}

void DockWidget::slot_parallel(int state){
    if (state > 1){
        Inp->dock_parallel = true;
        Inp->parallel_jobs = unsigned(this->parallel_jobs->value());
    }
    else {
        Inp->dock_parallel = false;
        Inp->parallel_jobs = 1;
    }
}

void DockWidget::slot_use_restraints(int state){
    switch(state){
    case Qt::Unchecked:
        Inp->use_Erestraints = false;
        this->docking_restraint_weight->setDisabled(true);
        break;

    case Qt::Checked:
        Inp->use_Erestraints = true;
        this->docking_restraint_weight->setEnabled(true);
        Inp->restraints_weight = this->docking_restraint_weight->value();
        break;
    }
}

void DockWidget::slot_use_grids(int state){
    switch (state){
    case Qt::Unchecked:
        Inp->use_grids = false;
        this->grid_file->setDisabled(true);
        this->grid_spacing->setDisabled(true);
        this->grid_spacing_lab->setDisabled(true);
        this->load_write_file->setDisabled(true);
        break;

    case Qt::Checked:
        Inp->use_grids = true;
        this->grid_file->setDisabled(false);
        this->grid_spacing->setDisabled(false);
        this->grid_spacing_lab->setDisabled(false);
        this->load_write_file->setDisabled(false);
        break;
    }
}

void DockWidget::slot_scale_ele(double value){
    Inp->scale_elec_energy = value;
}

void DockWidget::slot_scale_vdw(double value){
    Inp->scale_vdw_energy = value;
}

void DockWidget::Start(){
	this->set_parameters();
    RunEngine = new TEMP_SCHEME(Inp, Editor, this->progressbar);
    RunEngine->evaluation();
    delete RunEngine;
}

void DockWidget::set_parameters(){
	Inp->mode = "dock";
    Inp->dock_mode = true;
    Editor->appendPlainText("mode dock\n");

	if (mol2_aa->currentIndex() == 0){
		Inp->mol2_aa = false;
	}
	else if(mol2_aa->currentIndex() == 1){
		(Inp->mol2_aa = true);
	}

	Inp->scoring_function = scoring_function->currentIndex();

	Inp->diel = diel->value();
	Inp->sigma = sigma->value();
    Inp->solvation_alpha = this->alpha->value();
    Inp->solvation_beta = this->beta->value();

    switch (this->dielectric_model->currentIndex()){
    case 0:
        Inp->dielectric_model = "constant";
        break;
    case 1:
        Inp->dielectric_model = "r";
        break;
    }

    if (Inp->use_grids){
        Inp->grid_spacing = this->grid_spacing->value();
        Inp->x_dim = this->grid_box->value();
        Inp->y_dim = this->grid_box->value();
        Inp->z_dim = this->grid_box->value();
        switch(this->load_write_file->currentIndex()){
        case 1:
            Inp->load_grid_from_file = true;
            Inp->write_grids = false;
            break;
        case 2:
            Inp->load_grid_from_file = false;
            Inp->write_grids = true;
        }
    }

	Inp->deltaij6 = pow(deltaij_vdw->value(),6);
    Inp->deltaij_es3 = pow(deltaij_elec->value(), 3);
	Inp->deltaij_es6 = pow(deltaij_elec->value(), 6);

	Inp->elec_scale = elec_scale->value();
	Inp->vdw_scale = vdw_scale->value();

    Inp->search_box_x = box_size->value();
    Inp->search_box_y = box_size->value();
    Inp->search_box_z = box_size->value();

	Inp->min_delta = min_delta->value();
	Inp->min_tol = min_tol->value();

    Inp->min_timeout = 30;                              // WARNING
    Inp->dock_min_tol = this->dock_min_tol->value();
    Inp->min_timeout = this->dock_timeout->value();

	switch(overlay_optimizer->currentIndex()){
		case 0:
			Inp->overlay_optimizer = "cobyla";
			break;
		case 1:
            Inp->overlay_optimizer = "ln_auglag";
			break;
		case 2:
            Inp->overlay_optimizer = "ld_auglag";
			break;
		case 3:
            Inp->overlay_optimizer = "lbfgs";
			break;
		case 4:
            Inp->overlay_optimizer = "mma";
			break;
        case 5:
            Inp->overlay_optimizer = "subplex";
			break;
        case 6:
            Inp->overlay_optimizer = "direct";
            break;
	}

	switch(energy_optimizer->currentIndex()){
	case 0:
		Inp->energy_optimizer = "cobyla";
		break;
	case 1:
        Inp->energy_optimizer = "ln_auglag";
		break;
	case 2:
        Inp->energy_optimizer = "ld_auglag";
		break;
	case 3:
        Inp->energy_optimizer = "lbfgs";
		break;
	case 4:
        Inp->energy_optimizer = "mma";
		break;
	case 5:
        Inp->energy_optimizer = "subplex";
		break;
	case 6:
        Inp->energy_optimizer = "simplex";
		break;
    case 7:
        Inp->energy_optimizer = "direct";
		break;
    case 8:
        Inp->energy_optimizer = "crs";
        break;
    case 9:
        Inp->energy_optimizer = "stogo";
        break;
    case 10:
        Inp->energy_optimizer = "isres";
        break;
    case 11:
        Inp->energy_optimizer = "esch";
        break;
    case 12:
        Inp->energy_optimizer = "none";
        break;
	}

	Inp->output = output_prefix->text().toStdString();

	if (dock_no_h->isChecked()){
		Inp->dock_no_h = true;
	}
	else {
		Inp->dock_no_h = false;
	}
    if (Inp->generate_conformers){
        Inp->lig_conformers = this->conformers->value();
        Inp->conformers_to_evaluate = this->conformers_to_evaluate->value();
    }

    switch (ligand_energy_model->currentIndex()){
    case 0:
        Inp->ligand_energy_model = "GAFF";
        break;
    case 1:
        Inp->ligand_energy_model = "MMFF94";
        break;
    }

    if (this->use_cutoff_score->isChecked()){
        Inp->use_writeMol2_score_cutoff = true;
    }

    if (this->use_cutoff_energy->isChecked()){
        Inp->use_writeMol2_energy_cutoff = true;
    }

    Inp->writeMol2_score_cutoff = this->cutoff_score->value();
    Inp->writeMol2_energy_cutoff = this->cutoff_energy->value();

    Inp->scale_elec_energy = this->scale_ele_energy->value();
    Inp->scale_vdw_energy = this->scale_vdw_energy->value();

    // writing the list of molecules do dock

    FILE* multifile = fopen("multimol.dat", "w");
    for (int i=0; i< Inp->docking_molecules.size(); i++){
        fprintf(multifile, "%s\n", Inp->docking_molecules[i].toUtf8().toStdString().c_str());
    }
    fprintf(multifile, "%s\n", "EOF");
    fclose(multifile);
    Inp->multifile = "multimol.dat";
}

void DockWidget::write_parameters(){
	this->set_parameters();
	QString line;
	QFile input_params("iMcLiBELa.inp");

	if (!input_params.open(QIODevice::WriteOnly | QIODevice::Text)){
		printf("Failed to open iMcLiBELa.inp for writing.\n");
	}

	line = ("mode dock\n");
	input_params.write(line.toUtf8());

    /*
     * INPUT & OUTPUT FILES
    */

	line = ("rec_mol2 " + choose_rec_mol2->text() + "\n");
	input_params.write(line.toUtf8());

	line = ("lig_mol2 " + choose_reflig_mol2->text() + "\n");
	input_params.write(line.toUtf8());

	line = ("reflig_mol2 " + choose_reflig_mol2->text() + "\n");
	input_params.write(line.toUtf8());

    if (Inp->write_mol2){
        line = ("write_mol2 yes \n");
        input_params.write(line.toUtf8());
    }
    else {
        line = ("write_mol2 no \n");
        input_params.write(line.toUtf8());
    }

    line = ("output_prefix " + QString::fromStdString(Inp->output) + "\n");
    input_params.write(line.toUtf8());

    /*
     * SCORING FUNCTION
     */

    line = ("scoring_function " + QString::number(scoring_function->currentIndex()) + "\n");
    input_params.write(line.toUtf8());

    line = ("diel " + QString::number(diel->value()) + "\n");
    input_params.write(line.toUtf8());

    if (this->dielectric_model->currentIndex() == 0){
        line = ("dielectric_model constant \n");
        input_params.write(line.toUtf8());
    }
    else if (this->dielectric_model->currentIndex() == 1){
        line = ("dielectric_model r \n");
        input_params.write(line.toUtf8());
    }

    line = ("sigma " + QString::number(sigma->value()) + "\n");
    input_params.write(line.toUtf8());

    line = ("solvation_alpha " + QString::number(this->alpha->value()) + "\n");
    input_params.write(line.toUtf8());

    line = ("solvation_beta " + QString::number(this->beta->value()) + "\n");
    input_params.write(line.toUtf8());

	if (mol2_aa->currentIndex() == 0){
		line = ("mol2_aa no\n");
		input_params.write(line.toUtf8());
	}
	else if(mol2_aa->currentIndex() == 1){
		line = ("mol2_aa yes\n");
		input_params.write(line.toUtf8());
	}

	line = ("deltaij6 " + QString::number(pow(deltaij_vdw->value(), 6)) + "\n");
	input_params.write(line.toUtf8());

	line = ("deltaij_es6 " + QString::number(pow(deltaij_elec->value(), 6)) + "\n");
	input_params.write(line.toUtf8());

	line = ("elec_scale " + QString::number(elec_scale->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("vdw_scale " + QString::number(vdw_scale->value()) + "\n");
	input_params.write(line.toUtf8());

//	line = ("box " + QString::number(box_size->value()) + " " + QString::number(box_size->value()) + " " + QString::number(box_size->value()) + "\n");
//	input_params.write(line.toUtf8());

    /*
     * OPTIMIZATION
     */

	line = ("minimization_delta " + QString::number(min_delta->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("minimization_tolerance " + QString::number(min_tol->value()) + "\n");
	input_params.write(line.toUtf8());

    line = ("dock_min_tol " + QString::number(this->dock_min_tol->value()) + "\n");
    input_params.write(line.toUtf8());

    line = ("minimization_timeout " + QString::number(this->min_timeout->value()) + "\n");
    input_params.write(line.toUtf8());

	line = ("overlay_optimizer " + QString::fromStdString(Inp->overlay_optimizer) + "\n");
	input_params.write(line.toUtf8());

	line = ("energy_optimizer " + QString::fromStdString(Inp->energy_optimizer) + "\n");
	input_params.write(line.toUtf8());

	if(dock_no_h->isChecked()){
		line=("ignore_h yes\n");
		input_params.write(line.toUtf8());
	}
	else{
		line=("ignore_h no\n");
		input_params.write(line.toUtf8());
	}

    /*
     * CONFORMERS
     */

    if (this->generate_conformers->isChecked()){
        line=("generate_conformers yes\n");
        input_params.write(line.toUtf8());
    }
    else {
        line=("generate_conformers yes\n");
        input_params.write(line.toUtf8());
    }


    line = ("number_of_conformers " + QString::number(this->conformers->value()) + "\n");
    input_params.write(line.toUtf8());

    /*
     * GRIDS
     */

    if (this->use_grids->isChecked()){
        line=("use_grids yes\n");
        input_params.write(line.toUtf8());
    }
    else {
        line=("use_grids no\n");
        input_params.write(line.toUtf8());
    }

    line = ("grid_spacing " + QString::number(this->grid_spacing->value()) + "\n");
    input_params.write(line.toUtf8());

    switch (this->load_write_file->currentIndex()){
    case 1:
        line=(QString("load_grids %1\n").arg(this->grid_file->text()));
        input_params.write(line.toUtf8());
        break;
    case 2:
        line=(QString("write_grids %1\n").arg(this->grid_file->text()));
        input_params.write(line.toUtf8());
        break;
    }

    if (Inp->dock_parallel){
        line=(QString("dock_parallel %1\n").arg("yes"));
        input_params.write(line.toUtf8());
        line=(QString("parallel_jobs %1\n").arg(Inp->parallel_jobs));
        input_params.write(line.toUtf8());
    }
    else {
        line=(QString("dock_parallel %1\n").arg("no"));
        input_params.write(line.toUtf8());
    }


	QLabel* write_label = new QLabel(tr("Parameter file iMcLiBELa.inp written!"));
	editLayout->addWidget(write_label);

	input_params.close();
}

DockWidget::~DockWidget()
{
}

void DockWidget::hide_advanced_options(){
    deltaij_vdw_lab->hide();
    deltaij_vdw->hide();
    deltaij_elec_lab->hide();
    deltaij_elec->hide();

    diel_lab->hide();
    diel->hide();
    sigma_lab->hide();
    sigma->hide();

    elec_scale_lab->hide();
    elec_scale->hide();
    vdw_scale_lab->hide();
    vdw_scale->hide();
    box_size->hide();
    box_size_lab->hide();

    min_delta_lab->hide();
    min_delta->hide();
    min_tol_lab->hide();
    min_tol->hide();
    min_timeout->hide();
    min_timeout_lab->hide();

    overlay_optimizer_lab->hide();
    overlay_optimizer->hide();
    energy_optimizer_lab->hide();
    energy_optimizer->hide();

    dock_no_h->hide();
    write_mol2->hide();

    dielectric_model->hide();
    dielectric_model_lab->hide();

    alpha->hide();
    alpha_lab->hide();
    beta->hide();
    beta_lab->hide();

    conformers_to_evaluate->hide();
    conformers_to_evaluate_lab->hide();

    ligand_energy_model_lab->hide();
    ligand_energy_model->hide();

    dock_min_tol->hide();
    dock_min_tol_lab->hide();
    dock_timeout->hide();
    dock_timeout_lab->hide();

//    parallel->hide();
//    parallel_jobs->hide();
//    parallel_jobs_lab->hide();

    grid_box_lab->hide();
    grid_box->hide();

    opt_input->hide();

    use_docking_restraint->hide();
    docking_restraint_lab->hide();
    docking_restraint_weight->hide();

    scale_ele_energy_lab->hide();
    scale_vdw_energy_lab->hide();
    scale_ele_energy->hide();
    scale_vdw_energy->hide();

}

void DockWidget::show_advanced_options(){
//    deltaij_vdw_lab->show();
//    deltaij_vdw->show();
//    deltaij_elec_lab->show();
//    deltaij_elec->show();
    diel_lab->show();
    diel->show();
    sigma_lab->show();
    sigma->show();
    elec_scale_lab->show();
    elec_scale->show();
    vdw_scale_lab->show();
    vdw_scale->show();
    min_delta_lab->show();
    min_delta->show();
    min_tol_lab->show();
    min_tol->show();
    min_timeout->show();
    min_timeout_lab->show();
    overlay_optimizer_lab->show();
    overlay_optimizer->show();
    energy_optimizer_lab->show();
    energy_optimizer->show();
    dock_no_h->show();
    write_mol2->show();
    dielectric_model->show();
    dielectric_model_lab->show();
    alpha->show();
    alpha_lab->show();
    beta->show();
    beta_lab->show();

    box_size->show();
    box_size_lab->show();

    grid_box_lab->show();
    grid_box->show();

    generate_conformers->show();
    conformers->show();
    conformers_lab->show();
    conformers_to_evaluate->show();
    conformers_to_evaluate_lab->show();
    ligand_energy_model_lab->show();
    ligand_energy_model->show();
    dock_min_tol->show();
    dock_min_tol_lab->show();
    dock_timeout->show();
    dock_timeout_lab->show();
    parallel->show();
    parallel_jobs->show();
    parallel_jobs_lab->show();

    opt_input->show();

    use_docking_restraint->show();
    docking_restraint_lab->show();
    docking_restraint_weight->show();

    scale_ele_energy_lab->show();
    scale_vdw_energy_lab->show();
    scale_ele_energy->show();
    scale_vdw_energy->show();
}

void DockWidget::set_initial_parameters(){
    this->slot_use_grids(false);
    deltaij_vdw_lab->hide();
    deltaij_vdw->hide();
    deltaij_elec_lab->hide();
    deltaij_elec->hide();
}

#include "DockWidget.h"

DockWidget::DockWidget(PARSER *Input, QTextEdit* Ed)
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
    deltaij_vdw->setValue(2.75);
	deltaij_vdw->setSingleStep(0.25);
	deltaij_vdw->setDecimals(2);
	deltaij_vdw_lab = new QLabel(tr("Softcore delta for VDW term: "));

	deltaij_elec = new QDoubleSpinBox;
	deltaij_elec->setMaximum(10.0);
	deltaij_elec->setMinimum(0.0);
	deltaij_elec->setValue(1.75);
	deltaij_elec->setSingleStep(0.25);
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
    energy_optimizer->addItem("NONE");                      //9
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



    files_label = new QLabel(tr("<h2>Input Files</h2>"));
    editLayout->addWidget(files_label, 0, 0);

    editLayout->addWidget(rec_mol2, 1, 0);
    editLayout->addWidget(choose_rec_mol2, 1, 1);

    editLayout->addWidget(reflig_mol2, 2, 0);
    editLayout->addWidget(choose_reflig_mol2, 2, 1);

    editLayout->addWidget(docking_mol2, 3, 0);
    editLayout->addWidget(choose_docking_mol2, 3, 1);

    editLayout->addWidget(mol2_aa_lab, 4, 0);
    editLayout->addWidget(mol2_aa, 4, 1);

    output_label = new QLabel(tr("<h2>Output Files</h2>"));
    editLayout->addWidget(output_label, 5, 0);

    editLayout->addWidget(output_prefix_lab, 6, 0);
    editLayout->addWidget(output_prefix, 6, 1);

    editLayout->addWidget(write_mol2, 7,0);

    editLayout->addWidget(use_cutoff_score_lab, 8, 0);
    editLayout->addWidget(use_cutoff_score, 8, 1);
    editLayout->addWidget(cutoff_score, 8, 2);

    editLayout->addWidget(use_cutoff_energy_lab, 9, 0);
    editLayout->addWidget(use_cutoff_energy, 9, 1);
    editLayout->addWidget(cutoff_energy, 9, 2);

    sf_input = new QLabel(tr("<h2>Scoring options</h2>"));
    editLayout->addWidget(sf_input, 10, 0);

    editLayout->addWidget(scoring_function_lab, 11, 0);
    editLayout->addWidget(scoring_function, 11, 1);

    editLayout->addWidget(box_size_lab, 12, 0);
    editLayout->addWidget(box_size, 12, 1);

    editLayout->addWidget(diel_lab, 13, 0);
    editLayout->addWidget(diel, 13, 1);

    editLayout->addWidget(dielectric_model_lab, 14, 0);
    editLayout->addWidget(dielectric_model, 14, 1);

    editLayout->addWidget(sigma_lab, 15, 0);
    editLayout->addWidget(sigma, 15, 1);

    editLayout->addWidget(alpha_lab, 16, 0);
    editLayout->addWidget(alpha, 16, 1);

    editLayout->addWidget(beta_lab, 17, 0);
    editLayout->addWidget(beta, 17, 1);

    editLayout->addWidget(deltaij_vdw_lab, 18, 0);
    editLayout->addWidget(deltaij_vdw, 18, 1);

    editLayout->addWidget(deltaij_elec_lab, 19, 0);
    editLayout->addWidget(deltaij_elec, 19, 1);

    editLayout->addWidget(vdw_scale_lab, 20, 0);
    editLayout->addWidget(vdw_scale, 20, 1);

    editLayout->addWidget(elec_scale_lab, 21, 0);
    editLayout->addWidget(elec_scale, 21, 1);

    opt_input = new QLabel(tr("<h2>Optimization options</h2>"));
    editLayout->addWidget(opt_input, 22, 0);

    editLayout->addWidget(min_delta_lab, 23, 0);
    editLayout->addWidget(min_delta, 23, 1);

    editLayout->addWidget(min_tol_lab, 24, 0);
    editLayout->addWidget(min_tol, 24, 1);

    editLayout->addWidget(dock_min_tol_lab, 25, 0);
    editLayout->addWidget(dock_min_tol, 25, 1);

    editLayout->addWidget(min_timeout_lab, 26, 0);
    editLayout->addWidget(min_timeout, 26, 1);

    editLayout->addWidget(dock_timeout_lab, 27, 0);
    editLayout->addWidget(dock_timeout, 27, 1);

    editLayout->addWidget(overlay_optimizer_lab, 28, 0);
    editLayout->addWidget(overlay_optimizer, 28, 1);

    editLayout->addWidget(energy_optimizer_lab, 29, 0);
    editLayout->addWidget(energy_optimizer, 29, 1);

    editLayout->addWidget(box_size_lab, 30, 0);
    editLayout->addWidget(box_size, 30, 1);

    editLayout->addWidget(dock_no_h, 31, 0);

    conf_input = new QLabel(tr("<h2>Conformer options</h2>"));
    editLayout->addWidget(conf_input, 32, 0);

    editLayout->addWidget(generate_conformers, 33, 0);

    editLayout->addWidget(conformers_lab, 34, 0);
    editLayout->addWidget(conformers, 34, 1);

    editLayout->addWidget(conformers_to_evaluate_lab, 35, 0);
    editLayout->addWidget(conformers_to_evaluate, 35, 1);

    parallel_label = new QLabel(tr("<h2>Parallel execution</h2>"));
    editLayout->addWidget(parallel_label, 36, 0);

    editLayout->addWidget(parallel, 37, 0);

    editLayout->addWidget(parallel_jobs_lab, 38, 0);
    editLayout->addWidget(parallel_jobs, 38, 1);

    grid_label = new QLabel(tr("<h2>Grid Potentials</h2>"));
    editLayout->addWidget(grid_label, 39, 0);

    editLayout->addWidget(use_grids, 40,0);

    editLayout->addWidget(grid_spacing_lab, 41, 0);
    editLayout->addWidget(grid_spacing, 41, 1);

    editLayout->addWidget(load_write_file, 42, 0);
    editLayout->addWidget(grid_file, 42, 1);

    editLayout->addWidget(grid_box_lab, 43, 0);
    editLayout->addWidget(grid_box, 43, 1);

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
        Inp->parallel_jobs = state;
    }
    else {
        Inp->dock_parallel = false;
        Inp->parallel_jobs = 1;
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

void DockWidget::Start(){
	this->set_parameters();
    RunEngine = new TEMP_SCHEME(Inp, Editor, this->progressbar);
    RunEngine->evaluation();
    delete RunEngine;
}

void DockWidget::set_parameters(){
	Inp->mode = "dock";
    Inp->dock_mode = true;
    Editor->append("mode dock\n");

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

    if (this->use_cutoff_score->isChecked()){
        Inp->use_writeMol2_score_cutoff = true;
    }

    if (this->use_cutoff_energy->isChecked()){
        Inp->use_writeMol2_energy_cutoff = true;
    }

    Inp->writeMol2_score_cutoff = this->cutoff_score->value();
    Inp->writeMol2_energy_cutoff = this->cutoff_energy->value();

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
    dock_min_tol->show();
    dock_min_tol_lab->show();
    dock_timeout->show();
    dock_timeout_lab->show();
    parallel->show();
    parallel_jobs->show();
    parallel_jobs_lab->show();

    opt_input->show();
}

void DockWidget::set_initial_parameters(){
    this->slot_use_grids(false);
    deltaij_vdw_lab->hide();
    deltaij_vdw->hide();
    deltaij_elec_lab->hide();
    deltaij_elec->hide();
}

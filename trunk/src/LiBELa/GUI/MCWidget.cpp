#include "MCWidget.h"

MCWidget::MCWidget(PARSER* Input, QTextEdit* Editor)
{
	this->Inp = Input;
	this->Edit = Editor;

	mainLayout = new QVBoxLayout(this);
	editLayout = new QGridLayout;
	buttonLayout = new QHBoxLayout;
	messageLayout = new QHBoxLayout;
	progress_layout = new QVBoxLayout;
	mainLayout->addLayout(editLayout);
	mainLayout->addStretch();
	mainLayout->addLayout(buttonLayout);

	rec_mol2 = new QPushButton(tr("Receptor MOL2 file:"));
	if (Inp->rec_mol2 != ""){
		choose_rec_mol2 = new QLineEdit(QString::fromStdString(Inp->rec_mol2));
	}
	else {
		choose_rec_mol2 = new QLineEdit(tr("Choose file..."));
	}


	lig_mol2 = new QPushButton(tr("Ligand MOL2 file:"));
	if (Inp->lig_mol2 != ""){
		choose_lig_mol2 = new QLineEdit(QString::fromStdString(Inp->lig_mol2));
	}
	else {
		choose_lig_mol2 = new QLineEdit(tr("Choose file..."));
	}

	mol2_aa = new QComboBox;
	mol2_aa->addItem("False"); // 0
	mol2_aa->addItem("True");  // 1
	mol2_aa->setCurrentIndex(0);
	mol2_aa_lab = new QLabel(tr("MOL2 File has Amber Atom Types? "));

	temperature = new QDoubleSpinBox;
	temperature->setMinimum(1.0);
	temperature->setMaximum(1000000.0);
	temperature->setSingleStep(10.0);
	temperature->setValue(Inp->temp);
	temperature->setDecimals(1);
    temperature_lab = new QLabel(tr("MC Temperature: "));

	nsteps = new QSpinBox;
	nsteps->setMinimum(0);
	nsteps->setMaximum(10000000);
	nsteps->setValue(Inp->number_steps);
	nsteps_lab = new QLabel(tr("Number of steps: "));

	scoring_function = new QComboBox;
	scoring_function->addItem("Amber Softcore + Desolvation");
	scoring_function->addItem("Amber Softcore");
	scoring_function->addItem("Amber FF + Desolvation");
	scoring_function->addItem("Amber FF");
    scoring_function->addItem("Gaussian Smoothed Amber FF + Desolvation");
    scoring_function->addItem("Gaussian Smoothed Amber FF");
	scoring_function->setCurrentIndex(Inp->scoring_function);
	scoring_function_lab = new QLabel(tr("Scoring function: "));

	diel = new QDoubleSpinBox;
	diel->setMaximum(80.0);
	diel->setMinimum(1.0);
	diel->setDecimals(1);
	diel->setSingleStep(1.0);
	diel->setValue(Inp->diel);
	diel_lab = new QLabel(tr("Dieletric constant: "));

	sigma = new QDoubleSpinBox;
	sigma->setDecimals(1);
	sigma->setSingleStep(0.1);
	sigma->setValue(Inp->sigma);
	sigma->setMinimum(1.0);
	sigma_lab = new QLabel(tr("Gaussian sigma: "));


	deltaij_vdw = new QDoubleSpinBox;
	deltaij_vdw->setMaximum(10.0);
	deltaij_vdw->setMinimum(0.0);
	deltaij_vdw->setValue(pow(Inp->deltaij6, (1.0/6.0)));
	deltaij_vdw->setSingleStep(0.25);
	deltaij_vdw->setDecimals(2);
	deltaij_vdw_lab = new QLabel(tr("Softcore delta for VDW term: "));

	deltaij_elec = new QDoubleSpinBox;
	deltaij_elec->setMaximum(10.0);
	deltaij_elec->setMinimum(0.0);
	deltaij_elec->setValue(1.75);
	deltaij_elec->setSingleStep(pow(Inp->deltaij_es6, (1.0/6.0)));
	deltaij_elec->setDecimals(2);
	deltaij_elec_lab = new QLabel(tr("Softcore delta for electrostatic term: "));

	box_size = new QDoubleSpinBox;
	box_size->setMinimum(20.0);
	box_size->setMaximum(50.0);
	box_size->setSingleStep(1.0);
	box_size->setDecimals(1);
	box_size->setValue(Inp->x_dim);
	box_size_lab = new QLabel(tr("Box size (isotropic): "));

	cushion = new QDoubleSpinBox;
	cushion->setMinimum(0.0);
	cushion->setDecimals(1);
	cushion->setSingleStep(1.0);
	cushion->setValue(Inp->cushion);
	cushion_lab = new QLabel(tr("Translation step: "));

	rotation_step = new QDoubleSpinBox;
	rotation_step->setMinimum(1.0);
	rotation_step->setSingleStep(10.0);
	rotation_step->setValue(Inp->rotation_step);
	rotation_step_lab = new QLabel(tr("Rotation step (deg): "));

	timeout = new QSpinBox;
	timeout->setMinimum(60);
	timeout->setMaximum(1000000);
	timeout->setSingleStep(10);
	timeout->setValue(Inp->timeout);
	timeout_lab = new QLabel(tr("Timeout (sec): "));

	if (Inp->output != ""){
		output_prefix = new QLineEdit(QString::fromStdString(Inp->output));
	}
	else {
		output_prefix = new QLineEdit(tr("iMcLiBELa"));
	}

	output_prefix_lab = new QLabel(tr("Output prefix: "));

	progressbar = new QProgressBar;
	progressbar->setRange(0, 100);
	progressbar->setValue(0);
	progress_label = new QLabel(tr("<h2>Progress</h2>"));
	progress_label->setAlignment(Qt::AlignHCenter);

	sa = new QRadioButton(tr("Use &Simulated Annealing"));

	sa_start_temp = new QDoubleSpinBox;
	sa_start_temp->setMinimum(temperature->value());
	sa_start_temp->setMaximum(1000000.0);
	sa_start_temp->setSingleStep(10.0);
	sa_start_temp->setValue(Inp->sa_start_temp);
	sa_start_temp_lab = new QLabel(tr("Simulated annealing starting temperature: "));

	sa_steps = new QSpinBox;
	sa_steps->setMinimum(0);
	sa_steps->setMaximum(1000000);
	sa_steps->setSingleStep(10);
	sa_steps->setValue(Inp->sa_steps);
	sa_steps_lab = new QLabel(tr("Simulated annealing steps: "));

	sa_start_temp->hide();
	sa_start_temp_lab->hide();
	sa_steps->hide();
	sa_steps_lab->hide();

	input_files_lab = new QLabel(tr("<h3> <align = center> Input Files </h3>"));
	editLayout->addWidget(input_files_lab);
	editLayout->addWidget(rec_mol2, 1, 0);
	editLayout->addWidget(choose_rec_mol2, 1, 1);
	editLayout->addWidget(lig_mol2, 2, 0);
	editLayout->addWidget(choose_lig_mol2, 2, 1);
	editLayout->addWidget(mol2_aa_lab, 3, 0);
	editLayout->addWidget(mol2_aa, 3, 1);
	input_sf_lab = new QLabel(tr("<h3> <align = center> Scoring Function </h3>"));
	editLayout->addWidget(input_sf_lab);
	editLayout->addWidget(temperature_lab, 5, 0);
	editLayout->addWidget(temperature, 5, 1);
	editLayout->addWidget(scoring_function_lab, 6, 0);
	editLayout->addWidget(scoring_function, 6, 1);
	editLayout->addWidget(diel_lab, 7, 0);
	editLayout->addWidget(diel, 7, 1);
	editLayout->addWidget(sigma_lab, 8, 0);
	editLayout->addWidget(sigma, 8, 1);
	editLayout->addWidget(deltaij_vdw_lab, 9, 0);
	editLayout->addWidget(deltaij_vdw, 9, 1);
	editLayout->addWidget(deltaij_elec_lab, 10, 0);
	editLayout->addWidget(deltaij_elec, 10, 1);
	input_sampling_lab = new QLabel(tr("<h3> <align = center> Sampling </h3>"));
	editLayout->addWidget(input_sampling_lab);
	editLayout->addWidget(nsteps_lab, 12, 0);
	editLayout->addWidget(nsteps, 12, 1);
	editLayout->addWidget(box_size_lab, 13, 0);
	editLayout->addWidget(box_size, 13, 1);
	editLayout->addWidget(cushion_lab, 14, 0);
	editLayout->addWidget(cushion, 14, 1);
	editLayout->addWidget(box_size_lab, 15, 0);
	editLayout->addWidget(box_size, 15, 1);
	input_misc_lab = new QLabel(tr("<h3> <align = center > Miscellaneous </h3>"));
	editLayout->addWidget(input_misc_lab);
	editLayout->addWidget(timeout_lab, 17, 0);
	editLayout->addWidget(timeout, 17, 1);
	editLayout->addWidget(output_prefix_lab, 18, 0);
	editLayout->addWidget(output_prefix, 18, 1);
	editLayout->addWidget(sa, 19, 0);
	input_sa_lab = new QLabel(tr("<h3> <align = center > Simulated Annealing </h3>"));
	input_sa_lab->hide();
	editLayout->addWidget(input_sa_lab);
	editLayout->addWidget(sa_start_temp_lab, 20, 0);
	editLayout->addWidget(sa_start_temp, 20, 1);
	editLayout->addWidget(sa_steps_lab, 21, 0);
	editLayout->addWidget(sa_steps, 21, 1);

	Start = new QPushButton(tr("&Run"));
	WriteButton = new QPushButton(tr("&Write parameters"));
	buttonLayout->addWidget(WriteButton);
	buttonLayout->addWidget(Start);

	connect(sa, SIGNAL(toggled(bool)), this, SLOT(show_sa_parameters(bool)));
	connect(Start, SIGNAL(clicked()), this, SLOT(RunEng()));
	connect(WriteButton, SIGNAL(clicked()), this, SLOT(write_parameters()));
	connect(rec_mol2, SIGNAL(clicked()), this, SLOT(choose_rec_file()));
	connect(lig_mol2, SIGNAL(clicked()), this, SLOT(choose_lig_file()));
}

void MCWidget::show_sa_parameters(bool state){
	if (state){
		input_sa_lab->show();
		sa_start_temp_lab->show();
		sa_start_temp->show();
		sa_steps_lab->show();
		sa_steps->show();
		Inp->sa_scheme = true;
	}
	else {
		input_sa_lab->hide();
		sa_start_temp_lab->hide();
		sa_start_temp->hide();
		sa_steps_lab->hide();
		sa_steps->hide();
		Inp->sa_scheme = false;
	}
}

void MCWidget::choose_rec_file(){
	QString filename = QFileDialog::getOpenFileName(this, tr("Choose a File"), "", tr("MOL2 Files (*.mol2)"));
	choose_rec_mol2->setText(filename.toUtf8());
	Inp->rec_mol2 = filename.toStdString();
}

void MCWidget::choose_lig_file(){
	QString filename = QFileDialog::getOpenFileName(this, tr("Choose a File"), "", tr("MOL2 Files (*.mol2)"));
	choose_lig_mol2->setText(filename.toUtf8());
	Inp->reflig_mol2 = filename.toStdString();
	Inp->lig_mol2 = Inp->reflig_mol2;
}

void MCWidget::set_parameters(){
	Inp->mode = "mc";
	Inp->cushion = cushion->value();
	Inp->number_steps = nsteps->value();
	Inp->deltaij6 = pow(deltaij_vdw->value(), 6);
	Inp->deltaij_es6 = pow(deltaij_elec->value(), 6);
	Inp->diel = diel->value();
	Inp->lig_mol2 = (choose_lig_mol2->text().toStdString());
	Inp->rec_mol2 = (choose_rec_mol2->text().toStdString());
	if (mol2_aa->currentIndex() == 0){
		Inp->mol2_aa = false;
	}
	else if(mol2_aa->currentIndex() == 1){
		(Inp->mol2_aa = true);
	}
	Inp->temp = temperature->value();
	Inp->scoring_function = scoring_function->currentIndex();
	Inp->x_dim = box_size->value();
	Inp->y_dim = Inp->x_dim;
	Inp->z_dim = Inp->x_dim;
	Inp->rotation_step = rotation_step->value();
	Inp->timeout = timeout->value();
	Inp->output = (output_prefix->text().toStdString());
	if (Inp->sa_scheme){
		Inp->sa_start_temp = sa_start_temp->value();
		Inp->sa_steps = sa_steps->value();
	}
}

void MCWidget::write_parameters(){
	this->set_parameters();
	QString line;
	QFile input_params("iMcLiBELa.inp");

	if (!input_params.open(QIODevice::WriteOnly | QIODevice::Text)){
		printf("Failed to open iMcLiBELa.inp for writing.\n");
	}

	line = ("mode mc\n");
	input_params.write(line.toUtf8());

	line = ("cushion " + QString::number(cushion->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("temperature " + QString::number(temperature->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("rotation_step " + QString::number(rotation_step->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("timeout " + QString::number(timeout->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("rec_mol2 " + choose_rec_mol2->text() + "\n");
	input_params.write(line.toUtf8());

	line = ("lig_mol2 " + choose_lig_mol2->text() + "\n");
	input_params.write(line.toUtf8());

	if (mol2_aa->currentIndex() == 0){
		line = ("mol2_aa no\n");
		input_params.write(line.toUtf8());
	}
	else if(mol2_aa->currentIndex() == 1){
		line = ("mol2_aa yes\n");
		input_params.write(line.toUtf8());
	}

	line = ("scoring_function " + QString::number(scoring_function->currentIndex()) + "\n");
	input_params.write(line.toUtf8());

	line = ("diel " + QString::number(diel->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("sigma " + QString::number(sigma->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("deltaij6 " + QString::number(pow(deltaij_vdw->value(), 6)) + "\n");
	input_params.write(line.toUtf8());

	line = ("deltaij_es6 " + QString::number(pow(deltaij_elec->value(), 6)) + "\n");
	input_params.write(line.toUtf8());

	line = ("box_size " + QString::number(box_size->value()) + " " + QString::number(box_size->value()) + " " + QString::number(box_size->value()) + "\n");
	input_params.write(line.toUtf8());

	line = ("ouput_prefix " + QString::fromStdString(Inp->output) + "\n");
	input_params.write(line.toUtf8());

	input_params.close();
}

void MCWidget::RunEng(){
	this->set_parameters();
    RunEngine = new TEMP_SCHEME(Inp, Edit, progressbar);
    RunEngine->evaluation();
}

MCWidget::~MCWidget()
{

}

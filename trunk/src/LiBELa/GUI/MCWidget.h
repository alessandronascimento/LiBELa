#ifndef MCWIDGET_H
#define MCWIDGET_H

#include <QWidget>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QRadioButton>
#include "../PARSER.h"
#include <QProgressBar>
#include<QVBoxLayout>
#include<QHBoxLayout>
#include<QGridLayout>
#include<QComboBox>
#include<QFileDialog>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QLineEdit>
#include <QFile>
#include <QLineEdit>
#include <QString>
#include <math.h>
#include <stdio.h>
#include "../RunEngine.h"


class MCWidget : public QWidget
{
    Q_OBJECT

public:
    MCWidget(PARSER* Input, QPlainTextEdit* Editor);
    ~MCWidget();

    PARSER* Inp;
    QPlainTextEdit* Edit;

    QGridLayout* editLayout;
	QVBoxLayout* mainLayout;
	QHBoxLayout* buttonLayout;
	QHBoxLayout* messageLayout;

	QVBoxLayout* progress_layout;
	QProgressBar *progressbar;
	QLabel* progress_label;

    QSpinBox* nsteps;
    QLabel* nsteps_lab;

    QDoubleSpinBox* temperature;
    QLabel* temperature_lab;

    QRadioButton* sa;

    QDoubleSpinBox* cushion;
    QLabel* cushion_lab;

    QDoubleSpinBox* rotation_step;
	QLabel* rotation_step_lab;

    QComboBox* scoring_function;
    QLabel* scoring_function_lab;

    QDoubleSpinBox* diel;
    QLabel* diel_lab;

    QDoubleSpinBox* sigma;
    QLabel* sigma_lab;

    QDoubleSpinBox* deltaij_vdw;
    QLabel* deltaij_vdw_lab;

    QDoubleSpinBox* deltaij_elec;
    QLabel* deltaij_elec_lab;

    QPushButton* rec_mol2;
    QLineEdit* choose_rec_mol2;

    QPushButton* lig_mol2;
    QLineEdit* choose_lig_mol2;

    QComboBox* mol2_aa;
	QLabel* mol2_aa_lab;

	QLabel* output_prefix_lab;
	QLineEdit* output_prefix;

	QSpinBox* timeout;
	QLabel* timeout_lab;

	QDoubleSpinBox* box_size;
	QLabel* box_size_lab;

	QDoubleSpinBox* sa_start_temp;
	QLabel* sa_start_temp_lab;

	QSpinBox* sa_steps;
	QLabel* sa_steps_lab;

	QLabel* input_files_lab;
	QLabel* input_sf_lab;
	QLabel* input_sampling_lab;
	QLabel* input_misc_lab;
	QLabel* input_sa_lab;

	QPushButton* Start;
	QPushButton* WriteButton;
	TEMP_SCHEME* RunEngine;

private:

private slots:
	void show_sa_parameters(bool state);
	void set_parameters();
	void write_parameters();
	void RunEng();
	void choose_rec_file();
	void choose_lig_file();

};

#endif // MCWIDGET_H

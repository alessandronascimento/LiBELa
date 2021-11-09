#ifndef DOCKWIDGET_H
#define DOCKWIDGET_H

#include<cstdio>
#include<QWidget>
#include<QVBoxLayout>
#include<QHBoxLayout>
#include<QGridLayout>
#include<QLabel>
#include<QLineEdit>
#include<QPushButton>
#include<QComboBox>
#include<QDoubleSpinBox>
#include<QSpinBox>
#include<iostream>
#include<QFileDialog>
#include<QPlainTextEdit>
#include<QStringList>
#include<QProgressBar>
#include <QCheckBox>
#include <QScrollArea>
#include <QScroller>
#include"../PARSER.h"
#include"../RunEngine.h"
#include <math.h>
#include <QDebug>

class DockWidget : public QWidget
{
    Q_OBJECT

public:
//    DockWidget(QWidget *parent = 0);
    DockWidget(PARSER *Input, QPlainTextEdit* Ed);
    void hide_advanced_options();
    void show_advanced_options();
    void set_initial_parameters();
    ~DockWidget();

    QScrollArea* scrollarea;

    PARSER* Inp;
    QPlainTextEdit* Editor;

    QGridLayout* editLayout;
    QVBoxLayout* mainLayout;
    QHBoxLayout* buttonLayout;
    QHBoxLayout* messageLayout;

    QVBoxLayout* progress_layout;
    QProgressBar *progressbar;
    QLabel* progress_label;

    QLabel* done;

    QPushButton* rec_mol2;
    QLineEdit* choose_rec_mol2;

    QComboBox* mol2_aa;
    QLabel* mol2_aa_lab;

    QPushButton* lig_mol2;
    QLineEdit* choose_lig_mol2;

    QPushButton* reflig_mol2;
    QLineEdit* choose_reflig_mol2;

    QPushButton* docking_mol2;
    QLineEdit* choose_docking_mol2;

    QComboBox* scoring_function;
    QLabel* scoring_function_lab;

    QDoubleSpinBox* diel;
    QLabel* diel_lab;

    QComboBox* dielectric_model;
    QLabel* dielectric_model_lab;

    QDoubleSpinBox* sigma;
    QLabel* sigma_lab;

    QDoubleSpinBox* alpha;
    QLabel* alpha_lab;

    QDoubleSpinBox* beta;
    QLabel* beta_lab;

    QDoubleSpinBox* deltaij_vdw;
    QLabel* deltaij_vdw_lab;

    QDoubleSpinBox* deltaij_elec;
    QLabel* deltaij_elec_lab;

    QDoubleSpinBox* elec_scale;
    QLabel* elec_scale_lab;

    QDoubleSpinBox* vdw_scale;
	QLabel* vdw_scale_lab;

	QDoubleSpinBox* box_size;
	QLabel* box_size_lab;

	QDoubleSpinBox* min_delta;
	QLabel* min_delta_lab;

	QDoubleSpinBox* min_tol;
	QLabel* min_tol_lab;

	QSpinBox* min_timeout;
	QLabel* min_timeout_lab;

	QComboBox* overlay_optimizer;
	QLabel* overlay_optimizer_lab;
	QComboBox* energy_optimizer;
	QLabel* energy_optimizer_lab;

    QLabel* output_prefix_lab;
	QLineEdit* output_prefix;

    QPushButton* Run;
    QPushButton* WriteButton;

    QCheckBox* advanced_settings;
    QCheckBox* dock_no_h;

    QCheckBox* generate_conformers;
    QComboBox* ligand_energy_model;
    QLabel* ligand_energy_model_lab;

    QSpinBox* conformers;
    QLabel* conformers_lab;

    QSpinBox* conformers_to_evaluate;
    QLabel* conformers_to_evaluate_lab;

    QDoubleSpinBox* dock_min_tol;
    QLabel* dock_min_tol_lab;

    QSpinBox* dock_timeout;
    QLabel* dock_timeout_lab;

    QCheckBox* parallel;
    QLabel* parallel_lab;

    QSpinBox* parallel_jobs;
    QLabel* parallel_jobs_lab;

    QCheckBox* write_mol2;

    QCheckBox* use_grids;
    QDoubleSpinBox* grid_spacing;
    QLabel* grid_spacing_lab;
    QComboBox* load_write_file;
    QLineEdit* grid_file;
    QCheckBox* only_score;
    QDoubleSpinBox* grid_box;
    QLabel* grid_box_lab;

    QLabel* files_label;
    QLabel* output_label;
    QLabel* sf_input;
    QLabel* opt_input;
    QLabel* conf_input;
    QLabel* parallel_label;
    QLabel* grid_label;

    QLabel* search_box_lab;
    QDoubleSpinBox* search_box;

    QCheckBox* use_write_cutoff_score;
    QCheckBox* use_cutoff_energy;
    QDoubleSpinBox* write_cutoff_score;
    QDoubleSpinBox* cutoff_energy;
//    QLabel* use_write_cutoff_score_lab;
//    QLabel* use_cutoff_energy_lab;

    QLabel* use_cutoff_score_lab;
    QDoubleSpinBox* cutoff_score;

    QLabel* docking_restraint_lab;
    QDoubleSpinBox* docking_restraint_weight;
    QLabel* use_docking_restraint_lab;
    QCheckBox* use_docking_restraint;

    QLabel* scale_ele_energy_lab;
    QLabel* scale_vdw_energy_lab;
    QDoubleSpinBox* scale_ele_energy;
    QDoubleSpinBox* scale_vdw_energy;

    QCheckBox* sort_conformers;

    QLabel* delphi_gsize_label;
    QCheckBox* use_delphi_phimap;
    QSpinBox* delphi_gsize;
    QLabel* phimap_file_label;
    QLineEdit* phimap_file;
    QString delphi_phimap_filename;


    TEMP_SCHEME* RunEngine;

public slots:
	void choose_rec_file();
	void choose_reflig_file();
	void choose_docking_mol2_files();
    void choose_grid_file();
    void slot_grid_file(int i);
    void slot_use_grids(int state);
	void show_advanced_setttings(bool state);
    void slot_write_mol2(int state);
    void slot_ignore_h(int state);
    void slot_generate_conformers(int state);
    void slot_parallel(int state);
    void slot_use_cutoff_score(int state);
    void slot_use_cutoff_energy(int state);
	void Start();
	void set_parameters();
	void write_parameters();
    void show_scoring_function_options(int sf);
    void slot_use_restraints(int state);
    void slot_scale_ele(double value);
    void slot_scale_vdw(double value);
    void slot_use_delphi_phimap(int state);
};

#endif // DOCKWIDGET_H

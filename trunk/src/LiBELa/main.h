/*
 * main.h
 *
 *  Created on: 23/03/2010
 *      Author: Alessandro
 */

#ifndef MAIN_H_
#define MAIN_H_

#include<iostream>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<time.h>
#include"nlopt.h"
#include"nlopt.hpp"
#include "Optimizer.h"
#include "iMcLiBELa.h"
#include "RunEngine.h"
#ifdef HAS_GUI
#include "GUI/GUI.h"
#include <QTimer>
#include <QIcon>
#include <QSplashScreen>
#include <QtGui>
#include <QApplication>
#endif

Mol2* Optimizer::Rec;
Mol2* Optimizer::RefLig;
PARSER* Optimizer::Parser;
Grid* Optimizer::Grids;


#endif /* MAIN_H_ */

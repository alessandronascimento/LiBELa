#include "main.h"


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


int main(int argc, char *argv[])
{

#ifdef HAS_GUI

    if (argc > 1){
    	// using terminal mode
    	if (argc < 2){
            printf("Usage: %s input _file\n", argv[0]);
    		exit(1);
    	}
    	else {
            TEMP_SCHEME* RunEngine = new TEMP_SCHEME(argc, argv);
    		RunEngine->evaluation();
            delete RunEngine;
    		return(0);
    	}
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
    if (argc > 1){
        	// using terminal mode
        	if (argc < 2){
        		printf("Usage: %s input_file\n", argv[0]);
        		exit(1);
        	}
        	else {
                TEMP_SCHEME* RunEngine = new TEMP_SCHEME(argc, argv);
        		RunEngine->evaluation();
        		return(0);
        	}
        }
#endif
}

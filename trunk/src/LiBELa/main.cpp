#include "main.h"


/*********************************************************************************************
 *           iMcLiBELa - Monte Carlo - Based Ligand Biding Energy Landscape                  *
 *                                                                                           *
 *               This program was written by Alessandro Nascimento - March/2012              *
 * 					ABC Federal University - Santo Andre, SP. Brazil                         *
 *                    Engineering, Modeling and Social Sciences Center                       *
 *                          alessandro.nascimento@ufabc.edu.br                               *
 *                                                                                           *
 *           Please visit http://nascimento.ifsc.usp.br/  for detailed instructions          *
 *                  Distributed under the terms of GNU GPLv3 license                         *
 *********************************************************************************************/


int main(int argc, char *argv[])
{

#ifdef HAS_GUI

    if (argc > 1){
    	// using terminal mode
    	if (argc < 2){
            printf("Usage: McLiBELa.x input _file\n");
    		exit(1);
    	}
    	else {
            TEMP_SCHEME* RunEngine = new TEMP_SCHEME(argc, argv);
    		RunEngine->evaluation();
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

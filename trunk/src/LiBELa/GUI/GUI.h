#ifndef GUI_H
#define GUI_H

#include <QMainWindow>
#include <QAction>
#include <QDockWidget>
#include <QTextBrowser>
#include <iostream>
#include <QStringList>
#include <vector>
#include <QFile>
#include <QStatusBar>
#include <QTextEdit>
#include <QPlainTextEdit>
#include <QMenuBar>
#include <QTabWidget>
#include "plotter.h"
#include <QApplication>
#include <QMessageBox>
#include <QScrollArea>
#include "DockWidget.h"
#include "../PARSER.h"
#include "MCWidget.h"
#include <QKeySequence>


using namespace std;

class GUI : public QMainWindow
{
    Q_OBJECT

public:
    GUI();
    ~GUI();

    void createActions();
		PARSER *Input;

        QPlainTextEdit* Editor;

        QMenu *InputFile;
        QAction *DockAction;
        QAction* MCAction;
        QAction* LogAction;

        QMenu *Results;
        QAction *ResultsAction;

        QMenu *Graphics;
        QAction *GraphicsAction;

        QMenu *Help;
        QAction *About_iMcLiBELa;
        QMessageBox *ElSA_Message;
        QAction *About_Qt;

        QMenu *Exit;
        QAction *ExitAction;

        QTabWidget* tabwidget;

        QDockWidget* DockEditor;

        QTextBrowser* ResultsText;

        QWidget* childWildget;

        DockWidget* DockInput;
        MCWidget* MCInput;
        QKeySequence MCKey;
        QKeySequence DockKey;
        QKeySequence LogKey;

        QScrollArea* scrollarea;

    private:

    private slots:
        void openInputWidget();
        void openMCInputWidget();
        void ViewResults();
        void GraphResult();
        void RunEng();
        void closethistab(int tindex);
        void about_libela();
        void openLogArea();
        void show_status_bar_message(QString message);
};

#endif // GUI_H

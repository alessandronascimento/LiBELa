#include "GUI.h"
#include "../iMcLiBELa.h"

//using namespace std;

GUI::GUI()
{
	setWindowTitle(tr(LONG_NAME));
	createActions();
}

GUI::~GUI()
{
}

void GUI::createActions()
{
	Input = new PARSER;

	tabwidget = new QTabWidget;
	tabwidget->setTabsClosable(true);
    tabwidget->setMovable(true);
    tabwidget->setAutoFillBackground(true);
	openLogArea();

    setCentralWidget(tabwidget);

    statusBar()->show();

    DockKey = (tr("Ctrl+d"));
    MCKey = (tr("Ctrl+m"));
    LogKey = (tr("Ctrl+l"));

    DockAction = new QAction(tr("&Dock Mode"), this);
    DockAction->setStatusTip(tr("Start a Docking Session"));
    DockAction->setShortcut(DockKey);
    InputFile = menuBar()->addMenu(tr("Input Parameters"));

    MCAction = new QAction(tr("&MC Mode"), this);
    MCAction->setStatusTip("Start a Monte Carlo Session");
    MCAction->setShortcut(MCKey);
//    MCAction->setDisabled(true);

    LogAction = new QAction(tr("Open &Log Area"), this);
    LogAction->setStatusTip("Opens a new Log Area");
    LogAction->setShortcut(LogKey);
    LogAction->setDisabled(true);


    InputFile->addAction(DockAction);
    InputFile->addAction(MCAction);
    InputFile->addAction(LogAction);

    Results = menuBar()->addMenu(tr("&Results"));
    ResultsAction = new QAction(tr("&View Results File"), this);
    ResultsAction->setStatusTip(tr("Visualize Results"));
    Results->addAction(ResultsAction);


    Graphics = menuBar()->addMenu(tr("&Graphics"));
    GraphicsAction = new QAction(tr("&Graph Results"), this);
    GraphicsAction->setStatusTip(tr("Graph dat result file"));
    Graphics->addAction(GraphicsAction);


    Help = menuBar()->addMenu(tr("&Help"));
    About_iMcLiBELa = new QAction(tr("&About iMcLiBELa"), this);
    Help->addAction(About_iMcLiBELa);
    About_iMcLiBELa->setStatusTip(tr("Help on how to use the program"));
    About_Qt = new QAction(tr("About Qt"), this);
    Help->addAction(About_Qt);
    About_Qt->setStatusTip(tr("About Qt, the language used in this GUI"));

    Exit = menuBar()->addMenu((tr("&Exit")));
    ExitAction = new QAction(tr("&Exit"), this);
    Exit->addAction(ExitAction);

    connect(ExitAction, SIGNAL(triggered()), this, SLOT(close()));
    connect(tabwidget, SIGNAL(tabCloseRequested(int)), this, SLOT(closethistab(int)));
    connect(About_iMcLiBELa, SIGNAL(triggered()), this, SLOT(about_libela()));
    connect(About_Qt, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
    connect(GraphicsAction, SIGNAL(triggered()), this, SLOT(GraphResult()));
    connect(ResultsAction, SIGNAL(triggered()), this, SLOT(ViewResults()));
    connect(DockAction, SIGNAL(triggered()), this, SLOT(openInputWidget()));
    connect(MCAction, SIGNAL(triggered()), this, SLOT(openMCInputWidget()));
    connect(LogAction, SIGNAL(triggered()), this, SLOT(openLogArea()));


    this->openInputWidget();

}

void GUI::openInputWidget(){
    DockInput = new DockWidget(Input, Editor);
    scrollarea = new QScrollArea;
    scrollarea->resize(DockInput->size());
    scrollarea->setWidget(DockInput);
    tabwidget->addTab(scrollarea, tr("Dock Session"));
    tabwidget->setFixedWidth(scrollarea->width());
}

void GUI::openMCInputWidget(){
	MCInput = new MCWidget(Input, Editor);
	tabwidget->addTab(MCInput, tr("Monte Carlo Session"));
}

void GUI::openLogArea(){
	Editor = new QTextEdit;
	DockEditor = new QDockWidget;
	DockEditor->setAllowedAreas(Qt::RightDockWidgetArea);
	DockEditor->setObjectName("Log Area");
	DockEditor->setWindowTitle("Log Area");
	addDockWidget(Qt::RightDockWidgetArea, DockEditor);
	DockEditor->setWidget(Editor);
	Editor->setDocumentTitle("iMcLiBELa LOG AREA");
}

void GUI::ViewResults(){
    ResultsText = new QTextBrowser(0);
    QString filename = QFileDialog::getOpenFileName(this, tr("Choose a Results File"),
                                                    "", tr("LOG Files (*.log)"));
    QFile file(filename);
    if(file.open(QIODevice::ReadOnly|QIODevice::Text)){
        ResultsText->setPlainText(QString::fromUtf8(file.readAll()));
        statusBar()->showMessage(tr("File loaded successfully."), 5000);

        tabwidget->addTab(ResultsText, tr("Log File"));
        int i = tabwidget->indexOf(ResultsText);
        tabwidget->setCurrentIndex(i);
    }
}

void GUI::GraphResult(){
    plotter* Plotter = new plotter;
    tabwidget->addTab(Plotter, tr("Graphed Results"));
    int i = tabwidget->indexOf(Plotter);
    tabwidget->setCurrentIndex(i);
}

void GUI::RunEng(){
}

void GUI::closethistab(int tindex){
	tabwidget->removeTab(tindex);
}

void GUI::about_libela(){
	ElSA_Message = new QMessageBox;
	ElSA_Message->information(this, tr("McLiBELa"),
			tr("This program was written to dock and simulate ligand dynamics in a\n"
			    "protein active site. The program is distributed in the hope that it\n"
				"will be useful, but WITHOUT ANY WARRANTY; without even the implied \n"
				"warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n"
                "Detailed information can be found in http://nascimento.ifsc.usp.br"));

	ElSA_Message->setIcon(QMessageBox::Information);
}

void GUI::show_status_bar_message(QString message){
	statusBar()->showMessage(message, 3000);
}

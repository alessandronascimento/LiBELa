TEMPLATE = app
TARGET = iMcLiBELa
QT += core \
    gui \
    widgets \
    concurrent
DEFINES += HAS_GUI
DEFINES += STATIC
DEFINES += BUILD=2
HEADERS += ../LiBELa/GUI/DockWidget.h \
    ../LiBELa/GUI/GUI.h \
    ../LiBELa/GUI/QtWriter.h \
    ../LiBELa/GUI/plotter.h \
    ../LiBELa/GUI/MCWidget.h \
    ../LiBELa/RunEngine.h \
    ../LiBELa/Conformer.h \
    ../LiBELa/COORD_MC.h \
    ../LiBELa/Deal.h \
    ../LiBELa/Docker.h \
#    ../LiBELa/ENERGY.h \
    ../LiBELa/Gaussian.h \
    ../LiBELa/McLiBELa.h \
    ../LiBELa/Mol2.h \
    ../LiBELa/MC.h \
    ../LiBELa/Optimizer.h \
    ../LiBELa/PARSER.h \
    ../LiBELa/RAND.h \
    ../LiBELa/WRITER.h \
    ../LiBELa/main.h \
    ../LiBELa/Grid.h \
    ../LiBELa/SA.h \
    ../LiBELa/iMcLiBELa.h \
    ../LiBELa/Energy2.h
SOURCES += ../LiBELa/GUI/DockWidget.cpp \
    ../LiBELa/GUI/GUI.cpp \
    ../LiBELa/GUI/QtWriter.cpp \
    ../LiBELa/GUI/plotter.cpp \
    ../LiBELa/GUI/MCWidget.cpp \
    ../LiBELa/RunEngine.cpp \
    ../LiBELa/Conformer.cpp \
    ../LiBELa/COORD_MC.cpp \
    ../LiBELa/Deal.cpp \
    ../LiBELa/Docker.cpp \
#    ../LiBELa/ENERGY.cpp \
    ../LiBELa/Gaussian.cpp \
    ../LiBELa/Mol2.cpp \
    ../LiBELa/MC.cpp \
    ../LiBELa/Optimizer.cpp \
    ../LiBELa/PARSER.cpp \
    ../LiBELa/RAND.cpp \
    ../LiBELa/WRITER.cpp \
    ../LiBELa/Grid.cpp \
    ../LiBELa/SA.cpp \
    ../LiBELa/main.cpp \
    ../LiBELa/Energy2.cpp
RESOURCES += ../LiBELa/GUI/GUI.qrc
INCLUDEPATH += ../../lib/openbabel/include/openbabel-2.0 \
    ../lib/eigen3 \
    ../lib/nlopt/include \
    ../lib/gsl/include
LIBS += -lnlopt_cxx \
    -L../lib/nlopt/lib \
    -L../lib/openbabel/lib \
    -L../lib/gsl/lib \
    -lopenbabel \
    -lz \
    -lgsl \
    -lgslcblas \
    -lm \
    -lgomp
#    -liomp5
#QMAKE_CXX = /opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
QMAKE_CXXFLAGS += -fopenmp -w -ffast-math -O3 -static
OTHER_FILES += ../tools/binning.cpp \
    ../tools/FESummation.cpp

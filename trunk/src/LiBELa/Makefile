#!/bin/bash
CC=/usr/bin/g++
CFLAGS=-c -O3 -ffast-math -w -DBUILD=2 -fopenmp -I/share/apps/iMcLiBELa/trunk/lib/openbabel/include/openbabel-2.0 -I/share/apps/iMcLiBELa/trunk/src/lib/eigen3 -I/share/apps/iMcLiBELa/trunk/lib/nlopt/include -I/share/apps/iMcLiBELa/trunk/lib/gsl/include
LDFLAGS= -lgsl -lgslcblas -lm -lopenbabel -lz -lnlopt_cxx -lgomp -L/share/apps/iMcLiBELa/trunk/lib/openbabel/lib -L/share/apps/iMcLiBELa/trunk/lib/nlopt/lib -L/share/apps/iMcLiBELa/trunk/lib/gsl/lib
SOURCES=Conformer.cpp COORD_MC.cpp Docker.cpp Energy2.cpp ENERGY.cpp Gaussian.cpp Grid.cpp main.cpp Mol2.cpp Optimizer.cpp PARSER.cpp RAND.cpp RunEngine.cpp SA.cpp WRITER.cpp MC.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../../bin/McLiBELa.openmp

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS)

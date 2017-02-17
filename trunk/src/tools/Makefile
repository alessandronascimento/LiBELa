#!/bin/bash

CC=g++
CFLAGS=-c -O3 -ffast-math -fopenmp -I../../lib/openbabel/include/openbabel-2.0 -I../lib/eigen3 -DBUILD=0 -I../../lib/gsl/include -I../../lib/boost/include
LDFLAGS=-lz -L../../lib/openbabel/lib -lopenbabel -lgomp -lpthread -L../../lib/gsl/lib -lgsl -lgslcblas -lm
SOURCES=McConf.cpp McGrid.cpp MrPocket.cpp mol2_rmsd.cpp McColony.cpp FESummation.cpp McBootstrapper.cpp McCharge.cpp AutoCorrelation.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: McConf McGrid mol2_rmsd MrPocket McColony McSummation McBootstrapper McCharge AutoCorr
	mv McConf McGrid mol2_rmsd McPocket McColony McSummation McBootstrapper McCharge AutoCorr ../../bin/

McConf: McConf.o
	$(CC) McConf.o -o McConf $(LDFLAGS)

McGrid: McGrid.o
	$(CC) McGrid.o -o McGrid $(LDFLAGS)

mol2_rmsd: mol2_rmsd.o
	$(CC) mol2_rmsd.o -o mol2_rmsd $(LDFLAGS)

MrPocket: MrPocket.o
	$(CC) MrPocket.o -o McPocket $(LDFLAGS)

McColony: McColony.o 
	$(CC) McColony.o -o McColony $(LDFLAGS)

McSummation: FESummation.o
	$(CC) FESummation.o -o McSummation $(LDFLAGS)

McBootstrapper: McBootstrapper.o
	$(CC) McBootstrapper.o -o McBootstrapper $(LDFLAGS)

McCharge: McCharge.o
	$(CC) McCharge.o -o McCharge $(LDFLAGS)

AutoCorr: AutoCorrelation.o
	$(CC) Autocorrelation.o -o AutoCorr $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS) McConf McGrid mol2_rmsd McPocket

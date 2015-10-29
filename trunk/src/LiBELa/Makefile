#!/bin/bash
CC=/opt/intel/composer_xe_2015.3.187/bin/intel64/icpc
CFLAGS=-c -O3 -ffast-math -w -DBUILD=2 -fopenmp -I/home/asn/workspace/libela/trunk/lib/openbabel/include/openbabel-2.0 -I/home/asn/workspace/libela/trunk/src/lib/eigen3 -I/home/asn/workspace/libela/trunk/lib/nlopt/include -I/home/asn/workspace/libela/trunk/lib/gsl/include
LDFLAGS= -lgsl -lgslcblas -lm -lopenbabel -lz -lnlopt_cxx -liomp5 -L/home/asn/workspace/libela/trunk/lib/openbabel/lib -L/home/asn/workspace/libela/trunk/lib/nlopt/lib -L/home/asn/workspace/libela/trunk/lib/gsl/lib
SOURCES=Conformer.cpp COORD_MC.cpp Docker.cpp Energy2.cpp Gaussian.cpp Grid.cpp main.cpp Mol2.cpp Optimizer.cpp PARSER.cpp RAND.cpp RunEngine.cpp SA.cpp WRITER.cpp MC.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=../../bin/McLiBELa.openmp

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm $(OBJECTS)

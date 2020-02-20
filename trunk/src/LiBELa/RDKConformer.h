#ifndef RDKCONFORMER_H
#define RDKCONFORMER_H

#include "PARSER.h"
#include <fstream>
#include <iostream>
#include <string>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>



class RDKConformer
{
public:
    RDKConformer();
     PARSER* Input;
};

#endif // RDKCONFORMER_H

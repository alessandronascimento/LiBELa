#include "RDKConformer.h"

RDKConformer::RDKConformer(PARSER* _Input)
{
    this->Input = _Input;
}

RDKit::ROMOL_SPTR RDKConformer::openMol(string molfile){
    RDKit::ROMOL_SPTR mol;
    boost::iostreams::filtering_istream ins;
    ins.push( boost::iostreams::gzip_decompressor() );
    std::string comp_file = molfile;
    ins.push( boost::iostreams::file_source(comp_file));
    RDKit::ROMOL_SPTR mol(&ins);
    return(mol);
}


RDKConformer::genereate_conformers(RDKit::ROMOL_SPTR mol){
    RDKit::MMFF::MMFFOptimizeMolecule(*mol1);
    RDKit::INT_VECT conformers = RDKit::DGeomHelpers::EmbedMultipleConfs(*mol, Input->lig_conformers);


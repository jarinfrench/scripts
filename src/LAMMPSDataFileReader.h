#ifndef LAMMPSDATAFILEREADER_H
#define LAMMPSDATAFILEREADER_H

#include <fstream>
#include <vector>
#include <map>
#include <string>
#include "atom.h"
#include "Box.h"

class LAMMPSDataFileReader {
private:
  std::vector <Atom> _atoms;
  std::string _comment; // header line comment
  int _N; // number of atoms
  int _bonds; // number of bonds
  int _n_types; // number of atom types
  Box _box; // Box dimensions and tilt factors
  std::string _file; // filename
  std::vector <std::string> _header_keywords;
  std::vector <std::string> _section_keywords;
  std::vector <std::string> _atom_style_keywords;
  std::vector <std::string> _pair_style_keywords;
  std::vector <std::string> _bond_style_keywords;
  std::vector <std::string> _angle_style_keywords;
  std::vector <std::string> _dihedral_style_keywords;
  std::vector <std::string> _improper_style_keywords;
  void _initKeywords();

public:
  LAMMPSDataFileReader(std::string file);
  void readFile();
};

#endif // LAMMPSDATAFILEREADER_H

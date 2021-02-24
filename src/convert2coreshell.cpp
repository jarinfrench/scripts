#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream> // may not need this one
#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "Box.h"
#include "error_code_defines.h"

using namespace std;

struct Header {
  int N, n_types;
};

struct Bond {
  int id, type, atom1, atom2;
};

tuple <vector <Atom>, Box, Header> readData(const string& input) {
  string str, comment;
  int N, n_types, id, type;
  Box box;
  Header header;
  double x, y, z;
  vector <Atom> atoms;

  ifstream fin(input.c_str());
  if (fin.fail()) {
    cerr << "Error opening file \"" << input << "\"\n";
    exit(FILE_OPEN_ERROR);
  }

  getline(fin, comment); // get the comment line
  getline(fin, str); // blank line
  fin >> N >> str; // get the number of atoms
  fin >> n_types >> str; // get the number of atom types
  header.N = N;
  header.n_types = n_types;
  fin >> x >> y >> str >> str; box.setXLow(x); box.setXHigh(y); // get the box bounds
  fin >> x >> y >> str >> str; box.setYLow(x); box.setYHigh(y);
  fin >> x >> y >> str >> str; box.setZLow(x); box.setZHigh(y);

  atoms.resize(N, Atom());

  fin.ignore(); // get rid of extraneous characters in stream
  getline(fin, str); // blank line
  getline(fin, str); // Atoms
  getline(fin, str); // blank line

  while (getline(fin, str)) {
    stringstream ss(str);
    if (!(ss >> id >> type >> x >> y >> z)) {
      cerr << "File \"" << input << "\" corrupted: " << str << "\n";
      exit(FILE_FORMAT_ERROR);
    }

    atoms[id - 1] = Atom(id, type, 0, Position(x, y, z));
  }

  return make_tuple(atoms, box, header);
}

void writeData(const vector <Atom>& atoms, const Box& box, const Header& header, const string& output, const string& input) {
  int n_out = 0;
  vector <Bond> bonds (header.n_types);
  ofstream fout(output.c_str());
  if (fout.fail()) {
    cerr << "Error opening file \"" << output << "\n for writing\n";
    exit(FILE_OPEN_ERROR);
  }

  fout << "This file was converted to the core-shell format from file " << input << "\n\n";
  fout << header.N * 2 << " atoms\n"; // number of atoms is double (core and shell atoms)
  fout << header.N << "bonds\n"; // number of bonds is the original number of atoms
  fout << header.n_types * 2 << " atom types\n"; // assumes each type has a core-shell interaction
  fout << header.n_types << " bond types\n\n";

  //box.writeBounds(fout);
  fout << "\n\nAtoms\n\n";
  for (size_t i = 0; i < atoms.size(); ++i) {
    // atom_id molecule_id(old atom id) type charge x y z
    fout << ++n_out << " " << atoms[i].getId() << " " << atoms[i].getType()
         << " c" << atoms[i].getType() << " " << atoms[i].getWrapped().getX()
         << " " << atoms[i].getWrapped().getY() << " "
         << atoms[i].getWrapped().getZ() << "\n";

    fout << ++n_out << " " << atoms[i].getId() << " " << atoms[i].getType() + header.n_types
         << " c" << atoms[i].getType() + header.n_types << " " << atoms[i].getWrapped().getX()
         << " " << atoms[i].getWrapped().getY() << " "
         << atoms[i].getWrapped().getZ() << "\n";
    bonds[i].id = i + 1;
    bonds[i].type = i % header.n_types + 1;
    bonds[i].atom1 = n_out - 1;
    bonds[i].atom2 = n_out;
  }

  fout << "\nBonds\n\n";
  for (size_t i = 0; i < bonds.size(); ++i) {
    // bond_id, bond_type, atom1, atom2
    fout << bonds[i].id << " " << bonds[i].type << " " << bonds[i].atom1 << " " << bonds[i].atom2 << "\n";
  }

  fout.close();
}

int main(int argc, char** argv) {
  string file, output;
  vector <Atom> atoms;

  try {
    cxxopts::Options options(argv[0], "Convert a file in the atomic data format to the core-shell format");
    options
      .positional_help("file")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "File to convert", cxxopts::value<string>(file), "file")
        ("o,output", "Output file name", cxxopts::value<string>(output)->default_value("converted_file.dat"), "file")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !(result.count("file"))) {
      cout << options.help() << "\n";
      return EXIT_SUCCESS;
    }

    if (result.count("file")) {
      vector <Atom> atoms;
      Box box;
      Header header;
      tie(atoms, box, header) = readData(file);
      writeData(atoms, box, header, output, file);
    }

  } catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

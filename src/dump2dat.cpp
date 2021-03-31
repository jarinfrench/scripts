/******************************************************************************
* This program parses output from LAMMPS minimized structures and formats it
* so LAMMPS can run a simulation with the results.
*******************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cxxopts.hpp>
#include <map>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"

using namespace std;

struct Bond {
  int id1, id2;
  int type1, type2;

  Bond() : id1(-1), id2(-1), type1(-1), type2(-1) {}
};

struct dataValues {
  int N, n_types = 0, n_bonds = 0, n_bond_types = 0; // number of atoms, number of atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, xy, xz, yz; // box bounds and tilt factors
  vector <string> data_labels;
  vector <Atom> atoms;
  vector <Bond> bonds;
  vector <int> bond_types;
  bool using_charge = false;
  bool has_tilt = false;
  bool has_bonds = false;
};

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

dataValues readData(const string& infile) {
  string str; // junk string variable
  dataValues data;
  int id = 0, type = 0, n_total = 0, mol_id = 0;
  double charge = 0.0, x = 0.0, y = 0.0, z = 0.0;

  ifstream fin(infile.c_str());
  checkFileStream(fin, infile);

  for (int i = 0; i < 3; ++i) {getline(fin, str);} // gets some header info we don't need
  fin >> data.N;

  data.atoms.resize(data.N, Atom());
  data.bonds.resize(data.N / 2, Bond()); // At minimum, there will be n_atoms / 2 bonds

  fin.ignore();
  getline(fin, str); // gets the box bounds info line
  if (str.find("xy") == string::npos) {
    fin >> data.xlow >> data.xhigh >> data.ylow >> data.yhigh >> data.zlow >> data.zhigh;
  }
  else
  {
    data.has_tilt = true;
    fin >> data.xlow >> data.xhigh >> data.xy >> data.ylow >> data.yhigh >> data.xz >> data.zlow >> data.zhigh >> data.yz;
  }

  fin.ignore();
  getline(fin, str); // gets the line "ITEM: ATOMS <atom data labels>"
  stringstream ss(str);
  ss >> str >> str;
  while (ss >> str) {data.data_labels.push_back(str);} // get the labels

  // Now get the atom data
  while (getline(fin, str)) {
    stringstream ss2(str);
    for (int i = 0; i < data.data_labels.size(); ++i) {
      if (data.data_labels[i].compare("id") == 0) {ss2 >> id;}
      else if (data.data_labels[i].compare("mol") == 0) {
        ss2 >> mol_id;
        data.has_bonds = true;
      }
      else if (data.data_labels[i].compare("type") == 0) {
        ss2 >> type;
        if (type > data.n_types) {data.n_types = type;}
      }
      else if (data.data_labels[i].compare("q") == 0) {
        ss2 >> charge;
        data.using_charge = true;
      }
      else if (data.data_labels[i].compare("x") == 0) {ss2 >> x;}
      else if (data.data_labels[i].compare("y") == 0) {ss2 >> y;}
      else if (data.data_labels[i].compare("z") == 0) {ss2 >> z;}
      else {continue;}
    }

    Position p(x,y,z);
    data.atoms[id - 1] = Atom(id, type, charge, p);
    if (data.has_bonds) {
      data.atoms[id - 1].setExtraInfo(mol_id, 0);
      data.atoms[id - 1].setExtraInfoNames(0, "mol");
      if (data.bonds[mol_id - 1].id1 == -1) {
        ++data.n_bonds;
        data.bonds[mol_id - 1].id1 = id;
        data.bonds[mol_id - 1].type1 = type;
      } else if (data.bonds[mol_id -1 ].id2 == -1) {
        data.bonds[mol_id - 1].id2 = id;
        data.bonds[mol_id - 1].type2 = type;
      } else {
        cout << "Error determining bond information for bond " << mol_id << "\n"
             << "Atom1 = " << data.bonds[mol_id - 1].id1 << " Atom2 = " << data.bonds[mol_id - 1].id2
             << "Extra atom = " << id << "\n";
      }
    }
    ++n_total;
  }

  data.bonds.resize(data.n_bonds, Bond()); // gets rid of the unallocated bonds
  data.bond_types.resize(data.n_bonds, 0);
  map <pair <int, int>, int> bond_types;
  for (size_t i = 0; i < data.bonds.size(); ++i) {
    Bond bond = data.bonds[i];
    // Check to make sure each bond has a valid atom id
    if (bond.id2 == -1) {
      cout << "Missing bond for atom " << data.bonds[i].id1 << "\n";
    } else {
      pair <int, int> b = (bond.type1 < bond.type2) ? make_pair(bond.type1, bond.type2) : make_pair(bond.type2, bond.type1);
      auto ret = bond_types.emplace(make_pair(b, data.n_bond_types));
      if (ret.second) { // successful emplacement
        ++data.n_bond_types;
        data.bond_types[i] = data.n_bond_types - 1; // We index the bond types from 0 here
      } else { // already exists
        data.bond_types[i] = bond_types[b];
      }
    }
  }

  fin.close();

  if (n_total != data.N) {
    cerr << "Error: " << n_total << " atoms were read out of " << data.N << " atoms in the file.\n";
    exit(ATOM_COUNT_ERROR);
  }

  return data;
}

void writeData(const string& outfile, const string& chem_formula, const dataValues& data) {
  ofstream fout(outfile.c_str());
  checkFileStream(fout, outfile);

  fout << "These " << chem_formula << " coordinates are from a LAMMPS dump file: [id ";
  if (find(data.data_labels.begin(), data.data_labels.end(), "mol") != data.data_labels.end()) {
    fout << "mol ";
  }
  fout << "type ";
  if (find(data.data_labels.begin(), data.data_labels.end(), "q") != data.data_labels.end()) {
    fout << "charge ";
  }
  fout << "x y z]\n\n"
       << data.N << " atoms\n"
       << data.n_types << " atom types\n";

  if (data.has_bonds) {
    fout << data.n_bonds << " bonds\n"
         << data.n_bond_types << " bond types\n";
  }

  fout.precision(6);
  if (data.has_tilt) {
    fout << data.xlow << " " << data.xhigh << " xlo xhi\n"
         << data.ylow << " " << data.yhigh << " ylo yhi\n"
         << data.zlow << " " << data.zhigh << " zlo zhi\n"
         << data.xy << " " << data.xz << " " << data.yz << " xy xz yz\n\n";
  } else {
    fout << data.xlow << " " << data.xhigh << " xlo xhi\n"
         << data.ylow << " " << data.yhigh << " ylo yhi\n"
         << data.zlow << " " << data.zhigh << " zlo zhi\n\n";
  }

  fout << "Atoms\n\n";

  for (unsigned int i = 0; i < data.atoms.size(); ++i) {
    fout.precision(0);
    fout << fixed;
    fout << data.atoms[i].getId() << " ";
    if (data.has_bonds) {
      fout << data.atoms[i].getExtraInfo()[0] << " ";
    }
    fout << data.atoms[i].getType() << " ";
    if (data.using_charge) {
      fout.precision(6);
      fout << data.atoms[i].getCharge() << " ";
    }
    fout.precision(6);
    fout << data.atoms[i].getWrapped()[0] << " " << data.atoms[i].getWrapped()[1] << " "
         << data.atoms[i].getWrapped()[2] << endl;
  }

  if (data.has_bonds) {
    fout << "\nBonds\n\n";
    for (size_t i = 0; i < data.bonds.size(); ++i) {
      fout << i + 1 << " " << data.bond_types[i] + 1 << " " << data.bonds[i].id1 << " " << data.bonds[i].id2 << "\n";
    }
  }
  fout.close();
}

int main(int argc, char** argv) {
  string infile, outfile, chem_formula; // input file, output file, chemical formula
  dataValues data;

  try {
    cxxopts::Options options(argv[0], "Parse a LAMMPS dump file for use as a LAMMPS input data file");
    options
      .positional_help("infile outfile")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "File to convert", cxxopts::value<string>(infile), "file")
        ("o,output", "Name of output file", cxxopts::value<string>(outfile), "file")
        ("e,element", "The chemical formula of the system", cxxopts::value<string>(chem_formula)->default_value("<none specified>"), "element")
        ("h,help", "Show the help");

    options.parse_positional({"file", "output"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !(result.count("file"))) {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (!(result.count("output")))   {
      outfile = infile.substr(0, infile.rfind(".")) + ".dat"; // strips the last period off the original filename, and adds ".dat"
    }

    if (result.count("file")) {
      data = readData(infile);
      writeData(outfile, chem_formula, data);
    }

  }
  catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

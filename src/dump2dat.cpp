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

#define ID_MIN_DEFAULT 100000000

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
  Position min = Position(0,0,0); // minimum value in each direction
  unsigned int id_min = ID_MIN_DEFAULT; // minimum id value
};

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

dataValues readData(const string& infile, const bool& keep_bounds) {
  string str; // junk string variable
  dataValues data;
  int id = 0, type = 0, n_total = 0, mol_id = 0, n_mol_read = 0;
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
      if (data.data_labels[i].compare("id") == 0) {
        ss2 >> id;
        if (id < data.id_min) data.id_min = id;
      }
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

    Position p;
    if (!keep_bounds) {
      // adjust positions to start in the same place with the box dimensions shifted to start at 0
      p = Position(x - data.xlow, y - data.ylow, z - data.zlow);
    } else {
      p = Position(x, y, z);
    }
    // Track the minimum atom position values for possible shifting later
    for (unsigned int i = 0; i < 3; ++i) {
      if (p[i] < data.min[i]) data.min[i] = p[i];
    }
    data.atoms[n_total] = Atom(id, type, charge, p);
    if (data.has_bonds) {
      data.atoms[n_total].setExtraInfo(mol_id, 0);
      data.atoms[n_total].setExtraInfoNames(0, "mol");
      if (data.bonds[n_mol_read].id1 == -1) {
        ++data.n_bonds;
        data.bonds[n_mol_read].id1 = id;
        data.bonds[n_mol_read].type1 = type;
      } else if (data.bonds[n_mol_read].id2 == -1) {
        data.bonds[n_mol_read].id2 = id;
        data.bonds[n_mol_read].type2 = type;
      } else {
        cout << "Error determining bond information for bond " << mol_id << "\n"
             << "Atom1 = " << data.bonds[n_mol_read].id1 << "; Atom2 = " << data.bonds[n_mol_read].id2 << "\n"
             << "Extra atom = " << id << "\n";
        exit(1);
      }
    ++n_mol_read;
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

  if (!keep_bounds) {
    // Set minimum values to 0, adjust the rest accordingly
    data.xhigh -= data.xlow;
    data.xlow -= data.xlow;
    data.yhigh -= data.ylow;
    data.ylow -= data.ylow;
    data.zhigh -= data.zlow;
    data.zlow -= data.zlow;
  }

  return data;
}

void writeData(const string& infile, const string& outfile, const string& chem_formula, const dataValues& data, const bool& shift) {
  ofstream fout(outfile.c_str());
  checkFileStream(fout, outfile);

  fout << "These " << chem_formula << " coordinates are from the LAMMPS dump file " << infile << ": [id ";
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
         << data.n_bond_types << " bond types\n\n";
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

  fout << "Atoms";
  if (data.has_bonds) {fout << " # full";}
  else if (find(data.data_labels.begin(), data.data_labels.end(), "q") != data.data_labels.end()) {fout << " # charge";}
  else {fout << " # atomic";}
  fout << "\n\n";

  for (unsigned int i = 0; i < data.atoms.size(); ++i) {
    fout.precision(0);
    fout << fixed;
    if (data.id_min == ID_MIN_DEFAULT) {
      fout << data.atoms[i].getId() << " "; // occurs if there are over 100,000,000 atoms, or if the minimum id value is > 100,000,000
    } else {
      fout << data.atoms[i].getId() - data.id_min + 1 << " "; // otherwise we try to start the count at 1
    }

    if (data.has_bonds) {
      fout << data.atoms[i].getExtraInfo()[0] << " ";
    }
    fout << data.atoms[i].getType() << " ";
    if (data.using_charge) {
      fout.precision(6);
      fout << data.atoms[i].getCharge() << " ";
    }
    fout.precision(6);
    if (shift) { // Shifts the atoms to start at (0,0,0), which may not be the case otherwise
      fout << data.atoms[i].getWrapped()[0] - data.min.getX() << " " << data.atoms[i].getWrapped()[1] - data.min.getY() << " "
           << data.atoms[i].getWrapped()[2] - data.min.getZ() << endl;
    } else { // Leave the atom positions as calculated in read_data
      fout << data.atoms[i].getWrapped()[0] << " " << data.atoms[i].getWrapped()[1] << " "
           << data.atoms[i].getWrapped()[2] << endl;
    }
  }

  if (data.has_bonds) {
    fout << "\nBonds\n\n";
    if (data.id_min != ID_MIN_DEFAULT) {
      cerr << "WARNING: id's have been shifted, but bond ids have not been updated. (IDs shifted by " << data.id_min - 1 << ")\n";
    }
    for (size_t i = 0; i < data.bonds.size(); ++i) {
      fout << i + 1 << " " << data.bond_types[i] + 1 << " " << data.bonds[i].id1 << " " << data.bonds[i].id2 << "\n";
    }
  }
  fout.close();
}

int main(int argc, char** argv) {
  string infile, outfile, chem_formula; // input file, output file, chemical formula
  dataValues data;
  bool keep_bounds = false, shift = false;

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
        ("k,keep-bounds", "Flag to keep the bounds as set in the dump file", cxxopts::value<bool>(keep_bounds)->default_value("false"))
        ("s,shift", "Flag to shift the atoms to start at the origin", cxxopts::value<bool>(shift)->default_value("false"))
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
      data = readData(infile, keep_bounds);
      writeData(infile, outfile, chem_formula, data, shift);
    }

  }
  catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

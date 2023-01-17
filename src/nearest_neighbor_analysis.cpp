#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"

struct Box {
  double xlow, xhigh;
  double ylow, yhigh;
  double zlow, zhigh;
  double Lx, Ly, Lz;

  void calculateBoxLengths() {
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;
  }
} box;

struct DataIndices {
  int id = -1;
  int type = -1;
  int x = -1;
  int y = -1;
  int z = -1;
} data_indices;

// This is actually the halfway point between 2nd and 3rd nearest neighbor distances
const std::map<std::string, double> second_nn_distance = {
  {"fcc", 1.11237},
  {"bcc", 1.20711},
  {"sc", 1.70710678}
};

std::vector<Atom> parseDumpFile(std::string);
std::map<std::string, int> parseAtomPropertyIndices(std::string);
std::vector<std::vector<int> > generateNeighborList(std::vector<Atom>, std::string, double);
void analyzeNearestNeighbors(std::vector<Atom>, std::vector<std::vector <int> >);

int main(int argc, char **argv) {
  // The variables we need
  std::string file; // Input file
  std::string structure; // Crystal structure
  double a0; // Lattice parameter
  std::vector<std::string> allowed_structures = {"bcc", "fcc", "sc"};
  std::vector<Atom> atoms; // list of all atoms
  std::vector<std::vector<int> > iatom; // nearest neighbor list for each atom
  try {
    cxxopts::Options options(argv[0], "Analyze the 1st and 2nd nearest neighbors of a LAMMPS dump file");
    options
      .positional_help("file")
      .show_positional_help();
    
    options
      .add_options()
        ("f,file", "dump file", cxxopts::value<std::string>(file), "file")
        ("s,structure", "crystal structure", cxxopts::value<std::string>(structure), "structure")
        ("a,a0", "lattice parameter", cxxopts::value<double>(a0), "a0");
    
    options.parse_positional({"file", "structure", "a0"});
    auto result = options.parse(argc, argv);

    if (find(allowed_structures.begin(), allowed_structures.end(), structure) == allowed_structures.end()) {
      std::cerr << "The crystal structure must be one of 'fcc', 'bcc', or 'sc'\n";
      exit(OPTION_PARSING_ERROR);
    }
    atoms = parseDumpFile(file);
    iatom = generateNeighborList(atoms, structure, a0);
    analyzeNearestNeighbors(atoms, iatom);

  } catch (const cxxopts::OptionException& e) {
    std::cerr << "Error parsing options: " << e.what() << "\n";
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

std::vector<Atom> parseDumpFile(std::string file) {
  std::string str;
  unsigned int n_atoms;
  std::vector<Atom> atoms;
  std::map<std::string, int> property_index_map;
  std::ifstream fin(file.c_str());
  if (fin.fail()) {
    std::cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }

  // We ignore the first three lines
  for (unsigned int i = 0; i < 3; ++i) {
    getline(fin, str);
  }
  // The next line is the number of atoms, we keep this to double check at the end of parsing
  fin >> n_atoms;
  atoms.resize(n_atoms, Atom());
  fin.ignore(); // Clear the newline from the stream

  getline(fin, str); // ITEM: BOX BOUNDS
  // NOTE: This ignores the possibility of tilt factors
  fin >> box.xlow >> box.xhigh >> box.ylow >> box.yhigh >> box.zlow >> box.zhigh;
  box.calculateBoxLengths();

  fin.ignore();
  getline(fin, str);

  property_index_map = parseAtomPropertyIndices(str);

  unsigned int n_atoms_read = 0;
  while (getline(fin, str)) {
    std::stringstream ss(str);
    std::vector<double> data;
    double dummy;
    unsigned int i = 0;
    while (ss >> dummy) {data.push_back(dummy);}

    int id = (int)(data[data_indices.id]);
    int type = (int)(data[data_indices.type]);
    double x = data[data_indices.x];
    double y = data[data_indices.y];
    double z = data[data_indices.z];

    // Adjust everything to start at the origin
    x -= box.xlow;
    y -= box.ylow;
    z -= box.zlow;

    atoms[id - 1] = Atom(id, type, 0 /*assuming no charge*/, Position(x,y,z));
    ++n_atoms_read;
  }

  fin.close();

  if (n_atoms_read != atoms.size()) {
    std::cerr << "Error: number of atoms read does not match number of atoms in the file.\n"
              << "N = " << atoms.size() << " != n_atom_read = " << n_atoms_read << "\n";
    exit(ATOM_COUNT_ERROR);
  }

  return atoms;
}

std::map<std::string, int> parseAtomPropertyIndices(std::string line) {
  std::map<std::string, int> property_index_map;
  std::stringstream ss(line);
  std::string str;
  unsigned int i = 0;

  ss >> str >> str; // Skip the 'ITEM:' and 'ATOMS' fields
  while (ss >> str) {
    property_index_map.insert({str, i++});
  }
  data_indices.id = property_index_map.find("id")->second;
  data_indices.type = property_index_map.find("type")->second;

  std::vector<std::string> valid_x = {"x", "xu", "xs", "xsu"};
  for (std::vector<std::string>::iterator it = valid_x.begin(); it != valid_x.end(); ++it) {
    auto item_it = property_index_map.find(*it);
    if (item_it != property_index_map.end()) {
      data_indices.x = item_it->second;
      break;
    }
  }

  std::vector<std::string> valid_y = {"y", "yu", "ys", "ysu"};
  for (std::vector<std::string>::iterator it = valid_y.begin(); it != valid_y.end(); ++it) {
    auto item_it = property_index_map.find(*it);
    if (item_it != property_index_map.end()) {
      data_indices.y = item_it->second;
      break;
    }
  }

  std::vector<std::string> valid_z = {"z", "zu", "zs", "zsu"};
  for (std::vector<std::string>::iterator it = valid_z.begin(); it != valid_z.end(); ++it) {
    auto item_it = property_index_map.find(*it);
    if (item_it != property_index_map.end()) {
      data_indices.z = item_it->second;
      break;
    }
  }
  
  return property_index_map;
}

std::vector<std::vector<int> > generateNeighborList(std::vector<Atom> atoms, std::string structure, double a0) {
  int ncellx, ncelly, ncellz; // number of cells in each direction
  int idx, idy, idz; // cell number in each direction
  double lcellx, lcelly, lcellz; // length of cells in each direction
  int n_atoms_per_cell; // number of atoms allowed per cell
  double drij_sq, rxij, ryij, rzij; // square of distance, x, y, and z separation.
  std::vector <std::vector <int> > iatom; // cell-linked list
  std::vector <std::vector <std::vector <int> > > icell; // cell index
  std::vector <std::vector <std::vector <std::vector <int> > > > pcell; // atom index in each cell

  double rcut = a0 * second_nn_distance.at(structure);
  double rcut_sq = rcut * rcut;

  // First we generate the number of cells in each direction
  ncellx = (int)(box.Lx / rcut) + 1;
  ncelly = (int)(box.Ly / rcut) + 1;
  ncellz = (int)(box.Lz / rcut) + 1;

  // Length of cells in each direction
  lcellx = box.Lx / ncellx;
  lcelly = box.Ly / ncelly;
  lcellz = box.Lz / ncellz;

  // Minimum number of atoms allowed of 100
  n_atoms_per_cell = std::max((int)(atoms.size() / (double)(ncellx * ncelly * ncellz)), 100);

  // resize the vectors
  icell.resize(ncellx, std::vector <std::vector <int> > // x dimension
              (ncelly, std::vector <int> // y dimension
              (ncellz, 0))); // z dimension
  pcell.resize(ncellx, std::vector <std::vector <std::vector <int> > > // x dimension
              (ncelly, std::vector <std::vector <int> > // y dimension
              (ncellz, std::vector <int> // z dimension
              (n_atoms_per_cell, 0)))); // atom number in cell.
  iatom.resize(n_atoms_per_cell, std::vector <int> (atoms.size(),0));

  // generate the pcell and icell matrices.
  for (size_t i = 0; i < atoms.size(); ++i) {

    // Assign this atom to a cell
    // Rounds towards 0 with the type cast
    idx = (int)(atoms[i].getWrapped()[0] / lcellx); // assign the x cell
    idy = (int)(atoms[i].getWrapped()[1] / lcelly); // assign the y cell
    idz = (int)(atoms[i].getWrapped()[2] / lcellz); // assign the z cell
    // Check if we went out of bounds
    // C++ indexes from 0, so we have to subtract 1 from the maximum value to
    // stay within our memory bounds
    if (idx >= ncellx) idx = ncellx - 1;
    if (idy >= ncelly) idy = ncelly - 1;
    if (idz >= ncellz) idz = ncellz - 1;

    ++icell[idx][idy][idz]; // increase the number of atoms in this cell
    // assign the atom number to this index.
    if (icell[idx][idy][idz] > n_atoms_per_cell) {
      std::cout << "Too many atoms found in cell " << idx << ", " << idy << ", " << idz << "\n";
      exit(1);
    }
    pcell[idx][idy][idz][icell[idx][idy][idz] - 1] = i;
  }

  for (int i = 0; i < ncellx; ++i) { // For each x cell
    for (int j = 0; j < ncelly; ++j) { // For each y cell
      for (int k = 0; k < ncellz; ++k) { // For each z cell
        for (int l = 0; l < icell[i][j][k]; ++l) { // For each atom in this cell
          int id = pcell[i][j][k][l]; // store this atom id
          // Now we check each sub cell around the current one
          for (int ii = -1; ii < 2; ++ii) { // allowed values: -1, 0, and 1
            for (int jj = -1; jj < 2; ++jj) {
              for (int kk = -1; kk < 2; ++kk) {
                int ia = i + ii; // min value: -1.  Max value: number of cells in dimension
                int ja = j + jj;
                int ka = k + kk;
                // Check to make sure we are still in bounds
                // C++ indexes from 0, so we accomodate.
                if (ia >= ncellx) ia = 0;
                if (ja >= ncelly) ja = 0;
                if (ka >= ncellz) ka = 0;
                if (ia < 0) ia = ncellx - 1;
                if (ja < 0) ja = ncelly - 1;
                if (ka < 0) ka = ncellz - 1;

                // Now check each atom in this cell
                for (int m = 0; m < icell[ia][ja][ka]; ++m) {
                  int jd = pcell[ia][ja][ka][m];
                  // If jd <= id, we've already dealt with this interaction
                  if (jd <= id) {continue;}

                  // Now the actual calculations!
                  rxij = atoms[id].getWrapped()[0] - atoms[jd].getWrapped()[0];
                  ryij = atoms[id].getWrapped()[1] - atoms[jd].getWrapped()[1];
                  rzij = atoms[id].getWrapped()[2] - atoms[jd].getWrapped()[2];

                  // Apply PBCs
                  rxij = rxij - std::round(rxij / box.Lx) * box.Lx;
                  ryij = ryij - std::round(ryij / box.Ly) * box.Ly;
                  rzij = rzij - std::round(rzij / box.Lz) * box.Lz;

                  // Now calculate the distance
                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                  // move to the next atom if we're too far away
                  if (drij_sq > rcut_sq) {continue;}

                  // Create the neighbor list
                  iatom[0][id] += 1; //for atom id - number of neighbors
                  if (iatom[0][id] >= n_atoms_per_cell) {
                    n_atoms_per_cell += 100;
                    iatom.resize(n_atoms_per_cell, std::vector <int> (atoms.size(),0));
                  }
                  iatom[(iatom[0][id])][id] = jd; // point to the next atom
                  iatom[0][jd] += 1; // for atom jd
                  if (iatom[0][jd] >= n_atoms_per_cell) {
                    n_atoms_per_cell += 100;
                    iatom.resize(n_atoms_per_cell, std::vector <int> (atoms.size(),0));
                  }
                  iatom[(iatom[0][jd])][jd] = id;
                } // m
              } //kk
            } //jj
          } //ii
        } // l
      } // k
    } // j
  } // i

  return iatom;
}

void analyzeNearestNeighbors(std::vector<Atom> atoms, std::vector<std::vector<int> > iatom) {

  std::ofstream fout("nearest_neighbors.txt");

  double average_n_impurities_in_range = 0;
  for (size_t i = 0; i < atoms.size(); ++i) {
    Atom atom_focus = atoms[i];
    unsigned int n_impurities_in_range = (atom_focus.getType() == 2 ? 1 : 0);
    for (int l = 1; l <= iatom[0][i]; ++l) {
      Atom atom_neighbor = atoms[iatom[l][i]];
      n_impurities_in_range += (atom_neighbor.getType() == 2 ? 1 : 0);
    }
    fout << atom_focus.getId() << " " << n_impurities_in_range / (iatom[0][i] + 1.0) << "\n";
    // divide the number of impurities found by the total number of atoms in this region (including the focus atom)
    average_n_impurities_in_range += n_impurities_in_range / (iatom[0][i] + 1.0);
  }
  average_n_impurities_in_range /= atoms.size();
  
  std::cout << "An average of " << average_n_impurities_in_range * 100 << "% of atoms had an impurity within 2 nearest neighbors\n";
}

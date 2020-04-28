#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include "atom.h"
#include "error_code_defines.h"

#include <cxxopts.hpp>

#define NO_CLUSTER_FOUND 10000000

using namespace std;

struct Header {
  string filename, file_type;
  double xlo = 0.0, xhi = 0.0, ylo = 0.0, yhi = 0.0, zlo = 0.0, zhi = 0.0; // box bounds
  double Lx = 0.0, Ly = 0.0, Lz = 0.0; // box lengths
  double xy = 0.0, xz = 0.0, yz = 0.0; // triclinic factors

  int id_index = -1, type_index = -1, charge_index = -1, x_index = -1;
  int y_index = -1, z_index = -1;
  int xu_index = -1, yu_index = -1, zu_index = -1;

  vector <string> data_types;
  bool has_charge = false;

  int N = -1, n_types = -1; // number of atoms/atom types

  // map from the data_list index to the extra_info vector list index
  // the string is the data type, the first unsigned int is the index number
  // for original list (as it appears in the dump or dat file), and the second
  // unsigned int is the index for the extra_info vector (the member from the
  // Atom class)
  map <string, pair<unsigned int, unsigned int> > indexConverter;
};

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

// Calculate the rounded value of x
double anInt(double x) {
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) {x += 0.5;}
  if (x < 0.0) {x -= 0.5;}
  temp = (int)(x);
  return (double)(temp);
}

void showCommandList() {
  cout << "\nr,read        Specify a new file to analyze\n"
       << "H,header      Show the header information\n"
       << "s,select      Select atom by various criteria\n"
       << "p,print       Print selected atom info\n"
       << "c,cluster     Perform a cluster analysis of the snapshot\n"
       << "q,quit,exit   Exit the software\n"
       << "h,help        Show this list of accepted commands\n\n";
}

void welcome() {
  cout << "Welcome to atom playground!  This software allows for in-depth\n"
       << "analysis of a single snapshot of an atomic structure.\n"
       << "Here is the list of accepted commands: \n";

  showCommandList();
}

void printHeaderInfo(const Header& header) {
  cout << "Box lengths:\n"
       << "  x: " << header.xlo << " - " << header.xhi << " (" << header.Lx << ")\n"
       << "  y: " << header.ylo << " - " << header.yhi << " (" << header.Ly << ")\n"
       << "  z: " << header.zlo << " - " << header.zhi << " (" << header.Lz << ")\n"
       << "Triclinic factors: " << header.xy << " " << header.xz << " " << header.yz << " xy xz yz\n"
       << "N: " << header.N << "\n"
       << "Number of atom types: " << header.n_types << "\n"
       << "Data types: \n  ";
  for (unsigned int i = 0; i < header.data_types.size(); ++i) {
    cout << header.data_types[i] << " ";
  }
  cout << endl;
}

Header getHeaderData(const string& data_type, ifstream& fin) {
  Header header;
  string str; // junk variable
  if (data_type.compare("dump") == 0) {
    getline(fin, str); // ITEM: TIMESTEP
    getline(fin, str); // actual timestep
    getline(fin, str); // ITEM: NUMBER OF ATOMS
    fin >> header.N; // number of atoms
    fin.ignore();
    getline(fin, str); // ITEM: BOX BOUNDS
    if (str.find("xy") == string::npos) {
      fin >> header.xlo >> header.xhi;
      fin >> header.ylo >> header.yhi;
      fin >> header.zlo >> header.zhi;
    } else {
      fin >> header.xlo >> header.xhi >> header.xy;
      fin >> header.ylo >> header.yhi >> header.xz;
      fin >> header.zlo >> header.zhi >> header.yz;
    }
    fin.ignore();
    getline(fin, str); // ITEM: ATOMS <data types>
    stringstream ss(str);
    string dummy;
    ss >> dummy >> dummy; // ITEMS: and ATOMS
    while (ss >> dummy) {
      if (dummy.compare("charge") == 0 || dummy.compare("q") == 0) {header.has_charge = true;}
      header.data_types.push_back(dummy);
    }
  } else if (data_type.compare("dat") == 0) {
    getline(fin, str); // first (comment) line
    if (str.find("[") == string::npos || str.find("]") == string::npos) {
      cout << "Unable to identify data types.  The first line in the file should\n"
           << "have one set of [] with the data types listed inside, separated by spaces.\n";
      exit(FILE_FORMAT_ERROR);
    }
    stringstream ss(str.substr(str.find("[") + 1, str.find("]", str.find("[")) - str.find("[") - 1));
    string dummy;
    while (ss >> dummy) {
      if (dummy.compare("charge") == 0 || dummy.compare("q") == 0) {header.has_charge = true;}
      header.data_types.push_back(dummy);
    }
    fin >> header.N >> str; // number of atoms
    fin >> header.n_types >> str >> str; // number of atom types
    fin >> header.xlo >> header.xhi >> str >> str;
    fin >> header.ylo >> header.yhi >> str >> str;
    fin >> header.zlo >> header.zhi >> str >> str;
    fin.ignore();
    getline(fin, str); // blank line or triclinic tilt factors
    if (str.find("xy") != string::npos) {// Triclinic tilt factors
      string tmp_str;
      stringstream s2(str);
      ss >> header.xy >> header.xz >> header.yz >> tmp_str >> tmp_str;
      getline(fin, str); // blank line
    }
    getline(fin, str); // Atoms
    getline(fin, str); // blank line
  }

  header.Lx = header.xhi - header.xlo;
  header.Ly = header.yhi - header.ylo;
  header.Lz = header.zhi - header.zlo;

  return header;
}

vector <Atom> getAtomData(ifstream& fin, Header& header) {
  int n_atoms_read = 0;
  string str;
  vector <double> data;
  vector <Atom> atoms (header.N, Atom());
  int atom_id, type;
  double charge, x, y, z, xu, yu, zu;


  for (unsigned int i = 0; i < header.data_types.size(); ++i) {
    string name = header.data_types[i];
    if (name.compare("ID") == 0 || name.compare("id") == 0) {header.id_index = i;}
    else if (name.compare("type") == 0) {header.type_index = i;}
    else if (name.compare("charge") == 0 || name.compare("q") == 0) {header.charge_index = i;}
    else if (name.compare("x") == 0) {header.x_index = i;}
    else if (name.compare("y") == 0) {header.y_index = i;}
    else if (name.compare("z") == 0) {header.z_index = i;}
    else if (name.compare("xu") == 0) {header.xu_index = i;}
    else if (name.compare("yu") == 0) {header.yu_index = i;}
    else if (name.compare("zu") == 0) {header.zu_index = i;}
  }

  while (getline(fin, str)) {
    stringstream ss(str);
    double dummy;
    vector <double> additional_data;
    while (ss >> dummy) {data.push_back(dummy);}
    for (unsigned int i = 0; i < data.size(); ++i) {
      if (i == header.id_index) {atom_id = (int)(data[i]);}
      else if (i == header.type_index) {
        type = (int)(data[i]);
        if (header.file_type.compare("dump") == 0) {
          if (type > header.n_types) {header.n_types = type;}
        } else if (header.file_type.compare("dat") == 0) {
          if (type > header.n_types) {
            cout << "Error: incorrect number of atom types.\n"
                 << "This atom type: " << type << " specified atom types: " << header.n_types << "\n";
            exit(ATOM_TYPE_ERROR);
          }
        }
      }
      else if (header.has_charge && i == header.charge_index) {charge = data[i];}
      else if (i == header.x_index) {x = data[i] - header.xlo;}
      else if (i == header.y_index) {y = data[i] - header.ylo;}
      else if (i == header.z_index) {z = data[i] - header.zlo;}
      else if (i == header.xu_index) {xu = data[i] - header.xlo;}
      else if (i == header.yu_index) {yu = data[i] - header.ylo;}
      else if (i == header.zu_index) {zu = data[i] - header.zlo;}
      else {
        additional_data.push_back(data[i]);
        pair <unsigned int, unsigned int> value = make_pair(i, additional_data.size() - 1);
        pair <string, pair <unsigned int, unsigned int> > key_value_pair = make_pair(header.data_types[i], value);
        header.indexConverter.insert(key_value_pair);
      }
    }

    data.clear();

    Position p(x,y,z);
    atoms[atom_id - 1] = Atom(atom_id, type, charge, p);

    if (header.xu_index != -1) {p.setX(xu);}
    if (header.yu_index != -1) {p.setY(yu);}
    if (header.zu_index != -1) {p.setZ(zu);}
    atoms[atom_id - 1].setUnwrapped(p);
    atoms[atom_id - 1].setExtraInfo(additional_data);
    ++n_atoms_read;
  }

  if (n_atoms_read != atoms.size()) {
    cout << "Error reading file: number of atoms read versus number of atoms assigned do not match:\n"
         << n_atoms_read << " != " << atoms.size() << "\n";
  }
  return atoms;
}

string getFileName() {
  string filename;
  cout << "Please enter the filename: ";
  cin  >> filename;
  return filename;
}

pair <vector <Atom>, Header> readFile(const string& filename) {
  vector <Atom> atom_tmp;
  Header header_tmp;

  ifstream fin(filename);
  checkFileStream(fin, filename);

  if (filename.substr(filename.length() - 4) == "dump") {
    header_tmp = getHeaderData("dump", fin);
    header_tmp.file_type = "dump";
    atom_tmp = getAtomData(fin, header_tmp);
  }
  else if (filename.substr(filename.length() - 3) == "dat") {
    header_tmp = getHeaderData("dat", fin);
    header_tmp.file_type = "dat";
    atom_tmp = getAtomData(fin, header_tmp);
  } else {
    bool known = false;
    string file_type;
    while (!known) {
      cout << "Unrecognized file type.  Please specify dump or dat: ";
      cin  >> file_type;
      if (file_type.compare("dump") == 0) {
        header_tmp = getHeaderData("dump", fin);
        header_tmp.file_type = "dump";
        atom_tmp = getAtomData(fin, header_tmp);
        known = true;
      }
      else if (file_type.compare("dat") == 0) {
        header_tmp = getHeaderData("dat", fin);
        header_tmp.file_type = "dat";
        atom_tmp = getAtomData(fin, header_tmp);
        known = true;
      }
    }
  }

  fin.close();
  header_tmp.filename = filename;
  return make_pair(atom_tmp,header_tmp);
}

vector <Atom> getAtomsByType(const vector <Atom>& atoms) {
  string atom_types_string;
  vector <int> atom_types;
  vector <Atom> selected_atoms;

  cout << "Please enter the atom type(s) as a space-separated list: ";
  cin.ignore();
  getline(cin, atom_types_string);
  stringstream ss(atom_types_string);

  int type_tmp;
  while (ss >> type_tmp) {atom_types.push_back(type_tmp);}
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    if (find(atom_types.begin(), atom_types.end(), atoms[i].getType()) != atom_types.end()) {
      selected_atoms.push_back(atoms[i]);
    }
  }

  return selected_atoms;
}

vector <int> getAtomIdsFromFile() {
  string file, str;
  vector <int> atom_ids;

  cout << "Please enter the file name containing the list of atom IDs: ";
  cin  >> file;

  ifstream fin(file);
  checkFileStream(fin, file);

  while (getline(fin, str)) {
    int id;
    stringstream ss(str);
    ss >> id;
    atom_ids.push_back(id);
  }

  fin.close();

  return atom_ids;
}

vector <int> getAtomIdsByRange(unsigned int num_atoms) {
  string range_string, str;
  vector <string> ranges;
  vector <int> atom_ids;

  cout << "Please enter the range(s) for atom ids (e.g. BEGIN-10 20-25 27 100-END)\n"
       << "Note that BEGIN = 1, and END = the number of atoms: ";
  cin  >> range_string;

  stringstream ss(range_string);
  while (ss >> str) {ranges.push_back(str);}

  for (unsigned int i = 0; i < ranges.size(); ++i) {
    size_t minus_pos = ranges[i].find("-");
    if (minus_pos == string::npos) {
      int single_id;
      stringstream s2(ranges[i]);
      s2 >> single_id;
      atom_ids.push_back(single_id);
    } else {// we found a minus sign, indicating a range
      unsigned int first_index = 0, last_index = 0;
      string temp = ranges[i].substr(0,minus_pos);

      if (temp.compare("BEGIN") == 0) {first_index = 1;}
      else {
        stringstream s3(temp);
        s3 >> first_index;
      }

      temp = ranges[i].substr(minus_pos + 1);
      if (temp.compare("END") == 0) {last_index = num_atoms;}
      else {
        stringstream s3(temp);
        s3 >> last_index;
      }

      for (unsigned int j = first_index; j <= last_index; ++j) {
        atom_ids.push_back(j);
      }
    }
  }

  return atom_ids;
}

vector <Atom> getAtomsById(const vector <Atom>& atoms, const vector <int>& atom_ids) {
  vector <Atom> selected_atoms (atom_ids.size(), Atom());

  for (unsigned int i = 0; i < atom_ids.size(); ++i) {
    selected_atoms[i] = atoms[atom_ids[i] - 1];
  }

  return selected_atoms;
}

vector <Atom> getSelectedAtoms(const vector <Atom>& atoms) {
  char atom_specification;
  vector <Atom> selected_atoms;

  cout << "Would you like to specify atoms by type(t), file(f), or range(r)?: ";
  cin  >> atom_specification;

  while (atom_specification != 't' && atom_specification != 'f' && atom_specification != 'r') {
    // If we have a cin error, clear the error state, wipe the input, and try again
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(),'\n');
    cout << "Please enter either \'t\', \'f\', or \'r\': ";
    cin  >> atom_specification;
  }

  if (atom_specification == 't') {
    selected_atoms = getAtomsByType(atoms);
  } else {
    vector <int> atom_ids;
    if (atom_specification == 'f') {
      atom_ids = getAtomIdsFromFile();
    } else {
      atom_ids = getAtomIdsByRange(atoms.size());
    }

    selected_atoms = getAtomsById(atoms, atom_ids);
  }

  cout << "Selected " << selected_atoms.size() << " atoms (" << (double)(selected_atoms.size()) / (double)(atoms.size()) * 100.0 << "%).\n";

  return selected_atoms;
}

void printSelectedAtomInfo(const vector <Atom>& selected_atoms, Header& header) {
  string desired_data_string, filename;
  vector <bool> print_data(header.data_types.size(), false);

  cout << "The following data is available: \n  ";
  for (unsigned int i = 0; i < header.data_types.size(); ++i) {
    cout << header.data_types[i] << " ";
  }
  cout << endl;

  cout << "Please enter the data types to print as they appear above: ";
  cin.ignore();
  getline(cin, desired_data_string);

  stringstream ss(desired_data_string);
  string dummy;

  while (ss >> dummy) {
    if (dummy.compare("id") == 0 || dummy.compare("ID") == 0) {print_data[header.id_index] = true;}
    else if (dummy.compare("type") == 0) {print_data[header.type_index] = true;}
    else if (dummy.compare("charge") == 0 || dummy.compare("q") == 0) {print_data[header.charge_index] = true;}
    else if (dummy.compare("x") == 0) {print_data[header.x_index] = true;}
    else if (dummy.compare("y") == 0) {print_data[header.y_index] = true;}
    else if (dummy.compare("z") == 0) {print_data[header.z_index] = true;}
    else if (dummy.compare("xu") == 0) {print_data[header.xu_index] = true;}
    else if (dummy.compare("yu") == 0) {print_data[header.yu_index] = true;}
    else if (dummy.compare("zu") == 0) {print_data[header.zu_index] = true;} else {
      if (header.indexConverter.find(dummy) != header.indexConverter.end()) {
        print_data[header.indexConverter[dummy].first] = true;
      } else {
        cout << "Unrecognized data type specified.\n";
      }
    }
  }

  cout << "Please specify the output filename: ";
  cin  >> filename;

  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);

  fout << "# ";
  for (unsigned int i = 0; i < print_data.size(); ++i) {
    if (print_data[i]) {fout << header.data_types[i] << " ";}
  }

  fout << endl;

  for (unsigned int i = 0; i < selected_atoms.size(); ++i) {
    for (unsigned int j = 0; j < print_data.size(); ++j) {
      if (print_data[j]) {
        if (j == header.id_index) {fout << selected_atoms[i].getId() << " ";}
        else if (j == header.type_index) {fout << selected_atoms[i].getType() << " ";}
        else if (j == header.charge_index) {fout << selected_atoms[i].getCharge() << " ";}
        else if (j == header.x_index) {fout << selected_atoms[i].getWrapped().getX() + header.xlo << " ";}
        else if (j == header.y_index) {fout << selected_atoms[i].getWrapped().getY() + header.ylo << " ";}
        else if (j == header.z_index) {fout << selected_atoms[i].getWrapped().getZ() + header.zlo << " ";}
        else if (j == header.xu_index) {fout << selected_atoms[i].getUnwrapped().getX() + header.xlo << " ";}
        else if (j == header.yu_index) {fout << selected_atoms[i].getUnwrapped().getY() + header.ylo << " ";}
        else if (j == header.zu_index) {fout << selected_atoms[i].getUnwrapped().getZ() + header.zlo << " ";} else {
          fout << selected_atoms[i].getExtraInfo()[header.indexConverter[header.data_types[j]].second] << " ";
        }
      }
    }
    fout << endl;
  }

  fout.close();
}

vector <vector <int> > generateCellLinkedList(const vector <Atom>& atoms,
                                              const double& rcut,
                                              const double& lx,
                                              const double& ly,
                                              const double& lz,
                                              const vector <int> types_to_keep) {
  int ncellx, ncelly, ncellz;
  int idx, idy, idz;
  double lcellx, lcelly, lcellz;
  int n_atoms_per_cell;
  double drij_sq, rxij, ryij, rzij;
  vector <vector <int> > iatom;
  vector <vector <vector <int> > > icell;
  vector <vector <vector <vector <int> > > > pcell;

  double rcut_sq = rcut * rcut;

  ncellx = (int)(lx / rcut) + 1;
  ncelly = (int)(ly / rcut) + 1;
  ncellz = (int)(lz / rcut) + 1;

  lcellx = lx / ncellx;
  lcelly = ly / ncelly;
  lcellz = lz / ncellz;

  n_atoms_per_cell = max((int)(atoms.size() / (double)(ncellx * ncelly * ncellz)), 100);

  icell.resize(ncellx, vector <vector <int> > // x dimension
              (ncelly, vector <int> // y dimension
              (ncellz, 0))); // z dimension
  pcell.resize(ncellx, vector <vector <vector <int> > > // x dimension
              (ncelly, vector <vector <int> > // y dimension
              (ncellz, vector <int> // z dimension
              (n_atoms_per_cell, 0)))); // atom number in cell.
  iatom.resize(n_atoms_per_cell, vector <int> (atoms.size(),0));

  // generate the pcell and icell matrices.
  for (size_t i = 0; i < atoms.size(); ++i) {
    // If we aren't interested in this atom type, skip it.
    if (find(types_to_keep.begin(), types_to_keep.end(), atoms[i].getType()) == types_to_keep.end()) {continue;}
    // Assign this atom to a cell
    // Rounds towards 0 with a type cast
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
                  rxij = rxij - anInt(rxij / lx) * lx;
                  ryij = ryij - anInt(ryij / ly) * ly;
                  rzij = rzij - anInt(rzij / lz) * lz;

                  // Now calculate the distance
                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                  // move to the next atom if we're too far away
                  if (drij_sq > rcut_sq) {continue;}

                  // Create the neighbor list
                  iatom[0][id] += 1; //for atom id - number of neighbors
                  if (iatom[0][id] >= n_atoms_per_cell) {
                    n_atoms_per_cell += 100;
                    iatom.resize(n_atoms_per_cell, vector <int> (atoms.size(),0));
                  }
                  iatom[(iatom[0][id])][id] = jd; // point to the next atom
                  iatom[0][jd] += 1; // for atom jd
                  if (iatom[0][jd] >= n_atoms_per_cell) {
                    n_atoms_per_cell += 100;
                    iatom.resize(n_atoms_per_cell, vector <int> (atoms.size(),0));
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

void clusterAnalysis(vector <Atom>& atoms, const Header& header) {
  string str, structure;
  double a0;
  vector <int> cluster_types;
  map <string, double> first_nn_distances = {{"fcc", 0.853553}, {"bcc", 0.933013}, {"sc", 1.2071068}};

  cout << "There are " << header.n_types << " atom types. Enter the atom type number(s) to examine for clusters, separated by a space: ";
  cin.ignore();
  getline(cin, str);
  stringstream ss(str);
  int tmp;
  while (ss >> tmp) {cluster_types.push_back(tmp);}

  cout << "Enter the crystal structure (fcc, bcc, sc): ";
  cin  >> structure;

  while (structure.compare("fcc") != 0 && structure.compare("bcc") != 0 && structure.compare("sc") != 0) {
    cout << "Please enter fcc, bcc, or sc: ";
    cin  >> structure;
  }

  cout << "Enter the lattice parameter in angstroms: ";
  cin  >> a0;

  while (a0 < 0) {
    cout << "Please enter a number greater than 0: ";
    cin  >> a0;
  }


  // generate a neighbor list with atoms that are only within the desired cutoff distance from each other
  vector <vector <int> > iatom = generateCellLinkedList(atoms, first_nn_distances[structure] * a0, header.Lx, header.Ly, header.Lz, cluster_types);
  vector <set <Atom> > clusters;
  int num_atoms = 0;
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    // If the current atom type is not listed as one to do cluster analysis for, skip it
    if (find(cluster_types.begin(), cluster_types.end(), atoms[i].getType()) == cluster_types.end()) {continue;}
    ++num_atoms;
    int cluster_index = NO_CLUSTER_FOUND;
    if (atoms[i].getMark() == 0) {
      for (int l = 1; l <= iatom[0][i]; ++l) { // for each neighbor
        if (atoms[iatom[l][i]].getMark() != 0 && atoms[iatom[l][i]].getMark() < cluster_index) {
          cluster_index = atoms[iatom[l][i]].getMark();
        }
      }
      if (cluster_index == NO_CLUSTER_FOUND) {cluster_index = clusters.size() + 1;}
      set <Atom> tmp;
      atoms[i].setMark(cluster_index);
      tmp.insert(atoms[i]);
      for (int l = 1; l <= iatom[0][i]; ++l) {
        atoms[iatom[l][i]].setMark(cluster_index);
        tmp.insert(atoms[iatom[l][i]]);
      }
      if (cluster_index > clusters.size()) {clusters.push_back(tmp);}
      else {clusters[cluster_index - 1].insert(tmp.begin(), tmp.end());}
    } // atom has not already been assigned to cluster
  } // all atoms

  cout << clusters.size() << " clusters were found, containing " << num_atoms << " atoms\n";

  string outfile = header.filename.substr(0,header.filename.find("." + header.file_type)) + "_cluster_analysis.txt";
  ofstream fout (outfile.c_str());

  fout << "# Cluster_id num_atom_in_cluster\n";
  for (unsigned int i = 0; i < clusters.size(); ++i)
  {
    fout << i + 1 << " " << clusters[i].size() << "\n";
  }
  fout.close();

  ofstream fouttmp ("tmp.dat");
  for (unsigned int i = 0; i < atoms.size(); ++i) {
    if (find(cluster_types.begin(), cluster_types.end(), atoms[i].getType()) == cluster_types.end()) {continue;}
    fouttmp << atoms[i].getId() << " " << atoms[i].getType() << " " << atoms[i].getCharge() << " "
            << atoms[i].getWrapped()[0] << " " << atoms[i].getWrapped()[1] << " "
            << atoms[i].getWrapped()[2] << " " << atoms[i].getMark() << "\n";
  }
}

void playground(string filename = "None") {
  bool leave_playground = false;
  string command;
  vector <Atom> atoms, selected_atoms;
  Header header;

  if (filename.compare("None") != 0) {
    pair<vector <Atom>, Header> tmp = readFile(filename);
    atoms = tmp.first;
    header = tmp.second;
  }

  while (!leave_playground) {
    cout << "<< ";
    cin  >> command;

    if (command.compare("quit") == 0 || command.compare("exit") == 0 || command.compare("q") == 0) {
      leave_playground = true;
    } else if (command.compare("help") == 0 || command.compare("h") == 0) {
      showCommandList();
    } else if (command.compare("read") == 0 || command.compare("r") == 0) {
      filename = getFileName();
      pair <vector <Atom>, Header> tmp = readFile(filename);
      atoms = tmp.first;
      header = tmp.second;
    } else if (command.compare("header") == 0 || command.compare("H") == 0) {
      if (atoms.size() == 0) {
        cout << "No atomic information has been read yet!  Specify a file to read with r or read.\n";
      }
      else {printHeaderInfo(header);}
    } else if (command.compare("select") == 0 || command.compare("s") == 0) {
      if (atoms.size() == 0) {
        cout << "No atomic information has been read yet!  Specify a file to read with r or read.\n";
      } else {selected_atoms = getSelectedAtoms(atoms);}
    } else if (command.compare("print") == 0 || command.compare("p") == 0) {
      if (selected_atoms.size() == 0) {selected_atoms = atoms;}
      printSelectedAtomInfo(selected_atoms, header);
    } else if (command.compare("cluster") == 0 || command.compare("c") == 0) {
      clusterAnalysis(atoms, header);
    } else {
      cout << "Unrecognized command!\n";
    }
  }
}

int main(int argc, char** argv) {
  string file;

  try
  {
    cxxopts::Options options(argv[0], "Playground for extracting and examining atom properties");
    options
      .positional_help("file")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Data file to play with", cxxopts::value<string>(file), "data_file")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    welcome();

    if (result.count("file")) {
      playground(file);
    } else {
      playground();
    }
  }
  catch (const cxxopts::OptionException& e) {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <set>
#include <cmath> // for sin and cos
#include <numeric> // for accumulate
#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"

using namespace std;

#define PI 3.14159265358979

// Global values
static bool is_sphere = false;
static bool marked = true;
static bool rotated = true;
static bool removed = true;
static bool printedInputHelp = false;
static bool type_one_only = false;
static map <string, double> first_nn_distances = {
  {"fcc", 0.707106781}, // 1/sqrt(2)
  {"fluorite", 0.707106781}, // same as fcc for the larger sublattice
  {"bcc", 0.866025404}, // sqrt(3)/2
  {"sc", 1.0},
  {"hcp", 0.707106781}, // same as fcc
  {"diamond", 0.433012702} // sqrt(3)/4
};

bool mapValueCmp(pair<pair<int,int>,double> a, pair<pair<int,int>,double> b);

struct neighborData
{
  pair <int, int> ids; // the ids of the two atoms
  pair <int, int> types; // the types of the two atoms
  double distance; // the distance between the two atoms
};

struct boxData
{
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double xy = 0.0, xz = 0.0, yz = 0.0;
  double Lx, Ly, Lz;
  bool is_triclinic = false;

  void calculateBoxLengths()
  {
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;
  }
};

struct inputData
{
  // strictly from input file
  string data_file, output_basename, molecule, crystal_structure;
  double r_grain = 0, theta = 0, a0 = 0;
  map <pair <int, int>, double> rcut;

  // calculated from input
  double r_grain_sq = 0, costheta = 0, sintheta = 0, r_cut_max = 0, r_cut_max_sq = 0;
  int n_types = 0;
  map <pair <int, int>, double> rcut_sq;
  //vector <int> element_ratios;
  map <string, int> element_ratios; // map of element to its ratio
  map <int, string> number_to_element;
  bool has_charge = false;
  boxData box;

  void calculateRCutSq()
  {
    r_cut_max_sq = r_cut_max * r_cut_max;
    if (a0 * a0 > r_cut_max_sq) {r_cut_max_sq = a0 * a0;} // This is primarily for neighbor lists
  }
  void calculateRGrainSq() {r_grain_sq = r_grain * r_grain;}
  void calculateRCutMax()
  {
    r_cut_max = (*max_element(rcut.begin(), rcut.end(), mapValueCmp)).second;
    calculateRCutSq();
  }
  void calculateTrig() {costheta = cos(theta * PI / 180.0); sintheta = sin(theta * PI / 180.0);}
  void calculateNumElements()
  {
    set <string> elements;
    for (unsigned int i = 0; i < molecule.size(); ++i)
    {
      if (isdigit(molecule[i]) || islower(molecule[i])) {continue;}
      if (isupper(molecule[i]))
      {
        if (i + 1 >= molecule.size()) {elements.insert(molecule.substr(i,1));}
        else if (islower(molecule[i + 1])) {elements.insert(molecule.substr(i, 2));}
        else if (isupper(molecule[i + 1]) || isdigit(molecule[i + 1])) {elements.insert(molecule.substr(i, 1));}
        else {cout << "ERROR!\n";}
      }
    }
    n_types = elements.size();
    //element_ratios.resize(n_types, 0);
  }

  const bool validate() const
  {
    vector <string> structures = {"fcc", "bcc", "sc", "diamond", "fluorite"};
    if (data_file.empty()) {return false;}
    if (molecule.empty()) {return false;}
    if (crystal_structure.empty() || find(structures.begin(), structures.end(), crystal_structure) == structures.end()) {return false;}
    if (r_grain <= 0) {return false;}
    if (((n_types + 1) * n_types) / 2 != rcut.size()) {return false;}

    return true;
  }
};

bool mapValueCmp(pair<pair<int,int>,double> a, pair<pair<int,int>,double> b)
{
  return (a.second < b.second);
}

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

template <typename K, typename V>
int findByValue(map<K,V> mapOfElemen, V value)
{
  typename map<K,V>::iterator it = mapOfElemen.begin();
  while (it != mapOfElemen.end())
  {
    if (it->second == value)
    {
      return it->first;
    }
    ++it;
  }
  return -1;
}

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

bool compareNeighbors(const neighborData &a, const neighborData &b)
{
  return a.distance < b.distance; // return the distance that is shorter
}

int findElementSubscript(const string& molecule, const size_t& index)
{
  int value;
  size_t last_digit_pos = molecule.find_first_not_of("0123456789", index);
  stringstream ss(molecule.substr(index, last_digit_pos - index));
  ss >> value;

  return value;
}

void printInputFileHelp()
{
  cout << "The input file must contain the following items in order on one line:\n"
       << "1) The data file to process\n"
       << "2) The output file(s) basename\n"
       << "3) The desired grain radius\n"
       << "4) The desired misorientation angle\n"
       << "5) The chemical composition of the structure (i.e. Fe or FeO2)\n"
       << "6) The crystal structure of the data file (bcc, fcc, sc, fluorite)\n"
       << "7) The lattice parameter in Angstroms\n"
       << "8) Cutoff distance as a percentage of 1NN distance. Input as the minimum\n"
       << "allowed distance between nearest neighbor atoms.\n\n"
       << "Note that the number of cutoff radii is dependent on the number of atom\n"
       << "interactions. For example, if there are two atom types, there are the two\n"
       << "same-element interactions, as well as the interaction between the two\n"
       << "elements. Include all the relevant cutof radii one element at a time, i.e.\n"
       << "the cutoff radii for a 3-element system would be input as 1-1, 1-2, 1-3,\n"
       << "2-2, 2-3, 3-3. Maintain consistent numbering with the data file. If all cutoff\n"
       << "values are the same, just a single value will work.\n";
}

void determineElementRatios(inputData& input)
{
  int elem_num = 1;

  for (size_t i = 0; i < input.molecule.size(); ++i)
  {
    if (isdigit(input.molecule[i]) || islower(input.molecule[i])) {continue;}
    if (isupper(input.molecule[i])) // current character is upper case
    {
      if (islower(input.molecule[i + 1])) // next character is lowercase
      {
        string element = input.molecule.substr(i, 2);
        int element_index = findByValue(input.number_to_element, element);
        if (element_index < 0)
        {
          // second character after current is a number
          if (isdigit(input.molecule[i + 2])) {input.element_ratios[element] = findElementSubscript(input.molecule, i + 2);}
          else {input.element_ratios[element] = 1;} // second character after current is not a number
          if (input.number_to_element.insert(pair<int,string>(elem_num, element)).second) {++elem_num;}
          else {cout << "Error inserting element " << element << endl;}
        }
        else
        {
          if (isdigit(input.molecule[i + 2])) {input.element_ratios[input.number_to_element[element_index]] += findElementSubscript(input.molecule, i + 2);}
          else {++input.element_ratios[input.number_to_element[element_index]];}
        }
      }
      else // next character is either upper or a digit
      {
        string element = input.molecule.substr(i, 1);
        int element_index = findByValue(input.number_to_element, element);
        if (element_index < 0)
        {
          if (isdigit(input.molecule[i + 1])) {input.element_ratios[element] = findElementSubscript(input.molecule, i + 1);}
          else {input.element_ratios[element] = 1;}
          if (input.number_to_element.insert(pair<int,string>(elem_num, element)).second) {++elem_num;}
          else {cout << "Error inserting element " << element << endl;}
        }
        else
        {
          if (isdigit(input.molecule[i + 1])) {input.element_ratios[element] = findElementSubscript(input.molecule, i + 1);}
          else {++input.element_ratios[input.number_to_element[element_index]];}
        }
      }
    }
  }

  if (input.n_types != input.number_to_element.size())
  {
    cout << "Error determining element types. Input file n_types = " << input.n_types
         << " != found elements = " << input.number_to_element.size() << endl;
    exit(ELEMENT_COUNT_ERROR);
  }
}

void checkGrainSize(const double& r_grain, const double& val)
{
  if (r_grain * 2.0 > val)
  {
    cout << "Error: grain diameter = " << r_grain * 2.0
         << " >= boundary = " << val << endl;
    exit(GRAIN_TOO_LARGE_ERROR);
  }
}

vector <Atom> readDataFile(inputData& input)
{
  string str; // junk string variable
  int N, n_types, n_total = 0; // number of atoms, number of atom types, number of atoms read
  vector <Atom> atoms;
  int atom_id, type;
  double charge, x, y, z;

  ifstream fin(input.data_file.c_str());
  checkFileStream(fin, input.data_file);

  getline(fin, str);
  transform(str.begin(), str.end(), str.begin(), ::tolower);
  if (str.find("charge") != string::npos) {input.has_charge = true;}
  fin >> N >> str;
  fin >> n_types >> str >> str;
  if (n_types != input.n_types)
  {
    cout << "WARNING: Atom types calculated vs atom types in data file do not match.\n"
         << "Calculated = " << input.n_types << " != data file = " << n_types << endl;
  }

  fin >> input.box.xlow >> input.box.xhigh >> str >> str
      >> input.box.ylow >> input.box.yhigh >> str >> str
      >> input.box.zlow >> input.box.zhigh >> str >> str;

  input.box.calculateBoxLengths();

  checkGrainSize(input.r_grain, input.box.Lx);
  checkGrainSize(input.r_grain, input.box.Ly);
  if (is_sphere) {checkGrainSize(input.r_grain, input.box.Lz);}

  fin.ignore();
  getline(fin, str); // read the extra stuff
  if (str.find("xy") != string::npos) // triclinic system check
  {
    stringstream tmp(str);
    tmp >> input.box.xy >> input.box.xz >> input.box.yz >> str >> str >> str;
    fin >> str;
    input.box.is_triclinic = true;
    fin.ignore();
  }
  else {getline(fin, str);} // line up the input file structures

  atoms.resize(N, Atom());
  getline(fin, str); // blank line before the data
  while (getline(fin, str))
  {
    vector <double> data;
    stringstream ss(str);
    double dummy;
    while (ss >> dummy)
    {
      data.push_back(dummy);
    }
    atom_id = (int)(data[0]);
    type = (int)(data[1]);

    switch (data.size())
    {
      // atom has charge
      case 6: charge = data[2];
              x = data[3]; y = data[4]; z = data[5];
              break;
      // atom does not have charge
      case 5: charge = 0.0;
              x = data[2]; y = data[3]; z = data[4];
              break;
      default: cout << "Unrecognized file format.  Expected format: id type charge* x y z.  * = optional.\n";
               exit(FILE_FORMAT_ERROR);
    }

    data.clear();
    ++n_total;

    if (type > input.n_types)
    {
      cout << "Error: Atom type = " << type << " is greater than the number of types specified = "
           << input.n_types << endl;
      exit(ATOM_TYPE_ERROR);
    }

    Position p(x, y, z);
    atoms[atom_id - 1] = Atom(atom_id, type, charge, p);
  }
  fin.close();

  if (n_total != N)
  {
    cout << "Error: n_total = " << n_total << " != N = " << N << endl;
    exit(ATOM_COUNT_ERROR);
  }

  return atoms;
}

void rotateAtoms(const inputData& input, vector <Atom>& atoms)
{
  double x1, y1, z1, xtemp, ytemp;

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    // change the origin to the center of the simulation for rotating the atoms
    x1 = atoms[i].getWrapped().getX() - input.box.Lx / 2.0;
    y1 = atoms[i].getWrapped().getY() - input.box.Ly / 2.0;
    z1 = atoms[i].getWrapped().getZ() - input.box.Lz / 2.0;

    // Only rotate atoms that are within the radius r_grain from the center
    if (is_sphere)
    {
      if (x1 * x1 + y1 * y1 + z1 * z1 <= input.r_grain_sq)
      {
        xtemp = x1 * input.costheta - y1 * input.sintheta;
        ytemp = x1 * input.sintheta + y1 * input.costheta;
        x1 = xtemp;
        y1 = ytemp;
      }
    }
    else
    {
      if (x1 * x1 + y1 * y1 <= input.r_grain_sq)
      {
        xtemp = x1 * input.costheta - y1 * input.sintheta;
        ytemp = x1 * input.sintheta + y1 * input.costheta;
        x1 = xtemp;
        y1 = ytemp;
      }
    }

    x1 += input.box.Lx / 2.0;
    y1 += input.box.Ly / 2.0;
    z1 += input.box.Lz / 2.0;

    Position p(x1, y1, z1);
    atoms[i].setWrapped(p);
  }
}

vector <vector <int> > generateCellLinkedList(const inputData& input, const vector <Atom>& atoms)
{
  int ncellx, ncelly, ncellz; // number of cells in each direction
  int idx, idy, idz; // cell number in each direction
  double lcellx, lcelly, lcellz; // length of cells in each direction
  int n_atoms_per_cell; // number of atoms allowed per cell
  double drij_sq, rxij, ryij, rzij; // square of distance, x, y, and z separation.
  vector <vector <int> > iatom; // the neighbor list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell

  // First we generate the number of cells in each direction - dependent on the lattice parameter and crystal structure
  ncellx = (int)(input.box.Lx / (input.a0 * first_nn_distances[input.crystal_structure])) + 1;
  ncelly = (int)(input.box.Ly / (input.a0 * first_nn_distances[input.crystal_structure])) + 1;
  ncellz = (int)(input.box.Lz / (input.a0 * first_nn_distances[input.crystal_structure])) + 1;

  // Length of cells in each direction
  lcellx = input.box.Lx / ncellx;
  lcelly = input.box.Ly / ncelly;
  lcellz = input.box.Lz / ncellz;

  // Minimum number of atoms allowed of 100
  n_atoms_per_cell = max((int)(atoms.size() / (double)(3 * ncellx * ncelly * ncellz)), 100);

  // resize the vectors
  icell.resize(ncellx, vector <vector <int> > // x dimension
              (ncelly, vector <int> // y dimension
              (ncellz, 0))); // z dimension
  pcell.resize(ncellx, vector <vector <vector <int> > > // x dimension
              (ncelly, vector <vector <int> > // y dimension
              (ncellz, vector <int> // z dimension
              (n_atoms_per_cell, 0)))); // atom number in cell.
  iatom.resize(n_atoms_per_cell, vector <int> (atoms.size(),0));

  // generate the pcell and icell matrices.
  for (unsigned int i = 0; i < atoms.size(); ++i) // Look at each atom
  {
    // Note that this function differs from that in find_grains because we want
    // to keep track of _every_ atom type, not just the type 1.
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

  for (int i = 0; i < ncellx; ++i) // For each x cell
  {
    for (int j = 0; j < ncelly; ++j) // For each y cell
    {
      for (int k = 0; k < ncellz; ++k) // For each z cell
      {
        for (int l = 0; l < icell[i][j][k]; ++l) // For each atom in this cell
        {
          int id = pcell[i][j][k][l]; // store this atom id
          // Now we check each sub cell around the current one
          for (int ii = -1; ii < 2; ++ii) // allowed values: -1, 0, and 1
          {
            for (int jj = -1; jj < 2; ++jj)
            {
              for (int kk = -1; kk < 2; ++kk)
              {
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
                for (int m = 0; m < icell[ia][ja][ka]; ++m)
                {
                  int jd = pcell[ia][ja][ka][m];
                  // If jd <= id, we've already dealt with this interaction
                  if (jd <= id)
                  {
                    continue;
                  }

                  // Now the actual calculations!
                  rxij = atoms[id].getWrapped()[0] - atoms[jd].getWrapped()[0];
                  ryij = atoms[id].getWrapped()[1] - atoms[jd].getWrapped()[1];
                  rzij = atoms[id].getWrapped()[2] - atoms[jd].getWrapped()[2];

                  // Apply PBCs
                  rxij = rxij - anInt(rxij / input.box.Lx) * input.box.Lx;
                  ryij = ryij - anInt(ryij / input.box.Ly) * input.box.Ly;
                  rzij = rzij - anInt(rzij / input.box.Lz) * input.box.Lz;

                  // Now calculate the distance
                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                  if (drij_sq > input.r_cut_max_sq)
                  {
                    continue; // move to the next atom if we're too far away
                  }

                  // Create the neighbor list
                  iatom[0][id] += 1; //for atom id - number of neighbors
                  iatom[(iatom[0][id])][id] = jd; // point to the next atom
                  iatom[0][jd] += 1; // for atom jd
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

int removeAtoms(vector <Atom>& atoms, const vector <vector <int> >& iatom, const inputData& input)
{
  double rxij, ryij, rzij, drij_sq, x1, y1, z1;
  vector <neighborData> neighbor_dataset;
  vector <int> n_removed(input.n_types, 0);

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    int i_type = atoms[i].getType();
    if (type_one_only && i_type != 1) {continue;} // only focusing on the type 1 atoms
    if (atoms[i].getMark() != 0) {continue;} // ignoring atoms we have already marked

    // Store the x, y, and z positions of the atom for easy comparison
    x1 = atoms[i].getWrapped()[0];
    y1 = atoms[i].getWrapped()[1];
    z1 = atoms[i].getWrapped()[2];

    // Now check it's nearest neighbors
    neighbor_dataset.clear();
    bool removed = false;
    for (int l = 1; l <= iatom[0][i]; ++l)
    {
      // Note that we check the distances between the original atom and ALL of
      // its neighbors to facilitate maintaining the same elemental ratios
      int j = iatom[l][i];
      int j_type = atoms[j].getType();

      if (atoms[j].getMark() == 0) // If the atom is not marked, check it
      {

        // calculate the distance between the two
        rxij = x1 - atoms[j].getWrapped()[0];
        ryij = y1 - atoms[j].getWrapped()[1];
        rzij = z1 - atoms[j].getWrapped()[2];

        // Apply PBCs
        rxij = rxij - anInt(rxij / input.box.Lx) * input.box.Lx;
        ryij = ryij - anInt(ryij / input.box.Ly) * input.box.Ly;
        rzij = rzij - anInt(rzij / input.box.Lz) * input.box.Lz;

        drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
        if (input.n_types != 1)
        { // This is only needed if there is more than 1 type of atom.
          neighborData tmp;
          tmp.ids = make_pair(atoms[i].getId(),atoms[j].getId());
          tmp.types = make_pair(i_type,j_type);
          tmp.distance = drij_sq;
          neighbor_dataset.push_back(tmp);
        }

        if (j_type == i_type && !removed) // We are only interested in the primary atom interactions (1-1)
        {
          pair <int, int> key =  (i_type < j_type) ? make_pair(i_type,j_type) : make_pair(j_type,i_type); // make sure the lower index comes first
          if (drij_sq < input.rcut_sq.at(key))
          {
            atoms[i].setMark(1); // if we are below the specified cutoff, mark for removal
            ++n_removed[i_type - 1]; // count the number removed of this type
            removed = true;
          }
        }
      }
    }
    if (removed && input.n_types != 1)
    {
      sort(neighbor_dataset.begin(), neighbor_dataset.end(), compareNeighbors);
      // Now we need to remove atoms to maintain the original ratio
      for (int type = 1; type <= input.element_ratios.size(); ++type)
      { // for each ratio
        int ratio = input.element_ratios.at(input.number_to_element.at(type));
        if (type == i_type) {ratio -= 1;} // if this is the atom we already removed, we don't want to remove an extra atom

        for (int ii = 0; ii < ratio; ++ii) // remove atoms until this ratio is met
        {
          for (vector <neighborData>::iterator neigh_it = neighbor_dataset.begin();
               neigh_it != neighbor_dataset.end(); ++neigh_it)
          { // check each neighbor's type.  If it matches the type we are looking for, remove it, and break from the loop
            if ((*neigh_it).types.second == type && atoms[(*neigh_it).ids.second - 1].getMark() == 0)
            {
              atoms[(*neigh_it).ids.second - 1].setMark(1);
              ++n_removed[type - 1];
              break;
            }
          }
        }
      }
    }
  }

  for (unsigned int i = 0; i < n_removed.size(); ++i)
  {
    cout << "  " << n_removed[i] << " " << input.number_to_element.at(i + 1) << " atoms will be removed.\n";

    for (unsigned int j = 0; j < input.element_ratios.size(); ++j)
    {
      if (input.element_ratios.at(input.number_to_element.at(j + 1)) * n_removed[i] != input.element_ratios.at(input.number_to_element.at(i + 1)) * n_removed[j])
      {
        cout << "Error: the element ratio has not been kept.\n"
             << input.element_ratios.at(input.number_to_element.at(j + 1)) << "*" << n_removed[i] << " != "
             << input.element_ratios.at(input.number_to_element.at(i + 1)) << "*" << n_removed[j] << endl;
        exit(ATOM_COUNT_ERROR);
      }
    }
  }

  return accumulate(n_removed.begin(), n_removed.end(), 0);
}

void writeMarkedFile(const string& filename, const vector <Atom>& atoms, const inputData& input)
{
  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);
  fout << fixed;

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    fout << atoms[i].getId() << " " << atoms[i].getType() << " ";
    if (input.has_charge) {fout << atoms[i].getCharge() << " ";}
    atoms[i].getWrapped().print3DSpace(fout);
    fout << " " << atoms[i].getMark() << endl;
  }
  fout.close();
}

void writeRotatedFile(const string& filename, const vector <Atom>& atoms, const inputData& input)
{
  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);
  fout << fixed;

  fout << "These " << input.molecule << " coordinates are rotated [ID type ";
  if (input.has_charge) {fout << "charge ";}
  fout << "x y z]\n\n"
       << atoms.size() << " atoms\n"
       << input.n_types << " atom types\n";
  fout.precision(6);
  fout << input.box.xlow << " " << input.box.xhigh << " xlo xhi\n"
       << input.box.ylow << " " << input.box.yhigh << " ylo yhi\n"
       << input.box.zlow << " " << input.box.zhigh << " zlo zhi\n";
  if (input.box.is_triclinic)
  {
    fout << input.box.xy << " " << input.box.xz << " " << input.box.yz << "xy xz yz\n";
  }
  fout << "\nAtoms\n\n";

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    fout.precision(0);
    fout << atoms[i].getId() << " " << atoms[i].getType() << " ";
    if (input.has_charge)
    {
      fout.precision(1);
      fout << atoms[i].getCharge() << " ";
    }
    fout.precision(6);
    atoms[i].getWrapped().print3DSpace(fout);
    fout << endl;
  }
  fout.close();
}

void writeRemovedFile(const string& filename, const vector <Atom>& atoms, const int& n_removed, const inputData& input)
{
  int n_total = 0;

  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);
  fout << fixed;

  fout << "These " << input.molecule << " coordinates are shifted and have atoms removed: [ID type ";
  if (input.has_charge) {fout << "charge ";}
  fout << "x y z]\n"
       << "\n"
       << atoms.size() - n_removed << "  atoms\n"
       << input.n_types << "  atom types\n"
       << input.box.xlow << " " << input.box.xhigh << "   xlo xhi\n"
       << input.box.ylow << " " << input.box.yhigh << "   ylo yhi\n"
       << input.box.zlow << " " << input.box.zhigh << "   zlo zhi\n";
  if (input.box.is_triclinic)
  {
    fout << input.box.xy << " " << input.box.xz << " " << input.box.yz << "  xy xz yz\n";
  }

  fout << "\nAtoms\n\n";

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getMark() == 0)
    {
      ++n_total;
      fout << n_total << " " << atoms[i].getType() << " ";
      if (input.has_charge) {fout << atoms[i].getCharge() << " ";}
      atoms[i].getWrapped().print3DSpace(fout);
    }
  }

  if (n_total != atoms.size() - n_removed)
  {
    cout << "Error: The final number of removed atoms is not balanced.\n"
         << "n_total = " << n_total << " != N - n_removed = "
         << atoms.size() - n_removed << endl;
    exit(ATOM_COUNT_ERROR);
  }

  fout.close();
}

void rotateAndRemove(inputData& input, const unsigned int& index)
{
  string marked_filename, rotated_filename, removed_filename;
  vector <Atom> atoms;
  vector <vector <int> > iatom;
  int n_removed;

  stringstream ss;
  ss << "_" << index << ".dat";
  marked_filename = input.output_basename + "_marked" + ss.str();
  rotated_filename = input.output_basename + "_rotated" + ss.str();
  removed_filename = input.output_basename + "_removed" + ss.str();

  atoms = readDataFile(input);
  rotateAtoms(input, atoms);
  iatom = generateCellLinkedList(input, atoms);
  n_removed = removeAtoms(atoms, iatom, input);

  if (marked) {writeMarkedFile(marked_filename, atoms, input);}
  if (rotated) {writeRotatedFile(rotated_filename, atoms, input);}
  if (removed) {writeRemovedFile(removed_filename, atoms, n_removed, input);}

  cout << endl;
}

void printInputData(const inputData& input)
{
  cout << "Datafile: " << input.data_file
       << "\n  Output basename: " << input.output_basename
       << "\n  Grain radius: " << input.r_grain
       << "\n  Misorientation angle: " << input.theta
       << "\n  Molecule: " << input.molecule
       << "\n  Crystal structure: " << input.crystal_structure
       << "\n  Lattice parameter: " << input.a0
       << "\n\n  Elements found (" << input.n_types << "): ";
  for (map <int, string>::const_iterator it = input.number_to_element.begin(); it != input.number_to_element.end();)
  {
    cout << "    " << it->second;
    if (++it != input.number_to_element.end()) {cout << ", ";}
    else {cout << endl;}
  }

  if (input.n_types > 1)
  {
    cout << "  The ratio of elements was calculated as: ";
    for (map <int, string>::const_iterator it = input.number_to_element.begin(); it != input.number_to_element.end();)
    {
      cout << input.element_ratios.at(it->second) << " " << it->second;
      if (++it != input.number_to_element.end()) {cout << " : ";}
      else {cout << endl;}
    }
  }
}

void writeLogFile(const inputData& input, const unsigned int& index)
{
  ofstream fout("log.rotate_and_remove", fstream::app);
  checkFileStream(fout, "log.rotate_and_remove");

  fout << "File: " << input.data_file << endl
       << "Theta: " << input.theta << endl
       << "Cutoff value(s): \n";
  for (map <pair <int, int>, double>::const_iterator it = input.rcut.begin(); it != input.rcut.end(); ++it)
  {
    fout << "\t" << (*it).first.first << "-" << (*it).first.second << ": "
         << (*it).second << endl;
  }
  stringstream ss;
  ss << "_" << index << ".dat";
  string marked_filename = input.output_basename + "_marked" + ss.str();
  string rotated_filename = input.output_basename + "_rotated" + ss.str();
  string removed_filename = input.output_basename + "_removed" + ss.str();

  fout << "Outputs: \n";
  if (marked) {fout << "\tMarked: " << marked_filename << "\n";}
  if (rotated) {fout << "\tRotated: " << rotated_filename << "\n";}
  if (removed) {fout << "\tRemoved: " << removed_filename << "\n\n";}

  fout.close();
}

void runCommands(vector <inputData>& commands)
{
  map <int, string> elements;
  for (unsigned int i = 0; i < commands.size(); ++i)
  {
    printInputData(commands[i]);
    rotateAndRemove(commands[i], i);
    writeLogFile(commands[i], i);
  }
}

vector <inputData> parseInputFile(const string& input_file)
{
  string str; // holds the whole line
  vector <inputData> commands;
  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);
  bool good_stream = true;

  while (getline(fin,str))
  {
    bool single_cutoff = false;
    inputData input;
    stringstream ss(str); // for parsing the command string
    if (ss.str()[0] == '#') {continue;}
    if (!(ss >> input.data_file >> input.output_basename >> input.r_grain >> input.theta >> input.molecule >> input.crystal_structure >> input.a0))
    {
      cout << "Incorrect line formatting.\n" << "Skipping... \"" << str << "\"\n\n";
      continue;
    }
    input.calculateRGrainSq();
    input.calculateTrig();
    input.calculateNumElements();
    determineElementRatios(input);
    int combinations = ((input.n_types + 1) * input.n_types) / 2;
    double first_cutoff;

    for (int i = 1; i <= input.n_types; ++i)
    {
      if (!good_stream) {break;}
      for (int j = i; j <= input.n_types; ++j)
      {
        double tmp;
        double multiplier;
        if (!(ss >> tmp))
        {
          if (i == 1) {first_cutoff = tmp;}
          if (i == 1 && j == 2) {single_cutoff = true;} // looking at the 1-2 interaction (the second value in the list of cutoffs)
          if (!single_cutoff)
          {
            cout << "Incorrect number of cutoff values - requires " << combinations << " cutoff values.\n\n";
            good_stream = false;
            break;
          }
        }
        if (input.crystal_structure.compare("fcc") == 0) {multiplier = first_nn_distances["fcc"];}
        else if (input.crystal_structure.compare("bcc") == 0) {multiplier = first_nn_distances["bcc"];}
        else if (input.crystal_structure.compare("sc") == 0) {multiplier = first_nn_distances["sc"];}
        else if (input.crystal_structure.compare("diamond") == 0) {multiplier = first_nn_distances["diamond"];}
        else if (input.crystal_structure.compare("fluorite") == 0)
        {
          if (i == j) // comparing the same element
          {
            if (input.element_ratios[input.number_to_element[i]] == 1) {multiplier = first_nn_distances["fcc"];} // fcc lattice points
            else if (input.element_ratios[input.number_to_element[i]] == 2) {multiplier = first_nn_distances["sc"] / 2.0;} // tetrahedral sites
            else
            {
              cout << "Error determining nearest neighbor multiplier.\n";
              exit(CRYSTAL_STRUCTURE_ERROR);
            }
          }
          else {multiplier = first_nn_distances["diamond"];}
        }

        if (single_cutoff)
        {
          double val = first_cutoff * multiplier * input.a0;
          input.rcut.insert(make_pair(make_pair(i,j), val));
          input.rcut_sq.insert(make_pair(make_pair(i,j), val * val));
        }
        else
        {
          double val = tmp * multiplier * input.a0;
          input.rcut.insert(make_pair(make_pair(i,j),val));
          input.rcut_sq.insert(make_pair(make_pair(i,j), val * val));
          if (val > input.r_cut_max) {input.r_cut_max = val;}
        }
      }
    }
    if (!good_stream)
    {
      cout << "Error processing line \"" << str << "\"\n\n";
      if (!printedInputHelp)
      {
        printInputFileHelp();
        printedInputHelp = true;
      }
      good_stream = true;
      continue;
    }
    input.calculateRCutMax();

    if (input.validate()) {commands.push_back(input);} // store the line
  }

  fin.close();
  return commands;
}

int main (int argc, char **argv)
{
  string input_file, outputs;
  map <int, string> elements;
  vector <inputData> commands; // the list of commands for the script
  try
  {
    cxxopts::Options options(argv[0], "Rotate atoms within a cylinder or sphere, and remove atoms that are too close to each other");
    options
      .positional_help("file")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
      ("f,file", "Input file", cxxopts::value<string>(input_file), "file")
      ("o,output", "Output files: (m)arked data files, with atoms marked for removal, r(o)tated data files, with the atoms rotated, and (r)emoved data files, with the marked atoms removed",
        cxxopts::value<string>(outputs)->default_value("r"), "mor")
      ("s,sphere", "Flag to create a spherical grain", cxxopts::value<bool>(is_sphere)->implicit_value("true"))
      ("1", "Flag to remove only type 1 atoms using the cutoff", cxxopts::value<bool>(type_one_only)->implicit_value("true"))
      ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !result.count("file"))
    {
      cout << options.help() << endl << endl;
      printInputFileHelp();
      return EXIT_SUCCESS;
    }

    if (result.count("output"))
    {
      marked = false;
      rotated = false;
      removed = false;

      if (outputs.size() > 3)
      {
        cout << "Error: specify the desired outputs as a combination of m|o|r\n";
        return OUTPUT_SPECIFICATION_ERROR;
      }
      else
      {
        for (string::iterator it = outputs.begin(); it != outputs.end(); ++it)
        {
          if ((*it) == 'm') {marked = true;}
          else if ((*it) == 'o') {rotated = true;}
          else if ((*it) == 'r') {removed = true;}
          else
          {
            cout << "Unknown output file format \'" << (*it) << "\'\n";
            return OUTPUT_SPECIFICATION_ERROR;
          }
        }
      }
    }

    if (result.count("file"))
    {
      commands = parseInputFile(input_file);
      runCommands(commands);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

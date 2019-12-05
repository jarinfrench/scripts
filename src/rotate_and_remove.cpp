/******************************************************************************
* This script requires as input the LAMMPS-formatted file for a single crystal
* at 0K, the grain radius, and the rotation angle (in degrees).  Prompts will be
* given if all three are not specified at the command line.  The atoms within
* the radius are rotated, and then the distances between atoms are checked.
* If the atoms are too close (as specified by the #define terms), one is removed,
* being sure to maintain charge neutrality (one U for every 2 O removed).  Note
* that this script is specifically for UO2.  Changes can be made for other systems.
******************************************************************************/

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"

using namespace std;

#define PI 3.14159265358979 // easier and faster to simply store these values here.
/*#define UU_RNN_CUT 2.0 // Cutoff value for U-U atoms too close (Basak Potential)
#define UO_RNN_CUT 4.0 // Cutoff value for U-O atoms too close
#define OO_RNN_CUT 0.63801 // Cutoff value for O-O atoms too close
#define CU_RNN_CUT 1.0 // Cutoff value for Cu-Cu atoms too close (using Mishin potential)
#define AL_RNN_CUT 2.4 // Cutoff value for Al-Al atoms too close (using Ercolessi-Adams potential)
#define AG_RNN_CUT 1.0 // Cutoff value for Ag-Ag atoms too close (using Mishin potential)
#define AU_RNN_CUT 1.0 // Cutoff value for Au-Au atoms too close (using Foiles Au1 potential)
#define NI_RNN_CUT 1.0 // Cutoff value for Ni-Ni atoms too close (using Foiles-Hoyt potential)
*/

// Boolean flags
bool is_sphere = false;
bool has_charge = false;
bool marked = true;
bool rotated = true;
bool removed = true;

// forward declarations
bool mapValueCmp(pair<pair<int,int>,double> a, pair<pair<int,int>,double> b);

struct ratio
{
  string chem_formula;
  int n_types;
  vector <int> ratio;
} compound_ratio;

struct neighbor_data
{
  pair<int, int> ids;
  pair<int, int> types;
  double distance;
};

struct inputVars
{
  string data_file;
  double r_grain, r_grain_sq;
  double theta, costheta, sintheta;
  double r_cut_max, r_cut_max_sq;
  int axis;
  map <pair <int, int>, double> rcut, rcut_sq;

  void calculateRCutSq() {r_cut_max_sq = r_cut_max * r_cut_max;}
  void calculateRGrainSq() {r_grain_sq = r_grain * r_grain;}
  void calculateRCutMax() {r_cut_max = (*max_element(rcut.begin(), rcut.end(), mapValueCmp)).second;}

} input;

struct boxData
{
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double xy, xz, yz;
  double Lx, Ly, Lz;
  bool is_triclinic = false;

  void calculateBoxLengths()
  {
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;
  }
} box;

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

bool mapValueCmp(pair<pair<int,int>,double> a, pair<pair<int,int>,double> b)
{
  return (a.second < b.second);
}

bool compareNeighbors(const neighbor_data &a, const neighbor_data &b)
{
  if (a.types.first  == b.types.first)
  {
    if (a.types.second == b.types.second)
    {
      return a.distance < b.distance;
    }
    else
    {
      return a.types.second < b.types.second;
    }
  }
  else
  {
    return a.types.first < b.types.first;
  }
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

void parseInputFile(const string& file)
{
  unsigned int combinations;

  ifstream finput(file.c_str());
  checkFileStream(finput, file);
  finput >> input.data_file >> input.r_grain >> compound_ratio.n_types;
  input.calculateRGrainSq();

  compound_ratio.ratio.resize(compound_ratio.n_types, 1);

  combinations = ((compound_ratio.n_types + 1) * compound_ratio.n_types) / 2;

  for (int i = 1; i <= compound_ratio.n_types; ++i)
  {
    for (int j = i; j <= compound_ratio.n_types; ++j)
    {
      double temp;
      finput >> temp;
      input.rcut.insert(make_pair(make_pair(i,j), temp));
      input.rcut_sq.insert(make_pair(make_pair(i,j),temp * temp));
      if (temp > input.r_cut_max) {input.r_cut_max = temp;}
    }
  }
  if (input.rcut.size() != combinations)
  {
    cout << "Error: There are " << combinations << " interaction radii needed.  You provided " << input.rcut.size() << ".\n";
    exit(INTERACTION_RADII_ERROR);
  }
  input.calculateRCutMax();
  input.calculateRCutSq();
  if (input.r_cut_max < 0.5)
  {
    cout << "Error: not enough memory available.\n";
    exit(MEMORY_EXCEEDED_ERROR);
  }

  finput.close();

  cout << "Input:"
       << "\n\tData file: " << input.data_file
       << "\n\tGrain radius: " << input.r_grain
       << "\n\tMisorientation angle: " << input.theta
       << "\n\tBoundary type: ";
  (is_sphere) ? cout << "sphere" : cout << "cylinder";
  cout << "\n\tNumber of atom types: " << compound_ratio.n_types
       << "\n\tCutoff radii: ";
  for (map <pair <int, int>, double>::iterator it = input.rcut.begin();
       it != input.rcut.end();)
  {
    cout << (*it).second;
    if (++it != input.rcut.end()) {cout << "; ";}
    else {cout << endl;}
  }
}

void determineAxisFromFilename()
{
  // This requires that (if not specifically specified on the command line) the
  // axis appears immediately before the file extension, meaning the filename is
  // <precursor_to_axis>_<axis>.<extension>
  istringstream iss(input.data_file.substr(input.data_file.find(".") - 3, 3));
  if (!(iss >> input.axis))
  {
    cout << "Error determining rotation axis.\n";
    exit(FILE_NAME_ERROR);
  }
}

map <int, string> determineElementRatios()
{
  size_t element_pos_start, element_pos_end, temp;
  map <int, string> elements;
  int elem_num = 1;

  // This next step requires LAMMPS to be in the filename, and immediately
  // followed by the elements, i.e. LAMMPS_<element(s)>_[<extra_info>]_axis.extension
  element_pos_start = input.data_file.find("LAMMPS_") + 7; // TODO: check that this is valid!  If adding 7 to string:npos works, what does it return?
  if (element_pos_start == string::npos)
  {
    cout << "Error determining element(s).\n";
    exit(FILE_NAME_ERROR);
  }

  temp = input.data_file.find("_N");
  temp = input.data_file.find("_N", temp + 1);
  if (temp == string::npos) {element_pos_end = input.data_file.find("_N");}
  else {element_pos_end = temp;}

  for (unsigned int i = element_pos_start; i < element_pos_end; ++i)
  {
    if (islower(input.data_file[i])) // if the current character is lower-case
    {
      if (isdigit(input.data_file[i + 1]))
      {
        compound_ratio.ratio[elem_num - 1] = input.data_file[i + 1] - '0';
      }
      int test = findByValue(elements, input.data_file.substr(i-1, 2));
      if (test < 0)
      {
        if ((elements.insert(pair<int,string>(elem_num,input.data_file.substr(i-1, 2)))).second) // check to see if insertion was successful
        {
          ++elem_num;
        }
        else
        {
          cout << "Error inserting element " << input.data_file.substr(i-1, 2);
        }
      }
      else
      {
        ++compound_ratio.ratio[test - 1];
      }
    }
    else
    {
      if (islower(input.data_file[i+1])) // if the next character is lower-case
      {
        continue;
      }
      else
      {
        if (isdigit(input.data_file[i]))
        {
          continue;
        }
        else
        {
          if (isdigit(input.data_file[i + 1]))
          {
            compound_ratio.ratio[elem_num - 1] = input.data_file[i + 1] - '0';
          }
          int test = findByValue(elements, input.data_file.substr(i, 1));
          if (test < 0)
          {
            if ((elements.insert(pair<int,string>(elem_num,input.data_file.substr(i, 1)))).second) // check to see if insertion was successful
            {
              ++elem_num;
            }
            else
            {
              cout << "Error inserting element " << input.data_file.substr(i, 1);
            }
          }
          else
          {
            ++compound_ratio.ratio[test - 1];
          }
        }
      }
    }
  }

  cout << endl;
  compound_ratio.chem_formula = input.data_file.substr(element_pos_start, element_pos_end - element_pos_start);
  if (compound_ratio.n_types != elements.size())
  {
    cout << "Error determining element types.  Input file n_types = "
         << compound_ratio.n_types << " != found elements = " << elements.size() << endl;
    exit(ELEMENT_COUNT_ERROR);
  }
  else
  {
    cout << "Elements found (" << compound_ratio.n_types << "): ";
    for (map <int, string>::iterator it = elements.begin(); it != elements.end();)
    {
      cout << (*it).second;
      if (++it != elements.end()) {cout << ", ";}
      else {cout << endl;}
    }
  }

  if (compound_ratio.n_types > 1)
  {
    cout << "The ratio of elements is calculated as: ";
    map <int, string>::iterator elem_it = elements.begin();
    for (vector <int>::iterator it = compound_ratio.ratio.begin();
         it != compound_ratio.ratio.end();)
    {
      cout << (*it) << " " << (*elem_it).second;
      ++elem_it;
      if (++it != compound_ratio.ratio.end()) {cout << " : ";}
      else {cout << endl;}
    }
  }

  return elements;
}

void generateCellLinkedList(const vector <Atom>& atoms, vector <vector <int> >& iatom)
{
  int ncellx, ncelly, ncellz; // number of cells in each direction
  int idx, idy, idz; // cell number in each direction
  double lcellx, lcelly, lcellz; // length of cells in each direction
  int n_atoms_per_cell; // number of atoms allowed per cell
  double drij_sq, rxij, ryij, rzij; // square of distance, x, y, and z separation.
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell

  // First we generate the number of cells in each direction
  ncellx = (int)(box.Lx / input.r_cut_max) + 1;
  ncelly = (int)(box.Ly / input.r_cut_max) + 1;
  ncellz = (int)(box.Lz / input.r_cut_max) + 1;

  // Length of cells in each direction
  lcellx = box.Lx / ncellx;
  lcelly = box.Ly / ncelly;
  lcellz = box.Lz / ncellz;

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

    // This is for unwrapped coordinates
    while (idx < 0.0) idx = ncellx + idx; // Note that this keeps things within the bounds set by lcellx
    while (idy < 0.0) idy = ncelly + idy; // Note that this keeps things within the bounds set by lcelly
    while (idz < 0.0) idz = ncellz + idz; // Note that this keeps things within the bounds set by lcellz

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
                  rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
                  ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
                  rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

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
}

vector <string> determineFilename()
{
  vector <string> name (3, "NONE");
  ostringstream fn1, fn2;
  fn1 << input.data_file.substr(0,input.data_file.find(".")).c_str()
      << "_" << input.theta << "degree_r" << input.r_grain;
  name[0] = fn1.str() + "A_marked.dat";
  name[1] = fn1.str() + "A_rotated.dat";

  fn2 << input.data_file.substr(0, input.data_file.find("N")).c_str()
      << input.axis << "_";
  fn2.precision(2);
  fn2 << fixed << input.theta << "degree_r";
  fn2.precision(0);
  fn2 << input.r_grain << "A_removed.dat";
  name[2] = fn2.str();

  return name;
}

void checkGrainSize(const double& val)
{
  if (input.r_grain * 2.0 > val)
  {
    cout << "Error: grain diameter = " << input.r_grain * 2.0
         << " >= boundary = " << val << endl;
    exit(GRAIN_TOO_LARGE_ERROR);
  }
}

vector <Atom> readDataFile()
{
  string str; // junk string variable
  int N, n_types, n_total = 0; // number of atoms, number of atom types
  double scale_factor_a, scale_factor_b, scale_factor_c;
  vector <Atom> atoms;
  int atom_id, type;
  double charge, x, y, z, x1, y1, z1, xtemp, ytemp;

  ifstream fin(input.data_file.c_str());
  checkFileStream(fin, input.data_file);

  getline(fin, str);
  transform(str.begin(), str.end(), str.begin(), ::tolower);
  if (str.find("charge") != string::npos)
  {
    has_charge = true;
  }
  fin >> N >> str;
  fin >> n_types >> str >> str;
  if (n_types != compound_ratio.n_types)
  {
    cout << "Atom types incorrectly determined.\n";
    exit(ATOM_TYPE_ERROR);
  }

  fin >> box.xlow >> box.xhigh >> str >> str
      >> box.ylow >> box.yhigh >> str >> str
      >> box.zlow >> box.zhigh >> str >> str;

  // Scaling factor
  scale_factor_a = 1.0;
  scale_factor_b = 1.0;
  scale_factor_c = 1.0;

  box.calculateBoxLengths();

  checkGrainSize(box.Lx);
  checkGrainSize(box.Ly);
  if (is_sphere) {checkGrainSize(box.Lz);}

  // scale the dimensions
  box.xlow *= scale_factor_a;
  box.xhigh *= scale_factor_a;
  box.ylow *= scale_factor_b;
  box.yhigh *= scale_factor_b;
  box.zlow *= scale_factor_c;
  box.zhigh *= scale_factor_c;
  box.calculateBoxLengths();

  fin.ignore();
  getline(fin,str); // read the extra stuff
  if (str.find("xy") != string::npos) // Triclinic system check
  {
    stringstream tmp(str);
    tmp >> box.xy >> box.xz >> box.yz >> str >> str >> str;
    fin >> str;
    box.is_triclinic = true;
    fin.ignore();
  }
  else
  {
    getline(fin, str); // required in order to line up the different input file structures
  }

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

    if (type > compound_ratio.n_types)
    {
      cout << "Error: Atom type = " << type << " is greater than the number of types specified = " << compound_ratio.n_types << endl;
      exit(ATOM_TYPE_ERROR);
    }

    // change the origin to the center of the simulation for rotating the atoms
    x1 = x - box.Lx / 2.0;
    y1 = y - box.Ly / 2.0;
    z1 = z - box.Lz / 2.0;

    // If we are smaller than the radius, rotate the atom by theta around the
    // z axis.
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

    x1 += box.Lx / 2.0;
    y1 += box.Ly / 2.0;
    z1 += box.Lz / 2.0;

    x1 *= scale_factor_a;
    y1 *= scale_factor_b;
    z1 *= scale_factor_c;

    Position p(x1,y1,z1);
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

int removeAtoms(vector <Atom>& atoms, vector <vector <int> >& iatom,
                 map <int, string>& elements)
{
  double rxij, ryij, rzij, drij_sq, x1, y1, z1;
  vector <neighbor_data> neighbor_dataset;
  vector <int> n_removed(compound_ratio.n_types, 0);

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    int j = atoms[i].getType();
    if (atoms[i].getMark() == 0)
    {
      // Store the x, y, and z positions of the atom for easy comparison
      x1 = atoms[i].getWrapped()[0];
      y1 = atoms[i].getWrapped()[1];
      z1 = atoms[i].getWrapped()[2];

      // Now check it's nearest neighbors
      neighbor_dataset.clear();
      bool removed = false;
      for (int l = 1; l <= iatom[0][i]; ++l)
      {
        int id = iatom[l][i];
        if (atoms[id].getMark() == 0) // If the atom is not marked, check it
        {
          int k = atoms[id].getType();
          // calculate the distance between the two
          rxij = x1 - atoms[id].getWrapped()[0];
          ryij = y1 - atoms[id].getWrapped()[1];
          rzij = z1 - atoms[id].getWrapped()[2];

          // Apply PBCs
          rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
          ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
          rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
          if (compound_ratio.n_types != 1)
          { // This is only needed if there is more than 1 type of atom.
            neighbor_data tmp;
            tmp.ids = make_pair(i,id);
            tmp.types = make_pair(j,k);
            tmp.distance = drij_sq;
            neighbor_dataset.push_back(tmp);
          }

          if (atoms[id].getType() == j && !removed) // We are really only interested in the atoms that are too close to the same type of atom.
          {
            pair <int, int> key =  (j < k) ? make_pair(j,k) : make_pair(k,j); // make sure the lower index comes first
            if (drij_sq < input.rcut_sq[key])
            {
              atoms[i].setMark(1); // if we are below the specified cutoff, mark for removal
              ++n_removed[j - 1]; // count the number removed of this type
              removed = true;
            }
          }
        }
      }
      if (removed && compound_ratio.n_types != 1)
      {
        sort(neighbor_dataset.begin(), neighbor_dataset.end(), compareNeighbors);
        // Now we need to remove atoms to maintain the original ratio
        for (int it = 1; it <= compound_ratio.ratio.size(); ++it)
        { // for each ratio
          int ratio = compound_ratio.ratio[it - 1];
          if (it == j) // if this is the atom we already removed, we don't want to remove an extra atom
          {
            ratio -= 1;
          }
          for (int ii = 0; ii < ratio; ++ii) // remove atoms until this ratio is met
          {
            for (vector <neighbor_data>::iterator neigh_it = neighbor_dataset.begin();
                 neigh_it != neighbor_dataset.end(); ++neigh_it)
            { // check each neighbor's type.  If it matches the type we are looking for, remove it, and break from the loop
              if ((*neigh_it).types.second == it && atoms[(*neigh_it).ids.second - 1].getMark() == 0)
              {
                atoms[(*neigh_it).ids.second - 1].setMark(1);
                ++n_removed[it - 1];
                break;
              }
            }
          }
        }
      }
    }
  }

  for (unsigned int i = 0; i < n_removed.size(); ++i)
  {
    for (unsigned int j = 0; j < compound_ratio.ratio.size(); ++j)
    {
      if (compound_ratio.ratio[j] * n_removed[i] != compound_ratio.ratio[i] * n_removed[j])
      {
        cout << "Error: the element ratio has not been kept.\n"
             << compound_ratio.ratio[j] << "*" << n_removed[i] << " != "
             << compound_ratio.ratio[i] << "*" << n_removed[j] << endl;
        exit(ATOM_COUNT_ERROR);
      }
    }
    cout << n_removed[i] << " " << elements[i+1] << " atoms will be removed.\n";
  }
  return accumulate(n_removed.begin(), n_removed.end(), 0);
}

void writeMarkedFile(const string& filename, const vector <Atom>& atoms)
{
  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);
  fout << fixed;

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    fout << atoms[i].getId() << " " << atoms[i].getType() << " ";
    if (has_charge) {fout << atoms[i].getCharge() << " ";}
    fout << atoms[i].getWrapped()[0] << " " << atoms[i].getWrapped()[1] << " " << atoms[i].getWrapped()[2]
         << " " << atoms[i].getMark() << endl;
  }
  fout.close();
}

void writeRotatedFile(const string& filename, const vector <Atom>& atoms)
{
  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);
  fout << fixed;

  fout << "These " << compound_ratio.chem_formula << " coordinates are rotated [ID type ";
  if (has_charge) {fout << "charge ";}
  fout << "x y z]\n\n"
       << atoms.size() << " atoms\n"
       << compound_ratio.n_types << " atom types\n";
  fout.precision(6);
  fout << box.xlow << " " << box.xhigh << " xlo xhi\n"
       << box.ylow << " " << box.yhigh << " ylo yhi\n"
       << box.zlow << " " << box.zhigh << " zlo zhi\n";
  if (box.is_triclinic)
  {
    fout << box.xy << " " << box.xz << " " << box.yz << " xz xy yz\n";
  }
  fout << "\nAtoms\n\n";

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    fout.precision(0);
    fout << atoms[i].getId() << " " << atoms[i].getType() << " ";
    if (has_charge)
    {
      fout.precision(1);
      fout << atoms[i].getCharge() << " ";
    }
    fout.precision(6);
    fout << atoms[i].getWrapped()[0] << " " << atoms[i].getWrapped()[1] << " " << atoms[i].getWrapped()[2] << endl;
  }
  fout.close();
}

void writeRemovedFile(const string& filename, const vector <Atom>& atoms, const int& n_removed)
{
  int n_total = 0;

  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);
  fout << fixed;

  fout << "These " << compound_ratio.chem_formula << " coordinates are shifted and have atoms removed: [ID type ";
  if (has_charge) {fout << "charge ";}
  fout << "x y z]\n"
       << "\n"
       << atoms.size() - n_removed << "  atoms\n"
       << compound_ratio.n_types << "  atom types\n"
       << box.xlow << " " << box.xhigh << "   xlo xhi\n"
       << box.ylow << " " << box.yhigh << "   ylo yhi\n"
       << box.zlow << " " << box.zhigh << "   zlo zhi\n";
  if (box.is_triclinic) {fout << box.xy << " " << box.xz << " " << box.yz << "  xy xz yz\n";}

  fout << "\nAtoms\n\n";

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getMark() == 0)
    {
      ++n_total;
      fout << n_total << " " << atoms[i].getType() << " ";
      if (has_charge) {fout << atoms[i].getCharge() << " ";}
      fout << atoms[i].getWrapped()[0] << " " << atoms[i].getWrapped()[1] << " "
           << atoms[i].getWrapped()[2] << endl;
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

void rotateAndRemove(map <int, string>& elements)
{
  string marked_filename, rotated_filename, removed_filename;
  vector <Atom> atoms;
  vector <vector <int> > iatom;
  int n_removed;

  // Temporary namespace for determining filenames
  {
    vector <string> temp = determineFilename();
    marked_filename = temp[0];
    rotated_filename = temp[1];
    removed_filename = temp[2];
  }

  atoms = readDataFile(); // rotates the atoms as well
  generateCellLinkedList(atoms, iatom);
  n_removed = removeAtoms(atoms, iatom, elements);
  if (marked) {writeMarkedFile(marked_filename, atoms);}
  if (rotated) {writeRotatedFile(rotated_filename, atoms);}
  if (removed) {writeRemovedFile(removed_filename, atoms, n_removed);}
}

void writeLogFile()
{
  ofstream fout("rotate_and_remove.log", fstream::app);
  checkFileStream(fout, "rotate_and_remove.log");

  fout << "File: " << input.data_file << endl
       << "Theta: " << input.theta << endl
       << "Cutoff value(s): \n";
  for (map <pair <int, int>, double>::iterator it = input.rcut.begin(); it != input.rcut.end(); ++it)
  {
    fout << "\t" << (*it).first.first << "-" << (*it).first.second << ": "
         << (*it).second << endl;
  }
  vector <string> temp = determineFilename();
  fout << "Outputs: \n";
  if (marked) {fout << "\tMarked: " << temp[0] << "\n";}
  if (rotated) {fout << "\tRotated: " << temp[1] << "\n";}
  if (removed) {fout << "\tRemoved: " << temp[2] << "\n\n";}

  fout.close();
}

int main(int argc, char **argv)
{
  string input_file, outputs;
  map <int, string> elements;
  try
  {
    cxxopts::Options options(argv[0], "Rotate atoms with a cylinder or sphere, and remove atoms that are too close to each other.");
    options
      .positional_help("file angle")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
      ("f,file", "Input file", cxxopts::value<string>(input_file), "file")
      ("t,theta", "Misorientation angle", cxxopts::value<double>(), "angle") // try to make this a vector
      ("o,output", "Output files: (m)arked data files, with atoms marked for removal, (r)otated data files, with the atoms rotated, and r(e)moved data files, with the marked atoms removed.  If the flag is specified alone, no output files will be produced",
        cxxopts::value<string>(outputs)->default_value("mre"), "mre")
      ("a,axis", "Orientation axis.  If not identified before the file extension, must be included", cxxopts::value<int>(input.axis)) // TODO: This can be better generalized
      ("s,sphere", "Flag to create a spherical grain", cxxopts::value<bool>(is_sphere)->default_value("false")->implicit_value("true"))
      ("h,help", "Show the help");

    options.parse_positional({"file", "theta"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("file") == 0 || result.count("theta") == 0)
    {
      cout << options.help() << endl << endl
           << "Please enter a valid input file containing (in order):\n"
           << "1) the data file to be processed\n"
           << "2) the desired grain radius\n"
           << "3) the number of atom types\n"
           << "4) cutoff radii\n\n"
           << "Note that the number of cutoff radii is dependent on the\n"
           << "number of atom interactions, so, for example, if there are\n"
           << "two atom types, there are the two same-element interactions,\n"
           << "and there is also the interaction between the two different\n"
           << "elements.  Include all the relevant cutoff radii one element\n"
           << "at a time, i.e. the cutoff radii for a 3 element system would\n"
           << "be input as 1-1, 1-2, 1-3, 2-2, 2-3, 3-3. Maintain consistent\n"
           << "numbering with the data file.\n";
      return EXIT_SUCCESS;
    }

    if (result.count("output"))
    {
      marked = false;
      rotated = false;
      removed = false;

      if (outputs.size() > 3)
      {
        cout << "Error: please specify which outputs you want printed to a file as m|r|e\n";
        return OUTPUT_SPECIFICATION_ERROR;
      }
      else
      {
        for (string::iterator it = outputs.begin(); it != outputs.end(); ++it)
        {
          if ((*it) == 'm') {marked = true;}
          else if ((*it) == 'r') {rotated = true;}
          else if ((*it) == 'e') {removed = true;}
          else
          {
            cout << "Unknown output file format \'" << (*it) << "\'\n";
            return OUTPUT_SPECIFICATION_ERROR;
          }
        }
      }
    }

    if (result.count("theta"))
    {
      input.theta = result["theta"].as<double>();
      input.costheta = cos(input.theta * PI / 180.0);
      input.sintheta = sin(input.theta * PI / 180.0);
    }

    if (result.count("file"))
    {
      parseInputFile(input_file);
      if (!result.count("axis"))
      {
        determineAxisFromFilename();
      }
      elements = determineElementRatios();
      rotateAndRemove(elements);
      writeLogFile();
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

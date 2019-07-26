#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"

//#include <chrono>
//#include <thread>

// using Clock = std::chrono::steady_clock;
// using std::chrono::time_point;
// using std::chrono::duration_cast;
// using std::chrono::milliseconds;
// using std::this_thread::sleep_for;
using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::map;
using std::pair;
using std::make_pair;
using std::ifstream;
using std::ofstream;
using std::fixed;
using std::stringstream;
using std::accumulate;

bool weight_percent = false;
bool atomic_percent = false;
bool remove_whole_molecule = true;
bool substitutions = true;
bool vacancies = true;
bool marked = true;
bool has_charge = false;

#define AVOGADRO 6.022E23

struct impurityDetails
{
  vector <int> atom_id; // The atom specifically chosen to be substituted
  int num_substituted = 0; // number of atoms to substitute
  int atom_type_to_substitute = 0; // the atom type we are allowed to substitute
};

struct boxData
{
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double Lx, Ly, Lz;

  void calculateBoxLengths()
  {
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;
  }
};

struct ratio
{
  string chem_formula = "";
  int n_types = 0;
  vector <int> ratio;
};

struct neighbor_data
{
  pair<int, int> ids;
  pair<int, int> types;
  double distance;
};

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

template <typename K, typename V>
int findByValue(map<K,V> mapOfElement, V value)
{
  typename map<K,V>::iterator it = mapOfElement.begin();
  while (it != mapOfElement.end())
  {
    if (it->second == value)
    {
      return it->first;
    }
    ++it;
  }
  return -1;
}

// Comparison function for comparing only the second values in a pair
bool pairCmp(pair<int, double> &a, pair<int, double> &b)
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

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

void showInputFileHelp()
{
  cout << "The input file may consist of multiple lines that may be formatted one of two ways:\n"
       << "\t(1) datafile molecule_name impurity_value atom_type seed r_cut1 [r_cut2 r_cut3 ....]\n"
       << "\t(2) datafile molecule_name atom_id r_cut1 [r_cut2 r_cut3 ....]\n"
       << "\t(3) datafile molecule_name file_containing_atom_ids r_cut1 [r_cut2 r_cut3 ....]\n\n"
       << "In all of these, the cutoff value must be specified for determining nearest neighbor\n"
       << "interactions. Note that the number of cutoff radii is dependent on the number\n"
       << "of unique elements in the structure.  For example, a system of two unique atoms\n"
       << "would require three cutoff radii: one for each of the same-element interactions,\n"
       << "and one for the interaction between the two elements.  All the radii must be\n"
       << "specified one element at a time, i.e. if there are three unique elements, the\n"
       << "interaction radii would be specified by 1-1, 1-2, 1-3, 2-2, 2-3, and finally 3-3.\n"
       << "Maintain consitent numbering with the data file.\n\n"
       << "Format (1) description:\n"
       << "   the impurity value can either be a decimal (0 < impurity_value < 1) or an\n"
       << "      integer (>1)\n"
       << "   atom_type specifies where substitutional atoms will be placed\n"
       << "   seed is used for random number generation\n"
       << "Format (2) and (3) description:\n"
       << "   the only parameter required is the atom ID number being replaced.\n"
       << "   or alternatively the file containing a list of atom IDs."
       << "   This is useful for situations where only defects are desired at\n"
       << "   specific locations.\n\n";
}

void writeVacancyFile(const string& vac_file, const vector <Atom>& atoms,
                      const vector <Atom>& substituted_atoms,
                      const boxData& box, const ratio& compound_ratio)
{
  int atom_id = 0;
  int num_atoms, num_marked = 0;
  int atoms_per_molecule = accumulate(compound_ratio.ratio.begin(), compound_ratio.ratio.end(), 0);

  if (remove_whole_molecule)
  {
    num_atoms = atoms.size() - substituted_atoms.size() * atoms_per_molecule;
  }
  else
  {
    num_atoms = atoms.size() - substituted_atoms.size();
  }

  ofstream fout(vac_file.c_str());
  checkFileStream(fout, vac_file);
  fout << "These " << compound_ratio.chem_formula << " coordinates have "
       << substituted_atoms.size() << " vacancies: [ID type charge x y z]\n\n"
       << num_atoms << "  atoms\n"
       << compound_ratio.n_types << "  atom types\n";

  fout.precision(6);
  fout << fixed;
  fout << box.xlow << " " << box.xhigh << " xlo xhi\n"
       << box.ylow << " " << box.yhigh << " ylo yhi\n"
       << box.zlow << " " << box.zhigh << " zlo zhi\n"
       << "\nAtoms\n\n";

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getMark() == 0)
    {
      fout << ++atom_id << " " << atoms[i].getType() << " ";
      if (has_charge) {fout << atoms[i].getCharge() << " ";}
      fout << atoms[i].getWrapped()[0] << " "
           << atoms[i].getWrapped()[1] << " " << atoms[i].getWrapped()[2] << endl;
    }
    else {++num_marked;}
  }

  fout.close();

  if (atom_id != num_atoms)
  {
    cout << "Error in vacancy file: n_written = " << atom_id << " != num_atoms = " << num_atoms << endl;
    exit(ATOM_COUNT_ERROR);
  }
}

void writeSubstitutionFile(const string& sub_file, const vector <Atom>& atoms,
                           const vector <Atom>& substituted_atoms,
                           const boxData& box, const ratio& compound_ratio)
{
  // NOTE: this does not necessarily keep the atom id the same between the original
  // file and the new file(s).
  int atom_id = 0;
  int num_atoms;

  ofstream fout(sub_file.c_str());
  checkFileStream(fout, sub_file);

  int atoms_per_molecule = accumulate(compound_ratio.ratio.begin(), compound_ratio.ratio.end(), 0);

  double top = substituted_atoms.size() * compound_ratio.ratio[substituted_atoms[0].getType() - 1];
  double bottom = (double) (atoms.size()) / (double)(atoms_per_molecule);
  double atomic_percent = (top / bottom) * 100.0;
  fout << "These " << compound_ratio.chem_formula << " coordinates have "
       << substituted_atoms.size() << " substitutions (" << atomic_percent
       << "at%): [ID type ";
  if (has_charge) {fout << "charge ";}
  fout << "x y z]\n\n";

  if (remove_whole_molecule)
  {
    num_atoms = atoms.size() - substituted_atoms.size() * atoms_per_molecule + substituted_atoms.size();
  }
  else
  {
    num_atoms = atoms.size();
  }

  fout << num_atoms << "  atoms\n"
       << compound_ratio.n_types + 1 << "  atom types\n";

  fout.precision(6);
  fout << fixed;
  fout << box.xlow << " " << box.xhigh << " xlo xhi\n"
       << box.ylow << " " << box.yhigh << " ylo yhi\n"
       << box.zlow << " " << box.zhigh << " zlo zhi\n"
       << "\nAtoms\n\n";

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getMark() == 0) // Unmarked atoms are rewritten
    {
      fout << ++atom_id << " " << atoms[i].getType() << " ";
      if (has_charge) {fout << atoms[i].getCharge() << " ";}
      fout << atoms[i].getWrapped()[0] << " " << atoms[i].getWrapped()[1] << " " 
           << atoms[i].getWrapped()[2] << endl;
    }
    else if (atoms[i].getMark() == 1) // substituted atoms
    {
      fout << ++atom_id << " " << compound_ratio.n_types + 1 << " " ;
      if (has_charge) {fout << 0.0 << " ";}
      fout << atoms[i].getWrapped()[0] << " " << atoms[i].getWrapped()[1] << " "
           << atoms[i].getWrapped()[2] << endl;
    }
    // If neither of these conditions are met, we are looking at the neighbors
    // of the molecule, which will not be written.
  }

  fout.close();

  if (atom_id != num_atoms)
  {
    cout << "Error in substitution file: n_written = " << atom_id << " != num_atoms = " << num_atoms << endl;
    exit(ATOM_COUNT_ERROR);
  }
}

void writeMarkedFile(const string& marked_file, const vector <Atom>& atoms)
{

  ofstream fout(marked_file.c_str());
  checkFileStream(fout, marked_file);

  fout.precision(6);
  fout << fixed;
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    fout << atoms[i].getId() << " " << atoms[i].getType() << " "
         << atoms[i].getCharge() << " " << atoms[i].getWrapped()[0] << " "
         << atoms[i].getWrapped()[1] << " " << atoms[i].getWrapped()[2] << " "
         << atoms[i].getMark() << endl;
  }

  fout.close();
}

void writeAtomsToFiles(const vector <Atom>& atoms, const string& datafile,
                       const vector <Atom> substituted_atoms, const boxData& box,
                       const ratio& compound_ratio, const string & outfile = "none")
{
  string vac_file, sub_file, marked_file;
  size_t ending_pos = datafile.find(".dat");

  stringstream ss;
  if (substituted_atoms.size() == 1)
  {
    ss << datafile.substr(0, ending_pos) << "_single_vac.dat";
  }
  else
  {
    ss << datafile.substr(0, ending_pos) << "_" << substituted_atoms.size() << "vac.dat";
  }
  ss >> vac_file;
  sub_file = vac_file.substr(0,vac_file.find("vac.dat")) + "sub.dat";
  marked_file = datafile.substr(0, ending_pos) + "_marked.dat";

  if (outfile.compare("none") != 0)
  {
    if (substitutions) {sub_file = outfile;}
    else if (vacancies) {vac_file = outfile;}
    else if (marked) {marked_file = outfile;}
    else
    {
      cout << "Error: no output specification found.\n";
    }
  }

  if (vacancies) {writeVacancyFile(vac_file, atoms, substituted_atoms, box, compound_ratio);}
  if (substitutions) {writeSubstitutionFile(sub_file, atoms, substituted_atoms, box, compound_ratio);}
  if (marked) {writeMarkedFile(marked_file, atoms);}
}

void removeWholeMolecule(vector <Atom>& substituted_atoms, vector <Atom>& atoms,
                         const vector <vector <int> >& iatom, const boxData& box,
                         map <int, string>& elements, const ratio& compound_ratio,
                         map <pair <int, int>, double>& r_cut_sq_map)
{
  double x, y, z, rxij, ryij, rzij, drij_sq;
  vector <pair <int, double> > distances;
  vector <neighbor_data> neighbor_dataset;
  vector <int> n_removed (compound_ratio.n_types, 0);

  // For each atom that is being substituted, find its neighbors within the molecule,
  // and remove them as well.
  for (unsigned int i = 0; i < substituted_atoms.size(); ++i) // Remove the whole molecule
  {
    int j = substituted_atoms[i].getType(); // get the type of the substituted atom
    int n = substituted_atoms[i].getId() - 1; // store the atom id (minus one due to 0-based indexing)

    ++n_removed[j - 1]; // add the substituted atom to the removed list.

    // Store the values for easy comparison
    x = atoms[n].getWrapped()[0];
    y = atoms[n].getWrapped()[1];
    z = atoms[n].getWrapped()[2];

    neighbor_dataset.clear(); // reset the nearest neighbor dataset

    // for each neighbor of atom i
    for (int l = 1; l <= iatom[0][n]; ++l)
    {
      int id = iatom[l][n];
      if (atoms[id].getMark() == 0) // if the neighbor has not yet been removed
      {
        int k = atoms[id].getType(); // store this atom's type

        // calculate the distances
        rxij = x - atoms[id].getWrapped()[0];
        ryij = y - atoms[id].getWrapped()[1];
        rzij = z - atoms[id].getWrapped()[2];

        //Apply PBCs
        rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
        ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
        rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

        drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
        if (compound_ratio.n_types != 1)
        { // we're only worried about this case if there are > 1 n_types
          neighbor_data tmp;
          tmp.ids = make_pair(n, id);
          tmp.types = make_pair(j,k);
          tmp.distance = drij_sq;
          neighbor_dataset.push_back(tmp);
        }
      }
    }

    if (compound_ratio.n_types != 1)
    {
      sort(neighbor_dataset.begin(), neighbor_dataset.end(), compareNeighbors);

      // Maintain the original ratio
      for (int it = 1; it <= compound_ratio.ratio.size(); ++it)
      {
        int elem_ratio = compound_ratio.ratio[it - 1];
        if (it == j) {elem_ratio -= 1;} // we've already handled this case, don't remove another one!

        for (int ii = 0; ii < elem_ratio; ++ii)
        {
          for (vector <neighbor_data>::iterator neigh_it = neighbor_dataset.begin();
               neigh_it != neighbor_dataset.end(); ++neigh_it)
          {
            if ((*neigh_it).types.second == it && atoms[(*neigh_it).ids.second].getMark() == 0)
            {
              atoms[(*neigh_it).ids.second].setMark(2);
              ++n_removed[it - 1];
              break;
            }
          } // neighbor_dataset iterator
        } // ii
      } // it
    } // number of types > 1
  } // i

  for (unsigned int i = 0; i < n_removed.size(); ++i)
  {
    for (unsigned int j = 0; j < compound_ratio.ratio.size(); ++j)
    {
      if (compound_ratio.ratio[j] * n_removed[i] != compound_ratio.ratio[i] * n_removed[j] && !remove_whole_molecule)
      {
        cout << "Error: the element ratio has not been kept.\n"
             << compound_ratio.ratio[j] << "*" << n_removed[i] << " != "
             << compound_ratio.ratio[i] << "*" << n_removed[j] << endl;
        exit(ATOM_COUNT_ERROR);
      }
    }
    cout << n_removed[i] << " " << elements[i + 1] << " atoms will be removed.\n";
  }
  cout << "A total of " << accumulate(n_removed.begin(), n_removed.end(), 0) << " atoms will be removed.\n";
}

void determineAllowedSubstitionalAtoms(vector <Atom>& substitional_atoms, const vector <Atom>& atoms, const int& allowed_type)
{
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getType() == allowed_type)
    {
      substitional_atoms.push_back(atoms[i]);
    }
  }
}

impurityDetails
calculateImpurityDetails(const vector <Atom>& atoms, const vector <string>& commands,
                         const int& num_cutoffs, const ratio& compound_ratio)
{
  impurityDetails impurity_details;
  double impurity_value;
  int seed;

  if (commands.size() == 5 + num_cutoffs) // we either have a specific number to replace (randomly) or a percentage
  {
    stringstream ss(commands[2]); // This is the impurity value for input type (1)
    ss >> impurity_value;

    ss.clear();
    ss.str(commands[3]);
    ss >> impurity_details.atom_type_to_substitute;

    ss.clear();
    ss.str(commands[4]);
    ss >> seed;

    cout << "The seed for random number generation is " << seed << endl;
    srand(seed);

    if (impurity_value > 1) {impurity_details.num_substituted = (int)impurity_value;}
    else if (impurity_value > 0)
    {
      int numerator = compound_ratio.ratio[impurity_details.atom_type_to_substitute - 1] * atoms.size();
      int denominator = std::accumulate(compound_ratio.ratio.begin(), compound_ratio.ratio.end(), 0);

      impurity_details.num_substituted = (int)((double)(numerator) / (double)(denominator) * impurity_value);
    }
  }
  else if (commands.size() == 3 + num_cutoffs) // either one id, or a file containing a list of IDs has beeen specified
  {
    // Check to see if the value specified in commands[2] is a file or a number
    bool using_file = std::find_if(commands[2].begin(), commands[2].end(),
                      [](char c) { return !std::isdigit(c); }) != commands[2].end();

    if (using_file)
    {
      string line;
      int tmp; // temporary variable to hold the atom id
      int atom_type;
      int n = 0; // line number;
      bool same_atom_type = true;
      ifstream fin(commands[2].c_str());
      checkFileStream(fin, commands[2]);
      while (getline(fin, line))
      {
        stringstream line_data(line);
        line_data >> tmp;
        impurity_details.atom_id.push_back(tmp);
        if (n == 0) {atom_type = atoms[impurity_details.atom_id[n] - 1].getType();}
        if (atoms[impurity_details.atom_id[n] - 1].getType() != atom_type) {same_atom_type = false;}
        ++impurity_details.num_substituted;
      }
      if (same_atom_type)
      {
        impurity_value = (impurity_details.num_substituted * std::accumulate(compound_ratio.ratio.begin(), compound_ratio.ratio.end(), 0)) / \
        (atoms.size() * compound_ratio.ratio[atoms[impurity_details.atom_id[0] - 1].getType() - 1]);
      }
      else
      {
        impurity_value = impurity_details.num_substituted / atoms.size();
      }
    }
    else
    {
      int tmp;
      stringstream ss(commands[2]);
      ss >> tmp;
      impurity_details.atom_id.push_back(tmp);
      impurity_value = (1.0 * std::accumulate(compound_ratio.ratio.begin(), compound_ratio.ratio.end(), 0)) / \
      (atoms.size() * compound_ratio.ratio[atoms[impurity_details.atom_id[0] - 1].getType() - 1]);
      ++impurity_details.num_substituted;
    }
  }

  cout << "There will be " << impurity_details.num_substituted
       << " substitutions ("
       << impurity_value << "at%) made." << endl; // FIXME

  return impurity_details;
}

vector <Atom> generateImpurities(vector <Atom>& atoms, const vector <string>& commands, const int& num_cutoffs, const ratio& compound_ratio)
{
  int n_removed = 0;
  impurityDetails impurity_details;
  vector <Atom> substitutional_atoms;

  impurity_details = calculateImpurityDetails(atoms, commands, num_cutoffs, compound_ratio);


  if (impurity_details.atom_id.size() > 0) // if an atom id has been specified
  {
    substitutional_atoms.resize(impurity_details.num_substituted, Atom());
    for (unsigned int i = 0; i < impurity_details.atom_id.size(); ++i)
    {
      atoms[impurity_details.atom_id[i] - 1].setMark(1);
      substitutional_atoms[i] = atoms[impurity_details.atom_id[i] - 1];
      ++n_removed;
    }
  }
  else
  {
    determineAllowedSubstitionalAtoms(substitutional_atoms, atoms, impurity_details.atom_type_to_substitute);
    std::random_shuffle(substitutional_atoms.begin(), substitutional_atoms.end()); // randomize which atoms we keep.
    substitutional_atoms.resize(impurity_details.num_substituted); // Only keeping the number of substituted atoms
    for (unsigned int i = 0; i < substitutional_atoms.size(); ++i)
    {
      if (atoms[substitutional_atoms[i].getId() - 1].getType() != impurity_details.atom_type_to_substitute)
      {
        cout << "Error in removing atoms.\n";
        exit(ATOM_TYPE_ERROR);
      }
      atoms[substitutional_atoms[i].getId() - 1].setMark(1);
      ++n_removed;
    }
  }

  if (n_removed != impurity_details.num_substituted)
  {
    cout << "Error substituting atoms: n_removed = " << n_removed << " != num_substituted = " << impurity_details.num_substituted << endl;
    exit(ATOM_COUNT_ERROR);
  }

  return substitutional_atoms;
}

void generateCellLinkedList(const vector <Atom>& atoms, vector <vector <int> >& iatom,
                            const boxData& box, const double& r_cut_max)
{
  int ncellx, ncelly, ncellz; // number of cells in each direction
  int idx, idy, idz; // cell number in each direction
  double lcellx, lcelly, lcellz; // length of cells in each direction
  int n_atoms_per_cell; // number of atoms allowed per cell
  double drij_sq, rxij, ryij, rzij; // square of distance, x, y, and z separation.
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell
  double r_cut_max_sq = r_cut_max * r_cut_max;

  // First we generate the number of cells in each direction
  ncellx = (int)(box.Lx / r_cut_max) + 1;
  ncelly = (int)(box.Ly / r_cut_max) + 1;
  ncellz = (int)(box.Lz / r_cut_max) + 1;

  // Length of cells in each direction
  lcellx = box.Lx / ncellx;
  lcelly = box.Ly / ncelly;
  lcellz = box.Lz / ncellz;

  // Minimum number of atoms allowed of 100
  n_atoms_per_cell = std::max((int)(atoms.size() / (double)(ncellx * ncelly * ncellz)), 200);

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
    // Assign this atom to a cell
    // Rounds towards 0 with a type cast
    idx = (int)(atoms[i].getWrapped()[0] / lcellx) * 1; // assign the x cell
    idy = (int)(atoms[i].getWrapped()[1] / lcelly) * 1; // assign the y cell
    idz = (int)(atoms[i].getWrapped()[2] / lcellz) * 1; // assign the z cell
    // Check if we went out of bounds
    // C++ indexes from 0, so we have to subtract 1 from the maximum value to
    // stay within our memory bounds
    if (idx >= ncellx) idx = ncellx - 1;
    if (idy >= ncelly) idy = ncelly - 1;
    if (idz >= ncellz) idz = ncellz - 1;

    // This is for unwrapped coordinates
    // while (idx < 0.0) idx = ncellx + idx; // Note that this keeps things within the bounds set by lcellx
    // while (idy < 0.0) idy = ncelly + idy; // Note that this keeps things within the bounds set by lcelly
    // while (idz < 0.0) idz = ncellz + idz; // Note that this keeps things within the bounds set by lcellz

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

                  if (drij_sq > r_cut_max_sq)
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

double extractCutoffs(const vector <string>& commands, const int& num_cutoffs,
                        const ratio& compound_ratio,
                        map <pair <int, int>, double>& r_cut_map,
                        map <pair <int, int>, double>& r_cut_sq_map)
{
  double r_cut_max = 0, temp;
  int i1 = 1, i2 = 1; // atom type indices for the cutoff maps

  stringstream ss;
  for (unsigned int i = commands.size() - num_cutoffs; i < commands.size(); ++i)
  {
    ss << commands[i] << " ";
  }

  for (int i = 1; i <= compound_ratio.n_types; ++i)
  {
    for (int j = i; j <= compound_ratio.n_types; ++j)
    {
      ss >> temp;
      r_cut_map.insert(make_pair(make_pair(i,j),temp));
      r_cut_sq_map.insert(make_pair(make_pair(i,j), temp * temp));

      if (temp > r_cut_max) {r_cut_max = temp;}
    }
  }

  return r_cut_max;
}

void readFile(vector <Atom>& atoms, const string& datafile, const ratio& compound_ratio, boxData& box)
{
  string str, str2;
  int N, ntotal = 0; // number of atoms, atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, xy, xz, yz; // box bounds, tilt factors
  double Lx, Ly, Lz; // Box lengths
  int atom_id, atom_type, n_types;
  double atom_charge, x, y, z;

  ifstream fin(datafile.c_str());
  checkFileStream(fin, datafile);

  getline(fin,str); // header line

  if (str.find("charge") != string::npos || str.find(" q ") != string::npos)
  {
    has_charge = true;
  }

  fin >> N >> str; // number of atoms

  atoms.resize(N, Atom());

  fin >> n_types >> str >> str; // get the number of atom types
  if (compound_ratio.n_types != n_types)
  {
    cout << "Error: Number of element types calculated does not match number of element types in the data file:\n"
         << "calculated value = " << compound_ratio.n_types << " != data file value = " << n_types << endl;
    exit(ATOM_TYPE_ERROR);
  }

  // Get the bounds of the system
  fin >> box.xlow >> box.xhigh >> str >> str;
  fin >> box.ylow >> box.yhigh >> str >> str;
  fin >> box.zlow >> box.zhigh >> str >> str;

  fin.ignore();
  getline(fin,str);
  if (str.find("xy") != string::npos)
  {
    fin >> xy >> xz >> yz >> str >> str >> str;
    fin.ignore();
  }

  box.calculateBoxLengths();

  fin >> str; // Gets the next line (Atoms)
  fin.ignore();
  getline(fin, str); // get the empty line

  // Actually read the atoms now
  while (getline(fin, str))
  {
    if (fin.fail())
    {
      cout << "Read error\n";
      break;
    }

    stringstream ss2(str);
    ss2 >> atom_id >> atom_type;
    if (has_charge) {ss2 >> atom_charge;}
    ss2 >> x >> y >> z;

    ++ntotal;
    if (atom_type > compound_ratio.n_types)
    {
      cout << "Error: Atom type = " << atom_type << " > n_types = " << compound_ratio.n_types << endl;
      exit(ATOM_TYPE_ERROR);
    }

    Position p(x,y,z);
    atoms[atom_id - 1] = Atom(atom_id, atom_type, atom_charge, p);
  }

  fin.close();

  if (ntotal != N)
  {
    cout << "Error reading atoms: ntotal = " << ntotal << " != N = " << N << endl;
    exit(ATOM_COUNT_ERROR);
  }
}

bool validateNumElements(const int& n_types, const vector <string>& commands, int& num_cutoffs)
{
  num_cutoffs = ((n_types + 1) * n_types) / 2;

  // The 5 and 3 come from the number of parameters in the different input formats.
  // combinations represents the number of cutoff radii required.
  if ((5 + num_cutoffs != commands.size()) && (3 + num_cutoffs != commands.size()))
  {
    cout << "Number of commands in line\n" << "  ";
    for (unsigned int i = 0; i < commands.size(); ++i)
    {
      cout << commands[i] << " ";
    }
    cout << endl
         << "not correct.\nNumber of commands found: " << commands.size()
         << " != expected number of commands = " << 5 + num_cutoffs
         << " (input format 1) or " << 3 + num_cutoffs << " (input format 2)\n\n";

    showInputFileHelp();
    return false;
  }
  return true;
}

map <int, string> determineElementRatios(const string& molecule_name, ratio& compound_ratio)
{
  // map the element number to a specific element.  Note that this will not
  // necessarily match the input file numbering.
  map <int, string> elements;
  int elem_num = 1;

  for (unsigned int i = 0; i < molecule_name.size(); ++i)
  {
    if (islower(molecule_name[i])) // if the current character is lower-case
    {
      if (isdigit(molecule_name[i + 1]))
      { //NOTE: This does not work if there are more than 9 of one particular element!
        compound_ratio.ratio[elem_num - 1] = molecule_name[i + 1] - '0'; // If the character is a digit, subtracting '0' gives us the actual value of the number, i.e. '9'-'0' = 9
      }
      int test = findByValue(elements, molecule_name.substr(i-1, 2)); // See if this element has already been inserted
      if (test < 0)
      {
        if ((elements.insert(pair<int,string>(elem_num,molecule_name.substr(i-1, 2)))).second) // check to see if insertion was successful
        {
          ++elem_num;
        }
        else
        {
          cout << "Error inserting element " << molecule_name.substr(i-1, 2);
        }
      }
      else
      {
        ++compound_ratio.ratio[test - 1];
      }
    }
    else // current character is upper case
    {
      if (islower(molecule_name[i+1])) // if the next character is lower-case
      {
        continue;
      }
      else
      {
        if (isdigit(molecule_name[i]))
        {
          continue;
        }
        else
        {
          // In this case, we know we will be adding another element to the list
          compound_ratio.ratio.push_back(1); // So we increase the size by one
          ++compound_ratio.n_types;

          if (isdigit(molecule_name[i + 1]))
          {
            compound_ratio.ratio[elem_num - 1] = molecule_name[i + 1] - '0';
          }
          int test = findByValue(elements, molecule_name.substr(i, 1));
          if (test < 0)
          {
            if ((elements.insert(pair<int,string>(elem_num,molecule_name.substr(i, 1)))).second) // check to see if insertion was successful
            {
              ++elem_num;
            }
            else
            {
              cout << "Error inserting element " << molecule_name.substr(i, 1);
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
  compound_ratio.chem_formula = molecule_name.substr(0, molecule_name.size());

  cout << "Elements found (" << compound_ratio.n_types << "): ";
  for (map <int, string>::iterator it = elements.begin(); it != elements.end();)
  {
    cout << (*it).second;
    if (++it != elements.end()) {cout << ", ";}
    else {cout << endl;}
  }

  if (compound_ratio.n_types > 1 || compound_ratio.ratio[0] > 1)
  {
    cout << "The ratio of elements (per molecule) is calculated as: ";
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

vector <vector <string> > parseInputFile(const string& input_file)
{
  string str; // variable to hold the line
  vector <vector <string> > impurity_commands;
  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);

  while (getline(fin, str))
  {
    string command; // string variable holding each individual element of the command line
    stringstream ss(str); // For parsing the number of commands
    vector <string> temp; // For holding the whole line as individual commands
    while (ss >> command)
    {
      temp.push_back(command);
    }

    impurity_commands.push_back(temp); // Store the line as individual commands
  }

  return impurity_commands;
}

int main(int argc, char **argv)
{
  string input_file, outfile;
  vector <vector <string> > impurity_commands;
  vector <vector <int> > iatom; // Cell-linked list
  vector <Atom> atoms;
  map <int, string> elements;
  boxData box;
  bool valid = true;

  try
  {
    cxxopts::Options options(argv[0], "Generate impurities in a given structure.");
    options
      .positional_help("file")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Input file", cxxopts::value<string>(input_file), "file")
        ("d,defect-type", "Specify which output files are wanted: (m)arked atoms for removal, atoms with (s)ubstitutions,  or (v)acancies.  Interstitials not implemented.", cxxopts::value<string>()->default_value("msv"), "m, s and/or v")
        ("o,outfile", "In the case of only one output file, specify the file name. (Assumes \'-d s\' unless otherwise specified).", cxxopts::value<string>(outfile), "Output file")
        ("no-molecule-removal", "Only remove the atom(s) specified for removal, without removing the whole molecule")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !result.count("file"))
    {
      cout << options.help() << endl << endl;
      showInputFileHelp();
    }

    if (result.count("no-molecule-removal")) {remove_whole_molecule = false;}

    if (result.count("defect-type"))
    {
      string type = result["defect-type"].as<string>();
      marked = false;
      substitutions = false;
      vacancies = false;
      if (type.size() > 3)
      {
        cout << "Error: please specify which output files you want as m|s|v\n";
        return OUTPUT_SPECIFICATION_ERROR;
      }
      else
      {
        for (string::iterator it = type.begin(); it != type.end(); ++it)
        {
          if ((*it) == 'm') {marked = true;}
          else if ((*it) == 's') {substitutions = true;}
          else if ((*it) == 'v') {vacancies = true;}
          else
          {
            cout << "Unknown output file specification \'" << (*it) << "\'.\n";
            return OUTPUT_SPECIFICATION_ERROR;
          }
        }
      }
    }

    if (result.count("outfile"))
    {
      if (!result.count("defect-type"))
      {
        marked = false;
        substitutions = true;
        vacancies = false;
      }
      else
      {
        if (marked + substitutions + vacancies > 1)
        {
          cout << "Output file name option only allowed with one output type!\n";
          return INPUT_FORMAT_ERROR;
        }
      }
    }

    if (result.count("file"))
    {
      // time_point<Clock> start = Clock::now();
      impurity_commands = parseInputFile(input_file);
      // time_point<Clock> end = Clock::now();
      // milliseconds diff = duration_cast<milliseconds>(end - start);
      // cout << diff.count() << "ms for parseInputFile\n";

      for (unsigned int i = 0; i < impurity_commands.size(); ++i) // for each command set
      {
        if (impurity_commands[i].size() < 1)
        {
          impurity_commands.erase(impurity_commands.begin() + i); //remove empty lines
          --i; // fix the counter
          continue;
        }

        int num_cutoffs;
        ratio compound_ratio;
        map <pair <int, int>, double> r_cut_map, r_cut_sq_map;
        vector <Atom> substituted_atoms;

        elements = determineElementRatios(impurity_commands[i][1], compound_ratio);

        if (!(validateNumElements(compound_ratio.n_types, impurity_commands[i], num_cutoffs)))
        {
          valid = false;
          continue;
        }

        double r_cut_max = extractCutoffs(impurity_commands[i], num_cutoffs,
                                          compound_ratio, r_cut_map, r_cut_sq_map);

        // If it's our first set of commands, or if the data file has changed between
        // the previous and current set of commands, read the data, and create the
        // cell-linked list.
        if (i == 0)
        {
          readFile(atoms, impurity_commands[i][0], compound_ratio, box);
          generateCellLinkedList(atoms, iatom, box, r_cut_max);
          valid = true;
        }
        else if (impurity_commands[i][0] != impurity_commands[i - 1][0] || !valid)
        {
          atoms.clear();
          iatom.clear(); // segfault with this vector if we don't do this.
          readFile(atoms, impurity_commands[i][0], compound_ratio, box);
          generateCellLinkedList(atoms, iatom, box, r_cut_max);
          valid = true;
        }
        else
        {
           // clear the atoms marked state.
          for (unsigned int i = 0; i < atoms.size(); ++i) {atoms[i].setMark(0);}
        }

        substituted_atoms = generateImpurities(atoms, impurity_commands[i], num_cutoffs, compound_ratio);

        if (!(result.count("no-molecule-removal")))
        {
          if (compound_ratio.n_types > 1 || (compound_ratio.n_types == 1 && compound_ratio.ratio[0] > 1))
          {
            cout << "WARNING: Removing whole molecule\n";
            removeWholeMolecule(substituted_atoms, atoms, iatom, box, elements,
                                compound_ratio, r_cut_sq_map);
          }
        }

        if (result.count("outfile"))
        {
          writeAtomsToFiles(atoms, impurity_commands[i][0], substituted_atoms,
                          box, compound_ratio, outfile);
        }
        else
        {
          writeAtomsToFiles(atoms, impurity_commands[i][0], substituted_atoms,
                          box, compound_ratio);
        }

      }

    }
  }
  catch (cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

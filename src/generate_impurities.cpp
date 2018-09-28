#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
//#include <random> // See https://msdn.microsoft.com/en-us/library/bb982398.aspx <-- requires C++11
#include <cxxopts.hpp>
#include "atom.h"
#include "error_code_defines.h"
#include <chrono>
#include <thread>

// using Clock = std::chrono::steady_clock;
// using std::chrono::time_point;
// using std::chrono::duration_cast;
// using std::chrono::milliseconds;
// using std::this_thread::sleep_for;
using namespace std;

bool weight_percent = false;
bool atomic_percent = false;
bool substitutions = true;
bool vacancies = true;
bool marked = true;

#define AVOGADRO 6.022E23
#define UO_RNN_CUT 4.0 // Cutoff values for U-O atoms too close

double anInt(double x);

//TODO: Make this general for any system, not just UO2
struct inputVars
{
  string datafile, outfile; // datafile for input, output file name
  bool using_atom_id = false, outfile_name = false; // true if using id, false otherwise; using a different output file name than the default
  char impurity_type;
  double impurity;
  int id, atom_type, n_types; // atom id to be substituted/removed, atom type to be replaced.
  int N_vac, N_sub, n2; // Number of UO2 molecules left over, N_vac + number of substitutional atoms, number of substituted atoms
  double r_cut, r_cut_sq; // This will need to be specified.  Note that 4.0 is the U-O interaction cutoff that I use.

  inputVars()
  {
    setDefaults();
  }

  void setDefaults()
  {
    datafile = "none";
    impurity = -1.0;
    impurity_type = 'x';
    id = -1.0;
    atom_type = -1.0;
    r_cut = -1.0;
  }

  int validate()
  {
    if (r_cut_sq != r_cut * r_cut) {calculateRCutSq();}

    if (r_cut < 0)
    {
      cout << "r_cut must be positive.\n";
      return INPUT_FORMAT_ERROR;
    }

    if (using_atom_id)
    {
      if (id < 0)
      {
        cout << "id must be an integer greater than 0.\n";
        return INPUT_FORMAT_ERROR;
      }
    }
    else
    {
      if (impurity < 0) {return INPUT_FORMAT_ERROR;}
      if (impurity_type != 'a' && impurity_type != 'w') {return INPUT_FORMAT_ERROR;}
      if (atom_type < 0) {return INPUT_FORMAT_ERROR;}
    }
    return 0;
  }

  void calculateImpurityDetails(int N)
  {
    if (impurity > 1) {n2 = impurity;}
    else if (impurity > 0) {n2 = anInt(N / 3.0 * impurity);}
    else if (impurity < 0 && id > 0) {n2 = 1;}

    N_vac = N - 3 * n2;
    N_sub = N - 2 * n2;
  }

  void calculateRCutSq()
  {
    r_cut_sq = r_cut * r_cut;
  }
} input;

struct boxData
{
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double Lx, Ly, Lz;

  double calculateBoxLengths()
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

// Comparison function for comparing only the second values in a pair
bool pairCmp(pair<int, double> &a, pair<int, double> &b)
{
  return (a.second < b.second);
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

void showInputFileHelp()
{
  cout << "The input file may consist of multiple lines that may be formatted one of two ways:\n"
       << "\t(1) datafile r_cut impurity_value impurity_type atom_type seed\n"
       << "\t(2) datafile r_cut atom_id\n"
       << "Note that currently only the first line is used.\n\n"
       << "In both, the cutoff value must be specified for determining nearest neighbor\n"
       << "interactions. In (1), the impurity value can either be a decimal \n"
       << "(0 < impurity_value < 1) or an integer (>1), impurity_type must be one\n"
       << "character of either w(eight percent) or a(tomic percent), atom_type specifies\n"
       << "where substitutional atoms will be placed, and the seed is used for random\n"
       << "number generation.  In (2), the only parameter required is the atom id number\n"
       << "being replaced.  This is useful for situations where only one defect is\n"
       << "desired, at a specific location.\n\n";
}

void parseInputFile(const string& input_file)
{
  string str;
  int seed;

  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);

  // TODO: Make it so there can be multiple 'jobs' given in an input file.
  // while (getline(fin,str))
  // {
  getline(fin,str);
  stringstream ss(str);
  stringstream::pos_type pos = ss.tellg(); // get the beginning position

  if (!(ss >> input.datafile >> input.r_cut >> input.impurity >> input.impurity_type >> input.atom_type >> seed))
  {
    ss.clear(); // clears the error state of the stream
    ss.seekg(pos, ss.beg); // go back to the beginning of the stream
    input.setDefaults();
    ss >> input.datafile >> input.r_cut >> input.id;
    input.using_atom_id = true;
  }
  else
  {
    // seed was given, so we seed the random number generator
    cout << "The seed for random number generation is " << seed << endl;
    srand(seed);
    if (input.impurity_type == 'w') {weight_percent = true;}
    else if (input.impurity_type == 'a') {atomic_percent = true;}
    else
    {
      cout << "Error: percent-type must be either \'w\' or \'a\'.\n";
      exit(INPUT_FORMAT_ERROR);
    }
  }

  input.calculateRCutSq();

  if (input.validate() != 0)
  {
    showInputFileHelp();
    exit(INPUT_FORMAT_ERROR);
  }
  // }
  fin.close();
}

void readFile(vector <Atom>& substituted_atoms, vector <Atom>& atoms)
{
  string str;
  int N, ntotal = 0; // number of atoms, atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, xy, xz, yz; // box bounds, tilt factors
  double Lx, Ly, Lz; // Box lengths
  int atom_id, atom_type;
  double atom_charge, x, y, z;

  ifstream fin(input.datafile.c_str());
  checkFileStream(fin, input.datafile);

  getline(fin,str); // header line
  fin >> N >> str; // number of atoms

  input.calculateImpurityDetails(N);

  atoms.resize(N, Atom());

  cout << "There will be " << input.n2 << " Xe substitutions (" << (double)(input.n2) / (double)(N / 3.0) * 100.0 << "at%) made.\n";

  fin >> input.n_types >> str >> str; // get the number of atom types

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

  // Actually read the atoms now
  while (fin >> atom_id >> atom_type >> atom_charge >> x >> y >> z) // TODO: This will need to be generalized for charge neutral atoms
  {
    if (fin.fail())
    {
      cout << "Read error\n";
      break;
    }

    ++ntotal;
    if (atom_type > input.n_types)
    {
      cout << "Error: Atom type = " << atom_type << " > n_types = " << input.n_types << endl;
      exit(ATOM_TYPE_ERROR);
    }

    atoms[atom_id - 1] = Atom(atom_id, atom_type, atom_charge, x, y, z);
    if (atom_type == input.atom_type)
    {
      substituted_atoms.push_back(atoms[atom_id - 1]);
    }

  }

  fin.close();

  if (ntotal != N)
  {
    cout << "Error reading atoms: ntotal = " << ntotal << " != N = " << N << endl;
    exit(ATOM_COUNT_ERROR);
  }
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
  ncellx = (int)(box.Lx / input.r_cut) + 1;
  ncelly = (int)(box.Ly / input.r_cut) + 1;
  ncellz = (int)(box.Lz / input.r_cut) + 1;

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
    // Assign this atom to a cell
    // Rounds towards 0 with a type cast
    idx = (int)(atoms[i].getX() / lcellx); // assign the x cell
    idy = (int)(atoms[i].getY() / lcelly); // assign the y cell
    idz = (int)(atoms[i].getZ() / lcellz); // assign the z cell
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
                  rxij = atoms[id].getX() - atoms[jd].getX();
                  ryij = atoms[id].getY() - atoms[jd].getY();
                  rzij = atoms[id].getZ() - atoms[jd].getZ();

                  // Apply PBCs
                  rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
                  ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
                  rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

                  // Now calculate the distance
                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                  if (drij_sq > input.r_cut_sq)
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

void generateImpurities(vector <Atom>& substituted_atoms, vector <Atom>& atoms, const vector <vector <int> >& iatom)
{
  int n_removed = 0, n_2_removed = 0;
  double x, y, z, rxij, ryij, rzij, drij_sq;
  vector <pair <int, double> > distances;

  if (input.using_atom_id)
  {
    atoms[input.id - 1].setMark(1);
    substituted_atoms.resize(1, Atom());
    substituted_atoms[0] = atoms[input.id - 1];
    ++n_removed;
  }
  else
  {
    random_shuffle(substituted_atoms.begin(), substituted_atoms.end()); // randomize which atoms we keep.
    substituted_atoms.resize(input.n2); // Only keeping the first n2 atoms
    for (unsigned int i = 0; i < substituted_atoms.size(); ++i)
    {
      if (atoms[substituted_atoms[i].getId() - 1].getType() != input.atom_type)
      {
        cout << "Error in removing atoms.\n";
        exit(ATOM_TYPE_ERROR);
      }
      atoms[substituted_atoms[i].getId() - 1].setMark(1);
      ++n_removed;
    }
  }

  if (n_removed != input.n2)
  {
    cout << "Error substituting atoms: n_removed = " << n_removed << " != n2 = " << input.n2 << endl;
    exit(ATOM_COUNT_ERROR);
  }

  for (unsigned int i = 0; i < substituted_atoms.size(); ++i) // Remove the whole molecule
  {
    x = atoms[substituted_atoms[i].getId() - 1].getX();
    y = atoms[substituted_atoms[i].getId() - 1].getY();
    z = atoms[substituted_atoms[i].getId() - 1].getZ();

    distances.clear(); // Clear out the last values.
    for (int l = 1; l <= iatom[0][substituted_atoms[i].getId() - 1]; ++l)
    {
      int id = iatom[l][substituted_atoms[i].getId() - 1];
      if (atoms[id].getType() == 2 && atoms[id].getMark() == 0) // logic will need to change to generalize this
      {
        rxij = x - atoms[id].getX();
        ryij = y - atoms[id].getY();
        rzij = z - atoms[id].getZ();

        //Apply PBCs
        rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
        ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
        rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

        drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
        if (drij_sq < input.r_cut_sq)
        {
          distances.push_back(make_pair(id, drij_sq));
        }
      }
    }

    sort(distances.begin(), distances.end(), pairCmp);
    for (unsigned int k = 0; k < distances.size(); ++k)
    {
      int atom_id = distances[k].first;
      if (atoms[atom_id].getType() == 2 &&
          atoms[atom_id].getMark() == 0);
      {
        atoms[atom_id].setMark(2); // We will change those atoms marked 1 to a different type, so we mark these differently
        ++n_2_removed;

        if (n_2_removed % 2 == 0) {break;} // remove atoms in pairs
      }
    }
  }

  if (n_2_removed != 2 * n_removed)
  {
    cout << "Error maintaining charge neutrality!\n"
         << "n_2_removed = " << n_2_removed << " != 2 * n_removed = " << 2 * n_removed << endl; // TODO: This will need to be generalized.
    exit(ATOM_COUNT_ERROR);
  }
}

void writeVacancyFile(const string& vac_file, const vector <Atom>& atoms)
{
  int atom_id = 0;

  ofstream fout(vac_file.c_str());
  checkFileStream(fout, vac_file);
  fout << "These UO2 coordinates have " << input.n2 << " vacancies: [ID type charge x y z]\n\n"
       << input.N_vac << "  atoms\n"
       << input.n_types << "  atom types\n";

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
      fout << ++atom_id << " " << atoms[i].getType() << " "
           << atoms[i].getCharge() << " " << atoms[i].getX() << " "
           << atoms[i].getY() << " " << atoms[i].getZ() << endl;
    }
  }

  fout.close();

  if (atom_id != input.N_vac)
  {
    cout << "Error in vacancy file: n_written = " << atom_id << " != N_vac = " << input.N_vac << endl;
    exit(ATOM_COUNT_ERROR);
  }
}

void writeSubstitutionFile(const string& sub_file, const vector <Atom>& atoms)
{
  // NOTE: this does not necessarily keep the atom id the same between the original
  // file and the new file(s).
  int atom_id = 0;

  ofstream fout(sub_file.c_str());
  checkFileStream(fout, sub_file);
  double atomic_percent = (double)(input.n2) / (double)(atoms.size() / 3.0) * 100.0; // This is too specific.  TODO: generalize
  fout << "These UO2 coordinates have " << input.n2 << " Xe substitutions (" << atomic_percent << "at%Xe): [ID type charge x y z]\n\n"
       << input.N_sub << "  atoms\n"
       << input.n_types + 1 << "  atom types\n";

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
      fout << ++atom_id << " " << atoms[i].getType() << " "
           << atoms[i].getCharge() << " " << atoms[i].getX() << " "
           << atoms[i].getY() << " " << atoms[i].getZ() << endl;
    }
    else if (atoms[i].getMark() == 1) // Marked U atoms -> changed to Xe
    {
      fout << ++atom_id << " " << 3 << " " << 0.0 << " "
           << atoms[i].getX() << " " << atoms[i].getY() << " "
           << atoms[i].getZ() << endl;
    }
    // If neither of these conditions are met, we are looking at the neighbors
    // of the molecule, which will not be written.
  }

  fout.close();

  if (atom_id != input.N_sub)
  {
    cout << "Error in substitution file: n_written = " << atom_id << " != N_sub = " << input.N_sub << endl;
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
         << atoms[i].getCharge() << " " << atoms[i].getX() << " "
         << atoms[i].getY() << " " << atoms[i].getZ() << " "
         << atoms[i].getMark() << endl;
  }

  fout.close();
}

void writeAtomsToFiles(const vector <Atom>& atoms)
{
  string vac_file, sub_file, marked_file;

  if (!input.outfile_name)
  {
    stringstream ss;
    if (input.using_atom_id)
    {
      ss << input.datafile.substr(0,input.datafile.find(".dat")) << "_single_vac.dat";
    }
    else
    {
      ss << input.datafile.substr(0,input.datafile.find(".dat")) << "_" << input.impurity << "vac.dat";
    }
    ss >> vac_file;
    sub_file = vac_file.substr(0,vac_file.find("vac.dat")) + "sub.dat"; // TODO: May be able to generalize this
    marked_file = input.datafile.substr(0,input.datafile.find(".dat")) + "_marked.dat";
  }
  else
  {
    vac_file = input.outfile;
    sub_file = input.outfile;
    marked_file = input.outfile;
  }


  // time_point<Clock> start = Clock::now();
  if (vacancies) {writeVacancyFile(vac_file, atoms);}
  // time_point<Clock> end = Clock::now();
  // milliseconds diff = duration_cast<milliseconds>(end - start);
  // cout << "\t" << diff.count() << "ms for writeVacancyFile\n";
  if (substitutions) {writeSubstitutionFile(sub_file, atoms);}
  // start = Clock::now();
  // diff = duration_cast<milliseconds>(start - end);
  // cout << "\t" << diff.count() << "ms for writeSubstitutionFile\n";
  if (marked) {writeMarkedFile(marked_file, atoms);}
  // end = Clock::now();
  // diff = duration_cast<milliseconds>(end - start);
  // cout << "\t" << diff.count() << "ms for writeMarkedFile\n";

}

int main(int argc, char **argv)
{
  string input_file;
  double impurity; // maybe try to make this a vector so multiple structures can be generated?
  vector <Atom> substituted_atoms, atoms;
  vector <vector <int> > iatom; // Cell-linked list
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
        ("o,output", "Output file name (if only one defect output specified)", cxxopts::value<string>())
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !result.count("file"))
    {
      cout << options.help() << endl << endl;
      showInputFileHelp();
    }

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
        if (type.size() == 1 && result.count("output"))
        {
          input.outfile = result["output"].as<string>();
          input.outfile_name = true;
        }

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

    if (result.count("file"))
    {
      // time_point<Clock> start = Clock::now();
      parseInputFile(input_file);
      // time_point<Clock> end = Clock::now();
      // milliseconds diff = duration_cast<milliseconds>(end - start);
      // cout << diff.count() << "ms for parseInputFile\n";
      readFile(substituted_atoms, atoms);
      // start = Clock::now();
      // diff = duration_cast<milliseconds>(start - end);
      // cout << diff.count() << "ms for readFile\n";
      generateCellLinkedList(atoms, iatom);
      // end = Clock::now();
      // diff = duration_cast<milliseconds>(end - start);
      // cout << diff.count() << "ms for generateCellLinkedList\n";
      generateImpurities(substituted_atoms, atoms, iatom);
      // start = Clock::now();
      // diff = duration_cast<milliseconds>(start - end);
      // cout << diff.count() << "ms for generateImpurities\n";
      writeAtomsToFiles(atoms);
      // end = Clock::now();
      // diff = duration_cast<milliseconds>(end - start);
      // cout << diff.count() << "ms for writeAtomsToFiles\n";
    }
  }
  catch (cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

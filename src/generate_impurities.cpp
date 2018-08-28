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

using namespace std;

double anInt(double x);

//TODO: Make this general for any system, not just UO2
struct inputVars
{
  string datafile;
  double impurity;
  int id; // atom id to be substituted/removed
  int N_vac, N_sub, n2; // Number of UO2 molecules left over, N_vac + number of substitutional atoms, number of substituted atoms
  double r_cut = 4.0, r_cut_sq; // This will need to be specified.  Note that 4.0 is the U-O interaction cutoff that I use.

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

bool weight_percent = false;
bool atomic_percent = false;
bool substitutions = true;
bool vacancies = true;
bool marked = true;

#define AVOGADRO 6.022E23
#define UO_RNN_CUT 4.0 // Cutoff values for U-O atoms too close

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
// Useful later in the program.
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

void readFile()
{
  string str;
  int N, ntypes, ntotal = 0; // number of atoms, atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, xy, xz, yz; // box bounds, tilt factors
  double Lx, Ly, Lz; // Box lengths
  int atom_id, atom_type;
  double atom_charge, x, y, z;
  vector <Atom> atoms, u_atoms;

  ifstream fin(input.datafile.c_str());
  checkFileStream(fin, input.datafile);

  getline(fin,str); // header line
  fin >> N >> str; // number of atoms

  input.calculateImpurityDetails(N);

  atoms.resize(N, Atom());

  cout << "There will be " << input.n2 << " Xe substitutions (" << (double)(input.n2) / (double)(N / 3.0) * 100.0 << "at%) made.\n";

  fin >> ntypes >> str >> str; // get the number of atom types

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
    if (atom_type > ntypes)
    {
      cout << "Error: Atom type = " << atom_type << " > ntypes = " << ntypes << endl;
      exit(ATOM_TYPE_ERROR);
    }

    atoms[atom_id - 1] = Atom(atom_id, atom_type, atom_charge, x, y, z);
    if (atom_type = 1)
    {
      u_atoms.push_back(atoms[atom_id - 1]);
    }

    fin.close();
  }

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

int main(int argc, char **argv)
{
  double impurity; // maybe try to make this a vector so multiple structures can be generated?
  try
  {
    cxxopts::Options options(argv[0], "Generate impurities in a given structure.");
    options
      .positional_help("file seed impurity OR file --substitute-atom atom_id")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Data file to generate impurities in", cxxopts::value<string>(input.datafile), "file")
        ("s,seed", "Seed for the random number generator", cxxopts::value<int>())
        ("i,impurity", "Percent of impurities (0-1) or number of impurities (>=1)", cxxopts::value<double>(input.impurity)->default_value("-1"))
        ("substitute-atom", "ID of the atom to be substituted.  Will remove whole molecule.", cxxopts::value<int>(input.id)->default_value("-1"), "atom_id")
        ("d,defect-type", "(m)arked atoms for removal, atoms with (s)ubstitutions,  or (v)acancies.  Interstitials not implemented.", cxxopts::value<string>()->default_value("msv"), "m, s and/or v")
        ("p,percent-type", "Specify either (w)eight percent, or (a)tomic percent", cxxopts::value<char>(), "w or a")
        ("h,help", "Show the help");

    options.parse_positional({"file", "seed", "impurity"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !result.count("file"))
    {
      cout << options.help() << endl;
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

    if (result.count("impurity"))
    {
      if (impurity > 0 && impurity < 1)
      {
        if (!result.count("percent-type"))
        {
          cout << "Percent type must be specified when using a percentage.";
          return INPUT_FORMAT_ERROR;
        }
        else
        {
          char percent_type = result["percent-type"].as<char>();
          {
            if (percent_type == 'w') {weight_percent = true;}
            else if (percent_type = 'a') {atomic_percent = true;}
            else
            {
              cout << "Error: percent-type must be either \'w\' or \'a\'.\n";
              return INPUT_FORMAT_ERROR;
            }
          }
        }
      }
    }

    if (result.count("file") && result.count("seed") && result.count("impurity"))
    {
      cout << "The seed for random number generation is " << result["seed"].as<int>() << endl;
      srand(result["seed"].as<int>());
      readFile();
    }
    else if (result.count("file") && result.count("substitute-atom"))
    {
      readFile();
    }
  }
  catch (cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

/*int main(int argc, char** argv)
{
  string file1, file2, file3, str; // filename to read, output files, junk var
  int seed; // random number generator seed
  double impurity; // Number of vacancies to generate.
  int N, N_vac, N_sub, ntypes, ntotal = 0, n2; // Number of: atoms, vacancies/subs, atom types, atoms read, U to remove
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounding box / box dimensions
  double xy, xz, yz; // Triclinic tilt terms
  int atom_id, atom_id2, atom_type; // atom id number, type number
  double atom_charge; // atom charge.
  double x, y, z; // position of atom.
  double rxij, ryij, rzij, drij_sq; // positional differences
  double uo_rnn_cut_sq = UO_RNN_CUT * UO_RNN_CUT;
  int n_O_removed, n_U_removed; // self-explanatory

  // Variables used for the cell-linked list
  int n_atoms_per_cell; // self-explanatory
  vector <vector <int> > iatom; // Cell-linked list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell
  int ncellx, ncelly, ncellz, idx, idy, idz; // Number of sub cells in each direction, cell number in each direction
  double lcellx, lcelly, lcellz; // length of sub cells in each direction

  vector <Atom> atoms; // the atoms in our simulation
  vector <Atom> u_atoms; // U atoms in the simulation
  vector <pair <int, double> > distances; // vector of id and distance
  //vector <int> atoms_removed;





  stringstream ss3;
  ss3 << file1.substr(0,file1.find(".dat")) << "_" << impurity << "vac.dat";
  ss3 >> file2;
  file3 = file2.substr(0,file2.find("vac.dat")) + "Xe.dat";



  ofstream fout1(file2.c_str());
  if (fout1.fail())
  {
    cout << "Error opening file " << file2 << endl;
    return 1;
  }

  ofstream fout2(file3.c_str());
  if (fout2.fail())
  {
    cout << "Error opening file " << file3 << endl;
    return 1;
  }




  fout1 << "These UO2 coordinates have " << n2 << " vacancies: [ID type charge x y z]\n\n";
  fout2 << "These UO2 coordinates have " << n2 << " Xe interstitials (" << (double)(n2) / (double)(N) * 100.0 << "at%Xe): [ID type charge x y z]\n\n";
  fout1 << N_vac << "  atoms\n"; // remove a U atom with it's two O neighbors
  fout2 << N_sub << "  atoms\n"; // same as above, but replace U with Xe


  fout1 << ntypes << "   atom types\n";
  fout2 << ntypes + 1 << "   atom types\n"; // Extra atom type because of Xe



  fout1.precision(6);
  fout2.precision(6);

  fout1 << xlow << "\t" << xhigh << "\txlo xhi\n";
  fout1 << ylow << "\t" << yhigh << "\tylo yhi\n";
  fout1 << zlow << "\t" << zhigh << "\tzlo zhi\n";
  fout2 << xlow << "\t" << xhigh << "\txlo xhi\n";
  fout2 << ylow << "\t" << yhigh << "\tylo yhi\n";
  fout2 << zlow << "\t" << zhigh << "\tzlo zhi\n";



  fout1 << "\nAtoms\n\n";
  fout2 << "\nAtoms\n\n";



  // Generate the cell-linked list for fast calculations
  // First generate the number of cells in each direction (minimum is 1)
  ncellx = (int)(Lx / UO_RNN_CUT) + 1;
  ncelly = (int)(Ly / UO_RNN_CUT) + 1;
  ncellz = (int)(Lz / UO_RNN_CUT) + 1;
  lcellx = Lx / ncellx; // Length of the cells in each direction
  lcelly = Ly / ncelly;
  lcellz = Lz / ncellz;

  // Number of atoms per cell based on cell size, with a minimum allowed of 200
  n_atoms_per_cell = max((int)(N / (double)(3 * ncellx * ncelly * ncellz)), 200);

  // resizes the vectors to be the correct length. Saves on time.
  // Defaults all values to 0
  icell.resize(ncellx, vector <vector <int> > // x dimension
              (ncelly, vector <int> // y dimension
              (ncellz, 0))); // z dimension
  pcell.resize(ncellx, vector <vector <vector <int> > > // x dimension
              (ncelly, vector <vector <int> > // y dimension
              (ncellz, vector <int> // z dimension
              (n_atoms_per_cell, 0)))); // atom number in cell.
  iatom.resize(n_atoms_per_cell, vector <int> (N,0)); // the actual list.

  /****************************************************************************
  /**************************CREATE CELL-LINKED LIST***************************
  /****************************************************************************
  for (unsigned int i = 0; i < atoms.size(); ++i) // Look at each atom
  {
    //if (atoms[i].getType() != 1) continue; // Only want U atoms
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
                  rxij = rxij - anInt(rxij / Lx) * Lx;
                  ryij = ryij - anInt(ryij / Ly) * Ly;
                  rzij = rzij - anInt(rzij / Lz) * Lz;

                  // Now calculate the distance
                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                  if (drij_sq > uo_rnn_cut_sq)
                  {
                    continue; // move to the next atom if we're too far away
                  }

                  if (drij_sq == 0.0) // This should never be hit, but just in case
                  {
                    continue; // This is the same atom!
                  }

                  // Create the neighbor list
                  iatom[0][id] += 1; //for atom id
                  iatom[(iatom[0][id])][id] = jd;
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
  ****************************************************************************
  **********************END GENERATE CELL-LINKED LIST*************************
  ****************************************************************************

  random_shuffle(u_atoms.begin(), u_atoms.end()); // Shuffle the atoms around
  u_atoms.resize(n2); // Keep the first n2 atoms.
  // Mark these U atoms for removal/replacement
  for (unsigned int i = 0; i < u_atoms.size(); ++i)
  {
    if (atoms[u_atoms[i].getId() - 1].getType() != 1)
    {
      cout << "Error in removing U atoms.\n";
      return 8;
    }
    atoms[u_atoms[i].getId() - 1].setMark(1);
    ++n_U_removed;
  }

  if (n_U_removed != n2)
  {
    cout << "Error marking U atoms! n_U_removed = " << n_U_removed << " != n2 = " << n2 << endl;
    return 5;
  }

  for (unsigned int i = 0; i < u_atoms.size(); ++i) // We need to mark the O atoms for removal now.
  {
    // store this atom's position.
    x = atoms[u_atoms[i].getId() - 1].getX();
    y = atoms[u_atoms[i].getId() - 1].getY();
    z = atoms[u_atoms[i].getId() - 1].getZ();

    distances.clear(); // Clear out the last values.
    for (int l = 1; l <= iatom[0][u_atoms[i].getId() - 1]; ++l) // for each neighbor of this atom.
    {
      int id = iatom[l][u_atoms[i].getId() - 1];
      if (atoms[id].getType() == 2 && atoms[id].getMark() == 0)
      {
        rxij = x - atoms[id].getX();
        ryij = y - atoms[id].getY();
        rzij = z - atoms[id].getZ();

        //Apply PBCs
        rxij = rxij - anInt(rxij / Lx) * Lx;
        ryij = ryij - anInt(ryij / Ly) * Ly;
        rzij = rzij - anInt(rzij / Lz) * Lz;

        drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
        if (drij_sq < uo_rnn_cut_sq)
        {
          distances.push_back(make_pair(id, drij_sq));
        }
      }
    }
    // sort the distances using the comparison function written above.
    sort(distances.begin(), distances.end(), pairCmp);
    for (unsigned int k = 0; k < distances.size(); ++k)
    {
      atom_id = distances[k].first;
      // If the atom we are looking at is O and is unmarked
      if (atoms[atom_id].getType() == 2 &&
          atoms[atom_id].getMark() == 0)
      {
        atoms[atom_id].setMark(2); // Different mark because we will change the U to Xe later.
        ++n_O_removed;

        // make sure we only remove atoms in pairs of two!
        if (n_O_removed % 2 == 0)
        {
          break;
        }
      }
    }
  }

  if (n_O_removed != 2 * n_U_removed)
  {
    cout << "Error maintaining charge neutrality!\n"
         << "n_O_removed = " << n_O_removed << " != 2 * n_U_removed = " << 2 * n_U_removed << endl;;
    return 6;
  }

  // Now let's write to the respective files!
  atom_id = 0;
  atom_id2 = 0;
  ofstream fout3((file2.substr(0,file2.find(".dat")) + "_marked.dat").c_str());
  if (fout3.fail())
  {
    cout << "Error opening file " << file2.substr(0,file2.find(".dat")) << "_marked.dat\n";
    return 1;
  }
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    fout3 << atoms[i].getId() << " " << atoms[i].getType() << " "
          << atoms[i].getCharge() << " " << atoms[i].getX() << " "
          << atoms[i].getY() << " " << atoms[i].getZ() << " "
          << atoms[i].getMark() << endl;
    if (atoms[i].getMark() == 0)
    {
      ++atom_id;
      ++atom_id2;
      fout1 << atom_id << " " << atoms[i].getType() << " "
            << atoms[i].getCharge() << " " << atoms[i].getX() << " "
            << atoms[i].getY() << " " << atoms[i].getZ() << endl;
      fout2 << atom_id2 << " " << atoms[i].getType() << " "
            << atoms[i].getCharge() << " " << atoms[i].getX() << " "
            << atoms[i].getY() << " " << atoms[i].getZ() << endl;
    }
    else if (atoms[i].getMark() == 1) // These are the marked U atoms
    {
      ++atom_id2;
      fout2 << atom_id2 << " " << 3 << " " << 0.0 << " "
            << atoms[i].getX() << " " << atoms[i].getY() << " "
            << atoms[i].getZ() << endl;
    }
  }

  if (atom_id != N_vac)
  {
    cout << "Error in vacancies! n_written = " << atom_id << " != N_vac = " << N_vac << endl;
    return 7;
  }
  if (atom_id2 != N_sub)
  {
    cout << "Error in substitutions! n_written = " << atom_id << " != N_sub = " << N_sub << endl;
    return 7;
  }

  return 0;
}*/

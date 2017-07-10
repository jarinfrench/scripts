#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
//#include <random> // See https://msdn.microsoft.com/en-us/library/bb982398.aspx <-- requires C++11
#include "atom.h"

using namespace std;

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

int main(int argc, char** argv)
{
  string file1, file2, file3, str; // filename to read, output files, junk var
  int seed, impurity; // random number generator seed
  int N, N_vac, N_sub, ntypes, ntotal = 0, n2; // Number of: atoms, vacancies/subs, atom types, atoms read, U to remove
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounding box / box dimensions
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

  if (argc != 4)
  {
    cout << "Please enter the rotated grain file name: ";
    cin  >> file1;

    cout << "Please enter the percent of impurities to generate (0-1): ";
    cin  >> impurity;

    cout << "Please enter the seed for the random number generator: ";
    cin  >> seed;
  }
  else
  {
    file1 = argv[1];
    // Convert command line arguments to integers.
    istringstream ss(argv[2]);
    if (!(ss >> impurity))
    {
      cout << "Error converting argument 2 to integer.\n";
      return 2;
    }
    istringstream ss2(argv[3]);
    if (!(ss2 >> seed))
    {
      cout << "Invalid seed " << seed << endl;
      return 2;
    }
  }

  // Make sure we aren't doing anything weird with the impurity level
  if (impurity >= 100 || impurity <= 0)
  {
    cout << "Invalid value of impurity percentage: 0 < impurity < 100\n";
    return 3;
  }

  cout << "The seed for random number generation is " << seed << endl;
  srand(seed); // Seed the random number generator.
  //mt19937 gen(rd()); // Seed the mersenne twister (requires C++11)
  //uniform_int_distribution<> dist(1,6); // distribute results between 1 and 6 inclusive.

  stringstream ss3;
  ss3 << file1.substr(0,file1.find(".dat")) << "_" << impurity << "vac.dat";
  ss3 >> file2;
  file3 = file2.substr(0,file2.find("vac.dat")) + "Xe.dat";

  // open the file for reading.
  ifstream fin(file1.c_str());
  if (fin.fail())
  {
    cout << "Unable to read file " << file1 << endl;
    return 1;
  }

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


  // Read the data!
  getline(fin, str); // Get the header line
  fout1 << "These UO2 coordinates have vacancies: [ID type charge x y z]\n\n";
  fout2 << "These UO2 coordinates have Xe interstitials: [ID type charge x y z]\n\n";

  // Get the number of atoms
  fin >> N >> str;
  n2 = anInt(N / 3.0 * impurity / 100.0);
  N_vac = N - 3 * n2; // calculate the number of UO2 vacancies
  N_sub = N - 2 * n2; // Calculate how many atoms are left after replace UO2 with Xe
  atoms.resize(N, Atom());

  fout1 << N_vac << "  atoms\n"; // remove a U atom with it's two O neighbors
  fout2 << N_sub << "  atoms\n"; // same as above, but replace U with Xe

  // get the number of atom types
  fin >> ntypes >> str >> str;
  fout1 << ntypes << "   atom types\n";
  fout2 << ntypes + 1 << "   atom types\n"; // Extra atom type because of Xe

  // Get the bounds of the system
  fin >> xlow >> xhigh >> str >> str;
  fin >> ylow >> yhigh >> str >> str;
  fin >> zlow >> zhigh >> str >> str;

  fout1.precision(6);
  fout2.precision(6);

  fout1 << xlow << "\t" << xhigh << "\txlo xhi\n";
  fout1 << ylow << "\t" << yhigh << "\tylo yhi\n";
  fout1 << zlow << "\t" << zhigh << "\tzlo zhi\n";
  fout2 << xlow << "\t" << xhigh << "\txlo xhi\n";
  fout2 << ylow << "\t" << yhigh << "\tylo yhi\n";
  fout2 << zlow << "\t" << zhigh << "\tzlo zhi\n";

  Lx = xhigh - xlow;
  Ly = yhigh - ylow;
  Lz = zhigh - zlow;

  fin >> str; // Read the extra stuff.

  fout1 << "\nAtoms\n\n";
  fout2 << "\nAtoms\n\n";

  // Now lets read in the atoms.  If the current atom being read is a U atom,
  // run a test to see if we remove/replace it.
  while (fin >> atom_id >> atom_type >> atom_charge >> x >> y >> z)
  {
    if (fin.fail())
    {
      cout << "Read error\n";
      break;
    }

    ++ntotal;
    if (atom_type > ntypes)
    {
      cout << "Error! Atom_type = " << atom_type << " > " << ntypes << endl;
      return 4;
    }
    atoms[atom_id - 1] = Atom(atom_id, atom_type, atom_charge, x, y, z);
    if (atom_type == 1)
    {
      u_atoms.push_back(Atom(atom_id, atom_type, atom_charge, x, y, z));
    }
  }

  if (ntotal != N)
  {
    cout << "Error reading atoms! ntotal = " << ntotal << " != N = " << N << endl;
    return 4;
  }
  fin.close(); // done reading the file.

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

  /****************************************************************************/
  /**************************CREATE CELL-LINKED LIST***************************/
  /****************************************************************************/
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
  /****************************************************************************/
  /**********************END GENERATE CELL-LINKED LIST*************************/
  /****************************************************************************/

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
         << "n_O_removed = " << n_O_removed << " != 2 * n_U_removed = " << 2 * n_U_removed;
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
}

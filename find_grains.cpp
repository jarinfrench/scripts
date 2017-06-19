#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath> // for cos, sin
#include <cstdlib>
#include "atom.h"

using namespace std;

#define PI 3.141592653589793

// Cutoff distances
#define UU_CUT 0.866 // from Bai et al. Acta Materialia 85 (2015) 95-106: 0.866 * a0
#define a0 5.6 // typically 5.453
#define IDEAL 1.0 // 6.0 // This is dependent on the coefficients chosen for coeffs!

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

int main(int argc, char** argv)
{
  double coeffs [2] = {3,2}; // Coefficients of the symmetry parameter

  string filename1, filename2, str; // filenames read from and written to, junk variable
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounds variables
  int N, n_type, n_atoms_read = 0; // number of atoms, atom types, number of atoms read
  vector <Atom> atoms; // all of the atoms from the file
  vector <double> symm; // a vector to hold the calculated symmetry parameters
  int atom_id, type; // id and type number of atom; used to read in the data
  double charge, x, y, z; // charge and position of atom, used to read in data
  double rxij, ryij, rzij, drij_sq; // positions and distance squared
  double uu_cut_sq = UU_CUT * UU_CUT; // cutoff distance
  double theta, sintheta_sq, total = 0.0; // misorientation angle, sin^2 of the angle, symmetry parameter (it's a sum, so starts at 0)
  double xtemp, ytemp, sintheta, costheta, cutoff; // rotated x position, y position, sin theta, cos theta, cutoff for which grain an atom belongs to.
  bool dump; // boolean value to determine if the read file is a LAMMPS dump file or not.
  unsigned int n_grain_1, n_grain_2; // counter for number of atoms in each grain

  // Variables used for the cell-linked list
  int n_atoms_per_cell; // self-explanatory
  vector <vector <int> > iatom; // Cell-linked list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell
  int ncellx, ncelly, ncellz, idx, idy, idz; // Number of sub cells in each direction, cell number in each direction, atoms in cell i
  double lcellx, lcelly, lcellz; // length of sub cells in each direction

  // Parse command line input. Prompt for what we need.
  if (argc == 2)
  {
    filename1 = argv[1];
    cout << "Please enter the misorientation angle of the grain in degrees: ";
    cin  >> theta;
  }
  else if (argc == 3)
  {
    filename1 = argv[1];
    theta = strtod(argv[2], NULL);
  }
  else
  {
    cout << "Please enter the data file to be read: ";
    cin  >> filename1;
    cout << "Please enter the misorientation angle of the grain in degrees: ";
    cin  >> theta;
  }
  filename2 = filename1.substr(0,filename1.find(".dat")) + "_interface.dat";
  sintheta = sin(theta * PI / 180.0); // best to calculate this once
  costheta = cos(theta * PI / 180.0);

  // Calculate the cutoff value for grain identification
  vector <double> xx (12,0.0); // x positions in terms of a0 for nearest neighbors
  vector <double> yy (12,0.0); // y positions in terms of a0 for nearest neighbors
  xx[0] = 0.5; xx[1] = 0.5; xx[2] = -0.5; xx[3] = -0.5; xx[4] = 0.5; xx[5] = 0.5; xx[6] = -0.5; xx[7] = -0.5;
  yy[0] = 0.5; yy[1] = -0.5; yy[2] = 0.5; yy[3] = -0.5; yy[8] = 0.5; yy[9] = 0.5; yy[10] = -0.5; yy[11] = -0.5;

  for (unsigned int i = 0; i < xx.size(); ++i)
  {
    xtemp = costheta * xx[i] - sintheta * yy[i];
    ytemp = sintheta * xx[i] + costheta * yy[i];
    // Uses the idea that sin^2 = 1-cos^2
    // cos = A.B / (|A||B|), and since B is (1,0,0), this simplifies to
    // cos = A_x / |A|, meaning that cos^2 = A_x^2/(A_x^2+A_y^2)
    sintheta_sq = 1 - ((xtemp * xtemp) / (xtemp * xtemp + ytemp * ytemp));
    // Symmetry parameter as defined by Zhang. Coefficients are changed to get
    // better resolution between grains.
    total += (coeffs[0] - coeffs[1] * sintheta_sq) * (coeffs[0] - coeffs[1] * sintheta_sq) * sintheta_sq;
  }
  total /= xx.size();
  // Cutoff is the midpoint between the ideal value and the rotated ideal value.
  // Note that this assumes that the outside grain is oriented with it's x axis
  // aligned with the x axis of the lab frame.
  // FIXME: this needs to be generalized for outside grains that aren't aligned
  cutoff = (IDEAL + total) / 2.0;

  // Output the values for documentation.
  cout << "The value for grain 1 is " << IDEAL
       << " and the value for grain 2 is " << total << endl;

  // Open up the files for reading and writing.
  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error opening file " << filename1 << endl;
    return 1;
  }

  ofstream fout(filename2.c_str());
  if (fout.fail())
  {
    cout << "Error opening file " << filename2 << endl;
    return 1;
  }

  // Pull out the relevant information from the heading
  // Determine if the file is a dump file or an input file
    if (filename1.find("dump") == string::npos) // if "dump" is not in the filename
  {
    dump = false; // it isn't a dump file, assume the file is an input file.
  }
  else
  {
    dump = true; // if it is in the filename, assume that it's a LAMMPS dump file
  }

  if (dump)
  {
    // This is for a LAMMPS dump file
    getline(fin, str); // Gets ITEM: TIMESTEP
    getline(fin, str); // Gets the timestep number
    getline(fin, str); // Gets ITEM: NUMBER OF ATOMS
    fin >> N;
    fin.ignore();
    getline(fin, str); //get ITEM: BOX BOUNDS
    fin >> xlow >> xhigh;
    fin >> ylow >> yhigh;
    fin >> zlow >> zhigh;
    fin.ignore();
    getline(fin, str); // Gets ITEM: ATOMS <data types>
    n_type = 2; // Assumes only two types of atoms: U and O
  }
  else
  {
    // This is for a LAMMPS input file
    getline(fin, str); // gets comment line
    fin >> N >> str; // number of atoms
    fin >> n_type >> str >> str;
    fin >> xlow >> xhigh >> str >> str;
    fin >> ylow >> yhigh >> str >> str;
    fin >> zlow >> zhigh >> str >> str;
    fin >> str;
  }
  // Convert the bounds in terms of a0
  xlow /= a0;
  xhigh /= a0;
  ylow /= a0;
  yhigh /= a0;
  zlow /= a0;
  zhigh /= a0;
  Lx = xhigh - xlow;
  Ly = yhigh - ylow;
  Lz = zhigh - zlow;

  // Preallocate the information (saves time, and allows the atoms to be written
  // in order)
  atoms.resize(N, Atom());
  // Read the data
  while (fin >> atom_id >> type >> charge >> x >> y >> z)
  {
    if (type > n_type)
    {
      cout << "Error: unexpected atom type.\n"
           << "n_types = " << n_type << " < this atom's type = " << type << endl;
      return 2;
    }

    // Read the atoms line by line, and put them into a vector for analysis.
    // We adjust the positions so that we start at the origin so we can easily assign to cells
    x = x / a0 - xlow;
    y = y / a0 - ylow;
    z = z / a0 - zlow;
    // We make the atom id match (almost) the index.  There is a difference of 1
    // because C++ indexes from 0.
    atoms[atom_id - 1] = Atom(atom_id, type, charge, x, y, z);
    ++n_atoms_read;
  }
  fin.close(); // Close the data file, we're done with it.

  // Compare to N
  if (n_atoms_read != N)
  {
    cout << "Error: number of atoms read does not match number of atoms in the simulation.\n"
         << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
    return 3;
  }

  // Generate the cell-linked list for fast calculations
  // First generate the number of cells in each direction (minimum is 1)
  ncellx = (int)(Lx / UU_CUT) + 1;
  ncelly = (int)(Ly / UU_CUT) + 1;
  ncellz = (int)(Lz / UU_CUT) + 1;
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
    if (atoms[i].getType() != 1) continue; // Only want U atoms
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

                  if (drij_sq > uu_cut_sq)
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

  // Now that we have the atoms safely stored, we can process them.
  symm.resize(atoms.size(), 0); // Assign each atom a symmetry parameter
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getType() != 1)
    {
      continue; // Only focus on U
    }

    x = atoms[i].getX();
    y = atoms[i].getY();
    z = atoms[i].getZ();
    total = 0.0; // reset the total
    // We start at l = 1 because if we start at l = 0, we just re-use the same
    // atom over and over.
    for (int l = 1; l <= iatom[0][i]; ++l)
    {
      int id = iatom[l][i];

      // calculate the distances
      // We project onto the xy plane, effectively ignoring the z coordinate
      rxij = atoms[id].getX() - x;
      ryij = atoms[id].getY() - y;

      // Apply PBCs
      rxij = rxij - anInt(rxij / Lx) * Lx;
      ryij = ryij - anInt(ryij / Ly) * Ly;

      // Calculate the magnitude of the distance
      drij_sq = (rxij * rxij) + (ryij * ryij);

      // This uses the relation sin^2 = 1-cos^2, where cos = dot(A,B) / (|A|*|B|)
      sintheta_sq = 1 - (rxij * rxij) / drij_sq;
      total += (coeffs[0] - coeffs[1] * sintheta_sq) * (coeffs[0] - coeffs[1] * sintheta_sq) * sintheta_sq;
    }
    total /= iatom[0][i]; // This may not always be 12!
    symm[i] = total; // Store them for analysis

    if (atoms[i].getType() != 1)
    {
      continue;
    }
    if (symm[i] <= cutoff)
    {
      atoms[i].setMark(1);
      ++n_grain_1;
    }
    else // greater than the cutoff
    {
      atoms[i].setMark(2);
      ++n_grain_2;
    }
  }

  cout << "There are " << n_grain_1 << " U atoms in grain 1, and " << n_grain_2
       << " U atoms in grain 2.\n";
  // Make sure we write the entire set of atoms
  // This writes things in a tecplot-readable format.
  n_atoms_read = 0;
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getType() != 1)
    {
      continue;
    }
    fout << atoms[i].getId() << " "
         << atoms[i].getType() << " "
         << atoms[i].getCharge() << " "
         << (atoms[i].getX() + xlow) * a0 << " "
         << (atoms[i].getY() + ylow) * a0 << " "
         << (atoms[i].getZ() + zlow) * a0 << " "
         << atoms[i].getMark() << " "
         << symm[i] << endl;
    ++n_atoms_read;
  }
  fout.close(); // Close the output file

  if (n_atoms_read != N / 3.0)
  {
    cout << "Error: number of atoms written does not match number of atoms in the simulation.\n"
         << "N = " << N / 3.0 << " != n_atoms_read = " << n_atoms_read << endl;
    return 6;
  }

  return 0;
}

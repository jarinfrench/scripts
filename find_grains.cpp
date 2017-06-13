#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath> // for acos, sin
#include "atom.h"

using namespace std;

#define PI 3.141592653589793

// Cutoff distances
#define UU_CUT 4.8496 // from Bai et al. Acta Materialia 85 (2015) 95-106: 0.866 * a0
//#define UO_CUT 2.5 // actual value is 0.43301
//#define OO_CUT 3.7 // actual value is 2.7625

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

// Comparison functions for sorting a vector of pairs
bool pairCmp(pair<int, double> &a, pair<int, double> &b)
{
  return (a.second == b.second);
}

bool pairSort(pair<int, double> &a, pair <int, double> &b)
{
  return (a.second < b.second);
}

int main(int argc, char** argv)
{
  string input_file, filename, filename2, str; // input file, and file to read/write, junk variable
  int nfiles, N, n_atoms_read; // number of files to read, number of atoms, number of atom types, number of atoms read
  vector <Atom> atoms; // all of the atoms from the file
  double rxij, ryij, rzij, drij_sq; // positions
  double r_cut, a0, n_cut; // cutoff distance^2, lattice parameter, cutoff distance in terms of a0

  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounds variables
  int id, type; // atom id, atom type
  double charge, x, y, z; // charge and position of atom

  int n_nn, n_atoms_per_cell;  //number of nearest neighbors, atoms per sub cell
  vector <double> r_nnx (12), r_nny (12), r_nnz (12), rxijk, ryijk, rzijk; // nearest neighbor fcc positions
  vector <vector <double> > r_mkx (12, vector <double> (90)), r_mky (12, vector <double> (90)), r_mkz (12, vector <double> (90)); // ideal nearest neighbor positions
  double rij, theta, xtemp, ytemp, costheta, sintheta, sintheta_sq, ornt_prt, magnitude; // 2*distance^2, angle, rotated x and y, cos and sin of theta, orientation parameter
  vector <double> symm; // Symmetry parameter

  vector <vector <int> > iatom; // Cell-linked list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell
  int ncellx, ncelly, ncellz, idx, idy, idz, n_icell; // Number of sub cells in each direction, cell number in each direction, atoms in cell i
  double lcellx, lcelly, lcellz; // length of sub cells in each direction

  if (argc != 2)
  {
    cout << "Please enter the input file name: ";
    cin  >> input_file;
  }
  else
  {
    input_file = argv[1];
  }

  ifstream fin(input_file.c_str());
  if (fin.fail())
  {
    cout << "Error opening file " << input_file << endl;
    return 1;
  }

  // Get the important information from the input file
  fin >> nfiles >> n_cut >> n_nn >> a0; // gets number of files, cutoff distance, number of nearest neighbors and lattice parameter

  // Resize the vectors we will use
  //r_nnx.resize(n_nn);
  //r_nny.resize(n_nn);
  //r_nnz.resize(n_nn);
  //r_mkx.resize(n_nn);
  //r_mky.resize(n_nn);
  //r_mkz.resize(n_nn);

  // Get the nearest neighbor positions
  for (int i = 0; i < n_nn; ++i)
  {
    fin >> r_nnx[i] >> r_nny[i] >> r_nnz[i];
    rij = (r_nnx[i] * r_nnx[i] + r_nny[i] * r_nny[i] +r_nnz[i] * r_nnz[i]) * 2;
  }

  rij = sqrt(rij);
  // Distances in terms of a0
  for (unsigned int i = 0; i < r_nnx.size(); ++i)
  {
    r_nnx[i] /= rij;
    r_nny[i] /= rij;
    r_nnz[i] /= rij;
  }

  // Create an array for all the rotated positions for each angle between 0 and 89
  for (unsigned int i = 0; i < 90; ++i)
  {
    for (unsigned int j = 0; j < 12; ++j)
    {
      theta = i * PI / 180.0;
      costheta = cos(theta);
      sintheta = sin(theta);
      xtemp = costheta * r_nnx[j] - sintheta * r_nny[j]; // Rotated about the z axis
      ytemp = sintheta * r_nnx[j] + costheta * r_nny[j]; // Rotated about the z axis

      // Store the rotated positions.
      r_mkx[j][i] = xtemp;
      r_mky[j][i] = ytemp;
      r_mkz[j][i] = r_nnz[j];
    }
  }

  // Now we will read each file listed
  for (int file = 0; file < nfiles; ++file)
  {
    fin >> filename; // get the next filename to be read.
    filename2 = filename.substr(0,filename.find(".dat")) + "_interface.dat";
    ifstream fin2(filename.c_str());
    if (fin2.fail())
    {
      cout << "Error opening file " << filename << endl;
      return 2;
    }
    // Now we read the header of the dump file
    getline(fin2, str); // ITEM: TIMESTEP
    getline(fin2, str); // <time step number>
    getline(fin2, str); // ITEM: NUMBER OF ATOMS
    fin2 >> N; // Get the number of atoms in the simulation
    fin2.ignore(); // Clear the newline from the buffer
    getline(fin2, str); // ITEM: BOX BOUNDS pp pp pp
    fin2 >> xlow >> xhigh;
    fin2 >> ylow >> yhigh;
    fin2 >> zlow >> zhigh;
    fin2.ignore(); // Again, ignore the newline in the buffer
    getline(fin2, str); // ITEM: ATOMS id type q x y z

    // Set the dimensions in terms of a0
    xlow /= a0;
    xhigh /= a0;
    ylow /= a0;
    yhigh /= a0;
    zlow /= a0;
    zhigh /= a0;
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;

    // Make sure to reset everything for every file
    atoms.clear();
    n_atoms_read = 0;
    while (fin2 >> id >> type >> charge >> x >> y >> z)
    {
      Atom temp(id, type, charge, x / a0 - xlow, y / a0 - ylow, z / a0 - zlow);
      atoms.push_back(temp);
      ++n_atoms_read;
    }

    if (n_atoms_read != N)
    {
      cout << "Error: number of atoms read does not match number of atoms in the simulation.\n"
           << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
      return 2;
    }

    ncellx = (int)(Lx / n_cut) + 1;
    ncelly = (int)(Ly / n_cut) + 1;
    ncellz = (int)(Lz / n_cut) + 1;
    lcellx = Lx / ncellx;
    lcelly = Ly / ncelly;
    lcellz = Lz / ncellz;

    n_atoms_per_cell = (int)(N / (double)(ncellx * ncelly * ncellz));
    n_atoms_per_cell = max(n_atoms_per_cell, 200);

    // resizes the vectors to be the correct length for this dump file.
    icell.resize(ncellx, vector <vector <int> > (ncelly, vector <int> (ncellz)));
    pcell.resize(ncellx, vector <vector <vector <int> > > (ncelly, vector <vector <int> > (ncellz, vector <int> (n_atoms_per_cell))));
    iatom.resize(n_atoms_per_cell, vector <int> (N,0));

    // This creates the cell-linked list
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      // Only focusing on U atoms now
      if (atoms[i].getType() != 1) continue;

      // Assign this atom to a cell
      // Rounds to 0 with a type cast
      idx = (int)(atoms[i].getX() / lcellx);
      idy = (int)(atoms[i].getY() / lcelly);
      idz = (int)(atoms[i].getZ() / lcellz);
      // Check if we went out of bounds
      // C++ indexes from 0, so we have to subtract 1 from the maximum value to stay within our memory bounds
      if (idx >= ncellx) idx = ncellx - 1;
      if (idy >= ncelly) idy = ncelly - 1;
      if (idz >= ncellz) idz = ncellz - 1;

      icell[idx][idy][idz] += 1; // count the number of atoms in this cell
      n_icell = icell[idx][idy][idz]; // number of the atom in this cell
      pcell[idx][idy][idz][n_icell] = i;
    }
    r_cut = n_cut * n_cut;
    for (int i = 0; i < ncellx; ++i)
    {
      for (int j = 0; j < ncelly; ++j)
      {
        for (int k = 0; k < ncellz; ++k)
        {
          for (int l = 1; l < icell[i][j][k]; ++l)
          {
            int id = pcell[i][j][k][l];
            // Now we check each sub cell around the current one
            for (int ii = -1; ii < 2; ++ii)
            {
              for (int jj = -1; jj < 2; ++jj)
              {
                for (int kk = -1; kk < 2; ++kk)
                {
                  int ia = i + ii;
                  int ja = j + jj;
                  int ka = k + kk;
                  // Check to make sure we are still in bounds
                  if (ia >= ncellx) ia = 0;
                  if (ja >= ncelly) ja = 0;
                  if (ka >= ncellz) ka = 0;
                  if (ia < 0) ia = ncellx - 1;
                  if (ja < 0) ja = ncelly - 1;
                  if (ka < 0) ka = ncellz - 1;

                  for (int m = 0; m < icell[ia][ja][ka]; ++m)
                  {
                    int jd = pcell[ia][ja][ka][m];
                    // If jd < id, we've already dealt with this interaction
                    if (jd < id) continue;

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

                    if (drij_sq > r_cut) continue; // move to the next atom if we're too far away

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

    // Not quite sure why 0.32 is used.
    costheta = cos(0.32); // approximately 18 degrees
    sintheta = sin(0.32);

    symm.resize(N, 0.0);
    rxijk.resize(n_atoms_per_cell, 0.0);
    ryijk.resize(n_atoms_per_cell, 0.0);
    rzijk.resize(n_atoms_per_cell, 0.0);

    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      if (atoms[i].getType() != 1) continue; // Only focusing on U atoms

      for (int l = 0; l < iatom[0][i]; ++l) // Each atom in the neighbor list
      {
        rxijk.clear();
        ryijk.clear();
        rzijk.clear();
        int id = iatom[l][i];

        // Calculate the distance
        rxijk[l] = atoms[id].getX() - atoms[i].getX();
        ryijk[l] = atoms[id].getY() - atoms[i].getY();
        rzijk[l] = atoms[id].getZ() - atoms[i].getZ();

        // Apply PBCs
        rxijk[l] = rxijk[l] - anInt(rxijk[l] / Lx) * Lx;
        ryijk[l] = ryijk[l] - anInt(ryijk[l] / Ly) * Ly;
        rzijk[l] = rzijk[l] - anInt(rzijk[l] / Lz) * Lz;

        // Note that this is a rotation about the x axis.
        xtemp = costheta * ryijk[l] - sintheta * rzijk[l];
        ytemp = sintheta * ryijk[l] + costheta * rzijk[l];

        rxij = xtemp * xtemp;
        ryij = ytemp * ytemp;
        drij_sq = rxij + ryij;
        if (drij_sq == 0.0)
        {
          rxij = 1.0;
        }
        else
        {
          rxij /= drij_sq;
          ryij /= drij_sq;
        }

        // Bai code
        symm[i] += (3 - 2 * rxij) * (3 - 2 * rxij) * rxij;

        // My code
        /*magnitude = sqrt(rxijk[l] * rxijk[l] + ryijk[l] * ryijk[l] + rzijk[l] * rzijk[l]);
        theta = acos(rxijk[l] / magnitude);
        sintheta = sin(theta);
        sintheta_sq = sintheta * sintheta;
        symm[i] += (3 - 4 * sintheta_sq) * (3 - 4 * sintheta_sq) * sintheta_sq - sintheta_sq;*/
      }
      symm[i] /= iatom[0][i];

      // This is the Trautt formulation (Acta Materialia 60 (2012) 2407-2424, Eq. 33)
      /*int amax = -500000;
      for (unsigned int j = 0; j < 90; ++j)
      {
        ornt_prt = 0.0;
        for (int l = 0; l < iatom[0][i]; ++l)
        {
          for (int m = 0; m < n_nn; ++m)
          {
            rxij = rxijk[l] - r_mkx[m][j];
            ryij = ryijk[l] - r_mky[m][j];
            rzij = rzijk[l] - r_mkz[m][j];
            ornt_prt += exp(-rxij * rxij);
            ornt_prt += exp(-ryij * ryij);
            ornt_prt += exp(-rzij * rzij);
          }
        }
        if (ornt_prt >= amax) amax = ornt_prt;
        if (ornt_prt >= amax) symm[i] = j;
      }
      //End Trautt*/
    }

    /*vector <int> nalpha (90,0);
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      for (unsigned int j = 0; j < 90; ++j)
      {
        if (symm[i] == j) nalpha[j] += 1;
      }
    }

    ofstream fout1("ideal_angles.csv");
    if (fout1.fail())
    {
      cout << "Error opening file ideal_angles.csv\n";
      return 3;
    }
    for (unsigned int i = 0; i < 90; ++i)
    {
      fout1 << i << " " << nalpha[i] << endl;
    }
    fout1.close();*/

    // Now we write the results to a file.
    ofstream fout(filename2.c_str());
    if (fout.fail())
    {
      cout << "Error: Unable to open file " << filename2 << endl;
      return 4;
    }
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      if (atoms[i].getType() != 1) continue;
      int grain_num;
      if (symm[i] > 1.25)
      {
        grain_num = 1;
      }
      else
      {
        grain_num = 2;
      }
      // Write the results, converting distances back to the original lengths
      fout << atoms[i].getId() << " " << atoms[i].getType() << " "
           << atoms[i].getCharge() << " " << atoms[i].getX() * a0 << " "
           << atoms[i].getY() * a0 << " " << atoms[i].getZ() * a0 << " "
           << grain_num << " " << symm[i] << endl;;

    }
    fout.close();
    fin2.close();
  } // End reading file

  return 0;
}

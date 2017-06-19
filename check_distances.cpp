#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include "atom.h"

using namespace std;

#define OO_RNN_CUT 0.63801
#define UO_RNN_CUT 4.0 // nearest U-O neighbors.  Actual value is 2.3612
#define UU_RNN_CUT 3.0
#define NN_UU 5.0 //3.8559 // NN U-U atoms with a0 = 5.453

//Returns the rounded value of x
double anInt(double x)
{
  int temp; //temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

int main(int argc, char **argv)
{
  string filename, str; //file containing our atoms, and a string to contain lines
  int id, type; // atom id number and type number
  double r_grain, r_grain_m, r_grain_p; // rotated grain radius
  double r_grain_m_sq, r_grain_p_sq; // squared values for distances
  double charge, x, y, z; // atom charge value, and position
  double x1, y1, z1; // temp variables
  double rxij, ryij, rzij, drij_sq; // positional differences
  double xhigh, xlow, yhigh, ylow, zhigh, zlow; // bounds
  double Lx, Ly, Lz; // box size
  vector<Atom> atoms_checked, atoms; //contains all of our atoms
  vector<pair<double, int> > distances; // distances and id of close U-O atoms
  double oo_rnn_cut_sq = OO_RNN_CUT * OO_RNN_CUT;
  double uu_rnn_cut_sq = UU_RNN_CUT * UU_RNN_CUT;
  int n_O_too_close = 0, n_U_too_close = 0; //counter variable

  if (argc != 3)
  {
    cout << "Please enter the filename containing the atoms: ";
    cin  >> filename;

    cout << "Please enter the radius of the rotated grain: ";
  }
  else
  {
    filename = argv[1];
    r_grain = strtod(argv[2], NULL);
  }

  r_grain_m = r_grain - NN_UU;
  r_grain_p = r_grain + NN_UU;
  r_grain_m_sq = r_grain_m * r_grain_m; // we use the squared values a lot,
  r_grain_p_sq = r_grain_p * r_grain_p; // so just calculate once

  ifstream fin(filename.c_str());
  if (fin.fail())
  {
    cout << "Error! Cannot read file " << filename << endl;
    return -1;
  }

  getline(fin, str); // Header line
  fin >> str >> str; // number of atoms
  fin >> str >> str >> str; // number of atom types
  fin >> xlow >> xhigh >> str >> str; // x bounds
  fin >> ylow >> yhigh >> str >> str; // y bounds
  fin >> zlow >> zhigh >> str >> str; // z bounds
  fin >> str; // Atoms line

  // Bounds
  Lx = xhigh - xlow;
  Ly = yhigh - ylow;
  Lz = zhigh - zlow;

  // Get all the data
  while (fin >> id >> type >> charge >> x >> y >> z)
  {
    x1 = x - Lx / 2.0; // center the system at the origin
    y1 = y - Ly / 2.0;
    z1 = z - Lz / 2.0;
    drij_sq = x1 * x1 + y1 * y1 + z1 * z1;

    // if it's in the boundary range
    if (drij_sq < r_grain_p_sq && drij_sq > r_grain_m_sq)
    {
      // make a subset
      Atom a(id, type, charge, x, y, z);
      atoms_checked.push_back(a);
    }

    // put all atoms in this vector
    Atom b(id, type, charge, x, y, z);
    atoms.push_back(b);
  }

  // for each atom in the subset, check how close the closest O atoms are
  for (unsigned int i = 0; i < atoms_checked.size(); ++i)
  {
    x1 = atoms_checked[i].getX();
    y1 = atoms_checked[i].getY();
    z1 = atoms_checked[i].getZ();
    if (atoms_checked[i].getType() == 2 && atoms_checked[i].getMark() == 0) // looking at unmarked oxygen
    {
      for (unsigned int j = i + 1; j < atoms_checked.size(); ++j)
      {
        if (atoms_checked[j].getType() == 2 && atoms_checked[j].getMark() == 0) // looking at unmarked oxygen
        {
          // calculate the distance
          rxij = x1 - atoms_checked[j].getX();
          ryij = y1 - atoms_checked[j].getY();
          rzij = z1 - atoms_checked[j].getZ();

          // Apply PBCs;
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;
          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
          if (drij_sq < oo_rnn_cut_sq) // If the O atoms are too close, mark the second for removal
          {
            atoms_checked[j].setMark(1);
            atoms[atoms_checked[j].getId() - 1].setMark(1);
            ++n_O_too_close;
          }
        }
      }
    }

    else if (atoms_checked[i].getType() == 1 && atoms_checked[i].getMark() == 0) // looking at unmarked U
    {
      for (unsigned int j = i + 1; j < atoms_checked.size(); ++j)
      {
        if (atoms_checked[j].getType() == 1 && atoms_checked[j].getMark() == 0) // looking at unmarked U
        {
          // calculate the distance
          rxij = x1 - atoms_checked[j].getX();
          ryij = y1 - atoms_checked[j].getY();
          rzij = z1 - atoms_checked[j].getZ();

          // Apply PBCs;
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;
          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
          if (drij_sq < uu_rnn_cut_sq) // If the U atoms are too close, mark the second for removal
          {
            atoms_checked[j].setMark(1);
            atoms[atoms_checked[j].getId() - 1].setMark(1);
            ++n_U_too_close;
          }
        }
      }
    }
  }
  cout << "There are " << n_U_too_close << " U atoms that are too close to each other.\n";
  cout << "There are " << n_O_too_close << " O atoms that are too close to each other.\n";

  return 0;
}

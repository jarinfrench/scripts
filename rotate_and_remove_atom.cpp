#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include "atom.h"

using namespace std;
// Note: this conversion is for UO2 only!
// For UO2, charge style, the charge must be specified
// The sequence is atom-ID atom-type q x y z

#define PI 3.14159265
#define SKIN 16.0 //skin depth (just under 3a0, a0 = 5.453)
#define UU_RNN_CUT 2.5 // Cutoff value for U-U atoms too close
#define UO_RNN_CUT 4.0 // Cutoff value for U-O atoms too close
#define OO_RNN_CUT 0.63801 // Cutoff value for O-O atoms too close

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

// Comparison function for comparing the only the second values in a pair
bool pairCmp(pair<int, double> &a, pair<int, double> &b)
{
  return (a.second < b.second);
}

int main(int argc, char **argv)
{
  // External values
  string filename1, filename2, filename3, filename4, str; //filenames and line variable
  double r_grain, r_grain_m, r_grain_p; //radius of the grain, with buffer zone
  double r_grain_sq, r_grain_m_sq, r_grain_p_sq; // squared values for convenience
  double theta; // angle of rotation
  double costheta, sintheta; // better to calculate this once.
  double uu_rnn_cut_sq = UU_RNN_CUT * UU_RNN_CUT; //easier to do it once
  double uo_rnn_cut_sq = UO_RNN_CUT * UO_RNN_CUT;
  double oo_rnn_cut_sq = OO_RNN_CUT * OO_RNN_CUT;
  double scale_factor_a, scale_factor_b, scale_factor_c; // dimensional scaling values
  double Lx, Ly, Lz; // box size
  int ntotal, n_atom_id; // total number of atoms that have been read/written
  double rxij, ryij, rzij, drij_sq; //positional differences, total distance^2
  int n_U_removed = 0, n_O_removed = 0; // Counters for U and O removal.

  // Values from file
  int N, ntypes; // number of atoms, and number of atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh; // atom bounds
  int atom_id, atom_type; // id number and type number
  double atom_charge, x, y, z; // charge and position values
  double x1, y1, z1, temp_x, temp_y, x2, y2, z2; // Store the original value and manipulate!

  // Containers
  vector <Atom> atoms_checked, atoms; // contains the atoms we look at, and the entire set.
  vector <pair<int, double> > distances; // vector of id and distance.

  if (argc != 4) // check command line arguments
  {
    // filename
    cout << "Please input the filename in LAAMPS's format at 0K:\n";
    cin  >> filename1;

    // Grain radius
    cout << "Please input the radius of the new grain:\n";
    cin  >> r_grain;

    // Rotation angle
    cout << "Please input the rotation angle (in degrees) of the new grain:\n";
    cin  >> theta;
    if (abs(theta) > 180.0)
      cout << "Caution!  The rotation angle is greater than 180 degrees!\n";
  }
  else
  {
    filename1 = argv[1];
    r_grain = strtod(argv[2], NULL);
    theta = strtod(argv[3], NULL);
  }

  r_grain_m = r_grain - SKIN;
  r_grain_p = r_grain + SKIN;
  r_grain_m_sq = r_grain_m * r_grain_m; // we use the squared values a lot,
  r_grain_p_sq = r_grain_p * r_grain_p; // so just calculate once
  r_grain_sq = r_grain * r_grain;


  ostringstream fn2, fn3, fn4; // String streams for easy file naming
  fn2 << filename1.substr(0,filename1.find(".")).c_str() << "_" << theta
      //<< "degree_r" << r_grain << "A_rotated.dat";
      << "degree_r" << r_grain << "A_rotated_rcut" << UU_RNN_CUT << ".dat";
  filename2 = fn2.str();

  fn3 << filename1.substr(0,filename1.find(".")).c_str() << "_" << theta
      //<< "degree_r" << r_grain << "A_marked.dat";
      << "degree_r" << r_grain << "A_marked_rcut" << UU_RNN_CUT << ".dat";
  filename3 = fn3.str();

  theta *= PI / 180.0; // convert theta to radians
  costheta = cos(theta); // just calculate this once!
  sintheta = sin(theta);

  ifstream fin(filename1.c_str()); // only reading this file
  if (fin.fail())
  {
    cout << "Error opening the file " << filename1 << endl;
    return -1;
  } // End error check

  ofstream fout(filename2.c_str()); // only writing to this file
  fout << fixed; // makes sure to always use the precision specified.
  if (fout.fail()) //Error check
  {
    cout << "Error opening the file " << filename2 << endl;
    return -1;
  } // End error check

  getline(fin, str);
  fout << "These UO2 coordinates are shifted: [ID type charge x y z]\n";
  fout << "\n";

  //Get the number of atoms
  fin  >> N >> str;
  fout << N << "  atoms\n";

  //Get the number of atom types
  fin  >> ntypes >> str >> str;
  fout << ntypes << "   atom types\n";

  // Get the bounds of the system
  fin  >> xlow >> xhigh >> str >> str;
  fin  >> ylow >> yhigh >> str >> str;
  fin  >> zlow >> zhigh >> str >> str;

  // Sets the scaling factor.
  scale_factor_a = 1.0;
  scale_factor_b = 1.0;
  scale_factor_c = 1.0;

  // Calculate the bounding box size in each direction
  Lx = xhigh - xlow;
  Ly = yhigh - ylow;
  Lz = zhigh - zlow;

  // Check to make sure the diameter of the new grain is smaller than the box
  // dimensions
  if (r_grain * 2.0 >= Lx)
  {
    cout << "Error! Grain diameter = " << r_grain * 2.0 << " >= Lx = " << Lx << endl;
    return -2;
  }

  if (r_grain * 2.0 >= Ly)
  {
    cout << "Error! Grain diameter = " << r_grain * 2.0 << " >= Ly = " << Ly << endl;
    return -2;
  }

  // Scale the dimensions
  xlow  *= scale_factor_a;
  xhigh *= scale_factor_a;

  ylow  *= scale_factor_b;
  yhigh *= scale_factor_b;

  zlow  *= scale_factor_c;
  zhigh *= scale_factor_c;

  // Write the modified max and mins to the new file
  fout.precision(6);
  fout << xlow << "\t" << xhigh << "\txlo xhi\n";
  fout << ylow << "\t" << yhigh << "\tylo yhi\n";
  fout << zlow << "\t" << zhigh << "\tzlo zhi\n";

  fin  >> str; // Read the extra stuff

  fout << "\nAtoms\n\n"; // We now want to write the atoms down

  ntotal = 0; // Number of atoms so far
  n_atom_id = 0; // number written so far

  while (fin >> atom_id >> atom_type >> atom_charge >> x >> y >> z)
  {
    // Check for reading errors
    if (fin.fail())
    {
      cout << "Read error\n";
      break;
    }

    ++ntotal; // Increment the number of atoms

    // Make sure there aren't more than ntypes of atoms
    if (atom_type > ntypes)
    {
      cout << "Error! Atom_type = " << atom_type << " is greater than " << ntypes << endl;
      return -3;
    }
    // change the origin to the center of the simulation for rotating the atoms
    x1 = x - Lx / 2.0;
    y1 = y - Ly / 2.0;
    z1 = z - Lz / 2.0;

    // If we are smaller than the radius, rotate the atom by theta around the
    // z axis.
    // TODO: make this an option to do twist or tilt boundaries for 100, 110, and 111 axes
    if ((x1 * x1 + y1 * y1) <= (r_grain_sq))
    {
      // This is <100> Twist
      temp_x = x1 * costheta - y1 * sintheta;
      temp_y = x1 * sintheta + y1 * costheta;
      x1 = temp_x;
      y1 = temp_y;
    }

    if (x1 * x1 + y1 * y1 > r_grain_m_sq &&
        x1 * x1 + y1 * y1 < r_grain_p_sq) // If the atom is in our range of interest
    {
      Atom a(atom_id, atom_type, atom_charge, x1, y1, z1); // store checked atoms
      atoms_checked.push_back(a);
    }
    x1 += Lx / 2.0;
    y1 += Ly / 2.0;
    z1 += Lz / 2.0;

    x1 *= scale_factor_a;
    y1 *= scale_factor_b;
    z1 *= scale_factor_c;

    // Write the rotated atoms to their own file
    ++ n_atom_id; // increment atom ID
    fout.precision(0);
    fout << n_atom_id << " " << atom_type << " ";
    fout.precision(1);
    fout << atom_charge << " ";
    fout.precision(6);
    fout << x1 << " " << y1  << " " << z1 << endl;

    Atom b(atom_id, atom_type, atom_charge, x1, y1, z1); // store all atoms
    atoms.push_back(b);
  }

  // Make sure we read all of the atoms
  if (ntotal != N)
  {
    cout << "ntotal = " << ntotal << " != N = " << N << endl;
    return -4;
  } // End error check
  fin.close(); // Done reading the file

  /*****************************************************************************
  ****************************ATOM REMOVAL**************************************
  *****************************************************************************/
  // Compare the distances of each atom

  for (int i = 0; i < atoms_checked.size() - 1; ++i)
  {
    // Get the positions of atom i
    x1 = atoms_checked[i].getX();
    y1 = atoms_checked[i].getY();
    z1 = atoms_checked[i].getZ();
    // Checking U atoms for being too close
    if (atoms_checked[i].getType() == 1 && atoms_checked[i].getMark() == 0) // only checking unmarked U atoms
    {
      for (int j = i + 1; j < atoms_checked.size(); ++j) // Don't double count!
      {
        if (atoms_checked[j].getType() == 1 && atoms_checked[j].getMark() == 0) // Only unmarked U atoms
        {
          // Calculate the distance
          rxij = x1 - atoms_checked[j].getX();
          ryij = y1 - atoms_checked[j].getY();
          rzij = z1 - atoms_checked[j].getZ();

          // Apply PBCs
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;

          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
          if (drij_sq < uu_rnn_cut_sq) // if the U atoms are too close
          {
            atoms_checked[i].setMark(1);
            // The -1 is here because the id's start at 1, and c++ indices start at 0
            atoms[atoms_checked[i].getId() - 1].setMark(1);
            ++n_U_removed;
            break; // since we've removed this atom, move on.
          }
        }
      }
    }

    // If this atom is U and has been marked
    if (atoms_checked[i].getType() == 1 && atoms_checked[i].getMark() == 1)
    {
      distances.clear(); // Clear out the old values.
      // Check each O atom, and find the closest two and remove them
      for (int j = 0; j < atoms_checked.size(); ++j)
      {
        if (i == j)
          continue;
        // Only check unmarked O atoms
        if (atoms_checked[j].getType() == 2 && atoms_checked[j].getMark() == 0)
        {
          // Calculate the distance
          rxij = x1 - atoms_checked[j].getX();
          ryij = y1 - atoms_checked[j].getY();
          rzij = z1 - atoms_checked[j].getZ();

          // Apply PBCs
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;

          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
          if (drij_sq < uo_rnn_cut_sq) // if smaller than the cutoff distance
          {
            distances.push_back(make_pair(j, drij_sq));
          }
        }
      }
      // sort the distances using the comparison function written above.
      sort(distances.begin(), distances.end(), pairCmp);
      for (int k = 0; k < distances.size(); ++k)
      {
        atom_id = distances[k].first; // we use this a lot over the next lines
        // If the atom we are looking at is O and unmarked
        if (atoms_checked[atom_id].getType() == 2 &&
            atoms_checked[atom_id].getMark() == 0)
        {
          // mark this atom in BOTH lists
          atoms_checked[atom_id].setMark(1);
          // The -1 is here because the id's start at 1, and c++ indices start at 0
          atoms[atoms_checked[atom_id].getId() - 1].setMark(1);
          ++n_O_removed; // increase the counter for O removed

          // if we have removed enough O atoms to maintain charge neutrality,
          // exit the loop.  We can only remove 2 O atoms per U atom!
          if (n_O_removed == 2 * n_U_removed)
            break;
        }
      }
    }
  }

  // Now, go through the list again, and remove the O atoms that are too close
  for (int i = 0; i < atoms_checked.size() - 1; ++i)
  {
    if (atoms_checked[i].getType() == 2 && atoms_checked[i].getMark() == 0) // Looking at unmarked O atoms
    {
      // Get the position of atom i
      x1 = atoms_checked[i].getX();
      y1 = atoms_checked[i].getY();
      z1 = atoms_checked[i].getZ();

      for (int j = i + 1; j < atoms_checked.size(); ++j)
      {
        if (atoms_checked[j].getType() == 2 && atoms_checked[j].getMark() == 0) // Unmarked O atoms
        {
          // Calculate the distance
          rxij = x1 - atoms_checked[j].getX();
          ryij = y1 - atoms_checked[j].getY();
          rzij = z1 - atoms_checked[j].getZ();

          // Apply PBCs
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;

          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

          if (drij_sq < oo_rnn_cut_sq)
          {
            // We will remove both of these
            atoms_checked[i].setMark(1);
            atoms_checked[j].setMark(1);
            atoms[atoms_checked[i].getId() - 1].setMark(1);
            atoms[atoms_checked[j].getId() - 1].setMark(1);
            n_O_removed += 2;

            // Now go through and find the closest U atom to these two
            x2 = atoms_checked[j].getX();
            y2 = atoms_checked[j].getY();
            z2 = atoms_checked[j].getZ();
            for (int k = 0; k < atoms_checked.size(); ++k)
            {
              if (atoms_checked[k].getType() == 1 && atoms_checked[k].getMark() == 0) //unmarked U atoms
              {
                //calculate the distances between the U atom and both O atoms
                rxij = x1 - atoms_checked[k].getX();
                ryij = y1 - atoms_checked[k].getY();
                rzij = z1 - atoms_checked[k].getZ();

                // Apply PBCs
                rxij = rxij - anInt(rxij / Lx) * Lx;
                ryij = ryij - anInt(ryij / Ly) * Ly;
                rzij = rzij - anInt(rzij / Lz) * Lz;

                drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
                if (drij_sq < uo_rnn_cut_sq)
                {
                  rxij = x2 - atoms_checked[k].getX();
                  ryij = y2 - atoms_checked[k].getY();
                  rzij = z2 - atoms_checked[k].getZ();

                  // Apply PBCs
                  rxij = rxij - anInt(rxij / Lx) * Lx;
                  ryij = ryij - anInt(ryij / Ly) * Ly;
                  rzij = rzij - anInt(rzij / Lz) * Lz;

                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
                  if (drij_sq < uo_rnn_cut_sq)
                  {
                    // Only if both O atoms are close enough to the U atom is
                    // this U atom removed
                    atoms_checked[k].setMark(1);
                    atoms[atoms_checked[k].getId() - 1].setMark(1);
                    ++n_U_removed;
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  cout << n_U_removed << " U atoms will be removed.\n";
  cout << n_O_removed << " O atoms will be removed.\n";

  // Error checking
  if (n_U_removed * 2 != n_O_removed)
  {
    cout << "Error: the removed U:O ratio must be 1:2!\n";
    return -5;
  }

  ofstream fout2(filename3.c_str());
  fout2 << fixed;
  if (fout2.fail())
  {
    cout << "Error opening the file " << filename3 << endl;
    return -2;
  }

  fn4 << filename1.substr(0,filename1.find("N")).c_str() << anInt(theta*180 / PI)
      << "degree_r" << r_grain
      //<< "A_removed.dat";
      << "A_removed_rcut" << UU_RNN_CUT << ".dat";
  filename4 = fn4.str();

  ofstream fout3(filename4.c_str());
  fout3 << fixed;
  if (fout3.fail()) // error check
  {
    cout << "Error opening the file " << filename4 << endl;
    return -2;
  }

  // write the base data to the file
  fout3 << "These UO2 coordinates are shifted and have atoms removed:[ID type charge x y z]\n"
        << "\n"
        << N - n_U_removed-n_O_removed << "   atoms\n"
        << ntypes << "   atom types\n"
        << xlow << " " << xhigh << "   xlo xhi\n"
        << ylow << " " << yhigh << "   ylo yhi\n"
        << zlow << " " << zhigh << "   zlo zhi\n"
        << "\nAtoms\n\n";

  // Now write the atoms to the files.  filename3 has all the atoms including
  // the rotated ones and the tag. filename4 has the correct number of atoms.
  ntotal = 0;
  for (int i = 0; i < atoms.size(); ++i)
  {
    fout2 << atoms[i].getId() << " " << atoms[i].getType() << " "
          << atoms[i].getCharge() << " " << atoms[i].getX() << " "
          << atoms[i].getY() << " " << atoms[i].getZ() << " "
          << atoms[i].getMark() << endl;

    if (atoms[i].getMark() == 0) // We only write the atoms that are NOT marked for removal
    {
      ++ntotal;
      fout3 << ntotal << " " << atoms[i].getType() << " "
            << atoms[i].getCharge() << " " << atoms[i].getX() << " "
            << atoms[i].getY() << " " << atoms[i].getZ() << endl;
    }
  }

  if (ntotal != N - n_U_removed - n_O_removed) // One last check
  {
    cout << "Error! The final number of removed atoms is not balanced!\n"
         << "ntotal = " << ntotal << " != N - n_U_removed - n_O_removed = "
         << N - n_U_removed - n_O_removed << endl;
    return -6;
  }
  fout2.close();
  fout3.close();
  return 0;
}

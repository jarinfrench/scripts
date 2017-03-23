#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include "atom.h"

using namespace std;
//Note: this conversion is for UO2 only!
// For UO2, charge style, the charge must be specified
// The sequence is atom-ID atom-type q x y z

#define PI 3.14159265
#define NN_UU 10.0 //3.8559 // NN U-U atoms with a0 = 5.453
#define UU_RNN_CUT 3.0 // Nearest U-U neighbors
#define UO_RNN_CUT 4.0 // nearest U-O neighbors
#define OO_RNN_CUT 0.63801 // Nearest O-O neighbors

//Returns the rounded value of x
double anInt(double x)
{
  int temp; //temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

bool pairCmp(pair<int, double>& a, pair<int, double>& b)
{
  return (a.second < b.second);
}

int main(int argc, char **argv)
{
  string filename1, filename2, filename3, filename4, str; // filename of file to be read/written to
  double r_grain, r_grain_m, r_grain_p; // radius of the grain with buffer
  double theta; // angle of rotation
  vector<Atom> atoms_checked, atoms; // contains the atoms we look at
  vector<pair<int, double> > distances; //contains the id and distances between
  vector <pair<int, int> > u_vals; // used to removed u atoms
  int N, ntypes, ntotal, n_atom_id; // Number of atoms, number of atom types, number read, number written
  double xlow, xhigh, ylow, yhigh, zlow, zhigh; // Bounding box
  double scale_factor_a, scale_factor_b, scale_factor_c; //scaling factors for bounding box
  double Lx, Ly, Lz; // Box Size
  int atom_id, atom_type; // id and type numbers
  double atom_charge, x, y, z; // atom charge and position
  double x1, y1, z1, rxij, ryij, rzij, drij_sq; // shifted position, differences in position
  double temp_x, temp_y; // temp position values.
  double uu_rnn_cut_sq = UU_RNN_CUT * UU_RNN_CUT; //easier to do it once
  double uo_rnn_cut_sq = UO_RNN_CUT * UO_RNN_CUT;
  double oo_rnn_cut_sq = OO_RNN_CUT * OO_RNN_CUT;
  double costheta, sintheta; // we only want to use these functions once
  double r_grain_m_sq, r_grain_p_sq, r_grain_sq; // buffer radii squared
  int n_U_removed = 0, n_O_removed = 0, n_prev; // number of U and O atoms removed
  typedef map<int, unsigned int> CounterMap; // for counting
  CounterMap counts; // counter map

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

  r_grain_m = r_grain - NN_UU;
  r_grain_p = r_grain + NN_UU;
  r_grain_m_sq = r_grain_m * r_grain_m; // we use the squared values a lot,
  r_grain_p_sq = r_grain_p * r_grain_p; // so just calculate once
  r_grain_sq = r_grain * r_grain;

  // Note that the rotation axis is in the Z direction

  ostringstream fn2, fn3, fn4; // string streams to make naming files easy
  fn2 << filename1.substr(0,filename1.find(".")).c_str() << "_" << theta
      //<< "degree_r" << r_grain << "A_rotated.dat";
      << "degree_r" << r_grain << "A_rotated_rcut" << UU_RNN_CUT << ".dat";
  filename2 = fn2.str();

  fn3 << filename1.substr(0,filename1.find(".")).c_str() << "_" << theta
      //<< "degree_r" << r_grain << "A_marked.dat";
      << "degree_r" << r_grain << "A_marked_rcut" << UU_RNN_CUT << ".dat";
  filename3 = fn3.str();

  theta *= PI / 180.0; // convert theta to radians
  costheta = cos(theta);
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

  // Check to make sure the radius of the new grain is smaller than the box
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
  // Skip this if the scale_factors are 1
  if (scale_factor_a != 1.0 && scale_factor_b != 1.0 && scale_factor_c != 1.0)
  {
    xlow  *= scale_factor_a;
    xhigh *= scale_factor_a;

    ylow  *= scale_factor_b;
    yhigh *= scale_factor_b;

    zlow  *= scale_factor_c;
    zhigh *= scale_factor_c;
  }

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
  for (int i = 0; i < atoms_checked.size(); ++i)
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
          if (drij_sq < uu_rnn_cut_sq)
          {
            atoms_checked[i].setMark(1);
            atoms[atoms_checked[i].getId() - 1].setMark(1);
            ++n_U_removed;
            break;
          }
        }
      }
    }

    // Checking O atoms for being too close
    else if (atoms_checked[i].getType() == 2 && atoms_checked[i].getMark() == 0) // only checking unmarked O atoms
    {
      for (int j = i + 1; j < atoms_checked.size(); ++j)
      {
        if (atoms_checked[j].getType() == 2 && atoms_checked[j].getMark() == 0) // only unmarked O atoms
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
            atoms_checked[i].setMark(1);
            atoms[atoms_checked[i].getId() - 1].setMark(1);
            ++n_O_removed;
            break;
          }
        }
      }
    }
  }

  //Debugging
  //cout << "There are " << n_U_removed << " U atoms being removed after the first loop.\n"
  //     << "There are " << n_O_removed << " O atoms being removed after the first loop.\n";

  //cout << "So far " << n_U_removed << " U atoms will be removed.\n";
  //cout << "So far " << n_O_removed << " O atoms will be removed.\n";
  //num_marked = n_O_removed; // keep track of how many O atoms we've already removed

  // Now check the U atoms that are being removed, and remove their O neighbors
  for (int i = 0; i < atoms_checked.size(); ++i)
  {
    distances.clear();
    if (atoms_checked[i].getType() == 1 && atoms_checked[i].getMark() == 1) // only check marked U atoms
    {
      x1 = atoms_checked[i].getX();
      y1 = atoms_checked[i].getY();
      z1 = atoms_checked[i].getZ();
      for (int j = 0; j < atoms_checked.size(); ++j)
      {
        if (atoms_checked[j].getType() == 2 && atoms_checked[j].getMark() == 0) // only check O atoms that are unmarked
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
          if (drij_sq < uo_rnn_cut_sq)
          {
            distances.push_back(make_pair(j, drij_sq));
          }
        }
      }
      sort(distances.begin(), distances.end(), pairCmp); // sort the distances from least to greatest
      n_prev = n_O_removed;
      for (int k = 0; k < distances.size(); ++k)
      {
        if (atoms_checked[distances[k].first].getType() == 2 &&
            atoms_checked[distances[k].first].getMark() == 0)
        {
          atoms_checked[distances[k].first].setMark(1);
          atoms[atoms_checked[distances[k].first].getId() - 1].setMark(1);
          ++n_O_removed;
          if (n_O_removed == n_prev + 2) // only remove two atoms per U atom removed
            break;
        }
      }
    }
  }

  //Debugging
  //cout << "There are " << n_U_removed << " U atoms being removed after the second loop.\n"
  //     << "There are " << n_O_removed << " O atoms being removed after the second loop.\n";

  // Loop through the atoms again, this time checking the U atoms that are the
  // closest to two O atoms.
  distances.clear();
  for (int i = 0; i < atoms_checked.size(); ++i)
  {
    if (atoms_checked[i].getType() == 2 && atoms_checked[i].getMark() == 1) // Looking at marked O atoms
    {
      x1 = atoms_checked[i].getX();
      y1 = atoms_checked[i].getY();
      z1 = atoms_checked[i].getZ();
      for (int j = 0; j < atoms_checked.size(); ++j)
      {
        if (atoms_checked[j].getType() == 1 && atoms_checked[j].getMark() == 0) // Looking at unmarked U atoms
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
          if (drij_sq < uo_rnn_cut_sq)
          {
            distances.push_back(make_pair(j, drij_sq));
          }
        }
      }
    }
  }
  sort(distances.begin(), distances.end()); // we want the U atom closest to two O atoms

  // A bit of complicated coding that simply counts the occurrences of each
  // index from atoms_checked.
  for (int i = 0; i < distances.size(); ++i)
  {
    CounterMap::iterator k(counts.find(distances[i].first)); // iterate through the map and count the number of occurrences of each object
    // The keys for this map are the indices of the atoms_checked vector
    if (k != counts.end())
      ++k->second;
    else
      counts[distances[i].first] = 1; // if it's not found, the number of counts is 1
  }

  // Now iterate through that map, and append to the vector that keeps track of
  // which U atoms have the most O atoms too close
  for (CounterMap::iterator k = counts.begin(); k != counts.end(); ++k)
  {
    //k->first is atoms_checked index, k->second is number of occurrences
    u_vals.push_back(make_pair(k->second, k->first));
  }
  sort(u_vals.rbegin(), u_vals.rend()); // sorts from greatest to least

  int i = 0;
  for (int i = 0; i < u_vals.size(); ++i)
  {
    // check to make sure we don't remove more atoms than necessary
    if (2 * n_U_removed == n_O_removed)
      break;

    // make sure the atoms aren't already marked
    if (atoms_checked[u_vals[i].second].getMark() == 1)
    {
      cout << "Atom " << u_vals[i].second << " has already been marked!\n";
      continue;
    }
    // otherwise, mark for removal
    atoms_checked[u_vals[i].second].setMark(1);
    atoms[atoms_checked[u_vals[i].second].getId()].setMark(1);
    ++n_U_removed;
  }

  //Debugging
  //cout << "There are " << n_U_removed << " U atoms being removed after the third loop.\n"
  //     << "There are " << n_O_removed << " O atoms being removed after the third loop.\n";
  // One last check for O atoms too close
  /*for (int i = 0; i < atoms_checked.size(); ++i)
  {
    if (atoms_checked[i].getType() == 2 && atoms_checked[i].getMark() == 0) // only checking unmarked O atoms
    {
      x1 = atoms_checked[i].getX();
      y1 = atoms_checked[i].getY();
      z1 = atoms_checked[i].getZ();
      for (int j = i + 1; j < atoms_checked.size(); ++j)
      {
        if (atoms_checked[j].getType() == 2 && atoms_checked[j].getMark() == 0) // only unmarked O atoms
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
            atoms_checked[j].setMark(1);
            atoms[atoms_checked[j].getId()].setMark(1);
            ++n_O_removed;
          }
        }
      }
    }
  }*/
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

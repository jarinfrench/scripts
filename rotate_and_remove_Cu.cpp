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
#include <algorithm>
#include "atom.h"

using namespace std;
// Note: this conversion is for UO2 only!
// For UO2, charge style, the charge must be specified
// The sequence is atom-ID atom-type q x y z

#define PI 3.14159265358979 // easier and faster to simply store these values here.
#define SKIN 10.0 //skin depth (just under 3a0, a0 = 3.615)
#define RNN_CUT 1.0 // Cutoff value for Cu-Cu atoms too close (using Mishin potential)

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

int main(int argc, char **argv)
{
  // External values
  string filename1, filename2, filename3, filename4, str; //filenames and line variable
  int axis;
  double r_grain; //radius of the grain, with buffer zone
  double r_grain_sq; // squared values for convenience
  double theta, theta_conv; // angle of rotation
  double costheta, sintheta; // better to calculate this once.
  double rnn_cut_sq = RNN_CUT * RNN_CUT; //easier to do it once

  double scale_factor_a, scale_factor_b, scale_factor_c; // dimensional scaling values
  double Lx, Ly, Lz; // box size
  int ntotal, n_atom_id; // total number of atoms that have been read/written
  double rxij, ryij, rzij, drij_sq; //positional differences, total distance^2
  int n_Cu_removed = 0; // Counters for Cu removal.
  bool decimal; // boolean value to include decimal or not

  // Values from file
  int N, ntypes; // number of atoms, and number of atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh; // atom bounds
  int atom_id, atom_type; // id number and type number
  double x, y, z; // charge and position values
  double x1, y1, z1, temp_x, temp_y; // Store the original value and manipulate!

  // Containers
  vector <Atom> atoms; // contains the atoms we look at, and the entire set.
  vector <pair<int, double> > distances; // vector of id and distance.

  // Variables used for the cell-linked list
  int n_atoms_per_cell; // self-explanatory
  vector <vector <int> > iatom; // Cell-linked list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell
  int ncellx, ncelly, ncellz, idx, idy, idz; // Number of sub cells in each direction, cell number in each direction
  double lcellx, lcelly, lcellz; // length of sub cells in each direction

  if (argc != 4) // check command line arguments
  {
    // filename
    // Format of filename should be LAMMPS_UO2_SC_N######_{axis}.dat
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
    if (theta - (int)theta == 0)
      decimal = false;
    else
      decimal = true;
  }
  else
  {
    filename1 = argv[1];
    // the NULL is required to make this work.  see documentation for why.
    r_grain = strtod(argv[2], NULL);
    theta = strtod(argv[3], NULL);
    for (unsigned int i = 0; i < strlen(argv[3]); ++i) // check if a decimal point exists in argument 3
    {
      if (argv[3][i] == '.')
      {
        decimal = true; // it exists!
        break;
      }
      else
        decimal = false; // it doesn't exist
    }
  }

  // Calculate this once
  r_grain_sq = r_grain * r_grain;

  // String streams make things easy.  This particular line requires that the
  // filename be formatted such that before the file extension, the axis is
  // specified.  This means that the file must be of the format:
  // LAMMPS_UO2_SC_N######_{axis}.dat
  istringstream iss (filename1.substr(filename1.find(".") - 3, 3));
  if (!(iss >> axis))
  {
    cout << "Error determining rotation axis.\n";
    return 9;
  }

  ostringstream fn2, fn3, fn4; // String streams for easy file naming
  fn2 << filename1.substr(0,filename1.find(".")).c_str()
      << "_" << theta
      << "degree_r" << r_grain << "A_rotated.dat";
      // This next line was used in determining the ideal cutoff distance.
      //<< "degree_r" << r_grain << "A_rotated_rcut" << UU_RNN_CUT << ".dat";
  filename2 = fn2.str();

  fn3 << filename1.substr(0,filename1.find(".")).c_str()
      << "_" << theta
      << "degree_r" << r_grain << "A_marked.dat";
      //<< "degree_r" << r_grain << "A_marked_rcut" << UU_RNN_CUT << ".dat";
  filename3 = fn3.str();

  // We use theta_conv to do the calculations, theta is meant to just look pretty.
  theta_conv = theta * PI / 180.0; // convert theta_conv to radians
  costheta = cos(theta_conv); // just calculate this once!
  sintheta = sin(theta_conv);

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

  // Read and create the header of the dat file
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

  ntotal = 0; // Number of atoms read so far
  n_atom_id = 0; // number written so far
  atoms.resize(N, Atom());

  while (fin >> atom_id >> atom_type >> x >> y >> z) // read the data
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

    // If we are smaller than the radius, rotate the atom by theta_conv around the
    // z axis.
    // TODO: make this an option to do twist or tilt boundaries

    if ((x1 * x1 + y1 * y1) <= (r_grain_sq))
    {
      temp_x = x1 * costheta - y1 * sintheta;
      temp_y = x1 * sintheta + y1 * costheta;
      x1 = temp_x;
      y1 = temp_y;
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
    fout.precision(6);
    fout << x1 << " " << y1  << " " << z1 << endl;

    atoms[atom_id - 1] = Atom(atom_id, atom_type, 0.0, x1, y1, z1); // store all atoms
  }

  // Make sure we read all of the atoms
  if (ntotal != N)
  {
    cout << "ntotal = " << ntotal << " != N = " << N << endl;
    return -4;
  } // End error check
  fin.close(); // Done reading the file

  // Generate the cell-linked list for fast calculations
  // First generate the number of cells in each direction (minimum is 1)
  ncellx = (int)(Lx / RNN_CUT) + 1;
  ncelly = (int)(Ly / RNN_CUT) + 1;
  ncellz = (int)(Lz / RNN_CUT) + 1;
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

                  if (drij_sq > rnn_cut_sq)
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

  /*****************************************************************************
  ****************************ATOM REMOVAL**************************************
  *****************************************************************************/
  // Compare the distances of each atom
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getMark() == 0) // only checking unmarked Cu atoms
    {
      x1 = atoms[i].getX();
      y1 = atoms[i].getY();
      z1 = atoms[i].getZ();
      for (int l = 1; l <= iatom[0][i]; ++l)
      {
        int id = iatom[l][i];
        if (atoms[id].getMark() == 0) // Only unmarked Cu atoms
        {
          // Calculate the distance
          rxij = x1 - atoms[id].getX();
          ryij = y1 - atoms[id].getY();
          rzij = z1 - atoms[id].getZ();

          // Apply PBCs
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;

          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
          if (drij_sq < rnn_cut_sq) // If the Cu atoms are too close
          {
            atoms[i].setMark(1);
            ++n_Cu_removed;
            break; // since we've removed this atom, move on to the next one
          }
        }
      }
    }
  }

  cout << n_Cu_removed << " Cu atoms will be removed.\n";

  ofstream fout2(filename3.c_str());
  fout2 << fixed;
  if (fout2.fail())
  {
    cout << "Error opening the file " << filename3 << endl;
    return -2;
  }

  fn4 << filename1.substr(0,filename1.find("N")).c_str()
      << axis << "_";
  if (decimal) // Input is checked to determine if we print the decimal or not in the filename
    fn4.precision(2);
  else
    fn4.precision(0);
  fn4 << fixed << theta << "degree_r";
  fn4.precision(0);
  fn4 << r_grain
      << "A_removed.dat";
      //<< "A_removed_rcut" << UU_RNN_CUT << ".dat";
  filename4 = fn4.str();

  ofstream fout3(filename4.c_str());
  fout3 << fixed;
  if (fout3.fail()) // error check
  {
    cout << "Error opening the file " << filename4 << endl;
    return -2;
  }

  // write the base data to the file
  fout3 << "These Cu coordinates are shifted and have atoms removed:[ID type x y z]\n"
        << "\n"
        << N - n_Cu_removed << "   atoms\n"
        << ntypes << "   atom types\n"
        << xlow << " " << xhigh << "   xlo xhi\n"
        << ylow << " " << yhigh << "   ylo yhi\n"
        << zlow << " " << zhigh << "   zlo zhi\n"
        << "\nAtoms\n\n";

  // Now write the atoms to the files.  filename3 has all the atoms including
  // the rotated ones and the tag. filename4 has the correct number of atoms.
  ntotal = 0;
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    fout2 << atoms[i].getId() << " " << atoms[i].getType() << " "
          << atoms[i].getX() << " "
          << atoms[i].getY() << " " << atoms[i].getZ() << " "
          << atoms[i].getMark() << endl;

    if (atoms[i].getMark() == 0) // We only write the atoms that are NOT marked for removal
    {
      ++ntotal;
      fout3 << ntotal << " " << atoms[i].getType() << " "
            << atoms[i].getX() << " "
            << atoms[i].getY() << " " << atoms[i].getZ() << endl;
    }
  }

  // Close the file streams
  fout2.close();
  fout3.close();

  return 0;
}
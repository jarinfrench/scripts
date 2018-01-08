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
  string filename1, filename2, input_file, data_file, str; // filenames read from and written to, input file, data file, junk variable
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounds variables
  int N, n_type, n_atoms_read = 0; // number of atoms, atom types, number of atoms read
  vector <Atom> atoms; // all of the atoms from the file
  vector <double> symm; // a vector to hold the calculated symmetry parameters
  int atom_id, type; // id and type number of atom; used to read in the data
  double charge, x, y, z, xu, yu, zu; // charge and position of atom, used to read in data
  double rxij, ryij, rzij, drij_sq; // positions and distance squared
  double r_cut_sq; // cutoff distance
  bool dump; // boolean value to determine if the read file is a LAMMPS dump file or not.
  unsigned int n_grain_1, n_grain_2; // counter for number of atoms in each grain

  // Variables for the symmetry parameters
  double coeffs [2] = {3.0, 2.0}; // Coefficients of the symmetry parameter equation unrotated/rotated symmetry parameter values.
  double xtemp, ytemp, y2, cutoff, costheta_sq; // temp position variables, sin/cos of misorientation, cutoff value
  double total1 = 0.0; // symmetry parameters
  vector <double> xx (12,0.0); // x positions in terms of a0 for nearest neighbors
  vector <double> yy (12,0.0); // y positions in terms of a0 for nearest neighbors
  vector <double> zz (12,0.0); // z positions in terms of a0 for nearest neighbors
  vector <double> new_x_axis (3,0); // New x axis position
  vector <double> new_y_axis (3.0), new_rotated_y (3,0); // New y axis position, new rotated y
  vector <double> new_z_axis (3,0); // New z axis position
  double norm_z; // normalizing value for each axis
  int axis; // Miller indices of rotation axis
  double costheta, sintheta;

  // Input file parameters
  int n_files; // Number of files to be read, rotation axis
  double theta, r_cut, a0; // misorientation angle, cutoff distance (in terms of a0), lattice parameter

  // Variables used for the cell-linked list
  int n_atoms_per_cell; // self-explanatory
  vector <vector <int> > iatom; // Cell-linked list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell
  int ncellx, ncelly, ncellz, idx, idy, idz; // Number of sub cells in each direction, cell number in each direction, atoms in cell i
  double lcellx, lcelly, lcellz; // length of sub cells in each direction

  // Parse command line input. Prompt for what we need.
  if (argc != 2)
  {
    cout << "Please enter the input file to be read: ";
    cin  >> input_file;
  }
  else
  {
    input_file = argv[1];
  }

  // open the input file stream
  ifstream fin_input(input_file.c_str());
  if (fin_input.fail())
  {
    cout << "Error reading input file: " << input_file << endl;
    return 1;
  }

  // Open the output data file stream
  ofstream fout_data("data.txt");
  if (fout_data.fail())
  {
    cout << "Error opening file data.txt\n";
    return 1;
  }
  // The # is there for use with programs like gnuplot.
  fout_data << "# Data consists of: [timestep, atoms in grain 1, atoms in grain 2]\n";

  // Get the important information from the input file:
  // Number of files, misorientation angle, number of atom types, cutoff distance, lattice parameter
  // The way this is written, any values after a0 are ignored (I think)
  getline(fin_input, str);
  stringstream ss_input(str);
  if (!(ss_input >> n_files >> theta >> n_type >> r_cut >> a0))
  {
    cout << "Error reading the input file.  Did you forget a value?\n"
         << "Format of the first line of the input file is:\n"
         << "<number of files> <misorientation angle> <number of atom types> <cutoff distance in Angstroms> <lattice parameter in Angstroms>\n";
    return 9;
  }
  cout << "Input parameters:\n"
       << "\tn_files = " << n_files << endl
       << "\ttheta = " << theta << endl
       << "\tn_types = " << n_type << endl
       << "\tr_cut = " << r_cut << endl
       << "\ta0 = " << a0 << endl;
  r_cut_sq = r_cut * r_cut;

  // Now read the new axis locations
  getline(fin_input,str);
  ss_input.clear();
  ss_input.str(str);
  if (!(ss_input >> new_x_axis[0] >> new_x_axis[1] >> new_x_axis[2]))
  {
    cout << "Error reading new x axis.\n";
    return 9;
  }

  getline(fin_input,str);
  ss_input.clear();
  ss_input.str(str);
  if (!(ss_input >> new_y_axis[0] >> new_y_axis[1] >> new_y_axis[2]))
  {
    cout << "Error reading new y axis.\n";
    return 9;
  }

  getline(fin_input,str);
  ss_input.clear();
  ss_input.str(str);
  if (!(ss_input >> new_z_axis[0] >> new_z_axis[1] >> new_z_axis[2]))
  {
    cout << "Error reading new z axis.\n";
    return 9;
  }
  norm_z = sqrt(new_z_axis[0] * new_z_axis[0] + new_z_axis[1] * new_z_axis[1] + new_z_axis[2] * new_z_axis[2]);

  stringstream ss;
  ss << abs(new_z_axis[0]) << abs(new_z_axis[1]) << abs(new_z_axis[2]);
  ss >> axis;

  switch (axis)
  {
    // Note that these all assume a misorientation of 45 degrees.
    case 1:
    case 10:
    case 100:
      //theta = 0.0;
      costheta = 1.0;
      sintheta = 0.0;

      cutoff = 1.38;
      break;

    case 11:
    case 101:
    case 110:
      costheta = 0;
      sintheta = -1;

      cutoff = 1.45;
      break;

    case 111:
      costheta = sqrt(1.0 / 3.0);
      sintheta = -sqrt(2.0 / 3.0);

      cutoff = 1.34;
      break;

    default:
      cout << "This axis is not implemented explicitly. Using default values.\n";
      double theta_x = acos(new_z_axis[2] / norm_z);
      costheta = cos(theta_x);
      sintheta = -sin(theta_x);

      cutoff = 1.35;
      // Note that these cutoff values seem to work better for during the simulation.
      // For the initial crystal, values of 1.25, 1.2, 1.245, and 1.25 seem to
      // work better for the cutoff values (for 100, 110, 111, and default respectively)
  }

  cout << "\tRotated coordinate system:\n"
       << "\t  x = " << new_x_axis[0] << " " << new_x_axis[1] << " " << new_x_axis[2] << endl
       << "\t  y = " << new_y_axis[0] << " " << new_y_axis[1] << " " << new_y_axis[2] << endl
       << "\t  z = " << new_z_axis[0] << " " << new_z_axis[1] << " " << new_z_axis[2] << endl;

  int aa = 1;
  while (getline(fin_input, filename1) && aa <= n_files)
  {
    // Open up the files for reading and writing.
    ifstream fin(filename1.c_str());
    if (fin.fail())
    {
      cout << "Error opening file " << filename1 << endl;
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
      if (aa == 1)
      {
        if (str != "0")
        {
          cout << "Warning: first data file is not at timestep 0! "
               << "Ignore this warning if this is intentional.\n";
        }
      }
      getline(fin, str); // Gets ITEM: NUMBER OF ATOMS
      fin >> N;
      fin.ignore();
      getline(fin, str); //get ITEM: BOX BOUNDS
      fin >> xlow >> xhigh;
      fin >> ylow >> yhigh;
      fin >> zlow >> zhigh;
      fin.ignore();
      getline(fin, str); // Gets ITEM: ATOMS <data types>
      filename2 = filename1.substr(0,filename1.find(".dump")) + "_interface.dat";
    }
    else
    {
      // This is for a LAMMPS input file
      getline(fin, str); // gets comment line
      fin >> N >> str; // number of atoms
      fin >> str >> str >> str; // Number of types is specified in the input file
      fin >> xlow >> xhigh >> str >> str;
      fin >> ylow >> yhigh >> str >> str;
      fin >> zlow >> zhigh >> str >> str;
      fin >> str;
      fin.ignore();
      getline(fin, str);
      filename2 = filename1.substr(0,filename1.find(".dat")) + "_interface.dat";
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

    ofstream fout(filename2.c_str());
    if (fout.fail())
    {
      cout << "Error opening file " << filename2 << endl;
      return 1;
    }

    fout_data << filename2.substr(0,filename2.find("_interface")) << " ";

    // Preallocate the information (saves time, and allows the atoms to be written
    // in order)
    atoms.resize(N, Atom());
    n_atoms_read = 0;
    n_grain_1 = 0;
    n_grain_2 = 0;
    icell.clear();
    pcell.clear();
    iatom.clear();

    // Read the data
    while (getline(fin, str))
    {
      stringstream ss(str);
      if (n_type == 1)
      {
        if (!(ss >> atom_id >> type >> x >> y >> z >> xu >> yu >> zu))
        {
          xu = 0.0;
          yu = 0.0;
          zu = 0.0;
        }
        charge = 0.0;
      }
      else if (n_type == 2)
      {
        if (!(ss >> atom_id >> type >> charge >> x >> y >> z >> xu >> yu >> zu))
        {
          xu = 0.0;
          yu = 0.0;
          zu = 0.0;
        }
      }
      else
      {
        cout << "This case is not handled.  n_types = " << n_type << " != (1|2)\n";
        return 10;
      }

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
      atoms[atom_id - 1].setXu(xu);
      atoms[atom_id - 1].setYu(yu);
      atoms[atom_id - 1].setZu(zu);
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
    ncellx = (int)(Lx / r_cut) + 1;
    ncelly = (int)(Ly / r_cut) + 1;
    ncellz = (int)(Lz / r_cut) + 1;
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
                    rxij = rxij - anInt(rxij / Lx) * Lx;
                    ryij = ryij - anInt(ryij / Ly) * Ly;
                    rzij = rzij - anInt(rzij / Lz) * Lz;

                    // Now calculate the distance
                    drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                    if (drij_sq > r_cut_sq)
                    {
                      continue; // move to the next atom if we're too far away
                    }

                    if (drij_sq < 1.0E-8) // This should never be hit, but just in case
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
        continue; // Only focus on U (or the single atom type)
      }

      x = atoms[i].getX();
      y = atoms[i].getY();
      z = atoms[i].getZ();
      total1 = 0.0; // reset the total
      // We start at l = 1 because if we start at l = 0, we just re-use the same
      // atom over and over.
      for (int l = 1; l <= iatom[0][i]; ++l)
      {
        unsigned int id = iatom[l][i];
        if (id < i)
        {
          continue;
        }

        // calculate the distances
        rxij = atoms[id].getX() - x;
        ryij = atoms[id].getY() - y;
        rzij = atoms[id].getZ() - z;

        // Apply PBCs
        rxij = rxij - anInt(rxij / Lx) * Lx;
        ryij = ryij - anInt(ryij / Ly) * Ly;
        rzij = rzij - anInt(rzij / Lz) * Lz;

        // Rotation about the x axis
        ytemp = ryij * costheta - rzij * sintheta;

        // Rotation about the new z axis
        xtemp = rxij * costheta - ytemp * sintheta;
        y2 = rxij * sintheta + ytemp * costheta;

        // Calculate the magnitude of the distance, projected onto the xy plane
        drij_sq = (xtemp * xtemp) + (y2 * y2);

        if (drij_sq < 1.0E-8) // Handles the case where the projected position of the atom is right on top of the current atom.
        {
          total1 += 1;
          cout << "Note: drij_sq = 0.0\n";
          continue;
        }
        // cos = dot(A,B) / (|A|*|B|)
        costheta_sq = ((xtemp * xtemp) / drij_sq);
        symm[i] += (coeffs[0] - coeffs[1] * costheta_sq) * (coeffs[0] - coeffs[1] * costheta_sq) * costheta_sq;
        symm[id] += (coeffs[0] - coeffs[1] * costheta_sq) * (coeffs[0] - coeffs[1] * costheta_sq) * costheta_sq;
      }

    }
    for (unsigned int i = 0; i < symm.size(); ++i)
    {
      symm[i] /= iatom[0][i];

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

    fout_data << n_grain_1 << " " << n_grain_2 << endl;

    // Make sure we write the entire set of atoms
    // This writes things in a tecplot-readable format.
    n_atoms_read = 0;
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      if (atoms[i].getType() != 1)
      {
        continue;
      }
      if (n_type == 1)
      {
        if (i == 0)
        {
          fout << "VARIABLES = \"Atom ID\", \"Atom Type\", \"X\", \"Y\", \"Z\", \"Grain Number\", \"Orientation Parameter\",\"Xu\", \"Yu\", \"Zu\"\n";
        }
        fout << atoms[i].getId() << " "
             << atoms[i].getType() << " "
             << (atoms[i].getX() + xlow) * a0 << " "
             << (atoms[i].getY() + ylow) * a0 << " "
             << (atoms[i].getZ() + zlow) * a0 << " "
             << atoms[i].getMark() << " "
             << symm[i] << " " << atoms[i].getXu() << " " << atoms[i].getYu()
             << " " << atoms[i].getZu() << endl;
      }
      else if (n_type == 2)
      {
        if (i == 0)
        {
          fout << "VARIABLES = \"Atom ID\", \"Atom Type\", \"Atom Charge\",\"X\", \"Y\", \"Z\", \"Grain Number\", \"Orientation Parameter\",\"Xu\", \"Yu\", \"Zu\"\n";
        }
        fout << atoms[i].getId() << " "
             << atoms[i].getType() << " "
             << atoms[i].getCharge() << " "
             << (atoms[i].getX() + xlow) * a0 << " "
             << (atoms[i].getY() + ylow) * a0 << " "
             << (atoms[i].getZ() + zlow) * a0 << " "
             << atoms[i].getMark() << " "
             << symm[i] << " " << atoms[i].getXu() << " " << atoms[i].getYu()
             << " " << atoms[i].getZu() << endl;
      }
      else
      {
        cout << "Error: n_types != (1|2), n_types = " << n_type << endl;
        return 10;
      }

      ++n_atoms_read;
    }
    fout.close(); // Close the output file

    // Allows for a different number of atom types to be checked.
    if (n_type == 2)
    {
      n_atoms_read *= 3;
    }
    if (n_atoms_read != N)
    {
      cout << "Error: number of atoms written does not match number of atoms in the simulation.\n"
           << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
      return 6;
    }

    cout << "Processing of file \"" << filename1 << "\" completed.\n";
    ++aa;
  }

  fin_input.close();
  fout_data.close();
  return 0;
}

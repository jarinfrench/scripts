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
#include <map>
#include <algorithm>
#include "atom.h"

using namespace std;

#define PI 3.14159265358979 // easier and faster to simply store these values here.
#define UU_RNN_CUT 2.0 // Cutoff value for U-U atoms too close (Basak Potential)
#define UO_RNN_CUT 4.0 // Cutoff value for U-O atoms too close
#define OO_RNN_CUT 0.63801 // Cutoff value for O-O atoms too close
#define CU_RNN_CUT 1.0 // Cutoff value for Cu-Cu atoms too close (using Mishin potential)
#define AL_RNN_CUT 2.4 // Cutoff value for Al-Al atoms too close (using Ercolessi-Adams potential)
#define AG_RNN_CUT 1.0 // Cutoff value for Ag-Ag atoms too close (using Mishin potential)
#define AU_RNN_CUT 1.0 // Cutoff value for Au-Au atoms too close (using Foiles Au1 potential)
#define NI_RNN_CUT 1.0 // Cutoff value for Ni-Ni atoms too close (using Foiles-Hoyt potential)

struct ratio {
  int n_types;
  vector <int> ratio;
} compound_ratio;

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

template <typename K, typename V>
int findByValue(map<K,V> mapOfElemen, V value)
{
  typename map<K,V>::iterator it = mapOfElemen.begin();
  while (it != mapOfElemen.end())
  {
    if (it->second == value)
    {
      return it->first;
    }
    ++it;
  }
  return -1;
}

int main(int argc, char **argv)
{
  // External values
  string filename1, filename2, filename3, filename4, input_file, str; //filenames and line variable
  string element, chem_formula; // element type
  string boundary_type;
  bool is_sphere;
  int axis;
  double r_grain; //radius of the grain, with buffer zone
  double r_grain_sq; // squared values for convenience
  double theta, theta_conv; // angle of rotation
  double costheta, sintheta; // better to calculate this once.
  double uo_rnn_cut_sq = UO_RNN_CUT * UO_RNN_CUT;
  double oo_rnn_cut_sq = OO_RNN_CUT * OO_RNN_CUT;

  double scale_factor_a, scale_factor_b, scale_factor_c; // dimensional scaling values
  double Lx, Ly, Lz; // box size
  int ntotal, n_atom_id; // total number of atoms that have been read/written
  double rxij, ryij, rzij, drij_sq; //positional differences, total distance^2
  int n_1_removed = 0, n_2_removed = 0;

  // Values from file
  int N, ntypes, ntypes_test; // number of atoms, and number of atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh; // atom bounds
  int atom_id, atom_type; // id number and type number
  double atom_charge, x, y, z; // charge and position values
  double x1, y1, z1, temp_x, temp_y, x2, y2, z2; // Store the original value and manipulate!

  // Containers
  vector <double> rcut, rcut_sq; // cutoff radii for each interaction
  vector <int> n_removed; // number of atoms removed of each type
  vector <Atom> atoms; // contains the atoms we look at, and the entire set.
  vector <pair<int, double> > distances; // vector of id and distance.
  map <int, string> elements;

  // Variables used for the cell-linked list
  int n_atoms_per_cell; // self-explanatory
  vector <vector <int> > iatom; // Cell-linked list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell
  int ncellx, ncelly, ncellz, idx, idy, idz; // Number of sub cells in each direction, cell number in each direction
  double lcellx, lcelly, lcellz; // length of sub cells in each direction

  if (argc != 2) // check for input file
  {
    cout << "Please enter a valid input file containing (in order):\n"
         << "\tthe data file to be processed\n\tdesired grain radius\n\t"
         << "rotation angle\n\tboundary type (cylinder|sphere)\n\tnumber of atom types\n\tcutoff radius\n"
         << "Note that the number of cutoff radii is dependent on the number of atom interactions,\n"
         << "so, for example, if there are two atom types, there are the two same-element interactions,\n"
         << "and there is also the interaction between the two different elements.  Include all the relevant,\n"
         << "cutoff radii in the order of same-element interactions, followed by inter-element interactions.\n"
         << "Keep the numbering consistent with the data file.\n";
    return 1;
  }
  else
  {
    // Process the input file
    input_file = argv[1];
    ifstream finput(input_file.c_str());
    if (finput.fail())
    {
      cout << "Error reading file \"" << input_file << "\"\n";
      return 2;
    }
    else
    {
      finput >> filename1 >> r_grain >> theta >> boundary_type >> ntypes;
      if ((boundary_type.compare("cylinder") != 0) && (boundary_type.compare("sphere") != 0))
      {
        cout << "Error reading boundary type. Boundary type \"" <<  boundary_type << "\" not recognized.\n";
        return 3;
      }
      n_removed.resize(ntypes,0);
      compound_ratio.n_types = ntypes;
      compound_ratio.ratio.resize(compound_ratio.n_types,1);
      unsigned int combinations = ((ntypes + 1) * ntypes) / 2;
      double temp;
      while (finput >> temp)
      {
        rcut.push_back(temp);
        rcut_sq.push_back(temp * temp);
      }
      if (rcut.size() != combinations)
      {
        cout << "There are " << combinations << " total interaction radii needed.  You provided " << rcut.size() << ".\n";
        return 3;
      }

      (boundary_type.compare("sphere") == 0) ? is_sphere = true : is_sphere = false;

      cout << "Input information:"
           << "\n\tData file: " << filename1
           << "\n\tGrain radius: " << r_grain
           << "\n\tMisorientation angle: " << theta
           << "\n\tBoundary type: " << boundary_type
           << "\n\tNumber of atom types: " << ntypes << "\n\tCutoff radii: ";
      for (vector <double>::iterator it = rcut.begin(); it != rcut.end();)
      {
        cout << *it;
        if (++it != rcut.end())
        {
          cout << "; ";
        }
        else
        {
          cout << endl;
        }
      }

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
    return 4;
  }
  size_t element_pos = filename1.find("LAMMPS_") + 7;
  size_t element_pos_end;
  {
    size_t temp = filename1.find("_N");
    temp = filename1.find("_N", temp + 1);
    if (temp == string::npos)
    {
      element_pos_end = filename1.find("_N");
    }
    else
    {
      element_pos_end = temp;
    }
  }
  if (element_pos == string::npos)
  {
    cout << "Error determining element(s).\n";
    return 4;
  }
  else
  {
    int elem_num = 1;
    for (unsigned int i = element_pos; i < element_pos_end; ++i)
    {
      if (islower(filename1[i])) // if the current character is lower-case
      {
        if (isdigit(filename1[i + 1]))
        {
          compound_ratio.ratio[elem_num - 1] = filename1[i + 1] - '0';
        }
        int test = findByValue(elements, filename1.substr(i-1, 2));
        if (test < 0)
        {
          if ((elements.insert(pair<int,string>(elem_num,filename1.substr(i-1, 2)))).second) // check to see if insertion was successful
          {
            ++elem_num;
          }
          else
          {
            cout << "Error inserting element " << filename1.substr(i-1, 2);
          }
        }
        else
        {
          ++compound_ratio.ratio[test - 1];
        }
      }
      else
      {
        if (islower(filename1[i+1])) // if the next character is lower-case
        {
          continue;
        }
        else
        {
          if (isdigit(filename1[i]))
          {
            continue;
          }
          else
          {
            if (isdigit(filename1[i + 1]))
            {
              compound_ratio.ratio[elem_num - 1] = filename1[i + 1] - '0';
            }
            int test = findByValue(elements, filename1.substr(i, 1));
            if (test < 0)
            {
              if ((elements.insert(pair<int,string>(elem_num,filename1.substr(i, 1)))).second) // check to see if insertion was successful
              {
                ++elem_num;
              }
              else
              {
                cout << "Error inserting element " << filename1.substr(i, 1);
              }
            }
            else
            {
              ++compound_ratio.ratio[test - 1];
            }
          }
        }
      }
    }
  }

  cout << endl;
  chem_formula = filename1.substr(element_pos, element_pos_end - element_pos);
  if (ntypes != elements.size())
  {
    cout << "Error determining element types. Input file ntypes = " << ntypes << " != found elements = " << elements.size() << "\n";
    return 4;
  }
  else
  {
    cout << "Elements found (" << ntypes << "): ";
    for (map<int, string>::iterator it = elements.begin(); it != elements.end();)
    {
      cout << (*it).second;
      if (++it != elements.end())
      {
        cout << ", ";
      }
      else
      {
        cout << endl;
      }
    }
  }
  cout << "The ratio of elements is calculated as: ";
  map <int, string>::iterator elem_it = elements.begin();
  for (vector <int>::iterator it = compound_ratio.ratio.begin(); it != compound_ratio.ratio.end();)
  {

    cout << *it << " " << (*elem_it).second;
    ++elem_it;
    if (++it != compound_ratio.ratio.end())
    {
      cout << ":";
    }
    else
    {
      cout << endl;
    }
  }

  // This may be able to be removed after this update.
  if (ntypes != 1 && ntypes != 2)
  {
    cout << "This script can currently only handle 1 or 2 types of atoms.  You have asked for " << ntypes << " types of atoms.\n";
    return 8;
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
    return 1;
  } // End error check

  ofstream fout(filename2.c_str()); // only writing to this file
  fout << fixed; // makes sure to always use the precision specified.
  if (fout.fail()) //Error check
  {
    cout << "Error opening the file " << filename2 << endl;
    return 1;
  } // End error check

  // Read and create the header of the dat file
  getline(fin, str);

  fout << "These " << chem_formula << " coordinates are shifted [ID type charge* x y z] (*not included if always charge neutral)\n\n";

  //Get the number of atoms
  fin  >> N >> str;
  fout << N << "  atoms\n";

  //Get the number of atom types
  fin  >> ntypes_test >> str >> str;
  if (ntypes_test != ntypes)
  {
    cout << "Atom types incorrectly determined.\n";
    return 8;
  }
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
    return 2;
  }

  if (r_grain * 2.0 >= Ly)
  {
    cout << "Error! Grain diameter = " << r_grain * 2.0 << " >= Ly = " << Ly << endl;
    return 2;
  }

  // For spherical grain!
  if (is_sphere)
  {
    if (r_grain * 2.0 >= Lz)
    {
      cout << "Error! Grain diameter = " << r_grain * 2.0 << " >= Lz = " << Lz << endl;
      return 2;
    }
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

  fin.ignore();
  getline(fin,str); // gets the blank line before the data.
  while (getline(fin, str)) // read the data
  {
    stringstream ss(str);
    stringstream::pos_type pos = ss.tellg(); // Store the beginning position
    if (!(ss >> atom_id >> atom_type >> atom_charge >> x >> y >> z))
    {
      ss.clear();
      ss.seekg(pos,ss.beg);
      if (!(ss >> atom_id >> atom_type >> x >> y >> z))
      {
        cout << "Read error.\n";
        break;
      }
      atom_charge = 0.0;
    }

    ++ntotal; // Increment the number of atoms

    // Make sure there aren't more than ntypes of atoms
    if (atom_type > ntypes)
    {
      cout << "Error! Atom_type = " << atom_type << " is greater than " << ntypes << endl;
      return 3;
    }
    // change the origin to the center of the simulation for rotating the atoms
    x1 = x - Lx / 2.0;
    y1 = y - Ly / 2.0;
    z1 = z - Lz / 2.0;

    // If we are smaller than the radius, rotate the atom by theta_conv around the
    // z axis.
    // TODO: make this an option to do twist or tilt boundaries

    if (is_sphere)
    {
      if ((x1 * x1 + y1 * y1 + z1 * z1) <= (r_grain_sq))
      {
        temp_x = x1 * costheta - y1 * sintheta;
        temp_y = x1 * sintheta + y1 * costheta;
        x1 = temp_x;
        y1 = temp_y;
      }
    }
    else
    {
      if ((x1 * x1 + y1 * y1) <= (r_grain_sq))
      {
        temp_x = x1 * costheta - y1 * sintheta;
        temp_y = x1 * sintheta + y1 * costheta;
        x1 = temp_x;
        y1 = temp_y;
      }
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
    if (ntypes == 2) // Only configurations with 2 atoms have a charge
    {
      fout.precision(1);
      fout << atom_charge << " ";
    }
    fout.precision(6);
    fout << x1 << " " << y1  << " " << z1 << endl;

    atoms[atom_id - 1] = Atom(atom_id, atom_type, atom_charge, x1, y1, z1); // store all atoms
  }

  // Make sure we read all of the atoms
  if (ntotal != N)
  {
    cout << "ntotal = " << ntotal << " != N = " << N << endl;
    return 4;
  } // End error check
  fin.close(); // Done reading the file

  // Generate the cell-linked list for fast calculations
  // First generate the number of cells in each direction (minimum is 1)
  ncellx = (int)(Lx / *max_element(rcut.begin(), rcut.end())) + 1;
  ncelly = (int)(Ly / *max_element(rcut.begin(), rcut.end())) + 1;
  ncellz = (int)(Lz / *max_element(rcut.begin(), rcut.end())) + 1;

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

  /*****************************************************************************
  ****************************ATOM REMOVAL**************************************
  *****************************************************************************/
  // Compare the distances of each atom
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    int cur_type = atoms[i].getType();
    for (int j = 1; j <= ntypes; ++j)
    { // First check atoms of the same type
      if (atoms[i].getType() == j && atoms[i].getMark() == 0) // unmarked type j atoms
      {
        x1 = atoms[i].getX();
        y1 = atoms[i].getY();
        z1 = atoms[i].getZ();
        for (int l = 1; l <= iatom[0][i]; ++l)
        {
          int id = iatom[l][i];
          if (atoms[id].getType() == j && atoms[i].getMark() == 0) // unmarked type j atoms
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
            if (drij_sq < rcut_sq[j - 1])
            {
              atoms[i].setMark(1);
              ++n_removed[j - 1];
              break;
            }
          }
        }
      }

      if (ntypes == 1)
      {
        break; // we're done for the single element case.
      }

      if (atoms[i].getType() == 1 && atoms[i].getMark() == 1)
      {
        distances.clear(); // Clear out the old values.
        // Check each neighboring O atoms, and find the closest two and remove them
        for (int l = 1; l <= iatom[0][i]; ++l)
        {
          int id = iatom[l][i];
          if (atoms[id].getType() == 2 && atoms[id].getMark() == 0)
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
          atom_id = distances[k].first; // we use this a lot over the next lines
          // If the atom we are looking at is O and unmarked
          if (atoms[atom_id].getType() == 2 &&
          atoms[atom_id].getMark() == 0)
          {
            // mark this atom in BOTH lists
            atoms[atom_id].setMark(1);
            ++n_2_removed; // increase the counter for O removed

            // if we have removed enough O atoms to maintain charge neutrality,
            // exit the loop.  We can only remove 2 O atoms per U atom!
            if (n_2_removed == 2 * n_1_removed)
            break;
          }
        }
      }
    }

    // Now go through the list again and remove the O atoms that are too close
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      if (ntypes != 2)
      break;

      if (atoms[i].getType() == 2 && atoms[i].getMark() == 0)
      {
        // Get the position of atom i
        x1 = atoms[i].getX();
        y1 = atoms[i].getY();
        z1 = atoms[i].getZ();

        for (int l = 1; l < iatom[0][i]; ++l)
        {
          int id = iatom[l][i];

          if (atoms[id].getType() == 2 && atoms[id].getMark() == 0) // Unmarked O atoms
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
            if (drij_sq < oo_rnn_cut_sq)
            {
              atoms[i].setMark(1);
              atoms[id].setMark(1);
              n_2_removed += 2;

              // Now go through and find the closest U atom to these two
              x2 = atoms[id].getX();
              y2 = atoms[id].getY();
              z2 = atoms[id].getZ();

              for (int m = 1; m < iatom[0][i]; ++m)
              {
                int jd = iatom[m][i];
                if (atoms[jd].getType() == 1 && atoms[jd].getMark() == 0) // unmarked U atoms
                {
                  //calculate the distances between the U atom and both O atoms
                  rxij = x1 - atoms[jd].getX();
                  ryij = y1 - atoms[jd].getY();
                  rzij = z1 - atoms[jd].getZ();

                  // Apply PBCs
                  rxij = rxij - anInt(rxij / Lx) * Lx;
                  ryij = ryij - anInt(ryij / Ly) * Ly;
                  rzij = rzij - anInt(rzij / Lz) * Lz;

                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
                  if (drij_sq < uo_rnn_cut_sq) // This may not work because there may not be a U atom that is within the cutoff distance for both O atoms - needs to be rewritten
                  {
                    rxij = x2 - atoms[jd].getX();
                    ryij = y2 - atoms[jd].getY();
                    rzij = z2 - atoms[jd].getZ();

                    // Apply PBCs
                    rxij = rxij - anInt(rxij / Lx) * Lx;
                    ryij = ryij - anInt(ryij / Ly) * Ly;
                    rzij = rzij - anInt(rzij / Lz) * Lz;

                    drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
                    if (drij_sq < uo_rnn_cut_sq)
                    {
                      // Only if both O atoms are close enough to the U atom is
                      // this U atom removed.
                      atoms[jd].setMark(1);
                      ++n_1_removed;
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
  }
  if (ntypes == 2)
  {
    cout << n_1_removed << " U atoms will be removed.\n";
    cout << n_2_removed << " O atoms will be removed.\n";
    // Error checking
    if (n_1_removed * 2 != n_2_removed)
    {
      cout << "Error: the removed U:O ratio must be 1:2!\n";
      return 5;
    }
  }
  else // ntypes == 1
  {
    cout << n_1_removed << " " << chem_formula << " atoms will be removed.\n";
  }

  ofstream fout2(filename3.c_str());
  fout2 << fixed;
  if (fout2.fail())
  {
    cout << "Error opening the file " << filename3 << endl;
    return 2;
  }

  fn4 << filename1.substr(0,filename1.find("N")).c_str()
      << axis << "_";
  fn4.precision(2);
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
    return 2;
  }

  // write the base data to the file
  fout3 << "These " << element << " coordinates are shifted and have atoms removed:[ID type x y z]\n"
        << "\n"
        << N - n_1_removed - n_2_removed << "   atoms\n"
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
    if (ntypes == 2)
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
    else // ntypes == 1
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

  }

  if ((ntotal != N - n_1_removed - n_2_removed) && ntypes == 2) // One last check
  {
    cout << "Error! The final number of removed atoms is not balanced!\n"
         << "ntotal = " << ntotal << " != N - n_1_removed - n_2_removed = "
         << N - n_1_removed - n_2_removed << endl;
    return 6;
  }

  // Close the file streams
  fout2.close();
  fout3.close();

  return 0;
}

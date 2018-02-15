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
#include <numeric>
#include "atom.h"

using namespace std;

#define PI 3.14159265358979 // easier and faster to simply store these values here.
/*#define UU_RNN_CUT 2.0 // Cutoff value for U-U atoms too close (Basak Potential)
#define UO_RNN_CUT 4.0 // Cutoff value for U-O atoms too close
#define OO_RNN_CUT 0.63801 // Cutoff value for O-O atoms too close
#define CU_RNN_CUT 1.0 // Cutoff value for Cu-Cu atoms too close (using Mishin potential)
#define AL_RNN_CUT 2.4 // Cutoff value for Al-Al atoms too close (using Ercolessi-Adams potential)
#define AG_RNN_CUT 1.0 // Cutoff value for Ag-Ag atoms too close (using Mishin potential)
#define AU_RNN_CUT 1.0 // Cutoff value for Au-Au atoms too close (using Foiles Au1 potential)
#define NI_RNN_CUT 1.0 // Cutoff value for Ni-Ni atoms too close (using Foiles-Hoyt potential)
*/

struct ratio {
  int n_types;
  vector <int> ratio;
} compound_ratio;

struct neighbor_data {
  pair<int, int> ids;
  pair<int, int> types;
  double distance;
};

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

bool mapValueCmp(pair<pair<int,int>,double> a, pair<pair<int,int>,double> b)
{
  return (a.second < b.second);
}

bool compareNeighbors(const neighbor_data &a, const neighbor_data &b)
{
  if (a.types.first  == b.types.first)
  {
    if (a.types.second == b.types.second)
    {
      return a.distance < b.distance;
    }
    else
    {
      return a.types.second < b.types.second;
    }
  }
  else
  {
    return a.types.first < b.types.first;
  }
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

  double scale_factor_a, scale_factor_b, scale_factor_c; // dimensional scaling values
  double Lx, Ly, Lz; // box size
  int ntotal, n_atom_id; // total number of atoms that have been read/written
  double rxij, ryij, rzij, drij_sq; //positional differences, total distance^2

  // Values from file
  int N, ntypes, ntypes_test; // number of atoms, and number of atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh; // atom bounds
  int atom_id, atom_type; // id number and type number
  double atom_charge, x, y, z; // charge and position values
  double x1, y1, z1, temp_x, temp_y; // Store the original value and manipulate!

  // Containers
  vector <int> n_removed; // number of atoms removed of each type
  vector <Atom> atoms; // contains the atoms we look at, and the entire set.
  vector <neighbor_data> neighbor_dataset;
  map <int, string> elements;
  map <pair<int,int>, double> rcut, rcut_sq; // cutoff radii for each interaction

  // Variables used for the cell-linked list
  int n_atoms_per_cell; // self-explanatory
  vector <vector <int> > iatom; // Cell-linked list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell
  int ncellx, ncelly, ncellz, idx, idy, idz; // Number of sub cells in each direction, cell number in each direction
  double lcellx, lcelly, lcellz; // length of sub cells in each direction

  if (argc < 2) // check for input file
  {
    cout << "Please enter a valid input file containing (in order):\n"
         << "\tthe data file to be processed\n\tdesired grain radius\n\t"
         << "boundary type (cylinder|sphere)\n\tnumber of atom types\n\tcutoff radius\n"
         << "Note that the number of cutoff radii is dependent on the\n"
         << "number of atom interactions, so, for example, if there are\n"
         << "two atom types, there are the two same-element interactions,\n"
         << "and there is also the interaction between the two different\n"
         << "elements.  Include all the relevant, cutoff radii one element\n"
         << "at a time, i.e. the cutoff radii for a 3 element system would\n"
         << "be input as 1-1, 1-2, 1-3, 2-2, 2-3, 3-3. Maintain consistent\n"
         << "numbering with the data file.\n";
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
      finput >> filename1 >> r_grain >> boundary_type >> ntypes;
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
      for (int i = 1; i <= ntypes; ++i) // We start at 1 to be consistent with the atom type numbering
      {
        for (int j = i; j <= ntypes; ++j)
        {
          finput >> temp;
          rcut.insert(make_pair(make_pair(i,j),temp));
          rcut_sq.insert(make_pair(make_pair(i,j), temp * temp));
        }
      }
      if (rcut.size() != combinations)
      {
        cout << "There are " << combinations << " total interaction radii needed.  You provided " << rcut.size() << ".\n";
        return 3;
      }

      (boundary_type.compare("sphere") == 0) ? is_sphere = true : is_sphere = false;

      if (argc == 3)
      {
        theta = strtod(argv[2], NULL);
      }
      else
      {
        cout << "Please enter the misorientation angle between the two grains: ";
        cin  >> theta;
      }

      cout << "Input information:"
      << "\n\tData file: " << filename1
      << "\n\tGrain radius: " << r_grain
      << "\n\tMisorientation angle: " << theta
      << "\n\tBoundary type: " << boundary_type
      << "\n\tNumber of atom types: " << ntypes << "\n\tCutoff radii: ";
      for (map <pair<int,int>, double>::iterator it = rcut.begin(); it != rcut.end();)
      {
        cout << (*it).second;
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

  if ((*max_element(rcut.begin(), rcut.end(), mapValueCmp)).second < 1)
  {
    cout << "Error!  Not enough memory available!\n";
    return 15;
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
      cout << " : ";
    }
    else
    {
      cout << endl;
    }
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
  // We want the cells to be roughly the same size as a unit cell, so we use the largest
  // cutoff value given as a basis for that.
  ncellx = (int)(Lx / (*max_element(rcut.begin(), rcut.end(), mapValueCmp)).second) + 1;
  ncelly = (int)(Ly / (*max_element(rcut.begin(), rcut.end(), mapValueCmp)).second) + 1;
  ncellz = (int)(Lz / (*max_element(rcut.begin(), rcut.end(), mapValueCmp)).second) + 1;

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

                  if (drij_sq > (max_element(rcut_sq.begin(),rcut_sq.end(), mapValueCmp)->second))
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

  /*****************************************************************************
  ****************************ATOM REMOVAL**************************************
  *****************************************************************************/
  // Compare the distances of each atom
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    int j = atoms[i].getType();
    if (atoms[i].getMark() == 0)
    {
      // Store the x, y, and z positions of the atom for easy comparison
      x1 = atoms[i].getX();
      y1 = atoms[i].getY();
      z1 = atoms[i].getZ();

      // Now check it's nearest neighbors
      neighbor_dataset.clear();
      bool removed = false;
      for (int l = 1; l <= iatom[0][i]; ++l)
      {
        int id = iatom[l][i];
        if (atoms[id].getMark() == 0) // If the atom is not marked, check it
        {
          int k = atoms[id].getType();
          // calculate the distance between the two
          rxij = x1 - atoms[id].getX();
          ryij = y1 - atoms[id].getY();
          rzij = z1 - atoms[id].getZ();

          // Apply PBCs
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;

          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);
          if (ntypes != 1)
          { // This is only needed if there is more than 1 type of atom.
            neighbor_data tmp;
            tmp.ids = make_pair(i,id);
            tmp.types = make_pair(j,k);
            tmp.distance = drij_sq;
            neighbor_dataset.push_back(tmp);
          }

          if (atoms[id].getType() == j && !removed) // We are really only interested in the atoms that are too close to the same type of atom.
          {
            pair <int, int> key =  (j < k) ? make_pair(j,k) : make_pair(k,j); // make sure the lower index comes first
            if (drij_sq < rcut_sq[key])
            {
              atoms[i].setMark(1); // if we are below the specified cutoff, mark for removal
              ++n_removed[j - 1]; // count the number removed of this type
              removed = true;
            }
          }
        }
      }
      if (removed && ntypes != 1)
      {
        sort(neighbor_dataset.begin(), neighbor_dataset.end(), compareNeighbors);
        // Now we need to remove atoms to maintain the original ratio
        for (int it = 1; it <= compound_ratio.ratio.size(); ++it)
        { // for each ratio
          int ratio = compound_ratio.ratio[it - 1];
          if (it == j) // if this is the atom we already removed, we don't want to remove an extra atom
          {
            ratio -= 1;
          }
          for (int ii = 0; ii < ratio; ++ii) // remove atoms until this ratio is met
          {
            for (vector <neighbor_data>::iterator neigh_it = neighbor_dataset.begin();
                 neigh_it != neighbor_dataset.end(); ++neigh_it)
            { // check each neighbor's type.  If it matches the type we are looking for, remove it, and break from the loop
              if ((*neigh_it).types.second == it && atoms[(*neigh_it).ids.second - 1].getMark() == 0)
              {
                atoms[(*neigh_it).ids.second - 1].setMark(1);
                ++n_removed[it - 1];
                break;
              }
            }
          }
        }
      }
    }
  }

  for (unsigned int i = 0; i < n_removed.size(); ++i)
  {
    for (unsigned int j = 0; j < compound_ratio.ratio.size(); ++j)
    {
      if (compound_ratio.ratio[j] * n_removed[i] != compound_ratio.ratio[i] * n_removed[j])
      {
        cout << "Error: the element ratio has not been kept.\n";
        cout << compound_ratio.ratio[j] << "*" << n_removed[i] << "!=" << compound_ratio.ratio[i] << "*" << n_removed[j] << endl;
        return 5;
      }
    }
    cout << n_removed[i] << " " << elements[i+1] << " atoms will be removed.\n";
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
  int subtotal = accumulate(n_removed.begin(), n_removed.end(),0);
  fout3 << "These " << element << " coordinates are shifted and have atoms removed:[ID type x y z]\n"
        << "\n"
        << N - subtotal << "   atoms\n"
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

  if ((ntotal != N - subtotal)) // One last check
  {
    cout << "Error! The final number of removed atoms is not balanced!\n"
         << "ntotal = " << ntotal << " != N - subtotal = "
         << N - subtotal << endl;
    return 6;
  }

  // Close the file streams
  fout2.close();
  fout3.close();

  return 0;
}

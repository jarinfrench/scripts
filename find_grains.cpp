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

void extractAxis(int axis, vector <double> & v)
{
  if (axis > 999)
  {
    cout << "Axis indices are too high!  Must be <= 999.\n";
    return;
  }

  v[0] = (axis / 100) % 10;
  v[1] = (axis / 10) % 10;
  v[2] = axis % 10;
}

/* Note that this function does NOT correctly rotate the axes as given in GBstudio
 * however, it does allow the resulting simulations to be properly parsed.
 **/
vector <vector <double> > rotate2Axis(vector <double> axis1, vector <double> axis2)
{
  if (axis1.size() != 3 || axis2.size() != 3)
  {
    cout << "Error: size of vectors is incorrect.  Must be a 3D vector (3 elements)\n"
         << "\tSize of vector 1: " << axis1.size()
         << "\n\tSize of vector 2:" << axis2.size() << endl;
    exit(12);
  }

  vector <double> v (3,0); // cross product of axis 1 and axis 2
  vector <vector <double> > vx; // skew symmetric matrix given by the normalized axes
  vector <vector <double> > vx_sq; // dot product of vx with itself
  vector <vector <double> > rotation; // The resultant rotation matrix
  double c; // dot product of axis1 with axis2
  double norm1 = sqrt(axis1[0] * axis1[0] + axis1[1] * axis1[1] + axis1[2] * axis1[2]);
  double norm2 = sqrt(axis2[0] * axis2[0] + axis2[1] * axis2[1] + axis2[2] * axis2[2]);

  vx.resize(3, vector <double> (3, 0.0));
  vx_sq.resize(3, vector <double> (3, 0.0));
  rotation.resize(3, vector <double> (3, 0.0));

  if (axis1[0] == axis2[0] && axis1[1] == axis2[1] && axis1[2] == axis2[2])
  {
    v[0] = 0; v[1] = 0; v[2] = 0;
  }
  else
  {
    v[0] = axis1[1] / norm1 * axis2[2] / norm2 - axis1[2] / norm1 * axis2[1] / norm2;
    v[1] = -(axis1[0] / norm1 * axis2[2] / norm2 - axis1[2] / norm1 * axis2[0] / norm2);
    v[2] = axis1[0] / norm1 * axis2[1] / norm2 - axis1[1] / norm1 * axis2[0] / norm2;
  }

  c = axis1[0] / norm1 * axis2[0] / norm2 + axis1[1] / norm1 * axis2[1] / norm2 + axis1[2] / norm1 * axis2[2] / norm2;
  vx[0][1] = -v[2]; vx[0][2] = v[1];
  vx[1][0] = v[2]; vx[1][2] = -v[0];
  vx[2][0] = -v[1]; vx[2][1] = v[0];

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        vx_sq[i][j] += vx[i][k] * vx[k][j];
      }
    }
  }

  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      if (i == j)
      {
        rotation[i][j] += 1.0;
      }

      rotation[i][j] += vx[i][j] + vx_sq[i][j] / (1 + c);
    }
  }

  return rotation;
}

int main(int argc, char** argv)
{
  string filename1, filename2, input_file, data_file, str; // filenames read from and written to, input file, data file, junk variable
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounds variables
  int N, n_type, n_atoms_read = 0, smaller; // number of atoms, atom types, number of atoms read
  vector <Atom> atoms; // all of the atoms from the file
  vector <double> symm; // a vector to hold the calculated symmetry parameters
  int atom_id, type; // id and type number of atom; used to read in the data
  double charge, x, y, z; // charge and position of atom, used to read in data
  double rxij, ryij, rzij, drij_sq; // positions and distance squared
  double r_cut_sq; // cutoff distance
  double sintheta_sq, total1 = 0.0, total2 = 0.0; // sin^2 of the angle, symmetry parameter (it's a sum, so starts at 0)
  double xtemp, ytemp, sintheta, costheta, cutoff; // rotated x position, y position, sin theta, cos theta, cutoff for which grain an atom belongs to.
  bool dump; // boolean value to determine if the read file is a LAMMPS dump file or not.
  unsigned int n_grain_1, n_grain_2; // counter for number of atoms in each grain
  double coeffs [2] = {3,2}; // Coefficients of the symmetry parameter (default)
  vector <vector <double> > rotation; //rotation matrix

  // Input file parameters
  int n_files, rot_axis; // Number of files to be read, rotation axis
  double theta, r_cut, a0, ideal_symm; // misorientation angle, cutoff distance (in terms of a0), a0, ideal symmetry parameter.
  /* Note that the ideal symmetry parameter is calculated by taking the orientation
  * of the larger grain (or the outside grain) and calculating the orientation
  * parameter as defined by Bai et al.
  */

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
    cout << "Error reading input file " << input_file << endl;
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
  fout_data << "# Data consists of: [timestep, atoms in grain 2]\n";

  // Get the important information from the input file:
  // Number of files, misorientation angle, number of atom types, cutoff distance, lattice parameter, ideal symmetry parameter
  getline(fin_input, str);
  stringstream ss_input(str);
  if (!(ss_input >> n_files >> theta >> n_type >> r_cut >> a0 >> rot_axis ))
  {
    cout << "Error reading the input file.  Did you forget a value?\n"
         << "Format of the first line of the input file is:\n"
         << "<number of files> <misorientation angle> <number of atom types> <cutoff distance in Angstroms> <lattice parameter in Angstroms> <rotation axis>\n";
    return 9;
  }
  cout << "Input parameters:\n"
       << "\tn_files = " << n_files << endl
       << "\ttheta = " << theta << endl
       << "\tn_types = " << n_type << endl
       << "\tr_cut = " << r_cut << endl
       << "\ta0 = " << a0 << endl
       << "\trotation_axis = " << rot_axis << endl;
  r_cut_sq = r_cut * r_cut;

  sintheta = sin(theta * PI / 180.0); // best to calculate this once
  costheta = cos(theta * PI / 180.0);

  // Calculate the cutoff value for grain identification
  // Assumes FCC structure!
  vector <double> xx (12,0.0); // x positions in terms of a0 for nearest neighbors
  vector <double> yy (12,0.0); // y positions in terms of a0 for nearest neighbors
  vector <double> zz (12,0.0); // z positions in terms of a0 for nearest neighbors

  vector <double> rotation_axis (3,0);
  vector <double> z_axis (3,0);
  z_axis[2] = 1.0;
  // Ideally, a simple function could calculate the positions of the first nearest
  // neighbors just based on the rotation axis, thus eliminating the need for
  // this switch statement.  Note that GBStudio (the applet that I generally use
  // to create my structures) aligns the rotation axis with the z axis, and the
  // function rotate2Axis above seems to rotate things to align with the x axis
  // Perhaps I need an additional rotation about the Y axis to align x with z?
  /*switch (rot_axis)
  {
    case 100 :
      //ideal_symm = 1.0;
      xx[0]  =  0.0000; yy[0]  = -0.5000; zz[0]  = -0.5000; // (   0, -1/2, -1/2)
      xx[1]  =  0.0000; yy[1]  = -0.5000; zz[1]  =  0.5000; // (   0, -1/2,  1/2)
      xx[2]  =  0.0000; yy[2]  =  0.5000; zz[2]  = -0.5000; // (   0,  1/2, -1/2)
      xx[3]  =  0.0000; yy[3]  =  0.5000; zz[3]  =  0.5000; // (   0,  1/2,  1/2)
      xx[4]  = -0.5000; yy[4]  =  0.0000; zz[4]  = -0.5000; // (-1/2,    0, -1/2)
      xx[5]  = -0.5000; yy[5]  =  0.0000; zz[5]  =  0.5000; // (-1/2,    0,  1/2)
      xx[6]  =  0.5000; yy[6]  =  0.0000; zz[6]  = -0.5000; // ( 1/2,    0, -1/2)
      xx[7]  =  0.5000; yy[7]  =  0.0000; zz[7]  =  0.5000; // ( 1/2,    0,  1/2)
      xx[8]  = -0.5000; yy[8]  = -0.5000; zz[8]  =  0.0000; // (-1/2, -1/2,    0)
      xx[9]  = -0.5000; yy[9]  =  0.5000; zz[9]  =  0.0000; // (-1/2,  1/2,    0)
      xx[10] =  0.5000; yy[10] = -0.5000; zz[10] =  0.0000; // ( 1/2, -1/2,    0)
      xx[11] =  0.5000; yy[11] =  0.5000; zz[11] =  0.0000; // ( 1/2,  1/2,    0)
      break;

    case 110 :
      //ideal_symm = 1.54320928882;
      xx[0]  =  0.0000; yy[0]  = -0.7071; zz[0]  =  0.0000; // (   0,  -1/sqrt(2),           0)
      xx[1]  =  0.0000; yy[1]  =  0.7071; zz[1]  =  0.0000; // (   0,   1/sqrt(2),           0)
      xx[2]  = -0.5000; yy[2]  =  0.3536; zz[2]  = -0.3536; // (-1/2,  1/2sqrt(2), -1/2sqrt(2))
      xx[3]  =  0.5000; yy[3]  =  0.3536; zz[3]  = -0.3536; // ( 1/2,  1/2sqrt(2), -1/2sqrt(2))
      xx[4]  = -0.5000; yy[4]  = -0.3536; zz[4]  = -0.3536; // (-1/2, -1/2sqrt(2), -1/2sqrt(2))
      xx[5]  =  0.5000; yy[5]  = -0.3536; zz[5]  = -0.3536; // ( 1/2, -1/2sqrt(2), -1/2sqrt(2))
      xx[6]  =  0.0000; yy[6]  =  0.0000; zz[6]  = -0.7071; // (   0,           0,  -1/sqrt(2))
      xx[7]  =  0.0000; yy[7]  =  0.0000; zz[7]  =  0.7071; // (   0,           0,   1/sqrt(2))
      xx[8]  = -0.5000; yy[8]  = -0.3536; zz[8]  =  0.3536; // (-1/2, -1/2sqrt(2),  1/2sqrt(2))
      xx[9]  =  0.5000; yy[9]  = -0.3536; zz[9]  =  0.3536; // ( 1/2, -1/2sqrt(2),  1/2sqrt(2))
      xx[10] = -0.5000; yy[10] =  0.3536; zz[10] =  0.3536; // (-1/2,  1/2sqrt(2),  1/2sqrt(2))
      xx[11] =  0.5000; yy[11] =  0.3536; zz[11] =  0.3536; // ( 1/2,  1/2sqrt(2),  1/2sqrt(2))
      break;

    case 111 :
      //ideal_symm = 1.25000033391;
      xx[0]  =  0.7071; yy[0]  =  0.0000; zz[0]  =  0.0000;
      xx[1]  = -0.7071; yy[1]  =  0.0000; zz[1]  =  0.0000;
      xx[2]  =  0.0000; yy[2]  = -0.4082; zz[2]  = -0.5774;
      xx[3]  =  0.3536; yy[3]  =  0.6124; zz[3]  =  0.0000;
      xx[4]  = -0.3536; yy[4]  =  0.6124; zz[4]  =  0.0000;
      xx[5]  = -0.3536; yy[5]  =  0.2041; zz[5]  = -0.5774;
      xx[6]  =  0.3536; yy[6]  =  0.2041; zz[6]  = -0.5774;
      xx[7]  =  0.0000; yy[7]  =  0.4082; zz[7]  =  0.5774;
      xx[8]  =  0.3536; yy[8]  = -0.2041; zz[8]  =  0.5774;
      xx[9]  = -0.3536; yy[9]  = -0.2041; zz[9]  =  0.5774;
      xx[10] = -0.3536; yy[10] = -0.6124; zz[10] =  0.0000;
      xx[11] =  0.3536; yy[11] = -0.6124; zz[11] =  0.0000;
      break;

    case 112 :
      //ideal_symm = 1.04770825671;
      xx[0]  =  0.0000; yy[0]  = -0.7071; zz[0]  =  0.0000;
      xx[1]  =  0.0000; yy[1]  =  0.7071; zz[1]  =  0.0000;
      xx[2]  = -0.5774; yy[2]  =  0.0000; zz[2]  = -0.4082;
      xx[3]  =  0.0000; yy[3]  =  0.3536; zz[3]  =  0.6124;
      xx[4]  =  0.0000; yy[4]  = -0.3536; zz[4]  =  0.6124;
      xx[5]  = -0.5774; yy[5]  = -0.3536; zz[5]  =  0.2041;
      xx[6]  = -0.5774; yy[6]  =  0.3536; zz[6]  =  0.2041;
      xx[7]  =  0.5774; yy[7]  =  0.0000; zz[7]  =  0.4082;
      xx[8]  =  0.5774; yy[8]  =  0.3536; zz[8]  = -0.2041;
      xx[9]  =  0.5774; yy[9]  = -0.3536; zz[9]  = -0.2041;
      xx[10] =  0.0000; yy[10] = -0.3536; zz[10] = -0.6124;
      xx[11] =  0.0000; yy[11] =  0.3536; zz[11] = -0.6124;
      break;

    case 113 :
      //ideal_symm = 0.990510750888;
      xx[0]  =  0.7071; yy[0]  =  0.0000; zz[0]  =  0.0000;
      xx[1]  = -0.7071; yy[1]  =  0.0000; zz[1]  =  0.0000;
      xx[2]  =  0.3536; yy[2]  = -0.1066; zz[2]  =  0.6030;
      xx[3]  = -0.3536; yy[3]  = -0.1066; zz[3]  =  0.6030;
      xx[4]  = -0.3536; yy[4]  =  0.5330; zz[4]  =  0.3015;
      xx[5]  =  0.3536; yy[5]  =  0.5330; zz[5]  =  0.3015;
      xx[6]  =  0.0000; yy[6]  = -0.6396; zz[6]  =  0.3015;
      xx[7]  = -0.3536; yy[7]  =  0.1066; zz[7]  = -0.6030;
      xx[8]  =  0.3536; yy[8]  =  0.1066; zz[8]  = -0.6030;
      xx[9]  =  0.0000; yy[9]  =  0.6396; zz[9]  = -0.3015;
      xx[10] =  0.3536; yy[10] = -0.5330; zz[10] = -0.3015;
      xx[11] = -0.3536; yy[11] = -0.5330; zz[11] = -0.3015;
      break;

    case 135 :
      //ideal_symm = 1.41319314949;
      xx[0]  =  0.6325; yy[0]  = -0.2673; zz[0]  = -0.1690;
      xx[1]  = -0.3162; yy[1]  = -0.5345; zz[1]  = -0.3381;
      xx[2]  =  0.4743; yy[2]  =  0.4009; zz[2]  = -0.3381;
      xx[3]  = -0.1581; yy[3]  =  0.6682; zz[3]  = -0.1690;
      xx[4]  = -0.4743; yy[4]  =  0.1336; zz[4]  = -0.5071;
      xx[5]  =  0.1581; yy[5]  = -0.1336; zz[5]  = -0.6761;
      xx[6]  =  0.3162; yy[6]  =  0.5345; zz[6]  =  0.3381;
      xx[7]  =  0.4743; yy[7]  = -0.1336; zz[7]  =  0.5071;
      xx[8]  = -0.1581; yy[8]  =  0.1336; zz[8]  =  0.6761;
      xx[9]  = -0.6325; yy[9]  =  0.2673; zz[9]  =  0.1690;
      xx[10] = -0.4743; yy[10] = -0.4009; zz[10] =  0.3381;
      xx[11] =  0.1581; yy[11] = -0.6682; zz[11] =  0.1690;
      break;
    default:
      cout << "The " << rot_axis << " axis has not been implemented yet.\n";
      return 11;
  }*/

  extractAxis(rot_axis, rotation_axis);
  rotation = rotate2Axis(rotation_axis, z_axis)
  ideal_symm = 0;
  // Calculate the ideal rotation symmetry parameter

  xx[0]  =  0.0000; yy[0]  = -0.5000; zz[0]  = -0.5000; // (   0, -1/2, -1/2)
  xx[1]  =  0.0000; yy[1]  = -0.5000; zz[1]  =  0.5000; // (   0, -1/2,  1/2)
  xx[2]  =  0.0000; yy[2]  =  0.5000; zz[2]  = -0.5000; // (   0,  1/2, -1/2)
  xx[3]  =  0.0000; yy[3]  =  0.5000; zz[3]  =  0.5000; // (   0,  1/2,  1/2)
  xx[4]  = -0.5000; yy[4]  =  0.0000; zz[4]  = -0.5000; // (-1/2,    0, -1/2)
  xx[5]  = -0.5000; yy[5]  =  0.0000; zz[5]  =  0.5000; // (-1/2,    0,  1/2)
  xx[6]  =  0.5000; yy[6]  =  0.0000; zz[6]  = -0.5000; // ( 1/2,    0, -1/2)
  xx[7]  =  0.5000; yy[7]  =  0.0000; zz[7]  =  0.5000; // ( 1/2,    0,  1/2)
  xx[8]  = -0.5000; yy[8]  = -0.5000; zz[8]  =  0.0000; // (-1/2, -1/2,    0)
  xx[9]  = -0.5000; yy[9]  =  0.5000; zz[9]  =  0.0000; // (-1/2,  1/2,    0)
  xx[10] =  0.5000; yy[10] = -0.5000; zz[10] =  0.0000; // ( 1/2, -1/2,    0)
  xx[11] =  0.5000; yy[11] =  0.5000; zz[11] =  0.0000; // ( 1/2,  1/2,    0)

  for (unsigned int i = 0; i < xx.size(); ++i)
  {
    xtemp = rotation[0][0] * xx[i] + rotation[0][1] * yy[i] + rotation[0][2] * zz[i];
    ytemp = rotation[1][0] * xx[i] + rotation[1][1] * yy[i] + rotation[1][2] * zz[i];

    sintheta_sq = 1 - ((xtemp * xtemp) / (xtemp * xtemp + ytemp * ytemp));
    if (isnan(sintheta_sq))
    {
      // This is when the atoms lie on top of each other when projected onto the
      // xy plane.
      sintheta_sq = 1;
    }

    ideal_symm += (coeffs[0] - coeffs[1] * sintheta_sq) * (coeffs[0] - coeffs[1] * sintheta_sq) * sintheta_sq;
  }

  ideal_symm /= xx.size();


  for (unsigned int i = 0; i < xx.size(); ++i)
  {
    x = rotation[0][0] * xx[i] + rotation[0][1] * yy[i] + rotation[0][2] * zz[i];
    y = rotation[1][0] * xx[i] + rotation[1][1] * yy[i] + rotation[1][2] * zz[i];

    xtemp = costheta * x - sintheta * y;
    ytemp = sintheta * x + costheta * y;

    /* Uses the idea that sin^2 = 1-cos^2
    * Projection onto the XY plane means we ignore the z coordinates
    * cos = A.B / (|A||B|); this simplifies to
    * cos = (A_x + B_x) / |A||B|, meaning that
    * cos^2 = (A_x * B_x + A_y * B_y + A_z * B_z)^2 / (A_x^2 + A_y^2 + A_z^2)(B_x^2 + B_y^2 + B_z^2)
    * Here, A is the vector representing the distance between a central atom i
    * and a nearest neighbor atom j, and B is the unit vector in the X direction
    * or (100).  This simplifies the above equation to:
    * cos^2 = A_x^2 / (A_x^2 + A_y^2)
    */

    sintheta_sq = 1 - ((xtemp * xtemp) / (xtemp * xtemp + ytemp * ytemp));
    if (isnan(sintheta_sq)) // Handles the case of dividing by 0
    {
      // This is when the atoms lie on top of each other when projected onto the
      // xy plane.
      sintheta_sq = 1;
    }
    // Symmetry parameter as defined by Zhang. Coefficients are changed to get
    // better resolution between grains.
    total1 += (coeffs[0] - coeffs[1] * sintheta_sq) * (coeffs[0] - coeffs[1] * sintheta_sq) * sintheta_sq;
  }
  total1 /= xx.size();
  // Cutoff is the midpoint between the ideal value and the rotated ideal value.
  // Note that this assumes that the outside grain is oriented with it's x axis
  // aligned with the x axis of the lab frame.
  cutoff = (ideal_symm + total1) / 2.0;

  cout << "\tIdeal symmetry parameter: " << ideal_symm << endl
       << "\tRotated symmetry parameter: " << total1 << endl;

  // Now read through each set of files
  int j = 1;
  while (getline(fin_input, filename1))
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
      if (j == 1)
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
        ss >> atom_id >> type >> x >> y >> z;
        charge = 0.0;
      }
      else if (n_type == 2)
      {
        ss >> atom_id >> type >> charge >> x >> y >> z;
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
        continue; // Only focus on U (or the single atom type)
      }

      x = atoms[i].getX();
      y = atoms[i].getY();
      z = atoms[i].getZ();
      total2 = 0.0; // reset the total
      // We start at l = 1 because if we start at l = 0, we just re-use the same
      // atom over and over.
      for (int l = 1; l <= iatom[0][i]; ++l)
      {
        int id = iatom[l][i];

        // calculate the distances
        // We project onto the xy plane, effectively ignoring the z coordinate
        rxij = atoms[id].getX() - x;
        ryij = atoms[id].getY() - y;
        rzij = atoms[id].getZ() - z;

        // Apply PBCs
        rxij = rxij - anInt(rxij / Lx) * Lx;
        ryij = ryij - anInt(ryij / Ly) * Ly;
        rzij = rzij - anInt(rzij / Lz) * Lz;

        //Apply the rotation as above
        xtemp = rotation[0][0] * rxij + rotation[0][1] * ryij + rotation[0][2] * rzij;
        ytemp = rotation[1][0] * rxij + rotation[1][1] * ryij + rotation[1][2] * rzij;

        // Calculate the magnitude of the distance
        drij_sq = (xtemp * xtemp) + (ytemp * ytemp);

        if (drij_sq == 0) // Handles the case where the projected position of the atom is right on top of the current atom.
        {
          total2 += 1;
          continue;
        }
        // This uses the relation sin^2 = 1-cos^2, where cos = dot(A,B) / (|A|*|B|)
        sintheta_sq = 1 - ((xtemp * xtemp) / drij_sq);
        total2 += (coeffs[0] - coeffs[1] * sintheta_sq) * (coeffs[0] - coeffs[1] * sintheta_sq) * sintheta_sq;
      }
      total2 /= iatom[0][i]; // This may not always be 12!
      symm[i] = total2; // Store them for analysis

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

    // This is a check to make sure we are outputting the values for the shrinking grain.
    if (j == 1)
    {
      if (n_grain_1 < n_grain_2)
      {
        smaller = 1;
      }
      else //n_grain_1 >= n_grain_2
      {
        smaller = 2;
      }
    }
    else // j != 1
    {
      if (smaller == 1)
      {
        if (!(n_grain_1 < n_grain_2))
        {
          cout << "Grains were labeled incorrectly.\n";
          return 9;
        }
      }
      else // (smaller == 2)
      {
        if (!(n_grain_1 > n_grain_2)) // Ignores the possibility of n1 == n2
        {
          cout << "Grains were labeled incorrectly.\n";
          return 9;
        }
      }
    }
    // We want to make sure we ouput the smaller value
    fout_data << (n_grain_1 < n_grain_2 ? n_grain_1 : n_grain_2) << endl;

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
        fout << atoms[i].getId() << " "
             << atoms[i].getType() << " "
             << (atoms[i].getX() + xlow) * a0 << " "
             << (atoms[i].getY() + ylow) * a0 << " "
             << (atoms[i].getZ() + zlow) * a0 << " "
             << atoms[i].getMark() << " "
             << symm[i] << endl;
      }
      else if (n_type == 2)
      {
        fout << atoms[i].getId() << " "
             << atoms[i].getType() << " "
             << atoms[i].getCharge() << " "
             << (atoms[i].getX() + xlow) * a0 << " "
             << (atoms[i].getY() + ylow) * a0 << " "
             << (atoms[i].getZ() + zlow) * a0 << " "
             << atoms[i].getMark() << " "
             << symm[i] << endl;
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
    ++j;
  }

  fin_input.close();
  fout_data.close();
  return 0;
}

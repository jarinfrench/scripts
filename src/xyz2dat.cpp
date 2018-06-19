/******************************************************************************
* This script converts a file in the format .xyz to the .dat format.  Note that
* if copying data from GBStudio, the output format MUST be xmol, otherwise there
* will be extra information that this script will not handle properly.
******************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cmath> // for cos/sin
#include <algorithm> // for min/max_element
#include <complex>
#include <cxxopts.hpp>
#include "atom.h"

#define PI 3.14159265358979

using namespace std;

bool compareAtomsX(const Atom &a, const Atom &b) {return a.getX() < b.getX();}
bool compareAtomsY(const Atom &a, const Atom &b) {return a.getY() < b.getY();}
bool compareAtomsZ(const Atom &a, const Atom &b) {return a.getZ() < b.getZ();}

int main(int argc, char** argv)
{
  string filename1, filename2, filename3, str;
  int n_atoms, N, type, n_types;
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double x, y, z, charge, temp = 1.0;
  double Lx, Ly, Lz; // Box lengths in each direction
  double a, b, c, alpha, beta, gamma; // lattice parameters, lattice angles
  string atom_type, chem_formula;
  map <string, int> elements;
  vector <Atom> atoms;
  bool show_charge = false;
  // cxxopts::Options options(argv[0], "Convert a file in the xyz format to the dat format");
  // options.add_options()
  //   ("f,file", "File name", cxxopts::value<string>())
  //   ("s,show-charge", "Whether or not to write the charge values");

  if (argc == 1)
  {
    cout << "Please enter the file to be converted: ";
    cin >> filename1;
  }
  else
  {
    filename1 = argv[1];
  }

  // Determine where in the filename the element is specified.  Assumes a filename
  // format of "LAMMPS_<element(s)>_N<number_of_atoms>[_<extra_info>].dat"
  size_t element_pos = filename1.find("LAMMPS_")+7;
  size_t element_pos_end;
  { // Create a local namespace so we only store the results of this subroutine
    size_t temp = filename1.find("_N"); // Find where the '_N' is.
     // There are some elements that begin with N, so we check to make sure we
     // found the actual element list
    temp = filename1.find("_N", temp + 1);
    // if the temp variable does not find any other '_N' strings in the filename...
    if (temp == string::npos)
    {
      // The end of the elements position is where the '_N' string begins
      element_pos_end = filename1.find("_N");
    }
    else
    {
      // Otherwise, the end of the elements position is where the second _N string begins.
      element_pos_end = temp;
    }
  }

  if (element_pos == string::npos) // Error if we didn't find the 'LAMMPS_' string in the filename
  {
    cout << "Error determining element(s).  Make sure the filename is in the format\n"
         << "\t\"LAMMPS_<element(s)>_N<number_of_atoms>[_<extra_info>].dat\"\n"
         << "Note that elements should be properly capitalized, i.e. table salt would be listed as NaCl.\n";
    return 8;
  }
  else // Determine the number of elements.  Note that because we are using a map, duplicated elements will not appear more than once in the list
  {
    int elem_num = 1;
    for (unsigned int i = element_pos; i < element_pos_end; ++i)
    {
      if (islower(filename1[i])) // if the current letter is lowercase
      {
        // Our element is the character before plus this character
        if ((elements.insert(pair<string,int>(filename1.substr(i-1,2), elem_num))).second) // check to see if insertion was successful - fails if a duplicate
        {
          ++elem_num;
        }
      }
      else // current letter is uppercase
      {
        if (islower(filename1[i+1])) // if the next character is lowercase, we handle that in the next loop
        {
          continue;
        }
        else
        {
          // if the current character is a number, move ahead
          // we don't need to calculate the element ratios
          if (isdigit(filename1[i]))
          {
            continue;
          }
          else // Otherwise, we have a single-letter element (like O for oxygen)
          {
            if ((elements.insert(pair<string, int> (filename1.substr(i,1), elem_num))).second) // check to see if insertion was successful - fails if a duplicate
            {
              ++elem_num;
            }
          }
        }
      }
    }
  }
  chem_formula = filename1.substr(element_pos, element_pos_end - element_pos);

  n_types = elements.size();

  filename2 = filename1.substr(0, filename1.find(".")) + ".dat";

  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Unable to open file " << filename1 << endl;
    return 1;
  }

  ofstream fout(filename2.c_str());
  if (fout.fail())
  {
    cout << "Unable to open file " << filename2 << endl;
    return 1;
  }

  fout << "This bulk " << chem_formula << " coordinates format: [ID type charge* x y z] (* not included if always charge neutral)\n\n";

  fin >> N; // Number of atoms
  fout << N << "  atoms\n";
  atoms.resize(N); // pre-allocate the vector size to save time and space.

  fout << n_types << "   atom types\n";

  fin >> Lx >> Ly >> Lz >> alpha >> beta >> gamma; // Boundaries

  n_atoms = 0;
  while (fin >> atom_type >> x >> y >> z) // Read in the data
  {
    if (atom_type.compare("U") == 0) // Given a U atom, write the correct type and charge
    {
      type = 1;
      charge = 2.4;
      //show_charge = true;
    }
    else if (atom_type.compare("O") == 0) // Same with the O atom
    {
      type = 2;
      charge = -1.2;
      show_charge = true;
    }
    else
    {
      type = elements[atom_type];
      charge = 0.0; // This will have to be manually fixed later if there is a charge
    }
    ++n_atoms;
    Atom tmp = Atom(n_atoms, type, charge, x, y, z);
    atoms[n_atoms - 1] = tmp;
  }
  if (n_atoms != N)
  {
    cout << "Error counting atoms: n_atoms = " << n_atoms << " != N = " << N << endl;
    return 3;
  }

  xlow = (*min_element(atoms.begin(), atoms.end(), compareAtomsX)).getX();
  ylow = (*min_element(atoms.begin(), atoms.end(), compareAtomsY)).getY();
  zlow = (*min_element(atoms.begin(), atoms.end(), compareAtomsZ)).getZ();
  xhigh = xlow + Lx;
  yhigh = ylow + Ly;
  zhigh = zlow + Lz;

  fout << fixed;
  fout << setprecision(6) << xlow << " " << xhigh << " xlo xhi\n";
  fout << ylow << " " << yhigh << " ylo yhi\n";
  fout << zlow << " " <<  zhigh << " zlo zhi\n";

  if (alpha != 90.00 || beta != 90.00 || gamma != 90.00)
  {
    x = Lx;
    y = Ly;
    z = Lz;

    complex <double> xy, xz, yz, cos_a, cos_b, cos_g, sin_g; // Tilt parameters, sines and cosines
    complex <double> ii;
    ii.imag(1); // the square root of -1
    double ly;

    cout << "Triclinic box detected; calculating tilt parameters.\n";
    cos_a = cos(alpha * PI / 180.0);
    cos_b = cos(beta * PI / 180.0);
    cos_g = cos(gamma * PI / 180.0);
    sin_g = sin(gamma * PI / 180.0);

    // This needs some work... it depends on the orientation of the parallelepiped
    Lx = x + (y * sin(abs(90.00 - gamma) * PI / 180.0)) + z * (cos_b.real());
    Ly = 1.0;

    if (sin_g.real() < 1E-8 && sin_g.imag() < 1E-8)
    {
      cout << "Error: z cannot be coplanar with the xy plane!  Gamma must be > 0 and < 180\n";
    }
    else
    {
      if (gamma == 90.00)
      {
        xy = 0.0;
      }
      else
      {
        xy = Ly * cos_g / sqrt(sin_g * sin_g);
      }

      if (beta == 90.00)
      {
        xz = 0.0;
      }
      else
      {
        xz = (cos_b * ii * sqrt(Lz * Lz * sin_g * sin_g)) / \
            (sqrt(-1.0 + cos_a * cos_a + cos_b + cos_b + cos_g * cos_g - 2.0 * cos_a * cos_b * cos_g));
      }
      yz = ((cos_a - cos_b * cos_g) * ii * sqrt(Lz * sin_g * sin_g)) / \
            (sqrt(-1.0 + cos_a * cos_a + cos_b + cos_b + cos_g * cos_g - 2.0 * cos_a * cos_b * cos_g) * sqrt(sin_g * sin_g));
    }

    if (xy.imag() > 1E-8 || xz.imag() > 1E-8 || yz.imag() > 1E-8)
    {
      cout << "Error: imaginary value determined.\n";
      return 10;
    }
    else
    {
      fout << xy.real() << " " << xz.real() << " " << yz.real() << " xy xz yz\n";
    }
  }

  fout << "\nAtoms\n\n";

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    fout << setprecision(0) << atoms[i].getId() << " " << atoms[i].getType() << " ";

    if (show_charge)
    {
      fout << setprecision(1) << atoms[i].getCharge() << " ";
    }
    fout << setprecision(6) << atoms[i].getX() << " " << atoms[i].getY() << " " << atoms[i].getZ() << endl;
  }

  // Close the file streams
  fin.close();
  fout.close();


  return 0;
}

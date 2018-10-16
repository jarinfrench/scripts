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
#include "error_code_defines.h"

#define PI 3.14159265358979

using namespace std;

// various comparison functions for the atoms
bool compareAtomsX(const Atom &a, const Atom &b) {return a.getX() < b.getX();}
bool compareAtomsY(const Atom &a, const Atom &b) {return a.getY() < b.getY();}
bool compareAtomsZ(const Atom &a, const Atom &b) {return a.getZ() < b.getZ();}

vector <double> determineTiltParameters(const double& Lx, const double& Ly,
                                        const double& Lz, const double& alpha,
                                        const double& beta, const double& gamma)
{
  double x = Lx, y = Ly, z = Lz;
  complex <double> xy, xz, yz, cos_a, cos_b, cos_g, sin_g; // tilt parameters, (co)sines of the angles
  complex <double> ii;
  ii.imag(1); // the square root of -1
  double ly;

  cout << "Triclinic box detected; calculating tilt parameters.\n";
  cos_a = cos(alpha * PI / 180.0);
  cos_b = cos(beta * PI / 180.0);
  cos_g = cos(gamma * PI / 180.0);
  sin_g = sin(gamma * PI / 180.0);

  // This needs some work... it depends on the orientation of the parallelepiped
  x = x + (y * sin(abs(90.00 - gamma) * PI / 180.0)) + z * (cos_b.real());
  y = 1.0;

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
      xy = y * cos_g / sqrt(sin_g * sin_g);
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
    exit(IMAGINARY_TILT_ERROR);
  }
  else
  {
    vector <double> tilt_params(3,0.0);
    tilt_params[0] = xy.real();
    tilt_params[1] = xz.real();
    tilt_params[2] = yz.real();
    return tilt_params;
  }
}

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

string determineElement(const string& file)
{
  // Filename must be in the format "LAMMPS_<element(s)>_N<number_of_atoms>[_<extra_info>].dat"
  size_t element_pos_start = file.find("LAMMPS_") + 7;
  size_t element_pos_end = file.find("_", element_pos_start);

  if (element_pos_start == string::npos)
  {
    cout << "Error determining element(s).  Make sure the filename is in the format\n"
         << "\t\"LAMMPS_<element(s)>_N<number_of_atoms>[_<extra_info>].dat\"\n"
         << "Note that elements should be properly capitalized, i.e. table salt would be NaCl.\n";
    exit(FILE_FORMAT_ERROR);
  }
  else // Determine the number of elements.  Note that because we are using a map, duplicated elements will not appear more than once in the list
  {
    return file.substr(element_pos_start, element_pos_end - element_pos_start);
  }
}

void convertFile(const string& file, const bool show_charge, const bool no_element_listed)
{
  string file2, str; // second file (we write to this one)
  string chem_formula = "<element>"; // chemical formula
  map <string, int> elements; // map of element name to element number
  int n_types = 1, N, type, n_atoms; // number of atom types, number of atoms, atom type, number of atoms read
  vector <Atom> atoms; // container for the atoms
  double Lx, Ly, Lz, alpha, beta, gamma; // Box bounds, angles
  string atom_type; // atom type listed in the original file
  double x, y, z, charge; // atom x, y, and z positions in the file, charge of atom
  double xlow, ylow, zlow, xhigh, yhigh, zhigh; // min and max of atom positions

  if (!no_element_listed) // If an element is specified in the filename, find it
  {
    chem_formula = determineElement(file);
    int elem_num = 1;
    for (unsigned int i = 0; i < chem_formula.size(); ++i)
    {
      if (islower(chem_formula[i])) // if the current letter is lowercase
      {
        // Our element is the character before plus this character
        // check to see if insertion was successful - fails if a duplicate
        if ((elements.insert(pair<string,int>(chem_formula.substr(i-1,2), elem_num))).second)
        {
          ++elem_num;
        }
      }
      else // current letter is uppercase
      {
        if (islower(chem_formula[i+1])) // if the next character is lowercase, we handle that in the next loop
        {
          continue;
        }
        else
        {
          // if the current character is a number, move ahead
          // we don't need to calculate the element ratios
          if (isdigit(chem_formula[i]))
          {
            continue;
          }
          else // Otherwise, we have a single-letter element (like O for oxygen)
          {
            if ((elements.insert(pair<string, int> (chem_formula.substr(i,1), elem_num))).second)
            {
              ++elem_num;
            }
          }
        }
      }
    }
    n_types = elements.size();
  }

  file2 = file.substr(0,file.find(".")) + ".dat"; // Assumes only one "." symbol in file.

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  ofstream fout(file2.c_str());
  checkFileStream(fout, file2);

  if (!show_charge)
  {
    fout << "This bulk " << chem_formula << " coordinates format: [ID type x y z]\n\n";
  }
  else
  {
    fout << "This bulk " << chem_formula << " coordinates format: [ID type charge x y z]\n\n";
  }

  fin >> N; // Number of atoms
  getline(fin, str); // get the remainder of the line
  fout << N << "  atoms\n";
  atoms.resize(N); // pre-allocate the vector size to save time and space.

  fout << n_types << "   atom types\n";

  fin >> Lx >> Ly >> Lz >> alpha >> beta >> gamma; // boundaries and angles
  getline(fin,str); // get the remainder of the line

  n_atoms = 0;
  while (fin >> atom_type >> x >> y >> z) // Read in the data
  {
    if (atom_type.compare("U") == 0)
    {
      type = 1;
      charge = 2.4;
    }
    else if (atom_type.compare("O") == 0)
    {
      type = 2;
      charge = -1.2;
    }
    else
    {
      type = elements[atom_type];
      charge = 0.0; // Needs to be manually fixed later if there is a charge
    }
    ++n_atoms;
    atoms[n_atoms - 1] = Atom(n_atoms, type, charge, x, y, z);
  }

  if (n_atoms != N)
  {
    cout << "Error counting atoms: n_atoms = " << n_atoms << " != N = " << N << endl;
    exit(ATOM_COUNT_ERROR);
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
    vector <double> tilt_params(3,0.0);
    tilt_params = determineTiltParameters(Lx,Ly,Lz,alpha,beta,gamma);
    fout << tilt_params[0] << " " << tilt_params[1] << " " << tilt_params[2] << " xy xz yz\n";
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
}

int main(int argc, char** argv)
{
  bool show_charge = false;
  bool no_element_listed = false;

  try
  {
    cxxopts::Options options(argv[0], "Convert a file in the xyz format to the dat format");
    options
      .positional_help("File")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Input file", cxxopts::value<vector<string> >(), "file")
        ("s,show-charge", "Include charge in the converted file", cxxopts::value<bool>(show_charge)->default_value("false"))
        ("n,no-element", "Ignores presence of element names", cxxopts::value<bool>(no_element_listed)->default_value("false"))
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("file") == 0)
    {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (result.count("file"))
    {
      auto& filenames = result["file"].as<vector<string> >();
      for (const auto& file : filenames)
      {
        convertFile(file, show_charge, no_element_listed);
      }
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

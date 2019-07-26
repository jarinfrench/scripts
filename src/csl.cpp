// This code is based on the paper found here:
// https://journals.iucr.org/q/issues/1966/08/00/a05208/a05208.pdf
// Ranganathan, S. Acta Cryst. 21 (1966) 197.

#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <cmath> // for atan
#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

#define PI 3.14159265358979

double calculateAxialRatio(const vector<int>& vec)
{
  return sqrt(inner_product(vec.begin(), vec.end(), vec.begin(), 0));
}

int calculateAtomsPerCell(const vector<int>& vec)
{
  return vec[1] * vec[1] + vec[2] * vec[2] - 1;
}

int calculateSigmaNumber(const vector<int>& vec)
{
  int sumOfSquares = inner_product(vec.begin(),vec.end(),vec.begin(),0);
  if (sumOfSquares % 2 == 0) {return sumOfSquares / 2;}
  else {return sumOfSquares;}
}

double calculateMisorientation(const vector<int>& normal)
{
  if (abs(normal[0]) < 1.0e-8)
  {
    return 90.00;
  }
  else
    return 2.0 * atan((double)(normal[1]) / (double)(normal[0])) * 180.0 / PI;
}

void assignOrthogonalDirections(const vector<int>& normal, vector <int>& dir_x, vector<int>& dir_y)
{
  dir_y[0] = -(normal[1] * normal[1] + normal[2] * normal[2]); // - (k^2 + l^2)
  dir_y[1] = normal[0] * normal[1]; // hk
  dir_y[2] = normal[0] * normal[2]; // hl

  dir_x[0] = 0;
  dir_x[1] = normal[2]; // l
  dir_x[2] = -normal[1]; // -k
}

int main(int argc, char** argv)
{
  int h,k,l;
  // number of atoms in the boundary unit cell, number of atoms in the unit cell of the crystal structure
  int n_per_boundary_unit_cell, n_per_unit_cell;
  int sigma;
  double axial_ratio, misorientation;
  string structure;
  vector <int> normal (3,0), dir_x (3,0), dir_y (3,0);
  try
  {
    cxxopts::Options options(argv[0], "Determine CSL info based on any given hkl.");
    options
      .positional_help("h k l")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("h", "First index of grain boundary rotation axis", cxxopts::value<int>(h), "<integer>")
        ("k", "Second index of grain boundary rotation axis", cxxopts::value<int>(k), "<integer>")
        ("l", "Third index of grain boundary rotation axis", cxxopts::value<int>(l), "<integer>")
        ("s,structure", "Crystal structure.  Must be sc, fcc, or bcc", cxxopts::value<string>(structure), "sc|fcc|bcc")
        ("unordered", "Keep the order as input (h !< k !< l)")
        ("help", "Show the help");

    options.parse_positional({"h", "k", "l"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !(result.count("h") && result.count("k") && result.count("l")))
    {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (result.count("structure"))
    {
      if (structure.compare("sc") == 0) {n_per_unit_cell = 1;}
      else if (structure.compare("fcc") == 0) {n_per_unit_cell = 4;}
      else if (structure.compare("bcc") == 0) {n_per_unit_cell = 2;}
    }
    else {n_per_unit_cell = 1;}

    normal[0] = h; normal[1] = k; normal[2] = l;
    if (!result.count("unordered"))
    {
      sort(normal.begin(), normal.end());
      cout << "Sorted values of hkl are ";
      for (unsigned int i = 0; i < normal.size(); ++i) {cout << normal[i] << " ";}
      cout << "\nIf you would like to keep the order as input, use the --unordered flag.\n";
    }

    assignOrthogonalDirections(normal, dir_x, dir_y);
    cout << "Assuming that hkl is the z direction:\n"
         << "  x = " << dir_x[0] << " " << dir_x[1] << " " << dir_x[2] << "\n"
         << "  y = " << dir_y[0] << " " << dir_y[1] << " " << dir_y[2] << endl;

    misorientation = calculateMisorientation(normal);
    cout << "Calculated misorientation angle: " << misorientation << endl;

    sigma = calculateSigmaNumber(normal);
    cout << "Sigma number:" << sigma << endl;

    n_per_boundary_unit_cell = calculateAtomsPerCell(normal);
    cout << "Number of atoms per boundary unit cell: "
         << n_per_boundary_unit_cell * n_per_unit_cell << endl;

    axial_ratio = calculateAxialRatio(normal);
    cout << "Axial ratio: " << axial_ratio << endl;
  }
  catch (cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

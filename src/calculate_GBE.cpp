/******************************************************************************
* This script reads a file as generated by extract_all_energies.sh, and
* calculates the per atom energy at the grain boundary.  The filename can be
* passed in via command line.  Note that this program expects a filename in the
* format of 111Tilt_total_energies.csv
* IMPORTANT PARAMETERS
*   grain radius: the radius of the grain that was specified when rotating the grain.
*     note that this assumes a cylindrical grain.  Given in Angstroms.
*   grain thickness: the thickness of the grain (determined by Lz).  This is given
*     in Angstroms.
******************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

#define PI 3.14159265
#define eV2J 16.0218 // converts eV/angstrom^2 to J/m^2

int main(int argc, char **argv)
{
  string filename1, filename2, boundary_type;
  bool is_sphere = false;
  int N0, N; // number of atoms for single grain, rotated grain
  double eGB, e_single, theta; // energy for rotated grain, single grain, angle of rotation
  double r_grain; //radius of the grain in angstroms
  double Lz; // thickness of the grain in angstroms
  double gbe, gb_area; // grain boundary energy, denominator term to calculate gbe

  if (argc != 3)
  {
    // Prompt for filename if we don't already have it.
    cout << "Please enter the filename containing the angle and energies: ";
    cin  >> filename1;

    cout << "Please enter the boundary type (cylinder|sphere): ";
    cin >> boundary_type;
  }
  else
  {
    filename1 = argv[1];
    boundary_type = argv[2];
  }

  transform(boundary_type.begin(), boundary_type.end(), boundary_type.begin(), ::tolower);
  if (boundary_type.compare("cylinder") != 0 && boundary_type.compare("sphere") != 0)
  {
    cout << "Unable to process boundary type " << boundary_type << ", please enter either \"cylinder\" or \"sphere\"\n";
    return 20;
  }
  else
  {
    if (boundary_type.compare("sphere") == 0)
    {
      is_sphere = true;
    }
  }


  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error! Cannot open file " << filename1 << endl;
    return -1;
  }

  // Prompt for the grain radius
  cout << "Please enter the grain radius used for these simulations: ";
  cin  >> r_grain;

  //Prompt for the grain thickness
  if (!is_sphere)
  {
    cout << "Please enter the grain thickness (Lz) in angstroms: ";
    cin  >> Lz;
  }

  // This assumes a filename given above in the format 111Tilt_total_energies.csv,
  // and creates a second file called 111Tilt_individual_energy.csv
  filename2 = filename1.substr(0,filename1.find("_total")) + "_individual_energy.csv";
  ofstream fout(filename2.c_str());
  if (fout.fail())
  {
    cout << "Error! Cannot open file " << filename2 << endl;
    return -1;
  }

  if (is_sphere)
  {
    gb_area = 4 * PI * r_grain * r_grain;
  }
  else // cylinder
  {
    gb_area = 2 * PI * r_grain * Lz; // area of the grain boundary
  }

  // This assumes that the first line contains the original, single grain
  fin >> theta >> e_single >> N0;
  fout << theta << "," << setprecision(15) << 0.0 << endl;
  while (fin >> theta >> eGB >> N) // Now read each line of data
  {
    if (fin.fail())
      cout << "Read error.\n";

    gbe = (eGB - (e_single / N0) * N) / gb_area * eV2J; // Calculate the GBE
    fout << theta << "," << setprecision(15) << gbe << endl; // Write to the new file
  }

  // Close the file streams
  fin.close();
  fout.close();

  return 0;
}
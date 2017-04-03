#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

#define PI 3.14159265
#define eV2J 16.0218 // converts eV/angstrom^2 to J/m^2

int main(int argc, char **argv)
{
  string filename1, filename2;
  int N0, N; // number of atoms for single grain, rotated grain
  double eGB, e_single, theta; // energy for rotated grain, single grain
  double r_grain = 100; //radius of the grain in angstroms
  double Lz = 27.2651; // thickness of the grain in angstroms
  double gbe, bot;

  if (argc != 2)
  {
    cout << "Please enter the filename containing the angle and energies: ";
    cin  >> filename1;
  }
  else
    filename1 = argv[1];

  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error! Cannot open file " << filename1 << endl;
    return -1;
  }

  filename2 = filename1.substr(0,filename1.find("_total")) + "_individual_energy.csv";
  ofstream fout(filename2.c_str());
  if (fout.fail())
  {
    cout << "Error! Cannot open file " << filename2 << endl;
    return -1;
  }

  bot = 4 * PI * r_grain * Lz; // twice the area of the grain boundary
  fin >> theta >> e_single >> N0;
  fout << theta << "," << setprecision(15) << 0.0 << endl;
  while (fin >> theta >> eGB >> N)
  {
    if (fin.fail())
      cout << "Read error.\n";

    gbe = (eGB - (e_single / N0) * N) / bot * eV2J;
    fout << theta << "," << setprecision(15) << gbe << endl;
  }

  fin.close();
  fout.close();

  return 0;
}

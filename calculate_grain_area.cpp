#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;

#define PI 3.141592653589793

double latticeParam(double T)
{
  // A is the y intercept, B is the linear coefficient, C is the parabolic
  // coefficient
  double A, B, C;
  if (T >= 0.0 && T <= 1000.0)
  {
    A = 5.45297739;
    B = 5.89383E-5;
    C = 0.0;
  }
  else if (T > 1000.0 && T <= 3300.00)
  {
    A = 5.44364863;
    B = 6.18966E-5;
    C = 5.27784E-9;
  }
  else
  {
    cout << "Temperature out of fitted range (0 K - 3300 K).\n";
    exit(10); // We don't want to continue with execution if we're out of range.
  }

  return A + B * T + C * T * T;
}

int main(int argc, char **argv)
{
  string filename, filename2, str; // data file, junk variable.
  double T, a0, M, r_sq; // temperature, lattice parameter, mobility
  int N0, N1; // Number of atoms in grain 2 at time 1 and time 2
  double t0, t1, Lz; // time at 1, time at 2, height of cylindrical grain, GBE
  // NOTE: eGB is generally taken as 1.5.  I should use a value close to what
  // was calculated in the GBE calculations.

  if (argc != 4)
  {
    cout << "Please enter the filename containing the number of atoms in each grain: ";
    cin  >> filename;

    cout << "Please enter the temperature (in K): ";
    cin  >> T;

    cout << "Please enter the height of the original cylinder (in Angstroms): ";
    cin  >> Lz;
  }
  else
  {
    filename = argv[1]; // get the filename
    T = strtod(argv[2], NULL); // get the temperature.
    Lz = strtod(argv[3], NULL); // get the cylinder height
  }

  a0 = latticeParam(T);
  // open up the file for reading
  ifstream fin(filename.c_str());
  if (fin.fail())
  {
    cout << "Error: unable to open file " << filename << endl;
    return 1;
  }

  ofstream fout("mobility_data.csv");
  if (fout.fail())
  {
    cout << "Error: unable to open file mobility_data.csv.\n";
    return 1;
  }
  fout << "# This is the area data for T = " << T << " K [time, area (Angstroms^2)]\n";

  getline(fin, str); // get the comment line
  fin >> t0 >> N0;
  if (t0 != 0)
  {
    cout << "Something is wrong in the input file - Perhaps a missing entry?\n"
         << "This script requires that the second line in the data file contain the timestep 0 information!\n";
    return 2;
  }
  fout << t0 * 0.002 << "," << N0 * a0 * a0 * a0 / (4 * Lz) << endl;
  fin >> str >> str; // ignore the minimization step

  while (fin >> t1 >> N1)
  {
    // if t1 is less than or equal to t0, the data file is corrupted and needs
    // to be looked at.
    if (t1 <= t0)
    {
      cout << "Error: data file corrupted.  t0 = " << t0 << " >= t1 = " << t1 << endl;
      return 3;
    }
    // Now that we've multiplied everything together, convert it to the correct
    // units.  Units should be m^4 / J*s.
    r_sq = N1 * a0 * a0 * a0 / (4 * Lz);

    fout << t1 * 0.002 << "," << r_sq << endl;
  }
  // Close the file stream
  fin.close();
  fout.close();

  return 0;
}

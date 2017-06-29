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
  double T, a0, M; // temperature, lattice parameter, mobility
  int N0, N1; // Number of atoms in grain 2 at time 1 and time 2
  double t0, t1, Lz, eGB; // time at 1, time at 2, height of cylindrical grain, GBE

  if (argc != 5)
  {
    cout << "Please enter the filename containing the number of atoms in each grain: ";
    cin  >> filename;

    cout << "Please enter the temperature (in K): ";
    cin  >> T;

    cout << "Please enter the height of the original cylinder (in Angstroms): ";
    cin  >> Lz;

    cout << "Please enter the grain boundary energy (in J/m^2): ";
    cin  >> eGB;
  }
  else
  {
    filename = argv[1]; // get the filename
    T = strtod(argv[2], NULL); // get the temperature.
    Lz = strtod(argv[3], NULL); // get the cylinder height
    eGB = strtod(argv[4], NULL); // get the GB energy
  }

  a0 = latticeParam(T);
  // open up the file for reading
  ifstream fin(filename.c_str());
  if (fin.fail())
  {
    cout << "Error: unable to open file " << filename << endl;
    return 1;
  }

  getline(fin, str); // get the comment line
  getline(fin, str); // get the initial, unminimized structure.
  fin >> t0 >> str >> str >> str >> N0; // Get the initial values

  int i = 0;
  while (fin >> t1 >> str >> str >> str >> N1)
  {
    // Now that we've multiplied everything together, convert it to the correct
    // units.  Units should be m^4 / J*s.  Note that we have N0 before N1 to
    // calculate positive values.
    M += (N0 - N1) * a0 * a0 * a0 / ((t1 - t0) * 0.002 * 8 * PI * Lz * eGB) * 1.0e-8;
    //cout << M << endl;
    ++i;
  }
  // Close the file stream
  fin.close();

  cout << "The mobility is " << M / i << endl;

  stringstream ss;
  ss << "mobility_T" << (int)T << "K.csv";
  filename2 = ss.str();

  ofstream fout(filename2.c_str(), ios_base::app);
  if (fout.fail())
  {
    cout << "Error: unable to open file " << filename2 << endl;
    return 1;
  }

  fout << M / i << " " << T << endl;

  return 0;
}

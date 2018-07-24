#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

struct fit
{
  string name;
  double T1, T2;
  double lin_y_int, lin_slope;
  double parabolic_int, parabolic_lin, parabolic_term;
};

struct inputData
{
  string datafile;
  double temperature, height, a0;
} input;

#define PI 3.141592653589793

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

void showInputFileHelp()
{
  cout << "The input file must contain the following on one line, in order:\n"
       << "\t1) The filename containing the number of atoms in each grain, formatted as output by find_grains\n"
       << "\t2) The temperature in Kelvin\n"
       << "\t3) The height of the original cylinder in Angstroms\n"
       << "\t4) The 0 K lattice parameter in Angstroms\n";
}

vector <fit> setPotentials(string& database_file)
{
  string str;
  vector <fit> fits;
  ifstream fin(database_file.c_str());
  checkFileStream(fin, database_file);

  for (int i = 0; i < 6; ++i) {getline(fin,str);} // get the comment lines
  while (getline(fin, str))
  {
    fit temp;
    stringstream ss(str);
    if (!(ss >> temp.name >> temp.T1 >> temp.lin_y_int >> temp.lin_slope >> temp.T2
             >> temp.parabolic_int >> temp.parabolic_lin >> temp.parabolic_term))
    {
      cout << "Error: Corrupted database file.\n";
      exit(FILE_FORMAT_ERROR);
    }
    fits.push_back(temp);
  }

  return fits;
}

void listPotentials(const vector <fit>& fits)
{
  cout << "There are " << fits.size() << " fits:\n";
  for (unsigned int i = 0; i < fits.size(); ++i)
  {
    cout << "  " << i + 1 << " - " << fits[i].name << endl;
  }
}

double latticeParam(double T, int potential = 0)
{
  string db_file = "/home/jarinf/projects/scripts/lattice_params.db", str;
  vector <fit> fits;
  // A is the y intercept, B is the linear coefficient, C is the parabolic
  // coefficient
  double A, B, C;

  ifstream fin(db_file.c_str());
  if (fin.fail())
  {
    cout << "Error opening database file \"" << db_file << "\"\n";
    exit(11);
  }

  // first 6 lines are blank or comments
  for (unsigned int i = 0; i < 6; ++i)
  {
    getline(fin,str);
  }

  // Read each fitted set
  while (getline(fin,str))
  {
    fit temp;
    stringstream ss(str);
    if (!(ss >> temp.name >> temp.T1 >> temp.lin_y_int >> temp.lin_slope >> temp.T2
             >> temp.parabolic_int >> temp.parabolic_lin >> temp.parabolic_term))
    {
      cout << "Read error.\n";
      exit(9);
    }
    fits.push_back(temp);
  }

  if (potential == 0)
  {
    cout << "There are " << fits.size() << " fits:\n";
    for (unsigned int i = 0; i < fits.size(); ++i)
    {
      cout << "  " << i + 1 << " - " << fits[i].name << endl;
    }
    cout << "Please enter the number of the potential you want to use: ";
    cin >> potential;
  }

  --potential;
  while (potential > fits.size() || potential < 0)
  {
    cout << "Not a valid number.  If there are additional potentials you want to use, add them to the \"" << db_file << "\" database file.\n";
    cout << "Please enter the number of the potential you want to use: ";
    cin >> potential;
    --potential;
  }

  cout << "\tUsing the " << fits[potential].name << " potential.\n";

  if (T >= 0.0 && T <= fits[potential].T1)
  {
    A = fits[potential].lin_y_int;
    B = fits[potential].lin_slope;
    C = 0.0;
  }
  else if (T > fits[potential].T1 && T <= fits[potential].T2)
  {
    A = fits[potential].parabolic_int;
    B = fits[potential].parabolic_lin;
    C = fits[potential].parabolic_term;
  }
  else
  {
    cout << "Temperature out of fitted range (0 K - " << fits[potential].T2 << " K).\n";
    exit(10); // We don't want to continue with execution if we're out of range.
  }

  return A + B * T + C * T * T;
}

void parseInput(string& filename)
{
  string str;
  ifstream fin(filename.c_str());
  checkFileStream(fin, filename);

  getline(fin, str);
  stringstream ss(str);
  if (!(ss >> input.datafile >> input.temperature >> input.height >> input.a0))
  {
    showInputFileHelp();
  }

  cout << "Input parameters:"
       << "\n\tDatafile: " << input.datafile
       << "\n\tSimulation temperature: " << input.temperature
       << "\n\tCylinder thickness: " << input.height
       << "\n\t0 K lattice parameter: " << input.a0 << endl;
}

int main(int argc, char **argv)
{
  string input_file, output_file, database_file, structure;
  vector <fit> fits;
  try
  {
    cxxopts::Options options(argv[0], "Calculates the area of a grain.");
    options
      .positional_help("File")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Input file", cxxopts::value<string>(input_file), "file")
        ("o,output", "Output filename", cxxopts::value<string>(output_file)->default_value("area_data.txt"), "filename")
        ("p,potential", "The potential to use as specified in the database file", cxxopts::value<int>())
        ("c,crystal-structure", "Crystal structure of the material - simple cubic (sc), face-centered cubic (fcc), or body-centered cubic (bcc)", cxxopts::value<string>(structure)->default_value("fcc"))
        ("show-lattice-parameter", "Shows the calculated lattice parameter at the given temperature.  Can only be used with a database file and a potential.", cxxopts::value<double>(), "temperature")
        ("list-fits", "List the names of the fits given in the database file")
        ("potential-file", "Full path of the fitted lattice parameter database to use", cxxopts::value<string>(database_file)->default_value((string)(getenv("HOME")) + "/projects/scripts/lattice_params.db"), "db_file")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("file") == 0)
    {
      cout << options.help() << endl << endl;
      showInputFileHelp();
    }

    if (result.count("show-lattice-parameter"))
    {
      if (!result.count("potential"))
      {
        cout << "Must include the potential to use.\n\n";
        fits = setPotentials(database_file);
        listPotentials(fits);
        return EXIT_SUCCESS;
      }
      else
      {
        double temperature = result["show-lattice-parameter"].as<double>();
        int potential = result["potential"].as<int>();
        cout << "The lattice parameter at " << temperature << " K is:\n"
             << latticeParam(temperature, potential) << " Angstroms\n";
      }
    }

    if (result.count("list-fits"))
    {
      fits = setPotentials(database_file);
      listPotentials(fits);
      return EXIT_SUCCESS;
    }

    if (result.count("file"))
    {
      parseInput(input_file);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}
/*int main(int argc, char **argv)
{
  string filename, filename2, str; // data file, junk variable.
  double T, a0, a0_0, area, scale_factor; // temperature, lattice parameter, area
  int N1_0, N2_0, N1_next, N2_next; // Number of atoms in grains 1 and 2 at times 1 and 2
  int potential;
  double t0, t1, Lz; // time at 1, time at 2, height of cylindrical grain, GBE
  // NOTE: eGB is generally taken as 1.5.  I should use a value close to what
  // was calculated in the GBE calculations.



  cout << "Input parameters:"
       << "\n\tSimulation temperature: " << T
       << "\n\tCylinder thickness: " << Lz
       << "\n\t0 K lattice parameter: " << a0_0
       << "\n\tLattice parameter: " << a0 << endl;

  scale_factor = a0 / a0_0; // Important to account for expansion of the lattice
  Lz *= scale_factor;

  // open up the file for reading
  ifstream fin(filename.c_str());
  if (fin.fail())
  {
    cout << "Error: unable to open file " << filename << endl;
    return 1;
  }

  ofstream fout("area_data.txt");
  if (fout.fail())
  {
    cout << "Error: unable to open file area_data.txt.\n";
    return 1;
  }
  double n_extra;
  n_extra = exp(-1 / (8.617E-5 * T)) / (exp(-1 / (8.617E-5 * T)) + 1);
  fout << "# This is the area data for T = " << T << " K [time(ps) area(Angstroms^2)]\n";

  getline(fin, str); // get the comment line
  fin >> t0 >> N1_0 >> N2_0;
  if (t0 != 0)
  {
    cout << "Something is wrong in the input file - Perhaps a missing entry?\n"
         << "This script requires that the second line in the data file contain the timestep 0 information!\n";
    return 2;
  }
  fout << "# Note that at this temperature a possiblity of " << setprecision(0) << fixed << n_extra * (N1_0 + N2_0) << " atoms will be misassigned.\n";

  fin >> str >> str >> str; // ignore the minimization step

  // Output the converted values.  NOTE: This is only for fcc materials!
  fout << setprecision(0) << t0 * 0.002 << " " << setprecision(2) << N1_0 * a0 * a0 * a0 / (4 * Lz) << " " << N2_0 * a0 * a0 * a0 / (4 * Lz) << endl;
  while (fin >> t1 >> N1_next >> N2_next)
  {
    if (t1 < t0)
    {
      cout << "Error: data file corrupted.  t0 = " << t0 << " >= t1 = " << t1 << endl;
      return 3;
    }
    area = N1_next * a0 * a0 * a0 / (4 * Lz);

    fout << setprecision(0) << t1 * 0.002 << " " <<  setprecision(2) << area;

    area = N2_next * a0 * a0 * a0 / (4 * Lz);

    fout << " " << area << endl;
  }
  // Close the file stream
  fin.close();
  fout.close();

  return 0;
}*/

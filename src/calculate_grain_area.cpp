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
  string datafile, structure;
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
       << "\t1) The filename containing the number of atoms in each grain, formatted\n\t   as output by find_grains\n"
       << "\t2) The temperature in Kelvin\n"
       << "\t3) The height of the original cylinder in Angstroms\n"
       << "\t4) The 0 K lattice parameter in Angstroms\n"
       << "\t5) The crystal structure (must be one of sc (simple cubic), bcc\n\t   (body-centered cubic) or fcc (face-centered cubic)\n)";
}

fit promptForPotential(const vector <fit>& fits)
{
  int potential;

  cout << "Please enter the number of the potential you would like to use: ";
  cin  >> potential;

  --potential;
  while (potential > fits.size() || potential < 0)
  {
    cout << "Not a valid number.  If there are additional potentials you want to use, add them to the database file.\n";
    cout << "Please enter the number of the potential you want to use: ";
    cin >> potential;
    --potential;
  }
  cout << "Using the " << fits[potential].name << " potential.\n";

  return fits[potential];
}

double estimateAtomFluctuations(const double T)
{
  // Based on the Boltzmann probability distribution.  Only two allowed states,
  // where the "activated" state is a positive fluctuation (atom assigned to
  // shrinking grain), and the "non-activated" state is when the atom is assigned
  // to the exterior grain.
  return exp(-1 / (8.617E-5 * T)) / (exp(-1 / (8.617E-5 * T)) + 1);
}

vector <fit> setPotentials(const string& database_file)
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

double latticeParam(const double T, const fit& lattice_fit)
{
  // A is the y intercept, B is the linear coefficient, C is the parabolic
  // coefficient
  double A, B, C;

  if (T >= 0.0 && T <= lattice_fit.T1)
  {
    A = lattice_fit.lin_y_int;
    B = lattice_fit.lin_slope;
    C = 0.0;
  }
  else if (T > lattice_fit.T1 && T <= lattice_fit.T2)
  {
    A = lattice_fit.parabolic_int;
    B = lattice_fit.parabolic_lin;
    C = lattice_fit.parabolic_term;
  }
  else
  {
    cout << "Temperature out of fitted range (0 K - " << lattice_fit.T2 << " K).\n";
    exit(10); // We don't want to continue with execution if we're out of range.
  }

  return A + B * T + C * T * T;
}

void parseInput(const string& filename)
{
  string str;
  ifstream fin(filename.c_str());
  checkFileStream(fin, filename);

  getline(fin, str);
  stringstream ss(str);
  if (!(ss >> input.datafile >> input.temperature >> input.height >> input.a0 >> input.structure))
  {
    showInputFileHelp();
  }

  cout << "Input parameters:"
       << "\n\tDatafile: " << input.datafile
       << "\n\tSimulation temperature: " << input.temperature
       << "\n\tCylinder thickness: " << input.height
       << "\n\t0 K lattice parameter: " << input.a0
       << "\n\tCrystal structure: " << input.structure << endl;
}

void calculateGrainArea(const fit& lattice_fit, const string& output_file)
{
  string str;
  double t0, t1, structure_factor; // Times
  int N1_0, N2_0, N1_next, N2_next; // atom numbers
  double scale_factor = input.a0 / latticeParam(input.temperature, lattice_fit);
  double Lz = input.height * scale_factor;
  double n_extra = estimateAtomFluctuations(input.temperature);

  ifstream fin(input.datafile.c_str());
  checkFileStream(fin, input.datafile);

  ofstream fout(output_file.c_str());
  checkFileStream(fout, output_file);

  fout << "# This is the area data for T = " << input.temperature << " K [time (ps) area (Angstroms^2)]\n";

  getline(fin, str);
  fin >> t0 >> N1_0 >> N2_0;

  if (t0 != 0.0)
  {
    cout << "This script requires that the second line in the data file contain the timestep 0 information.\n";
    exit(FILE_FORMAT_ERROR);
  }

  fout << "# Note that at this temperature a possibility of " << setprecision(0) << fixed << n_extra * (N1_0 + N2_0) << " atoms may be misassigned.\n";

  fin >> str >> str >> str; // ignore the minimization step

  structure_factor = input.a0 * input.a0 * input.a0;
  if (input.structure.compare("sc") == 0) {structure_factor /= Lz;}
  if (input.structure.compare("fcc") == 0) {structure_factor /= 4 * Lz;}
  if (input.structure.compare("bcc") == 0) {structure_factor /= 2 * Lz;}

  fout << setprecision(0) << t0 * 0.002 << " " << setprecision(2) << N1_0 * structure_factor << " " << N2_0 * structure_factor << endl;
  while (fin >> t1 >> N1_next >> N2_next)
  {
    if (t1 < t0)
    {
      cout << "Error: data file corrupted.  t0 = " << t0 << " >= t1 = " << t1 << endl;
      exit(FILE_FORMAT_ERROR);
    }

    fout << setprecision(0) << t1 * 0.002 << " " << setprecision(2)
         << N1_next * structure_factor << " "
         << N2_next * structure_factor << endl;
  }

  fin.close();
  fout.close();
}

int main(int argc, char **argv)
{
  string input_file, output_file, database_file, structure;
  vector <fit> fits;
  fit lattice_fit;
  int potential;
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
        ("p,potential", "The potential to use as specified in the database file", cxxopts::value<int>(potential)->default_value("0"))
        ("show-lattice-parameter", "Shows the calculated lattice parameter at the given temperature.  Can only be used with a database file and a potential.", cxxopts::value<double>(), "temperature")
        ("estimate-fluctuations", "Estimates the percentage of misassigned atoms based on the temperature given using the Boltzmann probability", cxxopts::value<double>(), "temperature")
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

    fits = setPotentials(database_file);

    if (!result.count("potential"))
    {
      listPotentials(fits);
      lattice_fit = promptForPotential(fits);
    }
    else {lattice_fit = fits[potential - 1];}

    if (result.count("show-lattice-parameter"))
    {
      double temperature = result["show-lattice-parameter"].as<double>();
      int potential = result["potential"].as<int>();
      cout << "The lattice parameter at " << temperature << " K is:\n"
           << latticeParam(temperature, fits[--potential]) << " Angstroms\n";
      return EXIT_SUCCESS;
    }

    if (result.count("estimate-fluctuations"))
    {
      double temperature = result["estimate-fluctuatons"].as<double>();
      cout << "The estimated percentage of misassigned atoms at " << temperature
           << " K is:\n" << estimateAtomFluctuations(temperature) << endl;
      return EXIT_SUCCESS;
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
      calculateGrainArea(lattice_fit, output_file);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

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
#include "verifyNewFile.h"

using namespace std;

bool extrapolate = false;

struct fit {
  string name;
  double x3, x2, x1, y3, y2, y1, xy, xy2, x2y, z0;
  double min_conc, max_conc;
  double min_T, max_T;
};

struct inputData {
  string datafile, structure;
  double temperature, height, a0, concentration;
} input;

#define PI 3.141592653589793

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

void showInputFileHelp() {
  cout << "The input file must contain the following on one line, in order:\n"
       << "\t1) The filename containing the number of atoms in each grain, formatted\n\t   as output by find_grains\n"
       << "\t2) The temperature in Kelvin\n"
       << "\t3) The concentration of solute\n"
       << "\t4) The height of the original cylinder in Angstroms\n"
       << "\t5) The 0 K lattice parameter in Angstroms\n"
       << "\t6) The crystal structure (must be one of sc (simple cubic), bcc\n\t   (body-centered cubic) or fcc (face-centered cubic)\n";
}

fit promptForPotential(const vector <fit>& fits) {
  int potential;

  cout << "Please enter the number of the potential you would like to use: ";
  cin  >> potential;

  --potential;
  while (potential > fits.size() || potential < 0) {
    cout << "Not a valid number.  If there are additional potentials you want to use, add them to the database file.\n";
    cout << "Please enter the number of the potential you want to use: ";
    cin >> potential;
    --potential;
  }
  cout << "Using the " << fits[potential].name << " potential.\n";

  return fits[potential];
}

vector <fit> setPotentials(const string& database_file) {
  string str;
  vector <fit> fits;
  ifstream fin(database_file.c_str());
  checkFileStream(fin, database_file);

  while (getline(fin, str)) {
    if (str[0] == '#' || str.find_first_not_of("\t\n\v\f\r ") == std::string::npos) {continue;} // ignore commented lines and blank lines
    fit temp;
    stringstream ss(str);
    if (!(ss >> temp.name >> temp.min_conc >> temp.max_conc >> temp.min_T
             >> temp.max_T >> temp.x3 >> temp.x2 >> temp.x1 >> temp.y3 >> temp.y2
             >> temp.y1 >> temp.xy >> temp.xy2 >> temp.x2y >> temp.z0)) {
      cerr << "Error: Corrupted database file.\n";
      exit(FILE_FORMAT_ERROR);
    }
    if (temp.max_T <= temp.min_T) { // TODO: add negative value checks
      cerr << "Error: Corrupted database file: T_max = " << temp.max_T << " <= T_min = " << temp.min_T << ".\n";
      exit(FILE_FORMAT_ERROR);
    }
    if (temp.max_conc < temp.min_conc) {
      cerr << "Error: Corrupted database file: conc_max = " << temp.max_conc << " <= conc_min = " << temp.min_conc << ".\n";
      exit(FILE_FORMAT_ERROR);
    }
    fits.push_back(temp);
  }

  return fits;
}

void listPotentials(const vector <fit>& fits) {
  cout << "There are " << fits.size() << " fits:\n";
  for (unsigned int i = 0; i < fits.size(); ++i) {
    cout << "  " << i + 1 << " - " << fits[i].name << "\n";
  }
}

double latticeParam(const double T, const double conc, const fit& lattice_fit) {
  // A is the y intercept, B is the linear coefficient, C is the parabolic
  // coefficient
  double x3, x2, x1, y3, y2, y1, xy, xy2, x2y, z0;
  x3 = lattice_fit.x3; y3 = lattice_fit.y3;
  x2 = lattice_fit.x2; y2 = lattice_fit.y2;
  x1 = lattice_fit.x1; y1 = lattice_fit.y1;
  xy = lattice_fit.xy; xy2 = lattice_fit.xy2; x2y = lattice_fit.x2y;
  z0 = lattice_fit.z0;

  if (T > lattice_fit.max_T && extrapolate) {
    cout << "Warning: extrapolating beyond fitted temperatures.\n";
  }
  else if (T > lattice_fit.max_T || T < lattice_fit.min_T) {
    cerr << "Temperature out of fitted range (" << lattice_fit.min_T << " K - " << lattice_fit.max_T << " K).\n";
    exit(BOUNDS_ERROR); // We don't want to continue with execution if we're out of range.
  }

  if (conc > lattice_fit.max_conc || conc < lattice_fit.min_conc) {
    cerr << "Concentration out of fitted range (" << lattice_fit.min_conc << " - " << lattice_fit.max_conc << ").\n";
    exit(BOUNDS_ERROR);
  }

  double conc_poly = x3 * conc * conc * conc + x2 * conc * conc + x1 * conc;
  double T_poly = y3 * T * T * T + y2 * T * T + y1 * T;
  double cross_terms = xy * conc * T + xy2 * conc * T * T + x2y * conc * conc * T;

  return conc_poly + T_poly + cross_terms + z0;
}

void parseInput(const string& filename, const fit& lattice_fit) {
  string str;
  ifstream fin(filename.c_str());
  checkFileStream(fin, filename);

  getline(fin, str);
  stringstream ss(str);
  if (!(ss >> input.datafile >> input.temperature >> input.concentration >> input.height >> input.a0 >> input.structure)) {
    showInputFileHelp();
  }

  if (input.concentration > 1) {
    cout << "Assuming concentration given in percent, converting to fraction: " << input.concentration << " --> " << input.concentration / 100.0 << "\n";
    input.concentration /= 100.0;
  }

  cout << "Input parameters:"
       << "\n\tDatafile: " << input.datafile
       << "\n\tSimulation temperature: " << input.temperature
       << "\n\tConcentration of solute: " << input.concentration
       << "\n\tCylinder thickness: " << input.height
       << "\n\t0 K lattice parameter: " << input.a0
       << "\n\tCrystal structure: " << input.structure
       << "\n\tFitted to: " << lattice_fit.name << "\n";
}

void calculateGrainArea(const fit& lattice_fit, const string& output_file, const double& dt) {
  string str;
  double t0, t1, structure_factor; // Times
  int N1_0, N2_0, N1_next, N2_next, n_gb; // atom numbers
  double lattice_param = latticeParam(input.temperature, input.concentration, lattice_fit);
  double scale_factor = lattice_param / input.a0;
  double Lz = input.height * scale_factor;

  ifstream fin(input.datafile.c_str());
  checkFileStream(fin, input.datafile);

  VerifyNewFile verify(output_file);
  ofstream fout(verify.validNewFile().c_str());
  checkFileStream(fout, verify.validNewFile());

  fout << "# This is the area data for T = " << input.temperature << " K [time (ps) area (Angstroms^2)]\n";

  getline(fin, str);

  size_t first_bracket = str.find('[',0);
  size_t last_bracket = str.rfind(']');
  string elements = str.substr(first_bracket, last_bracket-first_bracket);

  string::difference_type n_elements = std::count(elements.begin(), elements.end(), ',') + 1;


  fin >> t0 >> N1_0 >> N2_0;
  if (n_elements == 4) {
    fin >> n_gb;
  }

  structure_factor = lattice_param * lattice_param * lattice_param; // a0^3
  if (input.structure.compare("sc") == 0) {structure_factor /= Lz;}
  else if (input.structure.compare("fcc") == 0) {structure_factor /= 4 * Lz;}
  else if (input.structure.compare("bcc") == 0) {structure_factor /= 2 * Lz;}
  else {
    cerr << "Error: unknown input structure: " << input.structure << "\n";
    exit(INPUT_FORMAT_ERROR);
  }

  fout << fixed;

  fout << setprecision(0) << t0 * dt << " " << setprecision(2) << N1_0 * structure_factor << " " << N2_0 * structure_factor << "\n";
  while (fin >> t1 >> N1_next >> N2_next) {
    if (n_elements == 4) {
      fin >> n_gb;
    }
    if (t1 < t0) {
      cerr << "Error: data file corrupted.  t0 = " << t0 << " >= t1 = " << t1 << "\n";
      exit(FILE_FORMAT_ERROR);
    }

    fout << setprecision(0) << t1 * dt << " " << setprecision(2)
         << N1_next * structure_factor << " "
         << N2_next * structure_factor << "\n";
  }

  fin.close();
  fout.close();
}

int main(int argc, char **argv) {
  string input_file, output_file, database_file, structure;
  vector <fit> fits;
  fit lattice_fit;
  int potential;
  double dt;
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
        ("show-lattice-parameter", "Shows the calculated lattice parameter at the given temperature.  Can only be used with a database file and a potential. Assumes concentration is 0", cxxopts::value<double>(), "temperature")
        ("list-fits", "List the names of the fits given in the database file")
        ("dt", "Specify the timestep used in the simulations", cxxopts::value<double>(dt)->default_value("0.002"))
        ("potential-file", "Full path of the fitted lattice parameter database to use", cxxopts::value<string>(database_file)->default_value((string)(getenv("HOME")) + "/projects/scripts/lattice_params.db"), "db_file")
        ("e,extrapolate", "Proceed with calculations even if temperature is above fitted range.")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("list-fits")) {
      fits = setPotentials(database_file);
      listPotentials(fits);
      return EXIT_SUCCESS;
    }

    if (result.count("help") || result.count("file") == 0) {
      cout << options.help() << "\n" << "\n";
      showInputFileHelp();
      return EXIT_SUCCESS;
    }

    if (result.count("extrapolate")) {extrapolate = true;}

    fits = setPotentials(database_file);

    if (!result.count("potential")) {
      listPotentials(fits);
      lattice_fit = promptForPotential(fits);
    }
    else
    {
      if (potential > fits.size()) {
        cout << "Invalid potential\n";
        return BOUNDS_ERROR;
      }
      else {lattice_fit = fits[potential - 1];}
    }

    if (result.count("show-lattice-parameter")) {
      double temperature = result["show-lattice-parameter"].as<double>();
      double concentration = 0.0;
      int potential = result["potential"].as<int>();
      cout << "The lattice parameter at " << temperature << " K is:\n"
           << latticeParam(temperature, concentration, fits[potential - 1]) << " Angstroms\n";
      return EXIT_SUCCESS;
    }

    if (result.count("file")) {
      parseInput(input_file, lattice_fit);
      calculateGrainArea(lattice_fit, output_file, dt);
    }
  }
  catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << "\n";
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

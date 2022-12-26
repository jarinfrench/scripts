#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <random> //
#include <unistd.h>
#include <numeric> // iota
#include <algorithm> // sort, all_of
#include <omp.h> // for omp_get_cancellation - requires the -fopenmp option during compiling

#include <cxxopts.hpp>
#include "error_code_defines.h"
#include "position.h"

using namespace std;

#define PI 3.141592653589793

struct InputData
{
  // Directly from the input file
  // number of grid points in x, y, RNG seed, number of order parameters,
  // number of timesteps, and number of timesteps to run before outputting a file
  unsigned int grid_x, grid_y, seed, n_eta, numsteps, nskip;
  double dx, dt; // spatial step, time step
  string mobility_file; // file containing matrix of mobility values
  string ic_type, ic_file;

  // Calculated from the input file values
  double max_x, max_y; // maximum x and y coordinate (dx * grid)
  unsigned int ncheck; // number of timesteps to run before giving the user a progress update

  // Constants
  unsigned int EQUILIBRIUM = 500; // number of timesteps for the sharp interface to equilibrate

  void calculateValues()
  {
    int denominator = 1000; // hard coded value indicating that we want updates every 1 tenth of a percent
    max_x = grid_x * dx;
    max_y = grid_y * dx;

    ncheck = numsteps / denominator;
    while (ncheck == 0) // will be 0 anytime numsteps < denominator
    {
      denominator /= 2;
      ncheck = numsteps / denominator;
    }
  }
};

struct Field
{
  vector <vector <double> > values;
  Position center; // original center point
  bool is_active; // whether or not the field is active

  // constructor
  Field(unsigned int grid_x, unsigned int grid_y)
  {
    values.assign(grid_x, vector <double> (grid_y, 0.0));
    is_active = true;
  }
};

void printInputFileHelp()
{
  cout << "The input file should contain the following items:\n"
       << "  Line 1 is a comment line, and is ignored. The following lines with\n"
       << "    a * are required\n"
       << "  *GRID_X = <integer number of grid points in the x direction>\n"
       << "  *GRID_Y = <integer number of grid points in the y direction>\n"
       << "  *DX = <floating point number indicating the step size in both\n"
       << "    directions>\n"
       << "  *SEED = <integer number specifying the seed for the random number\n"
       << "    generator>\n"
       << "  *N_ETA = <integer number of order parameters to use in the\n"
       << "    simulation>\n"
       << "  *NUMSTEPS = <integer number of steps to run the simulation for>\n"
       << "  *NSKIP = <integer number of when to print structure files>\n"
       << "  *DT = <floating point number indicating the timestep size>\n"
       << "  IC = <string of either 'random', 'fluid', 'hex', 'centroid', or\n"
       << "   'input'>\n"
       << "    if IC = input or centroid is specified, the next line _must_ be:\n"
       << "    IC_FILE = <structure input file>\n"
       << "    where the structure file for input contains the NxN grid of (integer)\n"
       << "    order parameters, and for centroid contains N_ETA lines containing\n"
       << "    the fractional coordinates of the centroids.\n"
       << "    if not specified, the random initial condition is used\n"
       << "  mobility_file = <txt file containing the matrix of boundary\n"
       << "    mobilities>\n"
       << "    if not specified, mobility is assumed isotropic with a value of 1.0\n\n"
       << "Note that these lines can be specified in any particular order.\n";
}

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  // round to the largest absolute distance from zero
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

InputData parseInputFile(const string& file)
{
  string str, str2;
  map <string, bool> required_params_found = {
    {"GRID_X", false},
    {"GRID_Y", false},
    {"DX", false},
    {"SEED", false},
    {"N_ETA", false},
    {"NUMSTEPS", false},
    {"NSKIP", false},
    {"DT", false}
  };
  InputData input;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  getline(fin, str); // first line is a comment line
  while (getline(fin, str))
  {
    stringstream ss(str);
    ss >> str2; // get the variable name

    if (str2.compare("GRID_X") == 0) {ss >> str2 >> input.grid_x; required_params_found["GRID_X"] = true;}
    else if (str2.compare("GRID_Y") == 0) {ss >> str2 >> input.grid_y; required_params_found["GRID_Y"] = true;}
    else if (str2.compare("DX") == 0) {ss >> str2 >> input.dx; required_params_found["DX"] = true;}
    else if (str2.compare("SEED") == 0) {ss >> str2 >> input.seed; required_params_found["SEED"] = true;}
    else if (str2.compare("N_ETA") == 0) {ss >> str2 >> input.n_eta; required_params_found["N_ETA"] = true;}
    else if (str2.compare("NUMSTEPS") == 0) {ss >> str2 >> input.numsteps; required_params_found["NUMSTEPS"] = true;}
    else if (str2.compare("NSKIP") == 0) {ss >> str2 >> input.nskip; required_params_found["NSKIP"] = true;}
    else if (str2.compare("DT") == 0) {ss >> str2 >> input.dt; required_params_found["DT"] = true;}
    else if (str2.compare("IC") == 0)
    {
      ss >> str2 >> input.ic_type;
      if (input.ic_type.compare("input") == 0 || input.ic_type.compare("centroid") == 0)
      {
        getline(fin, str);
        stringstream ss2(str);
        ss2 >> str2;

        if (str2.compare("IC_FILE") == 0) {ss2 >> str2 >> input.ic_file;}
        else
        {
          cerr << "Error: Initial condition specified as input, but input file not given on following line.\n";
          exit(INPUT_FORMAT_ERROR);
        }
      }
    }
    else if (str2.compare("mobility_file") == 0) {ss >> str2 >> input.mobility_file;}
    else {cout << "Unrecognized entry in input file: " << str2 << "\n";}
  }

  fin.close();

  if (!all_of(required_params_found.begin(), required_params_found.end(), [](pair <string, bool> a) {return a.second;}))
  {
    cerr << "Error: Missing required parameter in input file.\n"
         << "Missing the following parameters: \n";
    for (map <string, bool>::iterator it = required_params_found.begin(); it != required_params_found.end(); ++it)
    {
      if (!(it->second)) {cout << it->first << "\n";}
    }

    exit(INPUT_FORMAT_ERROR);
  }

  input.calculateValues();

  return input;
}

vector <vector <double> > getMobilityCoefficients(const string& file, const unsigned int& n_eta)
{
  vector <vector <double> > _L (n_eta, vector <double> (n_eta, 0.0));
  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  for (unsigned int i = 0; i < n_eta; ++i)
  {
    for (unsigned int j = 0; j < n_eta; ++j)
    {
      fin >> _L[i][j];
    }
  }
  fin.close();

  return _L;
}

void assignGridToEtas(vector <Field>& etas, const InputData& input)
{
  for (unsigned int i = 0; i < input.grid_x; ++i)
  {
    for (unsigned int j = 0; j < input.grid_y; ++j)
    {
      int min_id = -1; // invalid id number
      double min = 1e12; // some arbitrarily large number

      for (unsigned int k = 0; k < input.n_eta; ++k)
      {
        // find the distance to each point
        double x = i * input.dx - etas[k].center.getX();
        double y = j * input.dx - etas[k].center.getY();

        x = x - anInt(x / input.max_x) * input.max_x;
        y = y - anInt(y / input.max_y) * input.max_y;

        double distance = (x * x) + (y * y);

        if (distance < min)
        {
          min = distance;
          min_id = k;
        }
      }

      etas[min_id].values[i][j] = 1.0;
    }
  }
}

vector <Field> generateInitialCondition(const InputData& input)
{
  vector <Field> etas (input.n_eta, Field(input.grid_x, input.grid_y));

  mt19937 rng(input.seed); // seed the random number generator
  // generate the distribution
  if (input.ic_type.compare("random") == 0)
  {
    uniform_real_distribution<double> dist(0.0, 1.0);
    for (unsigned int i = 0; i < input.n_eta; ++i)
    {
      int x = dist(rng) * input.max_x;
      int y = dist(rng) * input.max_y;
      etas[i].center = Position(x, y);
    }
    assignGridToEtas(etas, input);
  }
  else if (input.ic_type.compare("hex") == 0)
  {
    unsigned int sqrt_num_centroids = sqrt(input.n_eta);
    if (sqrt_num_centroids != sqrt(input.n_eta))
    {
      cerr << "The number of order parameters must be a perfect square.\n";
      exit(INPUT_FORMAT_ERROR);
    }
    uniform_real_distribution<double> dist(-1.0, 1.0);
    unsigned int num_centroids = 0;

    for (unsigned int i = 0; i < sqrt_num_centroids; ++i)
    {
      for (unsigned int j = 0; j < sqrt_num_centroids; ++j)
      {
        for (unsigned int k = 0; k < 1; ++k)
        {
          double x = ((double)(i) / sqrt_num_centroids + (0.5 / sqrt_num_centroids * (j % 2)) + 0.5 / sqrt_num_centroids + dist(rng) * 0.1) * input.max_x;
          double y = ((double)(j) / sqrt_num_centroids + (0.5 / sqrt_num_centroids * (k % 2)) + dist(rng) * 0.1) * input.max_y;
          etas[num_centroids].center = Position(x,y);
          ++num_centroids;
        }
      }
    }
    assignGridToEtas(etas, input);
  }
  else if (input.ic_type.compare("centroid") == 0)
  {
    // get the centroids from the specified file
    string str;
    int num_grains = 0;
    ifstream fin(input.ic_file.c_str());
    checkFileStream(fin, input.ic_file);

    getline(fin, str); // first two lines are comment lines
    getline(fin, str);

    while (num_grains < input.n_eta)
    {
      double x, y;
      getline(fin, str);
      stringstream ss(str);
      ss >> x >> y; // get the x and y positions of the centroids

      etas[num_grains++].center = Position(input.max_x * x, input.max_y * y);;
    }
    fin.close();

    assignGridToEtas(etas, input);
  }
  else if (input.ic_type.compare("fluid") == 0)
  {
    uniform_real_distribution<double> dist(-1.0, 1.0);
    for (unsigned int i = 0; i < input.grid_x; ++i)
    {
      for (unsigned int j = 0; j < input.grid_y; ++j)
      {
        for (unsigned int k = 0; k < input.n_eta; ++k)
        {
          etas[k].values[i][j] = dist(rng) * 0.001; // we are creating a slightly perturbed fluid here
        }
      }
    }
  }
  else if (input.ic_type.compare("input") == 0)
  {
    // read the file into the etas
    string str;
    ifstream fin(input.ic_file.c_str());
    checkFileStream(fin, input.ic_file);

    unsigned int row = 0, col, eta_index;

    while (getline(fin, str))
    {
      col = 0;
      stringstream ss(str);
      while (ss >> eta_index)
      {
        etas[eta_index - 1].values[row][col] = 1.0;
        ++col;
      }
      ++row;
    }

    fin.close();
  }

  return etas;
}

void runSimulation(vector <Field> etas_old, vector <Field> etas,
                   vector <Field> laplacian, vector <Field> dEta_dt,
                   const vector <vector <double> >& _L, const vector <vector <double> >& _kappa,
                   const InputData& input, const string& gb_lengths_file,
                   const string& grain_sizes_file, const string& multi_junctions_file)
{
  ofstream fout_gb_length(gb_lengths_file.c_str());
  checkFileStream(fout_gb_length, gb_lengths_file);

  ofstream fout_grain_size(grain_sizes_file.c_str());
  checkFileStream(fout_grain_size, grain_sizes_file);

  ofstream fout_multi_junctions(multi_junctions_file.c_str());
  checkFileStream(fout_multi_junctions, multi_junctions_file);

  if (gb_lengths_file.compare("none") == 0) {fout_gb_length.close();}
  else
  {
    fout_gb_length << "# ";
    for (unsigned int i = 0; i < input.n_eta - 1; ++i)
    {
      for (unsigned int j = i + 1; j < input.n_eta; ++j)
      {
        fout_gb_length << "gr" << i << "gr" << j << "_gb ";
      }
    }
    fout_gb_length << "total (total2 if different)\n";
  }

  if (grain_sizes_file.compare("none") == 0) {fout_grain_size.close();}
  else
  {
    fout_grain_size << "# Step ";
    for (unsigned int i = 0; i < input.n_eta; ++i)
    {
      fout_grain_size << "gr" << i << " ";
    }
    fout_grain_size << "\n";
  }
  if (multi_junctions_file.compare("none)") == 0) {fout_multi_junctions.close();}
  else
  {
    fout_multi_junctions << "# Step junction_position active_etas junction_neighbors junction_angles\n";
  }

  /* fill in the rest of the simulation commands here*/




  if (gb_lengths_file.compare("none") != 0) {fout_gb_length.close();}
  if (grain_sizes_file.compare("none") != 0) {fout_grain_size.close();}
  if (multi_junctions_file.compare("none") != 0) {fout_multi_junctions.close();}
}

int main(int argc, char** argv)
{
  string input_file, gb_lengths_file, grain_sizes_file, multi_junctions_file;

  try
  {
    cxxopts::Options options(argv[0], "Phase field simulation script");
    options
      .positional_help("file")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Input file", cxxopts::value<string>(input_file), "file")
        ("l,gb-length", "Print the grain boundary lengths to a file", cxxopts::value<string>(gb_lengths_file)->implicit_value("gb_lengths.txt")->default_value("none"), "file")
        ("a,gb-area", "Print the grain boundary areas to a file", cxxopts::value<string>(grain_sizes_file)->implicit_value("grain_areas.txt")->default_value("none"), "file")
        ("m,multi-junction", "Print info about the multi-junctions in the system", cxxopts::value<string>(multi_junctions_file)->implicit_value("multi_junctions.txt")->default_value("none"), "file")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("file") == 0)
    {
      cout << options.help();
      printInputFileHelp();

      return EXIT_SUCCESS;
    }

    if (result.count("file"))
    {
      InputData input = parseInputFile(input_file);
      vector <vector <double> > _L = getMobilityCoefficients(input.mobility_file, input.n_eta);
      vector <vector <double> > _kappa (input.n_eta, vector <double> (input.n_eta, 2.0)); // constant values for energy
      vector <Field> etas_old = generateInitialCondition(input);
      // set up the other field parameters
      vector <Field> etas (input.n_eta, Field(input.grid_x, input.grid_y));
      vector <Field> laplacian (input.n_eta, Field(input.grid_x, input.grid_y));
      vector <Field> dEta_dt (input.n_eta, Field(input.grid_x, input.grid_y));

      runSimulation(etas_old, etas, laplacian, dEta_dt, _L, _kappa, input, gb_lengths_file, grain_sizes_file, multi_junctions_file);

    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cerr << "Error parsing options: " << e.what() << "\n";
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

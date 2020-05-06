#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <cxxopts.hpp>

#include "error_code_defines.h"
#include "position.h"

using namespace std;

#define PI 3.141592653589793

struct inputData
{
  string unit_cell_file;
  double a, b, c, alpha, beta, gamma;
  double cos_alpha, cos_beta, cos_gamma;
  int h, k, l;
};

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

struct compareMQRatio
{
  bool operator()(const pair <int, int>& p1, const pair <int, int>& p2) const
  {
    if (p1.first == p2.first && p1.second == p2.second) {return false;} // same value

    double p1_ratio = (double)(p1.first) / (double)(p1.second);
    double p2_ratio = (double)(p2.first) / (double)(p2.second);
    return p1_ratio < p2_ratio;
  }
};

void printInputFileHelp()
{
  cout << "The input file should contain the following items on one line, in order:\n"
       << "  1. The name of the file containing the list of atoms in a unit cell\n"
       << "  2. The a lattice constant\n"
       << "  3. The b lattice constant\n"
       << "  4. The c lattice constant\n"
       << "  5. The alpha angle (in degrees)\n"
       << "  6. The beta angle (in degrees)\n"
       << "  7. The gamma angle (in degrees)\n";
}

double calculateDeterminant(const vector <vector <double> >& m)
{
  if (m.size() != 3 || m[0].size() != 3 || m[1].size() != 3 || m[2].size() != 3)
  {
    cerr << "Error: only calculating the determinant of a 3x3 matrix!\n";
    exit(ERROR_CODE_NOT_DEFINED);
  }

  double a = m[0][0], b = m[0][1], c = m[0][2]; // cofactors
  double det1, det2, det3;

  det1 = m[1][1] * m[2][2] - m[1][2] * m[2][1];
  det2 = m[1][0] * m[2][2] - m[1][2] * m[2][0];
  det3 = m[1][0] * m[2][1] - m[1][1] * m[2][0];
  return a * det1 - b * det2 + c * det3;
}

pair <int, int> decimal2Fraction (const double& target)
{
  if (target >= 1.0)
  {
    cerr << "Error: passed in target is >= 1.0: target = " << target << "\n";
    exit(INPUT_FORMAT_ERROR);
  }

  // Taken from https://gist.github.com/mikeando/7073d62385a34a61a6f7, based on the
  // question given by https://stackoverflow.com/questions/26643695/converting-decimal-to-fraction-c
  int Aprev[2] = {1, 0};
  int Bprev[2] = {0, 1};
  int denominator, numerator;

  double x = target;
  double eps = 1.0;
  double approx;
  unsigned int num_iters = 0;

  while (eps > 1.0e-8)
  {
    // Original algorithm puts x between 0 and 1 here, we skip that part
    int n = floor(x);
    x -= n;
    x = 1.0 / x;

    // Note that the [AB]prev[1] values in the original algorithm are multiplied by n,
    // which is always zero in our case.
    denominator = Aprev[0] + n * Aprev[1];
    Aprev[0] = Aprev[1];
    Aprev[1] = denominator;

    numerator = Bprev[0] + n * Bprev[1];
    Bprev[0] = Bprev[1];
    Bprev[1] = numerator;

    approx = (double)(numerator) / (double)(denominator);
    eps = abs(approx - target);
    ++num_iters;

    if (num_iters > 10000)
    {
      cout << "Fraction not determined after 10000 iterations:\n"
           << "Currently value: " << approx
           << "\nError: " << eps << "\n";
      break;
    }
  }

  return make_pair(numerator, denominator);
}

void printCitationInfo()
{
  cout << "Fan, Q. Journal of Applied Crystallography 49(5) (2016) 1454-1458\n";
}

inputData parseInputFile(const string& file, const unsigned int& verbosity)
{
  inputData input;
  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  if (!(fin >> input.unit_cell_file >> input.a >> input.b >> input.c >> input.alpha >> input.beta >> input.gamma))
  {
    cerr << "Error: not enough input values.\n";
    printInputFileHelp();
  }

  if (verbosity >= 1)
  {
    cout << "Input data:"
         << "\n  Unit cell file: " << input.unit_cell_file
         << "\n  Lattice parameters: a = " << input.a << ", b = " << input.b << ", c = " << input.c
         << "\n  Lattice angles: alpha = " << input.alpha << ", beta = " << input.beta << ", gamma = " << input.gamma;
  }
  // convert to radians
  input.alpha *= PI / 180.0;
  input.beta *= PI / 180.0;
  input.gamma *= PI / 180.0;

  input.cos_alpha = cos(input.alpha);
  input.cos_beta = cos(input.beta);
  input.cos_gamma = cos(input.gamma);

  fin.close();
  return input;
}

vector <Position> parseUnitCellFile(const string& file, const unsigned int& verbosity = 0)
{
  double x, y, z, dummy; // positions
  vector <Position> unit_cell_positions;
  string str;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);


  while (getline(fin, str))
  {
    stringstream ss(str);
    if (!(ss >> x >> y >> z))
    {
      cerr << "Error: Not enough elements in " << file << " for line \""
           << str << "\"\n"
           << "Each line must contain the x, y, and z position of the atoms in the unit cell.\n";
      exit(FILE_FORMAT_ERROR);
    }
    if (!(ss >> dummy)) {unit_cell_positions.push_back(Position(x, y, z));}
    else
    {
      cerr << "Error: Too many elements in " << file << " for line \""
           << str << "\"\n";
      exit(FILE_FORMAT_ERROR);
    }
  }

  if (verbosity >= 1)
  {
    cout << "Unit cell atom positions:\n";
    for (size_t i = 0; i < unit_cell_positions.size(); ++i)
    {
      cout << i + 1 << ". " << unit_cell_positions[i] << "\n";
    }
    cout << "\n";
  }

  return unit_cell_positions;
}

double calculateUnitCellVolume(const inputData& input, const unsigned int& verbosity = 0)
{
  double abc = input.a * input.b * input.c;
  double cos_mult = 2 * input.cos_alpha * input.cos_beta * input.cos_gamma;
  double cos_alpha_sq = input.cos_alpha * input.cos_alpha;
  double cos_beta_sq = input.cos_beta * input.cos_beta;
  double cos_gamma_sq = input.cos_gamma * input.cos_gamma;
  // Eq. 3
  double unit_cell_vol = abc * sqrt(1 - cos_alpha_sq - cos_beta_sq - cos_gamma_sq + cos_mult);

  if (verbosity >= 1)
  {
    cout << "Volume of unit cell: " << unit_cell_vol << "\n";
  }

  return unit_cell_vol;
}

double calculateInterplanarSpacing(const inputData& input, const unsigned int& verbosity = 0)
{
  vector <vector <double> > m1 (3, vector <double> (3,0.0));
  vector <vector <double> > m2 (3, vector <double> (3,0.0));
  vector <vector <double> > m3 (3, vector <double> (3,0.0));
  vector <vector <double> > m4 (3, vector <double> (3,0.0));
  double a = input.h / input.a, b = input.k / input.b, c = input.l / input.c;

  // Eq. 4
  m1[0][0] = input.h / input.a; m1[0][1] = input.cos_gamma; m1[0][2] = input.cos_beta;
  m1[1][0] = input.k / input.b; m1[1][1] = 1.0;             m1[1][2] = input.cos_alpha;
  m1[2][0] = input.l / input.c; m1[2][1] = input.cos_alpha; m1[2][2] = 1.0;

  m2[0][0] = 1.0;             m2[0][1] = m1[0][0]; m2[0][2] = input.cos_beta;
  m2[1][0] = input.cos_gamma; m2[1][1] = m1[1][0]; m2[1][2] = input.cos_alpha;
  m2[2][0] = input.cos_beta;  m2[2][1] = m1[2][0]; m2[2][2] = 1.0;

  m3[0][0] = 1.0;             m3[0][1] = input.cos_gamma; m3[0][2] = m1[0][0];
  m3[1][0] = input.cos_gamma; m3[1][1] = 1.0;             m3[1][2] = m1[1][0];
  m3[2][0] = input.cos_beta;  m3[2][1] = input.cos_alpha; m3[2][2] = m1[2][0];

  m4[0][0] = 1.0;             m4[0][1] = input.cos_gamma; m4[0][2] = input.cos_beta;
  m4[1][0] = input.cos_gamma; m4[1][1] = 1.0;             m4[1][2] = input.cos_alpha;
  m4[2][0] = input.cos_beta;  m4[2][1] = input.cos_alpha; m4[2][2] = 1.0;

  double det1 = calculateDeterminant(m1);
  double det2 = calculateDeterminant(m2);
  double det3 = calculateDeterminant(m3);
  double det4 = calculateDeterminant(m4);

  // Eq. 4
  double interplanar_spacing = 1.0 / sqrt((a * det1 + b * det2 + c * det3) / det4);
  if (verbosity >= 1)
  {
    cout << "Interplanar spacing: " << interplanar_spacing << "\n";
  }
  return interplanar_spacing;

}

void calculatePlanarDensities(const inputData& input,
                              const vector <Position>& positions,
                              const double& dhkl,
                              const double& v_cell,
                              const unsigned int& verbosity = 0)
{
  // calculate the position factors
  vector <double> position_factors (positions.size(), 0);
  vector <int> s_vals (positions.size(), 0), m_vals (positions.size(), 0), q_vals (positions.size(), 0);
  pair <int, int> m_q_val;
  map <pair <int, int>, int, compareMQRatio> unique_m_q_pairs_count;
  map <pair <int, int>, vector <string>, compareMQRatio> unique_m_q_pairs_vector;
  pair <int, int> key_0 = make_pair(0, 1);
  unique_m_q_pairs_count[key_0] = 0;
  for (size_t i = 0; i < positions.size(); ++i)
  {
    stringstream ss;
    ss << "P" << i + 1;
    // Eq. 5
    position_factors[i] = positions[i].getX() * input.h + positions[i].getY() * input.k + positions[i].getZ() * input.l;
    s_vals[i] = (int)(position_factors[i]);
    if (!(position_factors[i] - s_vals[i] < 1.0e-8))
    {
      m_q_val = decimal2Fraction(position_factors[i] - s_vals[i]);
      m_vals[i] = m_q_val.first; q_vals[i] = m_q_val.second;

      if (!(unique_m_q_pairs_count.insert(make_pair(m_q_val, 1)).second))
      {
        ++unique_m_q_pairs_count[m_q_val];
      }
      unique_m_q_pairs_vector[m_q_val].push_back(ss.str());
    }
    else
    {
      ++unique_m_q_pairs_count[key_0];
      unique_m_q_pairs_vector[key_0].push_back(ss.str());
      m_vals[i] = 0;
      q_vals[i] = 1;
    }

    if (verbosity >= 2)
    {
      cout << "P" << i + 1 << " " << positions[i] << " for (" << input.h << " " << input.k << " " << input.l << "): " << position_factors[i]
      << "\n  s value: " << s_vals[i]
      << "\n  m/q value: " << position_factors[i] - s_vals[i]
      << "\n  m value: " << m_vals[i]
      << "\n  q value: " << q_vals[i]
      << "\n  Distance from plane to origin: " << position_factors[i] * dhkl
      << "\n\n";
    }
  }

  map <pair <int, int>, vector <string> >::const_iterator map_vec_it = unique_m_q_pairs_vector.begin();
  for (map <pair <int, int>, int>::const_iterator it = unique_m_q_pairs_count.begin(); it != unique_m_q_pairs_count.end(); ++it)
  {
    if ((it->first).first == 0) {cout << "m/q = 0  ";}
    else {cout << "m/q = " << (it->first).first << "/" << (it->first).second;}
    cout << "; N = " << it->second << " (";
    for (vector <string>::const_iterator vec_it = (map_vec_it->second).begin(); vec_it != (map_vec_it->second).end();)
    {
      cout << (*vec_it);
      if (++vec_it != (map_vec_it->second).end()) {cout << ", ";}
    }
    cout << ")\n";
    ++map_vec_it;
  }
  cout << "\n";

  map_vec_it = unique_m_q_pairs_vector.begin();
  for (map <pair <int, int>, int>::const_iterator it = unique_m_q_pairs_count.begin(); it != unique_m_q_pairs_count.end(); ++it)
  {
    cout << "\u03C1^("; // lowercase rho
    if ((it->first).first == 0) {cout << "0";}
    else {cout << (it->first).first << "/" << (it->first).second;}
    cout << ")_(" << input.h << " " << input.k << " " << input.l << ") = "
         << "N^(";
    if ((it->first).first == 0) {cout << "0";}
    else {cout << (it->first).first << "/" << (it->first).second;}
    cout << ")_(" << input.h << " " << input.k << " " << input.l << ") d'_("
         << input.h << " " << input.k << " " << input.l << ")/V_cell = "
         << (map_vec_it->second).size() * dhkl / v_cell << "\n";
    ++map_vec_it;
  }
}

int main(int argc, char** argv)
{
  string input_file;
  inputData input;
  vector <Position> unit_cell_positions;
  double unit_cell_volume, interplanar_spacing;
  unsigned int verbosity = 0;
  int h,k,l;

  try
  {
    cxxopts::Options options(argv[0], "Script to determine the planar density based on the Position-Duplication-Number (PDN) method");
    options
      .positional_help("file h k l")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Input file", cxxopts::value<string>(input_file), "file")
        ("h", "h Miller index", cxxopts::value<int>(h), "h")
        ("k", "k Miller index", cxxopts::value<int>(k), "k")
        ("l", "l Miller index", cxxopts::value<int>(l), "l")
        ("v,verbose", "Verbosity - if specified, shows the intermediate calulations", cxxopts::value<bool>())
        ("c,citation", "Show the citation of the paper where the method originated", cxxopts::value<bool>()->implicit_value("true"))
        ("help", "Show the help");

    options.parse_positional({"file", "h", "k", "l"});
    auto result = options.parse(argc, argv);
    verbosity = result.count("verbose");

    if (result.count("help") || !result.count("file") || !result.count("h") || !result.count("k") || !result.count("l"))
    {
      cout << options.help() << "\n";
      printInputFileHelp();
      return EXIT_SUCCESS;
    }

    if (result.count("citation")) {printCitationInfo();}

    if (result.count("file"))
    {
      input = parseInputFile(input_file, verbosity);
      input.h = h; input.k = k; input.l = l;
      if (verbosity >= 1)
      {
        cout << "\n  Plane indices: " << input.h << " " << input.k << " " << input.l
        << "\n\n";
      }
      unit_cell_positions = parseUnitCellFile(input.unit_cell_file, verbosity);
      unit_cell_volume = calculateUnitCellVolume(input, verbosity);
      interplanar_spacing = calculateInterplanarSpacing(input, verbosity);
      calculatePlanarDensities(input, unit_cell_positions, interplanar_spacing, unit_cell_volume, verbosity);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cerr << "Error parsing options: " << e.what() << "\n";
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

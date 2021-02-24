#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <numeric>

#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

const static map <string, double> atomic_weights = {
  {"H", 1.0079}, {"He", 4.0026}, {"Li", 6.941}, {"Be", 9.0122}, {"B", 10.811}, {"C", 12.0107},
  {"N", 14.0067}, {"O", 15.9994}, {"F", 18.9984}, {"Ne", 20.1797}, {"Na", 22.9897}, {"Mg", 24.305},
  {"Al", 26.9815}, {"Si", 28.0855}, {"P", 30.9738}, {"S", 32.065}, {"Cl", 35.453}, {"K", 39.0983},
  {"Ar", 39.948}, {"Ca", 40.078}, {"Sc", 44.9559}, {"Ti", 47.867}, {"V", 50.9415}, {"Cr", 51.9961},
  {"Mn", 54.938}, {"Fe", 55.845}, {"Ni", 58.6934}, {"Co", 58.9332}, {"Cu", 63.546}, {"Zn", 65.39},
  {"Ga", 69.723}, {"Ge", 72.64}, {"As", 74.9216}, {"Se", 78.96}, {"Br", 79.904}, {"Kr", 83.8},
  {"Rb", 85.4678}, {"Sr", 87.62}, {"Y", 88.9059}, {"Zr", 91.224}, {"Nb", 92.9064}, {"Mo", 95.94},
  {"Tc", 98.0}, {"Ru", 101.07}, {"Rh", 102.9055}, {"Pd", 106.42}, {"Ag", 107.8682}, {"Cd", 112.411},
  {"In", 114.818}, {"Sn", 118.71}, {"Sb", 121.76}, {"I", 126.9045}, {"Te", 127.6}, {"Xe", 131.293},
  {"Cs", 132.9055}, {"Ba", 137.327}, {"La", 138.9055}, {"Ce", 140.116}, {"Pr", 140.9077}, {"Nd", 144.24},
  {"Pm", 145.0}, {"Sm", 150.36}, {"Eu", 151.964}, {"Gd", 157.25}, {"Tb", 158.9253}, {"Dy", 162.5},
  {"Ho", 164.9303}, {"Er", 167.259}, {"Tm", 168.9342}, {"Yb", 173.04}, {"Lu", 174.967}, {"Hf", 178.49},
  {"Ta", 180.9479}, {"W", 183.84}, {"Re", 186.207}, {"Os", 190.23}, {"Ir", 192.217}, {"Pt", 195.078},
  {"Au", 196.9665}, {"Hg", 200.59}, {"Tl", 204.3833}, {"Pb", 207.2}, {"Bi", 208.9804}, {"Po", 209.0},
  {"At", 210.0}, {"Rn", 222.0}, {"Fr", 223.0}, {"Ra", 226.0}, {"Ac", 227.0}, {"Pa", 231.0359},
  {"Th", 232.0381}, {"Np", 237.0}, {"U", 238.0289}, {"Am", 243.0}, {"Pu", 244.0}, {"Cm", 247.0},
  {"Bk", 247.0}, {"Cf", 251.0}, {"Es", 252.0}, {"Fm", 257.0}, {"Md", 258.0}, {"No", 259.0},
  {"Rf", 261.0}, {"Lr", 262.0}, {"Db", 262.0}, {"Bh", 264.0}, {"Sg", 266.0}, {"Mt", 268.0},
  {"Rg", 272.0}, {"Hs", 277.0}, {"Ds", 0.0}, {"Cn", 0.0}, {"Nh", 0.0}, {"Fl", 0.0}, {"Mc", 0.0},
  {"Lv", 0.0}, {"Ts", 0.0}, {"Og", 0.0}
};

void printAtomicWeights(const vector <string>& elements)
{
  for (size_t i = 0; i < elements.size(); ++i)
  {
    cout << elements[i] << ": " << atomic_weights.at(elements[i]) << " amu\n";
  }
}

vector <double> promptForNumbers(const string& percent_type, const vector <string>& elements)
{
  vector <double> nums (elements.size(), 0.0);
  for (size_t i = 0; i < elements.size() - 1; ++i)
  {
    cout << "Enter the " << percent_type << " percent of element " << elements[i] << " : ";
    cin  >> nums[i];
    while (nums[i] >= 1 || nums[i] <= 0)
    {
      cout << "The number must be between 0 and 1: ";
      cin  >> nums[i];
    }
  }

  nums[elements.size() - 1] = 1 - accumulate(nums.begin(), nums.end() - 1, 0.0);
  if (nums.back() < 0) {
    cerr << "Error: cannot have negative concentration (" << elements.back() << " = " << nums.back() << ")\n";
    exit(INPUT_FORMAT_ERROR);
  }
  cout << "Assumed the " << percent_type << " percent of element " << elements.back() << " was " << nums.back() << "\n";

  if (abs(1.0 - accumulate(nums.begin(), nums.end(), 0.0)) > 1.0e-8)
  {
    cerr << "Error: total " << percent_type << " percent must be equal to 1!\n";
    exit(INPUT_FORMAT_ERROR);
  }

  return nums;
}

void convertToAtomic(const vector <string>& elements)
{
  vector <double> terms (elements.size(), 0.0);
  vector <double> weight_percents = promptForNumbers("weight", elements);
  for (size_t i = 0; i < terms.size(); ++i)
  {
    terms[i] = weight_percents[i] / atomic_weights.at(elements[i]);
  }

  double denominator = accumulate(terms.begin(), terms.end(), 0.0);

  if (abs(denominator) < 1.0e-8) {cerr << "Error: Division by zero!\n"; exit(DIVIDE_BY_ZERO);}

  for (size_t i = 0; i < terms.size(); ++i)
  {
    cout << "Atomic percent of " << elements[i] << ": " << terms[i] / denominator << "\n";
  }
}

void convertToWeight(const vector <string>& elements)
{
  vector <double> terms (elements.size(), 0.0);
  vector <double> atomic_percents = promptForNumbers("atomic", elements);
  for (size_t i = 0; i < terms.size(); ++i)
  {
    terms[i] = atomic_percents[i] * atomic_weights.at(elements[i]);
  }

  double denominator = accumulate(terms.begin(), terms.end(), 0.0);

  if (abs(denominator) < 1.0e-8) {cerr << "Error: Division by zero!\n"; exit(DIVIDE_BY_ZERO);}

  for (size_t i = 0; i < terms.size(); ++i)
  {
    cout << "Weight percent of " << elements[i] << ": " << terms[i] / denominator << "\n";
  }
}

int main(int argc, char** argv)
{
  char conversion;
  vector <string> elements;
  try
  {
    cxxopts::Options options(argv[0], "Atomic percent to weight percent (and vice versa) calculator");
    options
      .positional_help("a|w element1 element2 ... elementN")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("c,conversion", "What to convert to: (a)tomic percent, or (w)eight percent", cxxopts::value<char>(conversion), "a|w")
        ("e,elements", "Space separated list of the elements", cxxopts::value<vector <string> >(elements), "element1 element2 ... elementN")
        ("atomic-weight", "Output the atomic weights of the specified elements", cxxopts::value<bool>())
        ("h,help", "Show the help");

    options.parse_positional({"conversion", "elements"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("conversion") == 0 || result.count("elements") == 0)
    {
      cout << options.help();

      return EXIT_SUCCESS;
    }

    if (result.count("atomic-weight")) {printAtomicWeights(elements);}

    if (result.count("conversion"))
    {
      if (!(conversion == 'a' || conversion == 'w'))
      {
        cout << "Conversion must be specified as 'a' or 'w' for atomic or weight percent respectively.\n";
        return INPUT_FORMAT_ERROR;
      }
      else
      {
        if (conversion == 'a') {convertToAtomic(elements);}
        else if (conversion == 'w') {convertToWeight(elements);}
      }
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cerr << "Error parsing options: " << e.what() << "\n";
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

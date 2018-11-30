/******************************************************************************
* This script takes as input the LAMMPS output file that contains the energies
* from the simulations, and the destination filename.  If the file already exists,
* the data is simply appended to it.
******************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

int N; // Number of atoms <-- Pretty sure this is bad practice

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

double extractAngleFromInfile(const string& file)
{
  double theta;

  stringstream th;
  size_t degree_pos = file.find("degree");

  if (degree_pos == string::npos)
  {
    cout << "Error: filename must contain the the angle in the format \"*_<angle>degree*\" in order to be processed correctly.\n";
    exit(FILE_NAME_ERROR);
  }
  size_t angle_pos = file.rfind("_", degree_pos) + 1;
  th << file.substr(angle_pos, degree_pos - angle_pos);
  th >> theta;

  return theta;
}

double extractMinEnergyFromFile(const string& file)
{
  string lost_atoms = "ERROR: Lost atoms: ";
  string energy_indicator = "  Energy initial, next-to-last, final = ";
  string str; // contains the extraneous info
  double en;
  vector <double> energies;
  bool atoms_lost = false;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  //ignore the first 5 lines of the file
  for (int i = 0; i < 5; ++i)
  {
    getline(fin, str);
  }

  fin >> N >> str;
  fin.ignore();

  while (getline(fin,str))
  {
    if (str.compare(energy_indicator) == 0)
    {
      for (int i = 0; i < 3; ++i)
      {
        fin >> en;
        energies.push_back(en);
      }
      fin.ignore();
    }

    if (str.substr(0,19).compare(lost_atoms) == 0)
    {
      cout << "Simulation error: lost atoms.\n";
      cout << "Energies found: ";
      for (unsigned int i = 0; i < energies.size(); ++i)
      {
        cout << "  Energy " << i + 1 << ": " << energies[i] << endl;
      }
      atoms_lost = true;
      break;
    }
  }

  fin.close();

  if (energies.size() == 0)
  {
    cout << "Did not find any energies to return in file " << file << "\n\n";
    exit(FILE_FORMAT_ERROR);
  }

  if (atoms_lost)
  {
    exit(ATOMS_LOST_ERROR);
  }

  return *min_element(energies.begin(), energies.end());
}

void writeDataToFile(const string& file, const double& theta, const double& min_energy)
{
  ofstream fout(file.c_str(), ofstream::app);
  checkFileStream(fout, file);

  fout << theta << " " << setprecision(15) << min_energy << " " << N << endl;
  fout.close();
}

int main(int argc, char** argv)
{
  string infile, outfile;
  double theta, min_energy;

  try
  {
    cxxopts::Options options(argv[0], "Extract energies from a LAMMPS minimization simulation.");
    options
    .positional_help("infile outfile")
    .show_positional_help();

    options
    .allow_unrecognised_options()
    .add_options()
      ("i,infile", "LAMMPS output file containing the energy data", cxxopts::value<string>(infile), "infile")
      ("o,outfile", "Output file to place the extracted energy. Appends to the file by default", cxxopts::value<string>(outfile), "outfile")
      ("a,angle", "The angle for which the data is being extracted", cxxopts::value<double>(theta), "angle")
      // ("output-all", "Outputs all of the energy values, sorted from least to greatest")
      ("h,help", "Show the help");

    options.parse_positional({"infile", "outfile"});
    auto result = options.parse(argc, argv);

    if (result.count("help") ||
       !(result.count("infile") && result.count("outfile")) )
    {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (result.count("infile") && result.count("outfile"))
    {
      if (!(result.count("angle")))
      {
        theta = extractAngleFromInfile(infile);
      }

      cout << "Angle: " << theta << endl;

      min_energy = extractMinEnergyFromFile(infile);
      writeDataToFile(outfile, theta, min_energy);
    }

  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

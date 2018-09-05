#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <cxxopts.hpp>
#include "atom.h"
#include "error_code_defines.h"

using namespace std;

bool print_atom_ids = false;

struct inputVars
{
  string outfile;
  double x_left, x_right;
  double y_left, y_right;
  double z_left, z_right;
  vector <int> ignored_atoms;

  void boundsSanityCheck()
  {
    if (x_left > x_right) {exit(INPUT_FORMAT_ERROR);}
    if (y_left > y_right) {exit(INPUT_FORMAT_ERROR);}
    if (z_left > z_right) {exit(INPUT_FORMAT_ERROR);}
    if (x_left < 0.0) {exit(INPUT_FORMAT_ERROR);}
    if (y_left < 0.0) {exit(INPUT_FORMAT_ERROR);}
    if (z_left < 0.0) {exit(INPUT_FORMAT_ERROR);}
    if (x_right > 1.0) {exit(INPUT_FORMAT_ERROR);}
    if (y_right > 1.0) {exit(INPUT_FORMAT_ERROR);}
    if (z_right > 1.0) {exit(INPUT_FORMAT_ERROR);}
  }
} input;

struct boxData
{
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double xy, xz, yz;
  double Lx, Ly, Lz;

  double calculateBoxLengths()
  {
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;
  }
} box;

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

void printAtomIds(const vector <Atom>& atoms)
{
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getMark() == 1)
    {
      cout << atoms[i].getId() << endl;
    }
  }
}

vector <string> getReferenceData(const string& file, vector <Atom>& reference)
{
  string str;
  int N;
  vector <string> vars;
  int id, type;
  double charge, x, y, z;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  fin >> str >> str; // ITEM: TIMESTEP
  fin >> str; //<timestep value>
  fin.ignore();

  getline(fin, str); // ITEM: NUMBER OF ATOMS
  fin >> N;
  fin.ignore();

  getline(fin, str); // ITEM: BOX BOUNDS
  fin >> box.xlow >> box.xhigh;
  fin >> box.ylow >> box.yhigh;
  fin >> box.zlow >> box.zhigh;
  box.calculateBoxLengths();
  fin.ignore();

  getline(fin, str); // TIME: ATOMS <data types...>
  stringstream tmp(str); // Get the data types
  tmp >> str >> str; // moves the stream to the interesting stuff
  while (tmp >> str) {vars.push_back(str);}

  reference.resize(N, Atom());
  while (getline(fin, str))
  {
    double num;
    stringstream ss(str);
    for (unsigned int i = 0; i < vars.size(); ++i)
    {
      ss >> num;
      if (vars[i].compare("id") == 0) {id = (int)(num);}
      else if (vars[i].compare("type") == 0) {type = (int)(num);}
      else if (vars[i].compare("q") == 0) {charge = num;}
      else if (vars[i].compare("x") == 0) {x = num;}
      else if (vars[i].compare("y") == 0) {y = num;}
      else if (vars[i].compare("z") == 0) {z = num;}
      else {continue;}
    }

    reference[id - 1] = Atom(id, type, charge, x, y, z);
    if (x > input.x_left * box.Lx && x < input.x_right * box.Lx &&
        y > input.y_left * box.Ly && y < input.y_right * box.Ly &&
        z > input.z_left * box.Lz && z < input.z_right * box.Lz)
    {
      if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(),
               reference[id - 1].getType()) == input.ignored_atoms.end())
        {
          reference[id - 1].setMark(1);
        }
    }
  }
  fin.close();

  return vars;
}

void getCurrentData(const string& file, vector <Atom>& current, const vector <string>& vars)
{
  string str;
  int id, type;
  double charge, x, y, z;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  for (int i = 0; i < 9; ++i)
  {
    getline(fin, str); // we are ignoring the first 9 lines (header of the comparison file)
  }

  while (getline(fin, str))
  {
    double num;
    stringstream ss(str);
    for (unsigned int i = 0; i < vars.size(); ++i)
    {
      ss >> num;
      if (vars[i].compare("id") == 0) {id = (int)(num);}
      else if (vars[i].compare("type") == 0) {type = (int)(num);}
      else if (vars[i].compare("q") == 0) {charge = num;}
      else if (vars[i].compare("x") == 0) {x = num;}
      else if (vars[i].compare("y") == 0) {y = num;}
      else if (vars[i].compare("z") == 0) {z = num;}
      else {continue;}
    }

    current[id - 1] = Atom(id, type, charge, x, y, z);
  }

  fin.close();
}

void compareTimesteps(vector <Atom>& reference, vector <Atom>& current)
{
  ofstream fout(input.outfile.c_str());
  checkFileStream(fout, input.outfile);

  for (unsigned int i = 0; i < reference.size(); ++i)
  {
    if (reference[i].getMark() == 1)
    {
      fout << current[i].getId() << " " << current[i].getType() << " "
           << current[i].getCharge() << " " << current[i].getX() << " "
           << current[i].getY() << " " << current[i].getZ() << endl;
    }
  }

  fout.close();
}

int main(int argc, char **argv)
{
  string reference_file;
  vector <string> current_file, vars;
  vector <Atom> reference_atoms, current_atoms;

  try
  {
    cxxopts::Options options(argv[0], "Track the positions of atoms from one frame to another.");
    options
      .positional_help("reference_file current_file")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("r,reference-file", "Reference file to compare against", cxxopts::value<string>(reference_file), "file")
        ("c,current-file", "Current file", cxxopts::value<vector<string> >(current_file), "file [file file file...]")
        ("o,output", "Output file name", cxxopts::value<string>(input.outfile)->default_value("*_tracked.dat"), "file")
        ("xlo", "low x boundary of the tracked atoms", cxxopts::value<double>(input.x_left)->default_value("0.45"), "value")
        ("xhi", "high x boundary of the tracked atoms", cxxopts::value<double>(input.x_right)->default_value("0.55"), "value")
        ("ylo", "low y boundary of the tracked atoms", cxxopts::value<double>(input.y_left)->default_value("0.0"), "value")
        ("yhi", "high y boundary of the tracked atoms", cxxopts::value<double>(input.y_right)->default_value("1.0"), "value")
        ("zlo", "low z boundary of the tracked atoms", cxxopts::value<double>(input.z_left)->default_value("0.0"), "value")
        ("zhi", "high z boundary of the tracked atoms", cxxopts::value<double>(input.z_right)->default_value("1.0"), "value")
        ("i,ignore", "Atom type to ignore", cxxopts::value<vector <int> >(input.ignored_atoms), "atom_type")
        ("print-atom-ids", "Prints a list of the atom ids found within the specified boundary. Only needs a reference structure.")
        ("h,help", "Show the help");

    options.parse_positional({"reference-file","current-file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !(result.count("reference-file")))
    {
      cout << options.help() << endl << endl
           << "Note that the reference file and current file _must_ be LAMMPS dump files right now.\n";
      return EXIT_SUCCESS;
    }

    if (result.count("print-atom-ids")) {print_atom_ids = true;}

    if (result.count("xlo") || result.count("xhi") ||
        result.count("ylo") || result.count("yhi") ||
        result.count("zlo") || result.count("zhi"))
    {
      input.boundsSanityCheck();
    }

    vars = getReferenceData(reference_file, reference_atoms); // sets up the variable list, and gets the reference atom data.
    if (print_atom_ids)
    {
      printAtomIds(reference_atoms);
      return EXIT_SUCCESS;
    }

    if (!(result.count("current-file")))
    {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    current_atoms.resize(reference_atoms.size(), Atom());
    for (unsigned int i = 0; i < current_file.size(); ++i)
    {
      if (!(result.count("output")))
      {
        input.outfile = current_file[i].substr(0,current_file[i].find(".dump")) + "_tracked.dat";
      }
      getCurrentData(current_file[i], current_atoms, vars);
      compareTimesteps(reference_atoms, current_atoms);
      cout << "\rFile " << current_file[i] << " processed.";
    }
    cout << endl;

  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

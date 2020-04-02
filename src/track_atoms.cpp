#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"

using namespace std;

struct inputVars
{
  string outfile;
  double x_left, x_right;
  double y_left, y_right;
  double z_left, z_right;
  vector <int> ignored_atoms;
  vector <int> tracked_atoms;
  Position r_center = Position(0.5,0.5,0.5);
  double r_lower_bound;
  double r_upper_bound;
  string radial_system = "none";
  bool using_radial = false;

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

  void radialSanityCheck()
  {
    if (r_lower_bound < 0.0) {exit(INPUT_FORMAT_ERROR);}
    if (r_upper_bound < r_lower_bound) {exit(INPUT_FORMAT_ERROR);}
    for (unsigned int i = 0; i < 3; ++i)
    {
      if (r_center[i] < 0.0) {exit(INPUT_FORMAT_ERROR);}
      if (r_center[i] > 1.0) {exit(INPUT_FORMAT_ERROR);}
    }
  }
} input;

struct boxData
{
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double xy, xz, yz;
  double Lx, Ly, Lz;

  void calculateBoxLengths()
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

void printAtomInfo(const vector <Atom>& atoms)
{
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getMark() == 1)
    {
      cout << atoms[i].getId() << " "
           << atoms[i].getType() << " "
           << atoms[i].getCharge() << " "
           << atoms[i].getWrapped()[0] << " "
           << atoms[i].getWrapped()[1] << " "
           << atoms[i].getWrapped()[2] << endl;
    }
  }
}

pair <int, vector <string> > processLAMMPSDump(istream& fin)
{
  string str;
  int N;
  vector <string> vars;

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

  pair <int, vector <string> > data = make_pair(N, vars);
  return data;
}

pair <int, vector <string> > processLAMMPSInput(istream& fin)
{
  string str;
  int N;
  vector <string> vars;
  pair <int, vector <string> > data;

  getline(fin, str);
  size_t left_bracket = str.find("[");
  size_t right_bracket = str.find("]", left_bracket);
  stringstream ss(str.substr(left_bracket + 1, right_bracket - left_bracket - 1));
  while (ss >> str)
  {
    transform(str.begin(), str.end(), str.begin(), ::tolower);
    if (str.compare("id") == 0) {vars.push_back("id");}
    else if (str.compare("charge") == 0) {vars.push_back("q");}
    else {vars.push_back(str);}
  }

  fin >> N >> str;
  fin >> str >> str >> str; // n_types atom types

  data = make_pair(N, vars);

  fin >> box.xlow >> box.xhigh >> str >> str;
  fin >> box.ylow >> box.yhigh >> str >> str;
  fin >> box.zlow >> box.zhigh >> str >> str;
  box.calculateBoxLengths();

  fin >> str;
  fin.ignore();
  getline(fin,str);

  return data;
}

vector <string> getReferenceData(const string& file, vector <Atom>& reference)
{
  string str;
  int N, num_tracked = 0;
  vector <string> vars;
  int id, type;
  double charge, x, y, z;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  getline(fin, str);
  fin.seekg(0, ios::beg);

  // Process the header information based on the file format
  if (str.find("ITEM: TIMESTEP") != string::npos)
  {
    pair <int, vector <string> > data;
    data = processLAMMPSDump(fin);
    N = data.first;
    vars = data.second;
  }
  else // Assumes that the file is an input file
  {
    pair <int, vector <string> > data;
    data = processLAMMPSInput(fin);
    N = data.first;
    vars = data.second;
  }

  // Sanity check - make sure we aren't looking for things that don't exist.
  // here we check to see if any of the IDs provided are higher than the number of atoms
  for (vector <int>::iterator it = input.tracked_atoms.begin();
       it != input.tracked_atoms.end();)
  {
    if (*it > N) {input.tracked_atoms.erase(it);}
    else {++it;}
  }

  if (input.tracked_atoms.size() != 0) {reference.resize(input.tracked_atoms.size(), Atom());}
  else {reference.resize(N, Atom());}

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

    Position p (x,y,z);
    if (!(input.tracked_atoms.empty()))
    {
      if (find(input.tracked_atoms.begin(), input.tracked_atoms.end(), id) != input.tracked_atoms.end())
      {
        reference[num_tracked] = Atom(id, type, charge, p);
        reference[num_tracked++].setMark(1);
      }
    }
    else
    {
      reference[id - 1] = Atom(id, type, charge, p);
      if (input.using_radial)
      {
        double rlo_sq = input.r_lower_bound * input.r_lower_bound;
        double rhi_sq = input.r_upper_bound * input.r_upper_bound;
        double x_new = x - (input.r_center[0] * box.Lx);
        double y_new = y - (input.r_center[1] * box.Ly);
        double x_sq = x_new * x_new;
        double y_sq = y_new * y_new;
        double r_sq = x_sq + y_sq;

        if (input.radial_system.compare("s") == 0)
        {
          double z_new = z - (input.r_center[2] * box.Lz);
          r_sq += z_new * z_new;
        }

        if (r_sq >= rlo_sq && r_sq <= rhi_sq)
        {
          if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(),
              reference[id - 1].getType()) == input.ignored_atoms.end())
          {
            reference[id - 1].setMark(1);
          }
        }

      }
      else
      {
        if (x >= input.x_left * box.Lx && x <= input.x_right * box.Lx &&
          y >= input.y_left * box.Ly && y <= input.y_right * box.Ly &&
          z >= input.z_left * box.Lz && z <= input.z_right * box.Lz)
        {
          if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(),
              reference[id - 1].getType()) == input.ignored_atoms.end())
          {
            reference[id - 1].setMark(1);
          }
        }
      }
    }
  }

  if (!(input.tracked_atoms.empty()) && num_tracked != reference.size())
  {
    cout << "Error tracking specific atoms.\n";
    exit(ATOM_COUNT_ERROR);
  }
  fin.close();

  return vars;
}

void getCurrentData(const string& file, vector <Atom>& current, const vector <string>& vars)
{
  // NOTE: This function currently assumes a LAMMPS dump file!
  string str;
  int id, type, num_tracked = 0;
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

    Position p (x,y,z);
    if (!(input.tracked_atoms.empty()))
    {
      if (find(input.tracked_atoms.begin(), input.tracked_atoms.end(), id) != input.tracked_atoms.end())
      {
        current[num_tracked++] = Atom(id, type, charge, p);
      }
    }
    else {current[id - 1] = Atom(id, type, charge, p);}
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
           << current[i].getCharge() << " " << current[i].getWrapped()[0] << " "
           << current[i].getWrapped()[1] << " " << current[i].getWrapped()[2] << endl;
    }
  }

  fout.close();
}

int main(int argc, char **argv)
{
  string reference_file;
  vector <string> current_file, vars;
  vector <Atom> reference_atoms, current_atoms;
  double center_x, center_y, center_z;

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
        ("xlo", "low x boundary of the tracked atoms in fractional coordinates", cxxopts::value<double>(input.x_left)->default_value("0.45"), "value")
        ("xhi", "high x boundary of the tracked atoms in fractional coordinates", cxxopts::value<double>(input.x_right)->default_value("0.55"), "value")
        ("ylo", "low y boundary of the tracked atoms in fractional coordinates", cxxopts::value<double>(input.y_left)->default_value("0.0"), "value")
        ("yhi", "high y boundary of the tracked atoms in fractional coordinates", cxxopts::value<double>(input.y_right)->default_value("1.0"), "value")
        ("zlo", "low z boundary of the tracked atoms in fractional coordinates", cxxopts::value<double>(input.z_left)->default_value("0.0"), "value")
        ("zhi", "high z boundary of the tracked atoms in fractional coordinates", cxxopts::value<double>(input.z_right)->default_value("1.0"), "value")
        ("rlo", "low radius value (in Angstroms)", cxxopts::value<double>(input.r_lower_bound)->default_value("0.0"), "value")
        ("rhi", "high radius value (in Angstroms)", cxxopts::value<double>(input.r_upper_bound)->default_value("10,0"), "value")
        ("r-center-x", "x value of the origin (in fractional coordinates) for defining the radius", cxxopts::value<double>(center_x)->default_value("0.5"), "value")
        ("r-center-y", "y value of the origin (in fractional coordinates) for defining the radius", cxxopts::value<double>(center_y)->default_value("0.5"), "value")
        ("r-center-z", "z value of the origin (in fractional coordinates) for defining the radius", cxxopts::value<double>(center_z)->default_value("0.5"), "value")
        ("r-type", "using a (c)ylindrical or (s)pherical radial system", cxxopts::value<string>(input.radial_system), "c|s")
        ("i,ignore", "Atom type to ignore", cxxopts::value<vector <int> >(input.ignored_atoms), "atom_type")
        ("a,atom", "Atom to specifically track.  Tracks atoms even if the atom type is specified using --ignore", cxxopts::value<vector <int> >(input.tracked_atoms), "atom_id")
        ("print-atom-ids", "Prints a list of the atom ids found within the specified boundary. Only needs a reference structure.")
        ("print-atom-info", "Prints the atom information.  Only needs a reference structure")
        ("h,help", "Show the help");

    options.parse_positional({"reference-file","current-file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !(result.count("reference-file")))
    {
      cout << options.help() << endl << endl;
      return EXIT_SUCCESS;
    }

    if (result.count("xlo") || result.count("xhi") ||
    result.count("ylo") || result.count("yhi") ||
    result.count("zlo") || result.count("zhi"))
    {
      input.boundsSanityCheck();
    }

    if (result.count("r-type"))
    {
      input.using_radial = true;
      if (!(input.radial_system.compare("c") == 0 ||
          input.radial_system.compare("s") == 0))
      {
        cout << "radius-type must be either \'c\' or \'s\'\n";
        return INPUT_FORMAT_ERROR;
      }
    }

    if ((result.count("rlo") || result.count("rhi")) && !result.count("r-type"))
    {
      cout << "Radius type (c|s) must be specified for radial tracking.\n";
      return INPUT_FORMAT_ERROR;
    }

    if (result.count("r-center-x") || result.count("r-center-y") || result.count("r-center-z"))
    {
      input.r_center = Position(center_x, center_y, center_z);
    }

    input.radialSanityCheck(); // Make sure radius inputs are valid

    vars = getReferenceData(reference_file, reference_atoms); // sets up the variable list, and gets the reference atom data.

    if (result.count("print-atom-ids"))
    {
      printAtomIds(reference_atoms);
      return EXIT_SUCCESS;
    }

    if (result.count("print-atom-info"))
    {
      printAtomInfo(reference_atoms);
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
      cout << flush;
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

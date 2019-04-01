#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

vector <int> getAtomIds(const string& file)
{
  vector <int> atom_ids;
  string str;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  while (getline(fin, str))
  {
    stringstream ss(str);
    int tmp;
    if (!(ss >> tmp))
    {
      cout << "Error reading atom_ids file.\n";
      fin.close();
      exit(FILE_FORMAT_ERROR);
    }
    atom_ids.push_back(tmp);
  }

  fin.close();

  return atom_ids;
}

pair <vector<int>, vector<string> > parseInputFile(const string& infile)
{
  string atom_ids_file, str;
  vector <int> atom_ids;
  vector <string> grain_files;
  pair <vector <int>, vector <string> > result;

  ifstream fin(infile.c_str());
  checkFileStream(fin, infile);

  getline(fin, atom_ids_file); // first line is assumed to be the atom_ids file

  while (getline(fin, str))
  {
    grain_files.push_back(str); // create the grain_files vector
  }
  fin.close();

  atom_ids = getAtomIds(atom_ids_file);

  return make_pair(atom_ids, grain_files);
}

void makeZonedData(const vector<int>& atom_ids, const vector<string>& grain_files)
{
  //TODO: Allow for arbitrary assigning of atoms to zones based on an inequality, or a value
  string str;
  vector <vector <string> > zones (3, vector <string> (0, ""));
  bool has_charge = false;

  for (unsigned int i = 0; i < grain_files.size(); ++i)
  {
    string outfile = grain_files[i].substr(0, grain_files[i].find(".dat")) + "_zoned.dat";
    ofstream fout(outfile.c_str());
    checkFileStream(fout, outfile);

    ifstream fin(grain_files[i].c_str());
    checkFileStream(fin, grain_files[i]);

    getline(fin, str);
    if (str.find("VARIABLES") == string::npos)
    {
      cout << "Error reading file " << grain_files[i] << endl;
      fout.close();
      fin.close();
      exit(FILE_FORMAT_ERROR);
    }
    fout << str << endl; // This is the variables line in the data file
    if (str.find("\"Atom Charge\"") != string::npos) {has_charge = true;}

    while (getline(fin, str))
    {
      int id, grain_num, int_tmp;
      float float_tmp, charge;

      stringstream ss(str);
      // line format: id type charge* x y z grain_num orientation_param xu yu zu
      // charge may not be included.
      ss >> id >> int_tmp;
      if (has_charge) {ss >> charge;}
      ss >> float_tmp >> float_tmp >> float_tmp >> grain_num
         >> float_tmp >> float_tmp >> float_tmp >> float_tmp;

      if (find(atom_ids.begin(), atom_ids.end(), id) == atom_ids.end())
      { // atom id not in list
        if (grain_num == 1) {zones[0].push_back(str);}
        else if (grain_num == 2) {zones[1].push_back(str);}
      }
      else {zones[2].push_back(str);}
    }

    for (unsigned int j = 0; j < zones.size(); ++j)
    {
      fout << "ZONE\n";

      for (unsigned int k = 0; k < zones[j].size(); ++k)
      {
        fout << zones[j][k] << endl;
      }
      zones[j].clear();
    }
    fout.close();
    fin.close();
  }
}

int main(int argc, char** argv)
{
  string infile;
  vector <int> atom_ids;
  vector <string> grain_files;
  try
  {
    cxxopts::Options options(argv[0], "Converts the results from find_grains.cpp into three zones: grain 1, grain 2, and a tracked slice of atoms.");
    options
      .positional_help("input_file")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("i,infile", "Input file", cxxopts::value<string>(infile), "infile")
        //("c,column", "Column to assign zones by (must be integers)", cxxopts::value<int>(), "col_num")
        ("h,help", "Show the help");

    options.parse_positional({"infile"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("infile") == 0)
    {
      cout << options.help() << endl << endl
           << "The input file must contain the following, in order, one entry per line:\n"
           << "    1) The file containing the list of atom ids in the tracked slice\n"
           << "    2-N) The output files from find_grains for the desired system.\n";
      return EXIT_SUCCESS;
    }

    if (result.count("infile"))
    {
      pair <vector<int>, vector<string> > input;
      input = parseInputFile(infile);
      atom_ids = input.first;
      grain_files = input.second;

      makeZonedData(atom_ids, grain_files);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

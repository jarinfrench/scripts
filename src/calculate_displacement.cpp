#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "atom.h"
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

bool checkUnwrapped(const vector <Atom> & a)
{
  for (unsigned int i = 0; i < a.size(); ++i)
  {
    if (a[i].getXu() > 1.0E-8 || a[i].getYu() > 1.0E-8 || a[i].getZu() > 1.0E-8)
    {
      return false;
    }
  }
  return true;
}

pair <string, vector <string> > parseInputFile(const string& input_file)
{
  string files, reference, str;
  vector <string> compared;

  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);

  getline(fin, reference); // get's the first line: assumes it to be the file

  while (getline(fin, str))
  {
    compared.push_back(str);
  }

  return make_pair(reference, compared);
}

bool processFile(const string& file, vector <Atom>& atoms)
{
  string str;
  unsigned int atom_id, n_read = 0;
  int atom_type, grain_num;
  bool warned = false; // user has been warned of no wrapped coordinates
  bool dump = false; // file type is a dump file
  double charge, x, y, z, f_i, xu, yu, zu;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  if (file.find(".dump") != string::npos)
  {
    dump = true;
    for (unsigned int i =0; i < 9; ++i) {getline(fin,str);} //ignore the first 9 lines

    if (str.find("xu") == string::npos || str.find("yu") == string::npos || str.find("zu") == string::npos)
    {
      if (!warned)
      {
        cout << "WARNING: Unwrapped coordinates do not exist.\n";
        warned = true;
      }
    }
    while (getline(fin, str))
    {
      stringstream ss(str);
      stringstream::pos_type pos = ss.tellg();

      if (!warned)
      {
        if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z >> xu >> yu >> zu))
        {
          ss.clear();
          ss.seekg(pos, ss.beg);
          charge = 0.0; // assumed we did not find charge
          if (!(ss >> atom_id >> atom_type >> x >> y >> z >> xu >> yu >> zu))
          {
            cout << "Error: data corrupted.  Expected id type x y z xu yu zu\n"
                 << "Line: " << str;
            exit(FILE_FORMAT_ERROR);
          }
        }
      }
      else
      {
        if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z))
        {
          ss.clear();
          ss.seekg(pos, ss.beg);
          charge = 0.0; // assumed we did not find charge
          if (!(ss >> atom_id >> atom_type >> x >> y >> z))
          {
            cout << "Error: data corrupted.  Expected id type x y z\n"
                 << "Line: " << str;
            exit(FILE_FORMAT_ERROR);
          }
          xu = x;
          yu = y;
          zu = z;
        }
      }

      if (atom_id > atoms.size()) {atoms.resize(atom_id, Atom());}

      atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, x, y, z);
      atoms[atom_id - 1].setXu(xu);
      atoms[atom_id - 1].setYu(yu);
      atoms[atom_id - 1].setZu(zu);
      ++n_read;
    }
  }
  else
  {
    getline(fin, str); // we ignore the first line
    while (getline(fin, str))
    {
      stringstream ss(str);
      stringstream::pos_type pos = ss.tellg();
      if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z >> grain_num >> f_i >> xu >> yu >> zu))
      {
        ss.clear(); // clear the error state of the stream
        ss.seekg(pos, ss.beg);
        charge = 0.0; // Assumed to not find any charge
        if (!(ss >> atom_id >> atom_type >> x >> y >> z >> grain_num >> f_i >> xu >> yu >> zu))
        {
          if (!warned)
          {
            cout << "WARNING: Not enough entries per line.  Assuming unwrapped coordinates do not exist.\n";
            warned = true;
          }

          xu = x;
          yu = y;
          zu = z;
        }
      }

      if (atom_id > atoms.size()) {atoms.resize(atom_id, Atom());}

      atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, x, y, z);
      atoms[atom_id - 1].setXu(xu);
      atoms[atom_id - 1].setYu(yu);
      atoms[atom_id - 1].setZu(zu);
      atoms[atom_id - 1].setMark(grain_num);
      ++n_read;
    }

  }

  fin.close();

  cout << "File " << file << " processed.\r";
  cout.flush();

  if (n_read != atoms.size())
  {
    cout << "Error reading file.  n_read = " << n_read << " != atoms.size() = " << atoms.size() << endl;
    exit(ATOM_COUNT_ERROR);
  }

  if (checkUnwrapped(atoms))
  {
    cout << "WARNING: All unwrapped coordinate values are zero.  Proceeding using wrapped coordinates.\n";
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      atoms[i].setXu(atoms[i].getX());
      atoms[i].setYu(atoms[i].getY());
      atoms[i].setZu(atoms[i].getZ());
    }
    return false; // Coordinates are NOT unwrapped
  }
  return true; // Coordinates ARE unwrapped
}

void writeData(const vector<Atom>& reference_atoms, const vector<Atom>& compared_atoms,
               const string& reference, const string& compared, string output_file,
               const bool& first_is_unwrapped, const bool& second_is_unwrapped)
{
  double disp_x, disp_y, disp_z, disp_mag;
  bool include_change_grain = true;

  if (reference.find(".dump") != string::npos || compared.find(".dump") != string::npos)
  {
    include_change_grain = false;
  }

  if (output_file.compare("none") == 0)
  {
    string first, second;
    if (reference.find(".dump") == string::npos)
    {
      first = reference.substr(0, reference.find("_"));
    }
    else
    {
      first = reference.substr(0, reference.find(".dump"));
    }
    if (compared.find(".dump") == string::npos)
    {
      second = compared.substr(0,compared.find("_"));
    }
    else
    {
      second = compared.substr(0, compared.find(".dump"));
    }

    output_file = first + "to" + second + "_diplacement_data.dat";
  }

  ofstream fout(output_file);
  checkFileStream(fout, output_file);


  fout << "VARIABLES = \"Atom ID\", \"Atom Type\", \"Atom Charge\", ";
  if (include_change_grain) { fout << "\"Changes Grain\", ";}
  fout << "\"Xu\", \"Yu\", \"Zu\", \"X(K)\", \"Y(K)\", \"Z(K)\", \"Magnitude\"\n";

  for (unsigned int i = 0; i < reference_atoms.size(); ++i)
  {
    int grain_change = 0;

    //if (second_is_wrapped) // something here that will calculated the unwrapped coordinates...?
    disp_x = compared_atoms[i].getXu() - reference_atoms[i].getXu();
    disp_y = compared_atoms[i].getYu() - reference_atoms[i].getYu();
    disp_z = compared_atoms[i].getZu() - reference_atoms[i].getZu();

    if (include_change_grain)
    {
      if (reference_atoms[i].getMark() != compared_atoms[i].getMark()) {grain_change = 1;}
    }

    disp_mag = sqrt(disp_x * disp_x + disp_y * disp_y + disp_z * disp_z);
    fout << reference_atoms[i].getId() << " " << reference_atoms[i].getType() << " "
         << reference_atoms[i].getCharge() << " " << reference_atoms[i].getXu() << " "
         << reference_atoms[i].getYu() << " " << reference_atoms[i].getZu() << " ";
    if (include_change_grain) {fout << grain_change << " ";}
    fout << disp_x << " " << disp_y << " " << disp_z << " " << disp_mag << endl;
  }
}

int main(int argc, char** argv)
{
  string reference, input_file, output_file;
  vector <string> compared;
  vector <Atom> reference_atoms, compared_atoms;
  try
  {
    cxxopts::Options options(argv[0], "Calculate the displacement vectors between two snapshots.");
    options
      .positional_help("reference compared [compared2 compared3 ...]")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("r,reference", "Reference state of atomic system", cxxopts::value<string>(reference), "reference")
        ("c,compared", "Atomic state(s) to compare to the reference", cxxopts::value<vector<string> >(compared), "compared")
        ("o,outfile", "Name of the output file", cxxopts::value<string>(output_file)->default_value("none"), "output_file")
        ("i,input", "Option to compare a list of files.  Specify the file containing this list.", cxxopts::value<string>(input_file), "input_file")
        ("first-only", "Compares all files to the first")
        ("h,help", "Show the help");

    options.parse_positional({"reference", "compared"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || ((result.count("reference") == 0 && result.count("compared") == 0) && result.count("input") == 0))
    {
      cout << options.help() << endl << endl
           << "Note that " << argv[0] << " --input input_file can be used in place of "
           << argv[0] << " reference compared\n";
      return EXIT_SUCCESS;
    }

    if (result.count("reference") && result.count("compared") && result.count("input"))
    {
      cout << "Error: either use an input file, or a command line list, not both.\n";
      return OPTION_PARSING_ERROR;
    }

    if (result.count("input"))
    {
      pair <string, vector <string> > files;
      files = parseInputFile(input_file);
      reference = files.first;
      compared = files.second;
    }

    bool first_is_unwrapped = processFile(reference, reference_atoms);

    for (unsigned int i = 0; i < compared.size(); ++i)
    {
      bool second_is_unwrapped = processFile(compared[i], compared_atoms);

      if (compared_atoms.size() != reference_atoms.size())
      {
        cout << "Error reading atom data.  reference_atoms.size() = " << reference_atoms.size() << " != compared_atoms.size() = " << compared_atoms.size() << endl;
        return VECTOR_SIZE_ERROR;
      }

      writeData(reference_atoms, compared_atoms, reference, compared[i], output_file, first_is_unwrapped, second_is_unwrapped);

      if (!result.count("first-only"))
      {
        reference_atoms = compared_atoms;
        reference = compared[i];
      }
    }

  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

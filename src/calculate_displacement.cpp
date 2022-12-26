#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm> // for find

#include "atom.h"
#include "position.h"
#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

bool checkUnwrapped(const vector <Atom> & a) {
  for (unsigned int i = 0; i < a.size(); ++i) {
    if (a[i].getUnwrapped()[0] > 1.0E-8 || a[i].getUnwrapped()[1] > 1.0E-8 || a[i].getUnwrapped()[2] > 1.0E-8) {
      return false;
    }
  }
  return true;
}

pair <string, vector <string> > parseInputFile(const string& input_file) {
  string files, reference, str;
  vector <string> compared;

  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);

  getline(fin, reference); // get's the first line: assumes it to be the file

  while (getline(fin, str)) {compared.push_back(str);}

  return make_pair(reference, compared);
}

vector <int> getIdsFromFile(const string& file) {
  vector <int> ids;
  int id;
  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  while (fin >> id) {ids.push_back(id);}

  fin.close();
  return ids;
}

bool processFile(const string& file, vector <Atom>& atoms, const vector <int>& ids) {
  string str;
  unsigned int atom_id, n_read = 0;
  int atom_type, grain_num;
  bool warned = false; // user has been warned of no wrapped coordinates
  double charge, x, y, z, f_i, xu, yu, zu;

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  if (file.find(".dump") != string::npos) { // dump file
    for (unsigned int i =0; i < 9; ++i) {getline(fin,str);} //ignore the first 9 lines

    if (str.find("xu") == string::npos || str.find("yu") == string::npos || str.find("zu") == string::npos) {
      if (!warned) {
        cout << "WARNING: Unwrapped coordinates do not exist.\n";
        warned = true;
      }
    }
    while (getline(fin, str)) {
      stringstream ss(str);
      stringstream::pos_type pos = ss.tellg();

      if (!warned) {
        if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z >> xu >> yu >> zu)) {
          ss.clear();
          ss.seekg(pos, ss.beg);
          charge = 0.0; // assumed we did not find charge
          if (!(ss >> atom_id >> atom_type >> x >> y >> z >> xu >> yu >> zu)) {
            cerr << "Error: data corrupted.  Expected id type x y z xu yu zu\n"
                 << "Line: " << str;
            exit(FILE_FORMAT_ERROR);
          }
        }
      } else {
        if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z)) {
          ss.clear();
          ss.seekg(pos, ss.beg);
          charge = 0.0; // assumed we did not find charge
          if (!(ss >> atom_id >> atom_type >> x >> y >> z)) {
            cerr << "Error: data corrupted.  Expected id type x y z\n"
                 << "Line: " << str;
            exit(FILE_FORMAT_ERROR);
          }
          xu = x;
          yu = y;
          zu = z;
        }
      }

      if (atom_id > atoms.size()) {atoms.resize(atom_id, Atom());}
      ++n_read;

      if (ids.size() != 0 && find(ids.begin(), ids.end(), atom_id) == ids.end()) continue;

      Position p (x,y,z);
      atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, p);
      p = Position(xu,yu,zu);
      atoms[atom_id - 1].setUnwrapped(p);
    }
  } else if (file.find(".dat") != string::npos) {
    string output_type = "";
    string header = "";
    getline(fin, header);
    size_t left_bracket = header.find("[");
    size_t right_bracket = header.find("]", left_bracket);
    while (output_type.find("Atoms") == string::npos) {
      getline(fin, output_type);
    }

    getline(fin, str); // blank line

    if (left_bracket != string::npos && right_bracket != string::npos) {
      string data_string = header.substr(left_bracket + 1, right_bracket - left_bracket - 1);
      stringstream ss(data_string);
      vector <string> vars;
      int x_ind = -1, y_ind = -1, z_ind = -1, id_ind = -1, type_ind = -1;
      int xu_ind = -1, yu_ind = -1, zu_ind = -1, q_ind = -1;
      string str2;
      while (ss >> str2) {
        transform(str2.begin(), str2.end(), str2.begin(), ::tolower);
        if (str2.compare("charge") == 0) {vars.push_back("q"); q_ind = vars.size() - 1;}
        else {
          vars.push_back(str2);
          if (str2.compare("xu") == 0) {xu_ind = vars.size() - 1; continue;}
          if (str2.compare("yu") == 0) {yu_ind = vars.size() - 1; continue;}
          if (str2.compare("zu") == 0) {zu_ind = vars.size() - 1; continue;}
          if (str2.compare("x") == 0) {x_ind = vars.size() - 1; continue;}
          if (str2.compare("y") == 0) {y_ind = vars.size() - 1; continue;}
          if (str2.compare("z") == 0) {z_ind = vars.size() - 1; continue;}
          if (str2.compare("id") == 0) {id_ind = vars.size() - 1; continue;}
          if (str2.compare("type") == 0) {type_ind = vars.size() - 1;}
        }
      }
      if (xu_ind == -1 && x_ind == -1) {
        cout << "Could not find x coordinate in data (data string: " << data_string << "). Exiting...\n";
        exit(FILE_FORMAT_ERROR);
      }
      if (yu_ind == -1 && y_ind == -1) {
        cout << "Could not find y coordinate in data (data string: " << data_string << "). Exiting...\n";
        exit(FILE_FORMAT_ERROR);
      }
      if (zu_ind == -1 && z_ind == -1) {
        cout << "Could not find z coordinate in data (data string: " << data_string << "). Exiting...\n";
        exit(FILE_FORMAT_ERROR);
      }
      if (id_ind == -1) {
        cout << "Could not find atom id in data (data string: " << data_string << "). Exiting...\n";
        exit(FILE_FORMAT_ERROR);
      }
      if (type_ind == -1) {
        cout << "Could not find type in data (data string: " << data_string << "). Exiting...\n";
        exit(FILE_FORMAT_ERROR);
      }

      while (getline(fin, str)) {
        vector <double> data;
        double tmp;
        stringstream ss(str);
        while (ss >> tmp) {data.push_back(tmp);}
        atom_id = (int)(data[id_ind]);
        atom_type = (int)(data[type_ind]);
        if (q_ind == -1) {charge = 0.0;}
        else {charge = data[q_ind];}
        if (atom_id > atoms.size()) {atoms.resize(atom_id, Atom());}
        ++n_read;

        if (ids.size() != 0 && find(ids.begin(), ids.end(), atom_id) == ids.end()) continue;
        Position p(data[x_ind],data[y_ind],data[z_ind]);
        atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, p);
        if (xu_ind == -1 || yu_ind == -1 || zu_ind == -1) {
          p = Position(0,0,0);
        } else { //(xu_ind != -1 && yu_ind != -1 && zu_ind != -1)
          p = Position(data[xu_ind], data[yu_ind], data[zu_ind]);
        }
        atoms[atom_id - 1].setUnwrapped(p);
      }
      if (xu_ind == -1 && yu_ind == -1 && zu_ind == -1) {
        if (!warned) {
          cout << "WARNING: Unwrapped coordinates do not exist.\n";
          warned = true;
        }
      }
    } else if (str.find("#") != string::npos) {
      cout << "Warning: image flags not currently taken into account.\n";
      stringstream ss(str);
      ss >> output_type >> output_type >> output_type; // we only need the last bit (e.g. Atoms # _charge_)
      if (output_type.compare("atomic") == 0) {
        // id type x y z
        while (getline(fin, str)) {
          stringstream s2(str);
          if (!(s2 >> atom_id >> atom_type >> x >> y >> z)) {
            cout << "Error: data corrupted. Expected id type x y z\n"
                 << "Line: " << str << "\n";
            exit(FILE_FORMAT_ERROR);
          }
          xu = 0;
          yu = 0;
          zu = 0;
          if (atom_id > atoms.size()) {atoms.resize(atom_id, Atom());}
          ++n_read;
          if (ids.size() != 0 && find(ids.begin(), ids.end(), atom_id) == ids.end()) continue;
          Position p(x,y,z);
          atoms[atom_id - 1] = Atom(atom_id, atom_type, 0, p);
          p = Position(xu, yu, zu);
          atoms[atom_id - 1].setUnwrapped(p);
        }
      } else if (output_type.compare("charge") == 0) {
        // id type q x y z
        while (getline(fin, str)) {
          stringstream s2(str);
          if (!(s2 >> atom_id >> atom_type >> charge >> x >> y >> z)) {
            cout << "Error: data corrupted. Expected id type q x y z\n"
                 << "Line: " << str << "\n";
            exit(FILE_FORMAT_ERROR);
          }
          xu = x;
          yu = y;
          zu = z;
          if (atom_id > atoms.size()) {atoms.resize(atom_id, Atom());}
          ++n_read;
          if (ids.size() != 0 && find(ids.begin(), ids.end(), atom_id) == ids.end()) continue;
          Position p(x,y,z);
          atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, p);
          p = Position(xu, yu, zu);
          atoms[atom_id - 1].setUnwrapped(p);
        }
      }
    }
  } else { // file format from find_grains output
    getline(fin, str); // we ignore the first line
    while (getline(fin, str)) {
      stringstream ss(str);
      stringstream::pos_type pos = ss.tellg();
      if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z >> grain_num >> f_i >> xu >> yu >> zu)) {
        ss.clear(); // clear the error state of the stream
        ss.seekg(pos, ss.beg);
        charge = 0.0; // Assumed to not find any charge
        if (!(ss >> atom_id >> atom_type >> x >> y >> z >> grain_num >> f_i >> xu >> yu >> zu)) {
          if (!warned) {
            cout << "WARNING: Not enough entries per line.  Assuming unwrapped coordinates do not exist.\n";
            warned = true;
          }

          xu = x;
          yu = y;
          zu = z;
        }
      }

      if (atom_id > atoms.size()) {atoms.resize(atom_id, Atom());}
      ++n_read;

      if (ids.size() != 0 && find(ids.begin(), ids.end(), atom_id) == ids.end()) continue;

      Position p (x,y,z);
      atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, p);
      p = Position(xu,yu,zu);
      atoms[atom_id - 1].setUnwrapped(p);
      atoms[atom_id - 1].setMark(grain_num);
    }
  }

  fin.close();

  cout << "File " << file << " processed.\r";
  cout.flush();

  if (n_read != atoms.size()) {
    cerr << "Error reading file.  n_read = " << n_read << " != atoms.size() = " << atoms.size() << endl;
    exit(ATOM_COUNT_ERROR);
  }

  if (checkUnwrapped(atoms)) {
    cout << "WARNING: All unwrapped coordinate values are zero.  Proceeding using wrapped coordinates.\n";
    for (unsigned int i = 0; i < atoms.size(); ++i) {
      atoms[i].setUnwrapped(atoms[i].getWrapped());
    }
    return false; // Coordinates are NOT unwrapped
  }
  return true; // Coordinates ARE unwrapped
}

void writeData(const vector<Atom>& reference_atoms, const vector<Atom>& compared_atoms,
               const string& reference, const string& compared, const vector<int>& ids,
               string output_file, const bool& first_is_unwrapped, const bool& second_is_unwrapped) {
  double disp_x, disp_y, disp_z, disp_mag;
  bool include_change_grain = true;

  if (reference.find(".dump") != string::npos || compared.find(".dump") != string::npos) {
    include_change_grain = false;
  }

  if (output_file.compare("none") == 0) {
    string first, second;
    if (reference.find(".dump") == string::npos) {
      first = reference.substr(0, reference.find("_"));
    } else {
      first = reference.substr(0, reference.find(".dump"));
    }
    if (compared.find(".dump") == string::npos) {
      second = compared.substr(0,compared.find("_"));
    } else {
      second = compared.substr(0, compared.find(".dump"));
    }

    output_file = first + "to" + second + "_displacement_data.dat";
  }

  ofstream fout(output_file);
  checkFileStream(fout, output_file);


  fout << "VARIABLES = \"Atom ID\", \"Atom Type\", \"Atom Charge\", ";
  if (include_change_grain) { fout << "\"Changes Grain\", ";}
  fout << "\"Xu\", \"Yu\", \"Zu\", \"X(K)\", \"Y(K)\", \"Z(K)\", \"Magnitude\"\n";

  for (unsigned int i = 0; i < reference_atoms.size(); ++i) {
    if (ids.size() != 0 && find(ids.begin(), ids.end(), reference_atoms[i].getId()) == ids.end()) continue;
    int grain_change = 0;

    //if (second_is_wrapped) // something here that will calculated the unwrapped coordinates...?
    disp_x = compared_atoms[i].getUnwrapped()[0] - reference_atoms[i].getUnwrapped()[0];
    disp_y = compared_atoms[i].getUnwrapped()[1] - reference_atoms[i].getUnwrapped()[1];
    disp_z = compared_atoms[i].getUnwrapped()[2] - reference_atoms[i].getUnwrapped()[2];

    if (include_change_grain) {
      if (reference_atoms[i].getMark() != compared_atoms[i].getMark()) {grain_change = 1;}
    }

    disp_mag = sqrt(disp_x * disp_x + disp_y * disp_y + disp_z * disp_z);
    fout << reference_atoms[i].getId() << " " << reference_atoms[i].getType() << " "
         << reference_atoms[i].getCharge() << " " << reference_atoms[i].getUnwrapped()[0] << " "
         << reference_atoms[i].getUnwrapped()[1] << " " << reference_atoms[i].getUnwrapped()[2] << " ";
    if (include_change_grain) {fout << grain_change << " ";}
    fout << disp_x << " " << disp_y << " " << disp_z << " " << disp_mag << endl;
  }
}

int main(int argc, char** argv) {
  string reference, input_file, output_file, ids_file;
  vector <string> compared;
  vector <Atom> reference_atoms, compared_atoms;
  vector <int> tracked_ids;
  try {
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
        ("f,ids-file", "File name containing the specific atoms to examine", cxxopts::value<string>(ids_file)->default_value("none"), "file")
        ("first-only", "Compares all files to the first")
        ("h,help", "Show the help");

    options.parse_positional({"reference", "compared"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || ((result.count("reference") == 0 && result.count("compared") == 0) && result.count("input") == 0)) {
      cout << options.help() << endl << endl
           << "Note that " << argv[0] << " --input input_file can be used in place of "
           << argv[0] << " reference compared\n";
      return EXIT_SUCCESS;
    }

    if (result.count("reference") && result.count("compared") && result.count("input")) {
      cerr << "Error: either use an input file, or a command line list, not both.\n";
      return OPTION_PARSING_ERROR;
    }

    if (result.count("ids-file")) {
      tracked_ids = getIdsFromFile(ids_file);
    }

    if (result.count("input")) {
      pair <string, vector <string> > files;
      files = parseInputFile(input_file);
      reference = files.first;
      compared = files.second;
    }

    bool first_is_unwrapped = processFile(reference, reference_atoms, tracked_ids);

    for (unsigned int i = 0; i < compared.size(); ++i) {
      bool second_is_unwrapped = processFile(compared[i], compared_atoms, tracked_ids);

      if (compared_atoms.size() != reference_atoms.size()) {
        cerr << "Error reading atom data.  reference_atoms.size() = " << reference_atoms.size() << " != compared_atoms.size() = " << compared_atoms.size() << endl;
        return VECTOR_SIZE_ERROR;
      }

      writeData(reference_atoms, compared_atoms, reference, compared[i], tracked_ids, output_file, first_is_unwrapped, second_is_unwrapped);

      if (!result.count("first-only")) {
        reference_atoms = compared_atoms;
        reference = compared[i];
      }
    }

  } catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

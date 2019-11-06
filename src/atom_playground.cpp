#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include "atom.h"
#include "error_code_defines.h"

#include <cxxopts.hpp>

using namespace std;

struct Header
{
  double xlo = 0.0, xhi = 0.0, ylo = 0.0, yhi = 0.0, zlo = 0.0, zhi = 0.0; // box lengths
  double Lx = 0.0, Ly = 0.0, Lz = 0.0; // box lengths
  double xy = 0.0, xz = 0.0, yz = 0.0; // triclinic factors

  int id_index = -1, type_index = -1, charge_index = -1, x_index = -1;
  int y_index = -1, z_index = -1;
  int xu_index = -1, yu_index = -1, zu_index = -1;

  vector <string> data_types;
  bool has_charge = false;

  int N = 0, n_types = 0; // number of atoms/atom types

  // map from the data_list index to the extra_info vector list index
  // the string is the data type, the first unsigned int is the index number
  // for original list (as it appears in the dump or dat file), and the second
  // unsigned int is the index for the extra_info vector (the member from the
  // Atom class)
  map <string, pair<unsigned int, unsigned int> > indexConverter;
};

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
  }
}

void showCommandList()
{
  cout << "\nr,read        Specify a new file to analyze\n"
       << "H,header      Show the header information\n"
       << "s,select      Select atom by various criteria\n"
       << "p,print       Print selected atom info\n"
       << "q,quit,exit   Exit the software\n"
       << "h,help        Show this list of accepted commands\n\n";
}

void welcome()
{
  cout << "Welcome to atom playground!  This software allows for in-depth\n"
       << "analysis of a single snapshot of an atomic structure.\n"
       << "Here is the list of accepted commands: \n";

  showCommandList();
}

void printHeaderInfo(const Header& header)
{
  cout << "Box lengths:\n"
       << "  x: " << header.xlo << " - " << header.xhi << " (" << header.Lx << ")\n"
       << "  y: " << header.ylo << " - " << header.yhi << " (" << header.Ly << ")\n"
       << "  z: " << header.zlo << " - " << header.zhi << " (" << header.Lz << ")\n"
       << "Triclinic factors: " << header.xy << " " << header.xz << " " << header.yz << " xy xz yz\n"
       << "N: " << header.N << "\n"
       << "Number of atom types: " << header.n_types << "\n"
       << "Data types: \n  ";
  for (unsigned int i = 0; i < header.data_types.size(); ++i)
  {
    cout << header.data_types[i] << " ";
  }
  cout << endl;
}

Header getHeaderData(const string& data_type, ifstream& fin)
{
  Header header;
  string str; // junk variable
  if (data_type.compare("dump") == 0)
  {
    getline(fin, str); // ITEM: TIMESTEP
    getline(fin, str); // actual timestep
    getline(fin, str); // ITEM: NUMBER OF ATOMS
    fin >> header.N; // number of atoms
    fin.ignore();
    getline(fin, str); // ITEM: BOX BOUNDS
    if (str.find("xy") == string::npos)
    {
      fin >> header.xlo >> header.xhi;
      fin >> header.ylo >> header.yhi;
      fin >> header.zlo >> header.zhi;
    }
    else
    {
      fin >> header.xlo >> header.xhi >> header.xy;
      fin >> header.ylo >> header.yhi >> header.xz;
      fin >> header.zlo >> header.zhi >> header.yz;
    }
    fin.ignore();
    getline(fin, str); // ITEM: ATOMS <data types>
    stringstream ss(str);
    string dummy;
    ss >> dummy >> dummy; // ITEMS: and ATOMS
    while (ss >> dummy)
    {
      if (dummy.compare("charge") == 0 || dummy.compare("q") == 0)
      {
        header.has_charge = true;
      }
      header.data_types.push_back(dummy);
    }
  }
  else if (data_type.compare("dat") == 0)
  {
    getline(fin, str); // first (comment) line
    if (str.find("[") == string::npos || str.find("]") == string::npos)
    {
      cout << "Unable to identify data types.  The first line in the file should\n"
           << "have one set of [] with the data types listed inside, separated by spaces.\n";
      exit(FILE_FORMAT_ERROR);
    }
    stringstream ss(str.substr(str.find("[") + 1, str.find("]", str.find("[")) - str.find("[") - 1));
    string dummy;
    while (ss >> dummy)
    {
      if (dummy.compare("charge") == 0 || dummy.compare("q") == 0)
      {
        header.has_charge = true;
      }
      header.data_types.push_back(dummy);
    }
    fin >> header.N >> str; // number of atoms
    fin >> header.n_types >> str >> str; // number of atom types
    fin >> header.xlo >> header.xhi >> str >> str;
    fin >> header.ylo >> header.yhi >> str >> str;
    fin >> header.zlo >> header.zhi >> str >> str;
    fin.ignore();
    getline(fin, str); // blank line or triclinic tilt factors
    if (str.find("xy") != string::npos) // Triclinic tilt factors
    {
      string tmp_str;
      stringstream s2(str);
      ss >> header.xy >> header.xz >> header.yz >> tmp_str >> tmp_str;
      getline(fin, str); // blank line
    }
    getline(fin, str); // Atoms
    getline(fin, str); // blank line
  }

  header.Lx = header.xhi - header.xlo;
  header.Ly = header.yhi - header.ylo;
  header.Lz = header.zhi - header.zlo;

  return header;
}

vector <Atom> getAtomData(ifstream& fin, Header& header)
{
  int n_atoms_read = 0;
  string str;
  vector <double> data;
  vector <Atom> atoms (header.N, Atom());
  int atom_id, type;
  double charge, x, y, z, xu, yu, zu;


  for (unsigned int i = 0; i < header.data_types.size(); ++i)
  {
    string name = header.data_types[i];
    if (name.compare("ID") == 0 || name.compare("id") == 0) {header.id_index = i;}
    else if (name.compare("type") == 0) {header.type_index = i;}
    else if (name.compare("charge") == 0 || name.compare("q") == 0) {header.charge_index = i;}
    else if (name.compare("x") == 0) {header.x_index = i;}
    else if (name.compare("y") == 0) {header.y_index = i;}
    else if (name.compare("z") == 0) {header.z_index = i;}
    else if (name.compare("xu") == 0) {header.xu_index = i;}
    else if (name.compare("yu") == 0) {header.yu_index = i;}
    else if (name.compare("zu") == 0) {header.zu_index = i;}
  }

  while (getline(fin, str))
  {
    stringstream ss(str);
    double dummy;
    vector <double> additional_data;
    while (ss >> dummy) {data.push_back(dummy);}
    for (unsigned int i = 0; i < data.size(); ++i)
    {
      if (i == header.id_index) {atom_id = (int)(data[i]);}
      else if (i == header.type_index) {type = (int)(data[i]);}
      else if (header.has_charge && i == header.charge_index) {charge = data[i];}
      else if (i == header.x_index) {x = data[i];}
      else if (i == header.y_index) {y = data[i];}
      else if (i == header.z_index) {z = data[i];}
      else if (i == header.xu_index) {xu = data[i];}
      else if (i == header.yu_index) {yu = data[i];}
      else if (i == header.zu_index) {zu = data[i];}
      else
      {
        additional_data.push_back(data[i]);
        pair <unsigned int, unsigned int> value = make_pair(i, additional_data.size() - 1);
        pair <string, pair <unsigned int, unsigned int> > key_value_pair = make_pair(header.data_types[i], value);
        header.indexConverter.insert(key_value_pair);
      }

      if (type > header.n_types) {header.n_types = type;}
    }

    data.clear();

    Position p(x,y,z);
    atoms[atom_id - 1] = Atom(atom_id, type, charge, p);

    if (header.xu_index != -1) {p.setX(xu);}
    if (header.yu_index != -1) {p.setY(yu);}
    if (header.zu_index != -1) {p.setZ(zu);}
    atoms[atom_id - 1].setUnwrapped(p);
    atoms[atom_id - 1].setExtraInfo(additional_data);
  }

  return atoms;
}

string getFileName()
{
  string filename;
  cout << "Please enter the filename: ";
  cin  >> filename;
  return filename;
}

pair <vector <Atom>, Header> readFile(const string& filename)
{
  vector <Atom> atom_tmp;
  Header header_tmp;

  ifstream fin(filename);
  checkFileStream(fin, filename);

  if (filename.substr(filename.length() - 4) == "dump")
  {
    header_tmp = getHeaderData("dump", fin);
    atom_tmp = getAtomData(fin, header_tmp);
  }
  else if (filename.substr(filename.length() - 3) == "dat")
  {
    header_tmp = getHeaderData("dat", fin);
    atom_tmp = getAtomData(fin, header_tmp);
  }
  else
  {
    bool known = false;
    string file_type;
    while (!known)
    {
      cout << "Unrecognized file type.  Please specify dump or dat: ";
      cin  >> file_type;
      if (file_type.compare("dump") == 0)
      {
        header_tmp = getHeaderData("dump", fin);
        atom_tmp = getAtomData(fin, header_tmp);
        known = true;
      }
      else if (file_type.compare("dat") == 0)
      {
        header_tmp = getHeaderData("dat", fin);
        atom_tmp = getAtomData(fin, header_tmp);
        known = true;
      }
    }
  }

  fin.close();
  return make_pair(atom_tmp,header_tmp);
}

vector <Atom> getAtomsByType(const vector <Atom>& atoms)
{
  string atom_types_string;
  vector <int> atom_types;
  vector <Atom> selected_atoms;

  cout << "Please enter the atom type(s) as a space-separated list: ";
  cin.ignore();
  getline(cin, atom_types_string);
  stringstream ss(atom_types_string);

  int type_tmp;
  while (ss >> type_tmp) {atom_types.push_back(type_tmp);}
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (find(atom_types.begin(), atom_types.end(), atoms[i].getType()) != atom_types.end())
    {
      selected_atoms.push_back(atoms[i]);
    }
  }

  return selected_atoms;
}

vector <int> getAtomIdsFromFile()
{
  string file, str;
  vector <int> atom_ids;

  cout << "Please enter the file name containing the list of atom IDs: ";
  cin  >> file;

  ifstream fin(file);
  checkFileStream(fin, file);

  while (getline(fin, str))
  {
    int id;
    stringstream ss(str);
    ss >> id;
    atom_ids.push_back(id);
  }

  fin.close();

  return atom_ids;
}

vector <int> getAtomIdsByRange(unsigned int num_atoms)
{
  string range_string, str;
  vector <string> ranges;
  vector <int> atom_ids;

  cout << "Please enter the range(s) for atom ids (e.g. BEGIN-10 20-25 27 100-END)\n"
       << "Note that BEGIN = 1, and END = the number of atoms: ";
  cin  >> range_string;

  stringstream ss(range_string);
  while (ss >> str) {ranges.push_back(str);}

  for (unsigned int i = 0; i < ranges.size(); ++i)
  {
    size_t minus_pos = ranges[i].find("-");
    if (minus_pos == string::npos)
    {
      int single_id;
      stringstream s2(ranges[i]);
      s2 >> single_id;
      atom_ids.push_back(single_id);
    }
    else // we found a minus sign, indicating a range
    {
      unsigned int first_index = 0, last_index = 0;
      string temp = ranges[i].substr(0,minus_pos);

      if (temp.compare("BEGIN") == 0) {first_index = 1;}
      else
      {
        stringstream s3(temp);
        s3 >> first_index;
      }

      temp = ranges[i].substr(minus_pos + 1);
      if (temp.compare("END") == 0) {last_index = num_atoms;}
      else
      {
        stringstream s3(temp);
        s3 >> last_index;
      }

      for (unsigned int j = first_index; j <= last_index; ++j)
      {
        atom_ids.push_back(j);
      }
    }
  }

  return atom_ids;
}

vector <Atom> getAtomsById(const vector <Atom>& atoms, const vector <int>& atom_ids)
{
  vector <Atom> selected_atoms (atom_ids.size(), Atom());

  for (unsigned int i = 0; i < atom_ids.size(); ++i)
  {
    selected_atoms[i] = atoms[atom_ids[i] - 1];
  }

  return selected_atoms;
}

vector <Atom> getSelectedAtoms(const vector <Atom>& atoms)
{
  char atom_specification;
  vector <Atom> selected_atoms;

  cout << "Would you like to specify atoms by type(t), file(f), or range(r)?: ";
  cin  >> atom_specification;

  while (atom_specification != 't' && atom_specification != 'f' && atom_specification != 'r')
  {
    // If we have a cin error, clear the error state, wipe the input, and try again
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(),'\n');
    cout << "Please enter either \'t\', \'f\', or \'r\': ";
    cin  >> atom_specification;
  }

  if (atom_specification == 't')
  {
    selected_atoms = getAtomsByType(atoms);
  }
  else
  {
    vector <int> atom_ids;
    if (atom_specification == 'f')
    {
      atom_ids = getAtomIdsFromFile();
    }
    else
    {
      atom_ids = getAtomIdsByRange(atoms.size());
    }

    selected_atoms = getAtomsById(atoms, atom_ids);
  }

  cout << "Selected " << selected_atoms.size() << " atoms (" << (double)(selected_atoms.size()) / (double)(atoms.size()) * 100.0 << "%).\n";

  return selected_atoms;
}

void printSelectedAtomInfo(const vector <Atom>& selected_atoms, Header& header)
{
  string desired_data_string, filename;
  vector <bool> print_data(header.data_types.size(), false);

  cout << "The following data is available: \n  ";
  for (unsigned int i = 0; i < header.data_types.size(); ++i)
  {
    cout << header.data_types[i] << " ";
  }
  cout << endl;

  cout << "Please enter the data types to print as they appear above: ";
  cin.ignore();
  getline(cin, desired_data_string);

  stringstream ss(desired_data_string);
  string dummy;

  while (ss >> dummy)
  {
    if (dummy.compare("id") == 0 || dummy.compare("ID") == 0) {print_data[header.id_index] = true;}
    else if (dummy.compare("type") == 0) {print_data[header.type_index] = true;}
    else if (dummy.compare("charge") == 0 || dummy.compare("q") == 0) {print_data[header.charge_index] = true;}
    else if (dummy.compare("x") == 0) {print_data[header.x_index] = true;}
    else if (dummy.compare("y") == 0) {print_data[header.y_index] = true;}
    else if (dummy.compare("z") == 0) {print_data[header.z_index] = true;}
    else if (dummy.compare("xu") == 0) {print_data[header.xu_index] = true;}
    else if (dummy.compare("yu") == 0) {print_data[header.yu_index] = true;}
    else if (dummy.compare("zu") == 0) {print_data[header.zu_index] = true;}
    else
    {
      if (header.indexConverter.find(dummy) != header.indexConverter.end())
      {
        print_data[header.indexConverter[dummy].first] == true;
      }
      else
      {
        cout << "Unrecognized data type specified.\n";
      }
    }
  }

  cout << "Please specify the output filename: ";
  cin  >> filename;

  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);

  fout << "# ";
  for (unsigned int i = 0; i < print_data.size(); ++i)
  {
    if (print_data[i]) {fout << header.data_types[i] << " ";}
  }

  fout << endl;

  for (unsigned int i = 0; i < selected_atoms.size(); ++i)
  {
    for (unsigned int j = 0; j < print_data.size(); ++j)
    {
      if (print_data[j])
      {
        if (j == header.id_index) {fout << selected_atoms[i].getId() << " ";}
        else if (j == header.type_index) {fout << selected_atoms[i].getType() << " ";}
        else if (j == header.charge_index) {fout << selected_atoms[i].getCharge() << " ";}
        else if (j == header.x_index) {fout << selected_atoms[i].getWrapped().getX() << " ";}
        else if (j == header.y_index) {fout << selected_atoms[i].getWrapped().getY() << " ";}
        else if (j == header.z_index) {fout << selected_atoms[i].getWrapped().getZ() << " ";}
        else if (j == header.xu_index) {fout << selected_atoms[i].getUnwrapped().getX() << " ";}
        else if (j == header.yu_index) {fout << selected_atoms[i].getUnwrapped().getY() << " ";}
        else if (j == header.zu_index) {fout << selected_atoms[i].getUnwrapped().getZ() << " ";}
        else
        {
          fout << selected_atoms[i].getExtraInfo()[header.indexConverter[header.data_types[j]].second] << " ";
        }
      }
    }
    fout << endl;
  }

  fout.close();

}

void playground(string filename = "None")
{
  bool leave_playground = false;
  string command;
  vector <Atom> atoms, selected_atoms;
  Header header;

  if (filename.compare("None") != 0)
  {
    pair<vector <Atom>, Header> tmp = readFile(filename);
    atoms = tmp.first;
    header = tmp.second;
  }

  while (!leave_playground)
  {
    cout << "<< ";
    cin  >> command;

    if (command.compare("quit") == 0 || command.compare("exit") == 0 || command.compare("q") == 0)
    {
      leave_playground = true;
    }
    else if (command.compare("help") == 0 || command.compare("h") == 0)
    {
      showCommandList();
    }
    else if (command.compare("read") == 0 || command.compare("r") == 0)
    {
      filename = getFileName();
      pair <vector <Atom>, Header> tmp = readFile(filename);
      atoms = tmp.first;
      header = tmp.second;
    }
    else if (command.compare("header") == 0 || command.compare("H") == 0)
    {
      if (atoms.size() == 0)
      {
        cout << "No atomic information has been read yet!  Specify a file to read with r or read.\n";
      }
      else {printHeaderInfo(header);}
    }
    else if (command.compare("select") == 0 || command.compare("s") == 0)
    {
      if (atoms.size() == 0)
      {
        cout << "No atomic information has been read yet!  Specify a file to read with r or read.\n";
      }
      else {selected_atoms = getSelectedAtoms(atoms);}
    }
    else if (command.compare("print") == 0 || command.compare("p") == 0)
    {
      if (selected_atoms.size() == 0) {selected_atoms = atoms;}
      printSelectedAtomInfo(selected_atoms, header);
    }
    else
    {
      cout << "Unrecognized command!\n";
    }
  }
}

int main(int argc, char** argv)
{
  string file;

  try
  {
    cxxopts::Options options(argv[0], "Playground for extracting and examining atom properties");
    options
      .positional_help("file")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Data file to play with", cxxopts::value<string>(file), "data_file")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    welcome();

    if (result.count("file"))
    {
      playground(file);
    }
    else
    {
      playground();
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

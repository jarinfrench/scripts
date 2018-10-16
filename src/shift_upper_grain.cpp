#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cxxopts.hpp>
#include "atom.h"
#include "error_code_defines.h"

using namespace std;

struct boxData
{
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double xy, xz, yz;
  double Lx, Ly, Lz;
  bool is_triclinic = false;

  double calculateBoxLengths()
  {
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;
  }
} box;

// Calculate the rounded value of x
double anInt(double x)
{
  // NOTE: This function is different from other version of this same function.
  // Here, we DO NOT round the value!  We take the floor() value!
  int temp; // temporary variable to hold the integer value of x
  temp = (int)(x);
  return (double)(temp);
}

bool compareAtomZ(const Atom& a, const Atom& b) {return a.getZ() < b.getZ();}
bool compareAtomCharge(const Atom& a, const Atom& b) {return a.getCharge() < b.getCharge();}

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

pair <int, int> getAtomData(const string& filename, vector <Atom>& atoms)
{
  string str;
  int N, ntypes;
  int n_total = 0;

  ifstream fin(filename.c_str());
  checkFileStream(fin, filename);

  getline(fin, str); // Header comment line
  fin >> N >> str;
  fin >> ntypes >> str >> str;
  fin >> box.xlow >> box.xhigh >> str >> str
      >> box.ylow >> box.yhigh >> str >> str
      >> box.zlow >> box.zhigh >> str >> str;

  box.calculateBoxLengths();

  fin.ignore();
  getline(fin,str); // read the extra stuff
  if (str.find("xy") != string::npos) // Triclinic system check
  {
    stringstream tmp(str);
    tmp >> box.xy >> box.xz >> box.yz >> str >> str >> str;
    fin >> str;
    box.is_triclinic = true;
    fin.ignore();
  }
  else
  {
    getline(fin, str); // required in order to line up the different input file structures
  }

  atoms.resize(N, Atom());

  getline(fin, str); // blank line before the data
  while (getline(fin, str))
  {
    vector <double> data;
    int atom_id, type;
    double charge, x, y, z;

    stringstream ss(str);
    double dummy;
    while (ss >> dummy) {data.push_back(dummy);}
    atom_id = (int)(data[0]);
    type = (int)(data[1]);

    switch (data.size())
    {
      //atom has charge
      case 6: charge = data[2];
              x = data[3]; y = data[4]; z = data[5];
              break;
      case 5: charge = 0.0;
              x = data[2]; y = data[3]; z = data[4];
              break;
      default: cout << "Unrecognized file format.  Expected format: id type charge* x y z. * --> optional.\n";
               exit(FILE_FORMAT_ERROR);
    }

    data.clear();
    ++n_total;
    if (type > ntypes)
    {
      cout << "Error: Atom type = " << type << " is greater than the number of types = " << ntypes << endl;
      exit(ATOM_TYPE_ERROR);
    }

    atoms[atom_id - 1] = Atom(atom_id, type, charge, x, y, z);
  }
  fin.close();

  if (n_total != N)
  {
    cout << "Error: n_total = " << n_total << " != N = " << N << endl;
    exit(ATOM_COUNT_ERROR);
  }

  return make_pair(N, ntypes);
}

void createShiftedBoundaries(const string& filename, const pair <int, int>& ncell,
                             const pair <int, int>& ngrid, vector <Atom>& atoms,
                             const pair <int, int>& num_atoms_num_types)
{
  string newFileName;
  int ntotal, num_types;
  int N = num_atoms_num_types.first, ntypes = num_atoms_num_types.second;
  double xshift, yshift, lcellx, lcelly; // total shift in the x/y direction, length of cells in x/y direction
  bool has_charge = false;
  vector <Atom> shifted_atoms;

  lcellx = box.Lx / ncell.first;
  lcelly = box.Ly / ncell.second;

  if ((*max_element(atoms.begin(), atoms.end(), compareAtomCharge)).getCharge() != 0.0) {has_charge = true;}

  // The atoms need to be sorted (from greatest to least) by their Z value
  sort(shifted_atoms.begin(), shifted_atoms.end(), compareAtomZ);

  // for each grid unit, create a separate structure
  for (int ix = 1; ix <= ngrid.first; ++ix)
  {
    for (int iy = 1; iy <= ngrid.second; ++iy)
    {
      ntotal = 0;

      // create the new filename
      stringstream ss;
      ss << filename.substr(0, filename.find(".dat")) << "_" << ix << "X" << iy << "Y.dat";
      ss >> newFileName;

      ofstream fout(newFileName.c_str());
      checkFileStream(fout, newFileName);

      fout << "These coordinates have the top grain shifted by " << ix << "X" << iy << "Y: [ID type ";
      if (has_charge) {fout << "charge ";}
      fout << "x y z]\n\n"
           << N << " atoms\n"
           << ntypes << " atom types\n"
           << box.xlow << " " << box.xhigh << " xlo xhi\n"
           << box.ylow << " " << box.yhigh << " ylo yhi\n"
           << box.zlow << " " << box.zhigh << " zlo zhi\n";
      if (box.is_triclinic) // NOTE: This code may not work with triclinic systems
      {
        cout << "Warning: This code may not work correctly with triclinic systems!\n";
        fout << box.xy << " " << box.xz << " " << box.yz << " xz xy yz\n";
      }
      fout << "\nAtoms\n\n";

      shifted_atoms = atoms;
      for (unsigned int i = 0; i < shifted_atoms.size(); ++i)
      {
        ++ntotal;

        if (shifted_atoms[i].getType() > num_types) {num_types = shifted_atoms[i].getType();}
        if (shifted_atoms[i].getType() > ntypes)
        {
          cout << "Error: atom type is greater than expected. atom_type = "
               << shifted_atoms[i].getType() << " > ntypes = " << ntypes << endl;
          exit(ATOM_TYPE_ERROR);
        }

        if (ntotal >= N / 2.0)
        {
          xshift = ((double)(ix) - 1.0) / (double)(ngrid.first - 1) * lcellx;
          yshift = ((double)(iy) - 1.0) / (double)(ngrid.second - 1) * lcelly;
        }
        else
        {
          xshift = 0.0;
          yshift = 0.0;
        }

        shifted_atoms[i].setX(shifted_atoms[i].getX() + xshift);
        shifted_atoms[i].setY(shifted_atoms[i].getY() + yshift);

        // Apply PBCs
        shifted_atoms[i].setX(shifted_atoms[i].getX() - anInt(shifted_atoms[i].getX() / box.Lx) * box.Lx);
        shifted_atoms[i].setY(shifted_atoms[i].getY() - anInt(shifted_atoms[i].getY() / box.Ly) * box.Ly);
        // We're not shifting in the z direction, so we don't need to worry about PBCs.

        if (shifted_atoms[i].getX() < box.xlow || shifted_atoms[i].getX() > box.xhigh)
        {
          cout << "Error: atom " << ntotal << " is outside the bounds set by xlow ("
               << box.xlow << ") and xhigh (" << box.xhigh << "): "
               << shifted_atoms[i].getX() << endl;
          exit(BOUNDS_ERROR);
        }

        if (shifted_atoms[i].getY() < box.ylow || shifted_atoms[i].getY() > box.yhigh)
        {
          cout << "Error: atom " << ntotal << " is outside the bounds set by ylow ("
               << box.ylow << ") and yhigh (" << box.yhigh << "): "
               << shifted_atoms[i].getY() << endl;
          exit(BOUNDS_ERROR);
        }

        if (shifted_atoms[i].getZ() < box.zlow || shifted_atoms[i].getZ() > box.zhigh)
        {
          cout << "Error: atom " << ntotal << " is outside the bounds set by zlow ("
               << box.zlow << ") and zhigh (" << box.zhigh << "): "
               << shifted_atoms[i].getZ() << endl;
          exit(BOUNDS_ERROR);
        }

        fout << shifted_atoms[i].getId() << " " << shifted_atoms[i].getType() << " ";

        if (has_charge) {fout << shifted_atoms[i].getCharge() << " ";}

        fout << shifted_atoms[i].getX() << " " << shifted_atoms[i].getY() << " "
             << shifted_atoms[i].getZ() << endl;
      }

      fout.close();

    }
  }
}

int main(int argc, char** argv)
{
  string filename;
  pair <int, int> ncell, ngrid;
  vector <Atom> atoms;

  try
  {
    cxxopts::Options options(argv[0], "Creates a series of grain boundary structures for gamma surface mapping.");
    options
      .positional_help("file ncell_x ncell_y ngrid_x ngrid_y")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Original grain boundary structure file", cxxopts::value<string>(filename), "file")
        ("ncell_x", "Number of cells in the x direction", cxxopts::value<int>(ncell.first), "n")
        ("ncell_y", "Number of cells in the y direction", cxxopts::value<int>(ncell.second), "n")
        ("ngrid_x", "Number of grid lines in the x direction", cxxopts::value<int>(ngrid.first), "n")
        ("ngrid_y", "Number of grid lines in the y direction", cxxopts::value<int>(ngrid.second), "n")
        ("h,help", "Show the help");

    options.parse_positional({"file", "ncell_x", "ncell_y", "ngrid_x", "ngrid_y"});
    auto result = options.parse(argc, argv);

    bool has_required_params = false;;
    if (result.count("file") && result.count("ncell_x") && result.count("ncell_y") &&
        result.count("ngrid_x") && result.count("ngrid_y"))
    {
      has_required_params = true;
    }

    if (result.count("help") || !has_required_params)
    {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (has_required_params)
    {
      pair <int, int> num_atoms_num_types;
      num_atoms_num_types = getAtomData(filename, atoms);
      createShiftedBoundaries(filename, ncell, ngrid, atoms, num_atoms_num_types);
    }

  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

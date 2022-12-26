#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"

using namespace std;

static string progress = "|>                                                  |";

struct boxData {
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double xy, xz, yz;
  double Lx, Ly, Lz;
  bool is_triclinic = false;

  void calculateBoxLengths() {
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;
  }
} box;

void showProgress(const int& current_iter, const int& max_iters) {
  // We increase the number of equal signs every 2% of the max iterations

  int num_replace = current_iter / (max_iters * 0.02); // the number of = replacements we make
  for (int i = 0; i < num_replace; ++i) {
    progress.replace(i+1, 2,"=>");
  }

  double tmp = (double)(current_iter) / (double)(max_iters) * 100.0;
  cout << progress << " (" << setprecision(4) << tmp << "%)    \r" << flush;

  if (current_iter == max_iters) {
    cout << progress << " (100%)" << endl;
  }
}

// Calculate the rounded value of x
double anInt(double x) {
  int temp; // temporary variable to hold the integer value of x
  // if (x > 0.0) {x += 0.5;}
  // if (x < 0.0) {x -= 0.5;}
  temp = (int)(x);
  return (double)(temp);
}

bool compareAtomX(const Atom& a, const Atom& b) {return a.getWrapped()[0] < b.getWrapped()[0];}
bool compareAtomY(const Atom& a, const Atom& b) {return a.getWrapped()[1] < b.getWrapped()[1];}
bool compareAtomZ(const Atom& a, const Atom& b) {return a.getWrapped()[2] < b.getWrapped()[2];}
bool compareAtomCharge(const Atom& a, const Atom& b) {return a.getCharge() < b.getCharge();}

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

pair <int, int> getAtomData(const string& filename, vector <Atom>& atoms) {
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
  if (str.find("xy") != string::npos) { // Triclinic system check
    stringstream tmp(str);
    tmp >> box.xy >> box.xz >> box.yz >> str >> str >> str;
    fin >> str;
    box.is_triclinic = true;
    fin.ignore();
  } else {
    getline(fin, str); // required in order to line up the different input file structures
  }

  atoms.resize(N, Atom());

  getline(fin, str); // blank line before the data
  while (getline(fin, str)) {
    vector <double> data;
    int atom_id, type;
    double charge, x, y, z;

    stringstream ss(str);
    double dummy;
    while (ss >> dummy) data.push_back(dummy);
    atom_id = (int)(data[0]);
    type = (int)(data[1]);

    switch (data.size()) {
      //atom has charge
      case 6: charge = data[2];
              x = data[3]; y = data[4]; z = data[5];
              break;
      case 5: charge = 0.0;
              x = data[2]; y = data[3]; z = data[4];
              break;
      default: cerr << "Unrecognized file format.  Expected format: id type charge* x y z. * --> optional.\n";
               exit(FILE_FORMAT_ERROR);
    }

    data.clear();
    ++n_total;
    if (type > ntypes) {
      cerr << "Error: Atom type = " << type << " is greater than the number of types = " << ntypes << endl;
      exit(ATOM_TYPE_ERROR);
    }

    Position p(x,y,z);
    atoms[atom_id - 1] = Atom(atom_id, type, charge, p);
  }
  fin.close();

  if (n_total != N) {
    cerr << "Error: n_total = " << n_total << " != N = " << N << endl;
    exit(ATOM_COUNT_ERROR);
  }

  return make_pair(N, ntypes);
}

void createShiftedBoundaries(const string& filename, const vector <double>& shifts,
                             const vector <int>& maxes, vector <Atom>& atoms,
                             const pair <int, int>& num_atoms_num_types,
                             const bool& is_specific, const string gb_plane)
{
  string newFileName;
  int ntotal;
  int N = num_atoms_num_types.first, ntypes = num_atoms_num_types.second;
  bool has_charge = false;
  vector <Atom> shifted_atoms;

  cout << "Total shift distance: ";
  for (size_t i = 0; i < shifts.size(); ++i) {
    cout << shifts[i] * maxes[i];
    if (i == 0) cout << " in x, ";
    if (i == 1) cout << " in y, ";
    if (i == 2) cout << " in z\n";
  }

  if ((*max_element(atoms.begin(), atoms.end(), compareAtomCharge)).getCharge() != 0.0 ||
      (*min_element(atoms.begin(), atoms.end(), compareAtomCharge)).getCharge() != 0.0) has_charge = true;

  int iter = 0;
  int max_iters;
  int iz_start;
  int ix_start;
  int iy_start;
  if (gb_plane.compare("xy") == 0) {
    (is_specific) ? max_iters = 1 : max_iters = (maxes[0] + 1) * maxes[1] * (maxes[2] + 1); // the x and z iterations start from 0, so we add 1
    ix_start = 0;
    iy_start = 1;
    iz_start = 0;
  } else if (gb_plane.compare("xz") == 0) {
    (is_specific) ? max_iters = 1 : max_iters = maxes[0] * (maxes[1] + 1) * (maxes[2] + 1);// the y and z iterations start from 0, so we add 1
    ix_start = 1;
    iy_start = 0;
    iz_start = 0;
  } else { //gb_plane.compare("yz") == 0
    (is_specific) ? max_iters = 1 : max_iters = (maxes[0] + 1) * (maxes[1] + 1) * maxes[2];// the x and y iterations start from 0, so we add 1
    ix_start = 0;
    iy_start = 0;
    iz_start = 1;
  }


  // for each grid unit, create a separate structure
  for (int iz = iz_start; iz <= maxes[2]; ++iz) { // This should always run through the outer loop at least once
    for (int ix = ix_start; ix <= maxes[0]; ++ix) {
      for (int iy = iy_start; iy <= maxes[1]; ++iy) {
        if (is_specific) {
          ix = maxes[0];
          iy = maxes[1];
          iz = maxes[2];
        }
        ntotal = 0;

        // create the new filename
        stringstream ss;
        size_t ext_pos = filename.find(".dat") == string::npos ? filename.find(".lmp") : filename.find(".dat");

        ss << filename.substr(0, ext_pos) << "_" << ix << "X" << iy << "Y" << iz << "Z.dat";
        ss >> newFileName;

        ofstream fout(newFileName.c_str());
        checkFileStream(fout, newFileName);

        fout << "These coordinates have the top grain shifted by " << ix * shifts[0]
             << " in x, " << iy * shifts[1] << " in y, and " << iz * shifts[2] << " in z: [ID type ";

        if (has_charge) fout << "charge ";
        fout << "x y z]\n\n"
             << N << " atoms\n"
             << ntypes << " atom types\n";
        if (gb_plane.compare("xy") == 0) { // Accomodate adding space between the grains
          fout << box.xlow << " " << box.xhigh << " xlo xhi\n"
               << box.ylow << " " << box.yhigh << " ylo yhi\n"
               << box.zlow << " " << box.zhigh + (iz * shifts[2]) << " zlo zhi\n";
        } else if (gb_plane.compare("xz") == 0) {
          fout << box.xlow << " " << box.xhigh << " xlo xhi\n"
               << box.ylow << " " << box.yhigh + (iy  * shifts[1]) << " ylo yhi\n"
               << box.zlow << " " << box.zhigh << " zlo zhi\n";
        } else { // gb_plane.compare("yz") == 0
        fout << box.xlow << " " << box.xhigh + (ix * shifts[0]) << " xlo xhi\n"
             << box.ylow << " " << box.yhigh << " ylo yhi\n"
             << box.zlow << " " << box.zhigh << " zlo zhi\n";
        }
        if (box.is_triclinic) { // NOTE: This code may not work with triclinic systems
          cout << "Warning: This code may not work correctly with triclinic systems!\n";
          fout << box.xy << " " << box.xz << " " << box.yz << " xz xy yz\n";
        }
        fout << "\nAtoms\n\n";

        shifted_atoms = atoms;

        // The atoms need to be sorted (from greatest to least) by the coordinate normal to the GB plane
        if (gb_plane.compare("xy") == 0) {
          sort(shifted_atoms.begin(), shifted_atoms.end(), compareAtomZ);
        } else if (gb_plane.compare("xz") == 0) {
          sort(shifted_atoms.begin(), shifted_atoms.end(), compareAtomY);
        } else { // gb_plane.compare("yz") == 0
          sort(shifted_atoms.begin(), shifted_atoms.end(), compareAtomX);
        }

        for (unsigned int i = 0; i < shifted_atoms.size(); ++i) {
          ++ntotal;

          if (shifted_atoms[i].getType() > ntypes) {
            cerr << "Error: atom type is greater than expected. atom_type = "
                 << shifted_atoms[i].getType() << " > ntypes = " << ntypes << endl;
            exit(ATOM_TYPE_ERROR);
          }

          // Shift the atoms according to the specified GB plane
          // For any shifts normal to the GB plane, PBCs are not taken into account
          if (gb_plane.compare("xy") == 0) {
            if (shifted_atoms[i].getWrapped()[2] >= (box.zhigh - box.zlow) / 2.0) {
              Position p = shifted_atoms[i].getWrapped();
              p += Position(ix * shifts[0], iy * shifts[1], 0.0);
              shifted_atoms[i].setWrapped(p);

              // Apply PBCs
              p.setX(shifted_atoms[i].getWrapped()[0] - anInt(shifted_atoms[i].getWrapped()[0] / box.Lx) * box.Lx);
              p.setY(shifted_atoms[i].getWrapped()[1] - anInt(shifted_atoms[i].getWrapped()[1] / box.Ly) * box.Ly);

              shifted_atoms[i].setWrapped(p);

              if (shifted_atoms[i].getWrapped()[2] < box.zlow || shifted_atoms[i].getWrapped()[2] > box.zhigh + (iz * shifts[2])) {
                cerr << "Error: atom " << ntotal << " is outside the bounds set by zlow ("
                     << box.zlow << ") and zhigh (" << box.zhigh << "): "
                     << shifted_atoms[i].getWrapped()[2] << endl;
                exit(BOUNDS_ERROR);
              }
            }
          } else if (gb_plane.compare("xz") == 0) {
            if (shifted_atoms[i].getWrapped()[1] >= (box.yhigh - box.ylow) / 2.0) {
              Position p = shifted_atoms[i].getWrapped();
              p += Position(ix * shifts[0], 0.0, iz * shifts[2]);
              shifted_atoms[i].setWrapped(p);

              // Apply PBCs
              p.setX(shifted_atoms[i].getWrapped()[0] - anInt(shifted_atoms[i].getWrapped()[0] / box.Lx) * box.Lx);
              p.setZ(shifted_atoms[i].getWrapped()[2] - anInt(shifted_atoms[i].getWrapped()[2] / box.Lz) * box.Lz);

              shifted_atoms[i].setWrapped(p);

              if (shifted_atoms[i].getWrapped()[1] < box.ylow || shifted_atoms[i].getWrapped()[1] > box.yhigh) {
                cerr << "Error: atom " << ntotal << " is outside the bounds set by ylow ("
                     << box.ylow << ") and yhigh (" << box.yhigh << "): "
                     << shifted_atoms[i].getWrapped()[1] << endl;
                exit(BOUNDS_ERROR);
              }
            }
          } else { // gb_plane.compare("yz") == 0
            if (shifted_atoms[i].getWrapped()[0] >= (box.xhigh - box.xlow) / 2.0) {
              Position p = shifted_atoms[i].getWrapped();
              p += Position(0.0, iy * shifts[1], iz * shifts[2]);
              shifted_atoms[i].setWrapped(p);

              // Apply PBCs
              p.setY(shifted_atoms[i].getWrapped()[1] - anInt(shifted_atoms[i].getWrapped()[1] / box.Ly) * box.Ly);
              p.setZ(shifted_atoms[i].getWrapped()[2] - anInt(shifted_atoms[i].getWrapped()[2] / box.Lz) * box.Lz);

              shifted_atoms[i].setWrapped(p);

              if (shifted_atoms[i].getWrapped()[0] < box.xlow || shifted_atoms[i].getWrapped()[0] > box.xhigh) {
                cerr << "Error: atom " << ntotal << " is outside the bounds set by xlow ("
                     << box.xlow << ") and xhigh (" << box.xhigh << "): "
                     << shifted_atoms[i].getWrapped()[0] << endl;
                exit(BOUNDS_ERROR);
              }
            }
          }

          fout << shifted_atoms[i].getId() << " " << shifted_atoms[i].getType() << " ";

          if (has_charge) fout << shifted_atoms[i].getCharge() << " ";

          fout << shifted_atoms[i].getWrapped()[0] << " " << shifted_atoms[i].getWrapped()[1] << " "
               << shifted_atoms[i].getWrapped()[2] << endl;
        } // shift each atom
        fout.close();
        ++iter;
      } // yshifts
      showProgress(iter, max_iters);
    } // xshifts
  }// zshifts
}

int main(int argc, char** argv)
{
  string filename;
  string gb_plane = "xy";
  vector <double> shifts;
  vector <int> maxes;
  double shift_x, shift_y, shift_z, shift1, shift2;
  int max_x, max_y, max_z, max1, max2;
  vector <Atom> atoms;
  bool is_specific = false;

  try {
    cxxopts::Options options(argv[0], "Creates a series of grain boundary structures for gamma surface mapping.");
    options
      .positional_help("file shift1 shift2 max1 max2")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Original grain boundary structure file", cxxopts::value<string>(filename), "file")
        ("shift1", "Shift size in the first direction", cxxopts::value<double>(shift1), "value")
        ("shift2", "Shift size in the second direction", cxxopts::value<double>(shift2), "value")
        ("max1", "Maximum number of displacements in the first direction", cxxopts::value<int>(max1), "n")
        ("max2", "Maximum number of displacements in the second direction", cxxopts::value<int>(max2), "n")
        ("x,x-shift", "Shift size in the x direction", cxxopts::value<double>(shift_x), "value")
        ("y,y-shift", "Shift size in the y direction", cxxopts::value<double>(shift_y), "value")
        ("z,z-shift", "Shift size in the z direction", cxxopts::value<double>(shift_z), "value")
        ("x-max", "Maximum number of displacements in x direction", cxxopts::value<int>(max_x), "n")
        ("y-max", "Maximum number of displacements in y direction", cxxopts::value<int>(max_y), "n")
        ("z-max", "Maximum number of displacements in z direction", cxxopts::value<int>(max_z), "n")
        ("p,plane", "The plane the GB is in (default = 'xy')", cxxopts::value<string>(gb_plane)->default_value("xy"), "[xy, xz, yz]")
        ("specific", "Only create the specified shift", cxxopts::value<bool>(is_specific)->default_value("false")->implicit_value("true"))
        ("h,help", "Show the help");

    options.parse_positional({"file", "shift1", "shift2", "max1", "max2"});
    auto result = options.parse(argc, argv);

    bool has_required_params = false;
    if (result.count("shift1") && result.count("shifts2") && result.count("max1") && result.count("max2")) has_required_params = true;

    if (gb_plane.compare("xy") == 0) { // xy GB plane
      if (result.count("shift1")) shift_x = shift1;
      if (result.count("shift2")) shift_y = shift2;
      if (result.count("max1")) max_x = max1;
      if (result.count("max2")) max_y = max2;

      if ((result.count("z-shift") && !result.count("z-max")) || (!result.count("z-shift") && result.count("z-max"))) {
        cerr << "Error: both the z shift and the max number of shifts in z must be specified.\n";
        return INPUT_FORMAT_ERROR;
      }

      if (!result.count("z-shift")) shift_z = 0.0;
      if (!result.count("z-max")) max_z = 0;
    } else if (gb_plane.compare("xz") == 0) { // xz GB plane
      if (result.count("shift1")) shift_x = shift1;
      if (result.count("shift2")) shift_z = shift2;
      if (result.count("max1")) max_x = max1;
      if (result.count("max2")) max_z = max2;

      if ((result.count("y-shift") && !result.count("y-max")) || (!result.count("y-shift") && result.count("y-max"))) {
        cerr << "Error: both the y shift and the max number of shifts in y must be specified.\n";
        return INPUT_FORMAT_ERROR;
      }

      if (!result.count("y-shift")) shift_y = 0.0;
      if (!result.count("y-max")) max_y = 0;
    } else if (gb_plane.compare("yz") == 0) { // yz GB plane
      if (result.count("shift1")) shift_y = shift1;
      if (result.count("shift2")) shift_z = shift2;
      if (result.count("max1")) max_y = max1;
      if (result.count("max2")) max_z = max2;

      if ((result.count("x-shift") && !result.count("x-max")) || (!result.count("x-shift") && result.count("x-max"))) {
        cerr << "Error: both the x shift and the max number of shifts in x must be specified.\n";
        return INPUT_FORMAT_ERROR;
      }

      if (!result.count("x-shift")) shift_x = 0.0;
      if (!result.count("x-max")) max_x = 0;
    } else {
      cerr << "Specified boundary plane must be one of 'xy', 'xz', or 'yz'\n";
      exit(INPUT_FORMAT_ERROR);
    }

    if (result.count("file")) {
      has_required_params = true;
      shifts.push_back(shift_x);
      shifts.push_back(shift_y);
      shifts.push_back(shift_z);

      maxes.push_back(max_x);
      maxes.push_back(max_y);
      maxes.push_back(max_z);
    }

    if (result.count("help") || !has_required_params) {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (has_required_params) {
      pair <int, int> num_atoms_num_types;
      num_atoms_num_types = getAtomData(filename, atoms);
      createShiftedBoundaries(filename, shifts, maxes, atoms, num_atoms_num_types, is_specific, gb_plane);
    }

  } catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

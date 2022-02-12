#include <iostream>
#include <cstdlib> // for getenv
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <sys/stat.h> // for stat
#include <cmath> // for sin, cos
#include <numeric> // for iota
#include <algorithm> // for min, max_element, fill

#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"
#include "verifyNewFile.h"

using namespace std;

#define PI 3.141592653589793

static bool show_warnings = true;
static bool has_charge = false;
static bool print_nearest_neighbors = false;
static const map <string, vector <Position> > first_nn_list = {
  {"fcc", // 12 1st nearest neighbors
    {Position(0.0, 0.5, -0.5), Position(0.0, -0.5, 0.5), Position( 0.0, -0.5, -0.5), Position( 0.0,  0.5, 0.5),
     Position(0.5, 0.0, -0.5), Position(0.5,  0.0, 0.5), Position(-0.5,  0.0, -0.5), Position(-0.5,  0.0, 0.5),
     Position(0.5, 0.5,  0.0), Position(0.5, -0.5, 0.0), Position(-0.5,  0.5,  0.0), Position(-0.5, -0.5, 0.0)}
  },
  {"bcc", // 8 1st nearest neighbors
    {Position(-0.5, -0.5, -0.5), Position(-0.5, -0.5, 0.5), Position(0.5, -0.5, -0.5), Position(0.5, -0.5,  0.5),
     Position(-0.5,  0.5, -0.5), Position(-0.5,  0.5, 0.5), Position(0.5,  0.5, -0.5), Position(0.5,  0.5, -0.5)}
  },
  {"sc", // 6 1st nearest neighbors (incidentally, also works for the 2nd nearest neighbors for fcc and bcc)
    {Position(-1.0, 0.0, 0.0), Position(1.0, 0.0,  0.0), Position(0.0, -1.0, 0.0),
     Position( 0.0, 1.0, 0.0), Position(0.0, 0.0, -1.0), Position(0.0,  0.0, 1.0)}
  }
};
static const map <string, vector <vector <double>> > basis_vectors = {
  {"fcc",
    {{0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}}
  },
  {"bcc",
    {{-0.5, 0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, -0.5}}
  }
};
static map <int, vector <Position> > rotated_1nn;

bool neighborComparison(pair <int, double> a, pair <int, double> b) {
  return a.second < b.second;
}

double calculateAxisMagnitude(const vector <double>& axis) {
  return sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
}

vector <double> normalize(const vector <double>& axis) {
  vector <double> result (3, 0.0);
  if (axis.size() != 3) {
    cerr << "Vector size is incorrect (" << axis.size() << " != 3). Cannot normalize.\n";
    exit(VECTOR_SIZE_ERROR);
  }

  double magnitude = calculateAxisMagnitude(axis);

  for (size_t i = 0; i < axis.size(); ++i) {result[i] = axis[i] / magnitude;}

  return result;
}

struct inputData {
  // directly from the input file
  int n_files, n_types;
  double theta, r_cut, fi_cut, a0;
  string crystal_structure;
  vector <double> new_x_axis {0.0,0.0,0.0}, new_y_axis{0.0,0.0,0.0}, new_z_axis{0.0,0.0,0.0};
  vector <double> old_x_axis {0.0,0.0,0.0}, old_y_axis{0.0,0.0,0.0}, old_z_axis{0.0,0.0,0.0};
  vector <string> files; // the data files to examine
  string r_axis; // the string representation of the orientation matrix, given as zzz

  // from the command line
  vector <int> ignored_atoms, atom_ids_list;
  char algorithm; // which algorithm to use
  string outfile, nn_filebase;
  unsigned int print_files_every = 0;

  // calculated from input file
  double sinin, cosin, janssens_symm; // sin and cos of the angle, the janssens symmetry parameter for ideal lattices.
  vector <Position> perfect_matrix_neighbors, perfect_grain_neighbors;
  vector <vector <double> > basis_1 {{0.0, 0.0, 0.0}, {0.0,0.0,0.0}, {0.0, 0.0, 0.0}};
  vector <vector <double> > basis_2 = basis_1, reciprocal_1 = basis_1, reciprocal_2 = basis_1;
  double scalenorm = 0;

  void calculateTrig() {
    sinin = sin(theta * PI / 180.0);
    cosin = cos(theta * PI / 180.0);
  }

  void calculateOldAxes() {
    old_x_axis[0] = new_x_axis[0];
    old_x_axis[1] = new_y_axis[0];
    old_x_axis[2] = new_z_axis[0];

    old_y_axis[0] = new_x_axis[1];
    old_y_axis[1] = new_y_axis[1];
    old_y_axis[2] = new_z_axis[1];

    old_z_axis[0] = new_x_axis[2];
    old_z_axis[1] = new_y_axis[2];
    old_z_axis[2] = new_z_axis[2];
  }

  // for Janssens orientation parameter
  void calculatePerfectNeighbors() {
    double x, y, z, rxij, ryij, rzij;
    perfect_matrix_neighbors = first_nn_list.at(crystal_structure);
    perfect_grain_neighbors = first_nn_list.at(crystal_structure);

    // expand the neighbor vectors to be the same size as the crystal structure,
    // and rotate the second matrix

    for (size_t i = 0; i < perfect_matrix_neighbors.size(); ++i) {
      x = perfect_matrix_neighbors[i][0] * a0;
      y = perfect_matrix_neighbors[i][1] * a0;
      z = perfect_matrix_neighbors[i][2] * a0;

      // rotate to same orientation
      perfect_matrix_neighbors[i][0] = new_x_axis[0] * x + new_x_axis[1] * y + new_x_axis[2] * z;
      perfect_matrix_neighbors[i][1] = new_y_axis[0] * x + new_y_axis[1] * y + new_y_axis[2] * z;
      perfect_matrix_neighbors[i][2] = new_z_axis[0] * x + new_z_axis[1] * y + new_z_axis[2] * z;

      perfect_grain_neighbors[i] = perfect_matrix_neighbors[i];

      // now rotate by the misorientation
      x = perfect_grain_neighbors[i][0];
      y = perfect_grain_neighbors[i][1];

      perfect_grain_neighbors[i][0] = x * cosin - y * sinin;
      perfect_grain_neighbors[i][1] = y * cosin + x * sinin;
    }

    vector <double> distances;
    janssens_symm = 0.0;
    for (size_t i = 0; i < perfect_matrix_neighbors.size(); ++i) {
      x = perfect_matrix_neighbors[i][0];
      y = perfect_matrix_neighbors[i][1];
      z = perfect_matrix_neighbors[i][2];
      distances.clear();
      for (size_t j = 0; j < perfect_grain_neighbors.size(); ++j) {
        rxij = perfect_grain_neighbors[i][0] - x;
        ryij = perfect_grain_neighbors[i][1] - y;
        rzij = perfect_grain_neighbors[i][2] - z;

        distances.push_back(rxij * rxij + ryij * ryij + rzij * rzij);
      }

      janssens_symm += sqrt(*min_element(distances.begin(), distances.end()));
    }
    janssens_symm /= perfect_matrix_neighbors.size();

  }

  // for Ulomek orientation parameter
  void calculateReciprocalLattices() {
    // hard coded numbers because these sizes do not change - they form the 3x3
    // orientation matrix of the two grains
    vector <vector <double> > rot_mat = {
      {cosin, -sinin, 0.0},
      {sinin, cosin, 0.0},
      {0.0, 0.0, 1.0}
    };

    vector <vector <double> > orient_mat_1 = {
      {normalize(new_x_axis)},
      {normalize(new_y_axis)},
      {normalize(new_z_axis)}
    };

    vector <vector <double> > orient_mat_2 (3, vector <double> (3, 0.0));
    for (size_t i = 0; i < orient_mat_2.size(); ++i) {
      for (size_t j = 0; j < orient_mat_2[i].size(); ++j) {
        orient_mat_2[i][j] = rot_mat[i][0] * orient_mat_1[0][j] + rot_mat[i][1] * orient_mat_1[1][j] + rot_mat[i][2] * orient_mat_1[2][j];
      }
    }

    for (size_t row = 0; row < 3; ++row) {
      for (size_t col = 0; col < 3; ++col) {
        // The correct result is the transpose of the matrix multiplication, hence the [col][row] indexing for basis_1
        basis_1[col][row] = orient_mat_1[row][0] * basis_vectors.at(crystal_structure)[0][col] * a0 +
                            orient_mat_1[row][1] * basis_vectors.at(crystal_structure)[1][col] * a0 +
                            orient_mat_1[row][2] * basis_vectors.at(crystal_structure)[2][col] * a0;
      }
    }

    for (size_t row = 0; row < 3; ++row) {
      for (size_t col = 0; col < 3; ++col) {
        // Need to transpose the basis_1 vector before multiplying by the rotation matrix
        // then transposing that result as the answer for basis_2, hence the [col][row] indexing
        // for both basis matrices.
        basis_2[col][row] = rot_mat[row][0] * basis_1[col][0] +
                            rot_mat[row][1] * basis_1[col][1] +
                            rot_mat[row][2] * basis_1[col][2];
      }
    }

    // This code is adapted from LAMMPS fix_eco_force.cpp.
    double b1_vol = 0.5 / PI *
              (basis_1[0][0] * (basis_1[1][1] * basis_1[2][2] - basis_1[1][2] * basis_1[2][1]) +
              (basis_1[0][1] * (basis_1[1][2] * basis_1[2][0] - basis_1[1][0] * basis_1[2][2])) +
              (basis_1[0][2] * (basis_1[1][0] * basis_1[2][1] - basis_1[1][1] * basis_1[2][0])));
    double b2_vol = 0.5 / PI *
              (basis_2[0][0] * (basis_2[1][1] * basis_2[2][2] - basis_2[1][2] * basis_2[2][1]) +
              (basis_2[0][1] * (basis_2[1][2] * basis_2[2][0] - basis_2[1][0] * basis_2[2][2])) +
              (basis_2[0][2] * (basis_2[1][0] * basis_2[2][1] - basis_2[1][1] * basis_2[2][0])));

    reciprocal_1[0][0] = (basis_1[1][1] * basis_1[2][2] - basis_1[1][2] * basis_1[2][1]) / b1_vol;
    reciprocal_1[0][1] = (basis_1[1][2] * basis_1[2][0] - basis_1[1][0] * basis_1[2][2]) / b1_vol;
    reciprocal_1[0][2] = (basis_1[1][0] * basis_1[2][1] - basis_1[1][1] * basis_1[2][0]) / b1_vol;
    reciprocal_1[1][0] = (basis_1[2][1] * basis_1[0][2] - basis_1[2][2] * basis_1[0][1]) / b1_vol;
    reciprocal_1[1][1] = (basis_1[2][2] * basis_1[0][0] - basis_1[2][0] * basis_1[0][2]) / b1_vol;
    reciprocal_1[1][2] = (basis_1[2][0] * basis_1[0][1] - basis_1[2][1] * basis_1[0][0]) / b1_vol;
    reciprocal_1[2][0] = (basis_1[0][1] * basis_1[1][2] - basis_1[0][2] * basis_1[1][1]) / b1_vol;
    reciprocal_1[2][1] = (basis_1[0][2] * basis_1[1][0] - basis_1[0][0] * basis_1[1][2]) / b1_vol;
    reciprocal_1[2][2] = (basis_1[0][0] * basis_1[1][1] - basis_1[0][1] * basis_1[1][0]) / b1_vol;

    reciprocal_2[0][0] = (basis_2[1][1] * basis_2[2][2] - basis_2[1][2] * basis_2[2][1]) / b2_vol;
    reciprocal_2[0][1] = (basis_2[1][2] * basis_2[2][0] - basis_2[1][0] * basis_2[2][2]) / b2_vol;
    reciprocal_2[0][2] = (basis_2[1][0] * basis_2[2][1] - basis_2[1][1] * basis_2[2][0]) / b2_vol;
    reciprocal_2[1][0] = (basis_2[2][1] * basis_2[0][2] - basis_2[2][2] * basis_2[0][1]) / b2_vol;
    reciprocal_2[1][1] = (basis_2[2][2] * basis_2[0][0] - basis_2[2][0] * basis_2[0][2]) / b2_vol;
    reciprocal_2[1][2] = (basis_2[2][0] * basis_2[0][1] - basis_2[2][1] * basis_2[0][0]) / b2_vol;
    reciprocal_2[2][0] = (basis_2[0][1] * basis_2[1][2] - basis_2[0][2] * basis_2[1][1]) / b2_vol;
    reciprocal_2[2][1] = (basis_2[0][2] * basis_2[1][0] - basis_2[0][0] * basis_2[1][2]) / b2_vol;
    reciprocal_2[2][2] = (basis_2[0][0] * basis_2[1][1] - basis_2[0][1] * basis_2[1][0]) / b2_vol;

    double dijx, dijy, dijz, dikx, diky, dikz, rsqij, rsqik, scalarprod, wij, wik;
    int layer = 4; // Not sure where this comes from
    double cut_sq = r_cut * r_cut * a0 * a0;
    int neigh = 0;

    for (int ix = -layer; ix <= layer; ++ix) {
      for (int iy = -layer; iy <= layer; ++iy) {
        for (int iz = -layer; iz <= layer; ++iz) {
          dijx = ix * orient_mat_2[0][0] + iy * orient_mat_2[1][0] + iz * orient_mat_2[2][0];
          dijy = ix * orient_mat_2[0][1] + iy * orient_mat_2[1][1] + iz * orient_mat_2[2][1];
          dijz = ix * orient_mat_2[0][2] + iy * orient_mat_2[1][2] + iz * orient_mat_2[2][2];
          rsqij = (dijx * dijx + dijy * dijy + dijz * dijz) / cut_sq;
          if (rsqij > 1e-8 && rsqij < 1.0) {
            ++neigh;
            for (int ia = -layer; ia <= layer; ++ia) {
              for (int ib = -layer; ib <= layer; ++ib) {
                for (int ic = -layer; ic <= layer; ++ic) {
                  dikx = ia * orient_mat_2[0][0] + ib * orient_mat_2[1][0] + ic * orient_mat_2[2][0];
                  diky = ia * orient_mat_2[0][1] + ib * orient_mat_2[1][1] + ic * orient_mat_2[2][1];
                  dikz = ia * orient_mat_2[0][2] + ib * orient_mat_2[1][2] + ic * orient_mat_2[2][2];
                  rsqik = (dikx * dikx + diky * diky + dikz * dikz) / cut_sq;
                  if (rsqik > 1e-8 && rsqik < 1.0) {
                    for (int alpha = 0; alpha < 3; ++alpha) {
                      scalarprod = reciprocal_1[alpha][0] * (dijx - dikx) + reciprocal_1[alpha][1] * (dijy - diky) + reciprocal_1[alpha][2] * (dijz - dikz);
                      wij = rsqij * (rsqij - 2) + 1;
                      wik = rsqik * (rsqik - 2) + 1;
                      scalenorm += wij * wik * (1 - cos(scalarprod));
                    } // alpha
                  } // rsqik
                } // ic
              } // ib
            } // ia
          } // rsqij
        } // iz
      } // iy
    } // ix

    cout << "  Number of neighbors: " << neigh << "\n";
  }

  template <typename T>
  string generateNNFilename(const T& id) const {
    if (nn_filebase.find("*") == string::npos) {return nn_filebase;}

    // This is limited in the sense that we only replace the first occurence
    // of "*" in the string
    string first = nn_filebase.substr(0,nn_filebase.find("*"));
    string last = nn_filebase.substr(nn_filebase.find("*") + 1);

    string tmp;
    stringstream ss;
    ss << id; ss >> tmp; // convert the id to a string

    return first + tmp + last;
  }
};

struct boxData {
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double xy, xz, yz;
  double Lx, Ly, Lz;

  void calculateBoxLengths() {
    Lx = xhigh - xlow;
    Ly = yhigh - ylow;
    Lz = zhigh - zlow;
  }

  void reset() {
    xlow = 0.0; xhigh = 0.0; ylow = 0.0; yhigh = 0.0; zlow = 0.0; zhigh = 0.0;
    xy = 0.0; xz = 0.0; yz = 0.0;
    Lx = 0.0; Ly = 0.0; Lz = 0.0;
  }
};

struct dataIndices {
  int id = -1;
  int type = -1;
  int x = -1;
  int y = -1;
  int z = -1;
  int q = -1; // charge
} data_indices;

void printInputFileHelp() {
  cout << "\n\nThe first line of the input file must contain the following items in order:\n"
       << "  1. The number of files to be processed.\n"
       << "  2. The additional rotation to apply to the system (for Zhang method) OR\n"
       << "     the misorientation of the embedded grain (for the other methods).\n"
       << "  3. The number of atom types in the simulation.\n"
       << "  4. The cutoff distance (r_cut) for generating the nearest neighbor list\n"
       << "     in terms of the lattice parameter.\n"
       << "  5. The cutoff value (fi_cut) for assigning orientation parameter values\n"
       << "     to a specific grain (for Zhang method) OR the cutoff value (< 0.5) for\n"
       << "     the Janssens method. Ignored by Ulomek method.\n"
       << "  6. The lattice parameter in Angstroms.\n"
       << "  7. The crystal structure of the lattice (fcc, bcc, sc).\n\n"
       << "These lines are then followed by the orientation matrix of the system, i.e.\n"
       << "  1 0 0 (orientation of the x axis)\n"
       << "  0 1 0 (orientation of the y axis)\n"
       << "  0 0 1 (orientation of the z axis)\n\n"
       << "The list of files (one file per line) should be placed afterwards.\n\n";
}

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

// Calculate the rounded value of x
double anInt(double x) {
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) {x += 0.5;}
  if (x < 0.0) {x -= 0.5;}
  temp = (int)(x);
  return (double)(temp);
}

double dotp(const vector <double>& v1, const vector <double>& v2) {
  if (v1.size() != 3 || v2.size() != 3) {
    cerr << "Error: vector size incorrect - should have three elements!\n";
    exit(VECTOR_SIZE_ERROR);
  }

  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

inline bool exists (const string& file) {
  struct stat buffer;
  return (stat (file.c_str(), &buffer) == 0);
}

vector <double> getAxisData(istream& fin) {
  string str;
  vector <double> axis (3, 0.0);
  getline(fin, str);
  stringstream ss(str);
  if (!(ss >> axis[0] >> axis[1] >> axis[2])) {
    cerr << "Error reading axis: " << str << "\n";
    exit(AXIS_READ_ERROR);
  }
  return axis;
}

void initializeRotatedPerfectNeighborList(const string& structure, const vector <Position>& perfect_neighbors, const string& r_axis, const double& a0) {
  string path, str;
  int theta;
  double x, y, z;
  bool rebuild_cache = false;

  if (const char* db_dir = getenv("TRAUTT_DB_DIR")) {
    path = string(db_dir) + "/" + structure + "_" + r_axis + "_rotated_1nn.db";
  } else {
    path = "./" + structure + "_" + r_axis + "_rotated_1nn.db";
  }

  if (exists(path)) {
    ifstream fin(path.c_str());
    checkFileStream(fin, path);
    int line_no = 0;

    while (getline(fin, str)) {
      ++line_no;
      int num = 0;
      vector <Position> ps;
      stringstream ss(str);
      ss >> theta;
      while (ss >> x >> y >> z) {
        ++num;
        ps.push_back(Position(x * a0, y * a0, z * a0));
      }

      switch (num) {
        case 12: if (structure.compare("fcc") != 0) {
            cerr << "Database corrupted at line " << line_no << "\n";
            rebuild_cache = true;
          }
          break;
        case 8: if (structure.compare("bcc") != 0) {
            cerr << "Database corrupted at line " << line_no << "\n";
            rebuild_cache = true;
          }
          break;
        case 6: if (structure.compare("sc") != 0) {
            cerr << "Database corrupted at line " << line_no << "\n";
            rebuild_cache = true;
          }
          break;
        default: cerr << "Database corrupted at line " << line_no << "\n";
          rebuild_cache = true;
      }
      rotated_1nn[theta] = ps;
    }
    fin.close();
  } else {rebuild_cache = true;}

  for (unsigned int i = 0; i < 90; ++i) {
    if (rebuild_cache) {break;} // no need to go through this if we already know we will rebuild
    if (rotated_1nn.find(i) == rotated_1nn.end()) {
      rebuild_cache = true;
      break;
    }
  }

  if (rebuild_cache) {
    // generate the rotation matrix, multiply it onto the 1nn positions, save and write the result
    ofstream fout(path.c_str());
    checkFileStream(fout, path);
    double costheta, sintheta, degrees = PI / 180.0, xtmp, ytmp;
    rotated_1nn.clear(); // empty the map since we're rebuilding

    // The non-rotated case
    fout << 0;
    for (const auto& pos: perfect_neighbors) {
      fout << " " <<  pos.getX() << " " << pos.getY() << " " << pos.getZ();
    }
    fout << "\n";

    for (unsigned int i = 0; i < 90; ++i) {
      vector <Position> positions;
      costheta = cos(i * degrees);
      sintheta = sin(i * degrees);

      fout << i; // angle;
      for (const auto& pos: perfect_neighbors) {
        xtmp = pos.getX() * costheta - pos.getY() * sintheta;
        ytmp = pos.getX() * sintheta + pos.getY() * costheta;
        positions.push_back(Position(xtmp * a0, ytmp * a0, pos.getZ() * a0));

        fout << " " << xtmp << " " << ytmp << " " << pos.getZ(); // positions. z position does not change when rotating about the z axis.
      }
      fout << "\n";
      rotated_1nn[i] = positions;
    }
  }

}

void parseInputFile(inputData& input, const string& input_file) {
  string str; // junk string variable

  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);

  getline(fin, str);
  stringstream ss(str);
  if (!(ss >> input.n_files >> input.theta >> input.n_types >> input.r_cut >> input.fi_cut >> input.a0 >> input.crystal_structure)) {
    cerr << "Error reading the input file. Did you forget a value?";
    printInputFileHelp();
    exit(INPUT_FORMAT_ERROR);
  }

  input.calculateTrig();

  vector <int> n(input.n_types);
  iota(n.begin(), n.end(), 1); // makes a list of numbers from 1 to n_types

  if (all_of(n.begin(), n.end(), [&](const int & i) {
      return find(input.ignored_atoms.begin(), input.ignored_atoms.end(), i) != input.ignored_atoms.end();})) {
    cout << "All atoms ignored. Exiting...\n";
    exit(EXIT_SUCCESS);
  }

  // Get the rotated coordinate system
  input.new_x_axis = getAxisData(fin);
  input.new_y_axis = getAxisData(fin);
  input.new_z_axis = getAxisData(fin);

  // store the integer values of the rotation axis into a string
  stringstream r_axis;
  r_axis << (int)(input.new_z_axis[0]) << (int)(input.new_z_axis[1]) << (int)(input.new_z_axis[2]);
  r_axis >> input.r_axis;

  cout << "Input parameters:"
       << "\n  n_files = " << input.n_files
       << "\n  theta = " << input.theta
       << "\n  n_types = " << input.n_types
       << "\n  r_cut = " << input.r_cut
       << "\n  fi_cut = " << input.fi_cut
       << "\n  a0 = " << input.a0
       << "\n  crystal structure = " << input.crystal_structure
       << "\n  Rotated coordinate system:"
       << "\n    x = " << input.new_x_axis[0] << " " << input.new_x_axis[1] << " " << input.new_x_axis[2]
       << "\n    y = " << input.new_y_axis[0] << " " << input.new_y_axis[1] << " " << input.new_y_axis[2]
       << "\n    z = " << input.new_z_axis[0] << " " << input.new_z_axis[1] << " " << input.new_z_axis[2] << "\n";

  input.new_x_axis = normalize(input.new_x_axis);
  input.new_y_axis = normalize(input.new_y_axis);
  input.new_z_axis = normalize(input.new_z_axis);

  // Calculate the transformation matrix to go from new_axis to the <100> orientation (identity matrix)
  input.calculateOldAxes();

  input.old_x_axis = normalize(input.old_x_axis);
  input.old_y_axis = normalize(input.old_y_axis);
  input.old_z_axis = normalize(input.old_z_axis);

  input.calculatePerfectNeighbors();
  input.calculateReciprocalLattices();

  // get the desired files
  int aa = 0;
  while (getline(fin, str) && aa <= input.n_files) {
    if (!str.empty()) {
      input.files.push_back(str);
      ++aa;
    }
  }

  if (show_warnings && input.files.size() != input.n_files) {
    cout << "Warning: " << input.n_files << " files specified in input file, but only " << input.files.size() << " found.\n";
  }
}

vector <Atom> getAtomData(ifstream& fin, const inputData& input, const boxData& box) {
  int n_atoms_read = 0; // number of atoms read
  string str; // holds the info
  vector <double> data; // holds the number of elements
  vector <Atom> atoms; // the set of atoms
  int atom_id, type; // id number and type of the atom
  double charge, x, y, z; // charge, and wrapped position of atom

  while (getline(fin, str)) {
    stringstream ss(str);
    double dummy;
    while (ss >> dummy) {data.push_back(dummy);} // counts the number of elements

    atom_id = (int)(data[data_indices.id]);
    type = (int)(data[data_indices.type]);
    x = data[data_indices.x];
    y = data[data_indices.y];
    z = data[data_indices.z];
    if (has_charge) {charge = data[data_indices.q];}
    else {charge = 0.0;}

    data.clear();

    if (type > input.n_types) {
      cerr << "Error: unexpected atom type.\n"
           << "n_types = " << input.n_types << " < this atom's type = " << type << "\n";
      exit(ATOM_TYPE_ERROR);
    }

    // We adjust the positions so that we start at the origin so we can easily
    // assign to cells
    x = x - box.xlow;
    y = y - box.ylow;
    z = z - box.zlow;

    Position p(x,y,z);
    if (atom_id > atoms.size()) {atoms.resize(atom_id, Atom());}
    atoms[atom_id - 1] = Atom(atom_id, type, charge, p);
    ++n_atoms_read;
  }
  fin.close();

  if (n_atoms_read != atoms.size()) {
    cerr << "Error: number of atoms read does not match the number of atoms in the simulation.\n"
         << "N = " << atoms.size() << " != n_atoms_read = " << atoms.size() << "\n";
    exit(ATOM_COUNT_ERROR);
  }
  return atoms;
}

vector <vector <int> > generateCellLinkedList(const vector <Atom>& atoms,
                                              const inputData& input,
                                              const boxData& box) {
  int ncellx, ncelly, ncellz; // number of cells in each direction
  int idx, idy, idz; // cell number in each direction
  double lcellx, lcelly, lcellz; // length of cells in each direction
  int n_atoms_per_cell; // number of atoms allowed per cell
  double drij_sq, rxij, ryij, rzij; // square of distance, x, y, and z separation.
  vector <vector <int> > iatom; // cell-linked list
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell

  double rcut = input.a0 * input.r_cut;
  double rcut_sq = rcut * rcut;

  // First we generate the number of cells in each direction
  ncellx = (int)(box.Lx / rcut) + 1;
  ncelly = (int)(box.Ly / rcut) + 1;
  ncellz = (int)(box.Lz / rcut) + 1;

  // Length of cells in each direction
  lcellx = box.Lx / ncellx;
  lcelly = box.Ly / ncelly;
  lcellz = box.Lz / ncellz;

  // Minimum number of atoms allowed of 100
  n_atoms_per_cell = max((int)(atoms.size() / (double)(ncellx * ncelly * ncellz)), 100);

  // resize the vectors
  icell.resize(ncellx, vector <vector <int> > // x dimension
              (ncelly, vector <int> // y dimension
              (ncellz, 0))); // z dimension
  pcell.resize(ncellx, vector <vector <vector <int> > > // x dimension
              (ncelly, vector <vector <int> > // y dimension
              (ncellz, vector <int> // z dimension
              (n_atoms_per_cell, 0)))); // atom number in cell.
  iatom.resize(n_atoms_per_cell, vector <int> (atoms.size(),0));

  // generate the pcell and icell matrices.
  for (size_t i = 0; i < atoms.size(); ++i) {
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(),
        atoms[i].getType()) != input.ignored_atoms.end()) {continue;}

    // Assign this atom to a cell
    // Rounds towards 0 with a type cast
    idx = (int)(atoms[i].getWrapped()[0] / lcellx); // assign the x cell
    idy = (int)(atoms[i].getWrapped()[1] / lcelly); // assign the y cell
    idz = (int)(atoms[i].getWrapped()[2] / lcellz); // assign the z cell
    // Check if we went out of bounds
    // C++ indexes from 0, so we have to subtract 1 from the maximum value to
    // stay within our memory bounds
    if (idx >= ncellx) idx = ncellx - 1;
    if (idy >= ncelly) idy = ncelly - 1;
    if (idz >= ncellz) idz = ncellz - 1;

    ++icell[idx][idy][idz]; // increase the number of atoms in this cell
    // assign the atom number to this index.
    pcell[idx][idy][idz][icell[idx][idy][idz] - 1] = i;
  }

  for (int i = 0; i < ncellx; ++i) { // For each x cell
    for (int j = 0; j < ncelly; ++j) { // For each y cell
      for (int k = 0; k < ncellz; ++k) { // For each z cell
        for (int l = 0; l < icell[i][j][k]; ++l) { // For each atom in this cell
          int id = pcell[i][j][k][l]; // store this atom id
          // Now we check each sub cell around the current one
          for (int ii = -1; ii < 2; ++ii) { // allowed values: -1, 0, and 1
            for (int jj = -1; jj < 2; ++jj) {
              for (int kk = -1; kk < 2; ++kk) {
                int ia = i + ii; // min value: -1.  Max value: number of cells in dimension
                int ja = j + jj;
                int ka = k + kk;
                // Check to make sure we are still in bounds
                // C++ indexes from 0, so we accomodate.
                if (ia >= ncellx) ia = 0;
                if (ja >= ncelly) ja = 0;
                if (ka >= ncellz) ka = 0;
                if (ia < 0) ia = ncellx - 1;
                if (ja < 0) ja = ncelly - 1;
                if (ka < 0) ka = ncellz - 1;

                // Now check each atom in this cell
                for (int m = 0; m < icell[ia][ja][ka]; ++m) {
                  int jd = pcell[ia][ja][ka][m];
                  // If jd <= id, we've already dealt with this interaction
                  if (jd <= id) {continue;}

                  // Now the actual calculations!
                  rxij = atoms[id].getWrapped()[0] - atoms[jd].getWrapped()[0];
                  ryij = atoms[id].getWrapped()[1] - atoms[jd].getWrapped()[1];
                  rzij = atoms[id].getWrapped()[2] - atoms[jd].getWrapped()[2];

                  // Apply PBCs
                  rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
                  ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
                  rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

                  // Now calculate the distance
                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                  // move to the next atom if we're too far away
                  if (drij_sq > rcut_sq) {continue;}

                  // Create the neighbor list
                  iatom[0][id] += 1; //for atom id - number of neighbors
                  if (iatom[0][id] >= n_atoms_per_cell) {
                    n_atoms_per_cell += 100;
                    iatom.resize(n_atoms_per_cell, vector <int> (atoms.size(),0));
                  }
                  iatom[(iatom[0][id])][id] = jd; // point to the next atom
                  iatom[0][jd] += 1; // for atom jd
                  if (iatom[0][jd] >= n_atoms_per_cell) {
                    n_atoms_per_cell += 100;
                    iatom.resize(n_atoms_per_cell, vector <int> (atoms.size(),0));
                  }
                  iatom[(iatom[0][jd])][jd] = id;
                } // m
              } //kk
            } //jj
          } //ii
        } // l
      } // k
    } // j
  } // i

  return iatom;
}

vector <double> zhangSymmetryParameter(vector <Atom>& atoms,
                                       const vector <vector <int> >& iatom,
                                       const vector <bool> allowed_atoms,
                                       const inputData& input,
                                       const boxData& box) {
  vector <double> symm (atoms.size(), 0.0);
  // positions, distances, changed positions, square of distance, cos^2 of angle, orientation parameter
  double x, y, z, rxij, ryij, rzij, xtemp, ytemp, drij_sq, costheta_sq, val;
  double coeffs [2] = {3.0, 2.0}; // coefficients of the orientation parameter calculation

  for (size_t i = 0; i < atoms.size(); ++i) {
    // Check if we are ignoring this atom type
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}

    // Check if we are calculating the symmetry parameter for this particular atom
    if (!allowed_atoms[i]) {continue;}

    x = atoms[i].getWrapped()[0];
    y = atoms[i].getWrapped()[1];
    z = atoms[i].getWrapped()[2];

    // We start at l = 1 because l = 0 indicates the number of neighbors for atom i
    for (int l = 1; l <= iatom[0][i]; ++l) {
      unsigned int id = iatom[l][i];

      // Position difference vector
      rxij = atoms[id].getWrapped()[0] - x;
      ryij = atoms[id].getWrapped()[1] - y;
      rzij = atoms[id].getWrapped()[2] - z;

      // Apply PBCs. Note that applying PBCs with the positions projected in the
      // <100> reference frame messes up the calculations!
      rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
      ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
      rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

      // Additional rotation
      xtemp = rxij * input.cosin - ryij * input.sinin;
      ytemp = ryij * input.cosin + rxij * input.sinin;

      rxij = xtemp;
      ryij = ytemp;

      // Rotate the vector into the <100> reference frame
      // NOTE: in order for this method to work well, the correct cutoff distance needs to be used!
      xtemp = rxij * input.new_x_axis[0] + ryij * input.new_y_axis[0] + rzij * input.new_z_axis[0];
      ytemp = rxij * input.new_x_axis[1] + ryij * input.new_y_axis[1] + rzij * input.new_z_axis[1];

      // Project onto the xy plane
      drij_sq = xtemp * xtemp + ytemp * ytemp;

      // Handles the case there the projected position of the atom is right on top of the current atom
      if (drij_sq < 1.0E-8) {
        symm[i] += 1.0;
        continue;
      }
      // cos = dot(A,B) / (|A| * |B|)
      costheta_sq = (xtemp * xtemp) / drij_sq;
      val = (coeffs[0] - coeffs[1] * costheta_sq) * (coeffs[0] - coeffs[1] * costheta_sq) * costheta_sq;
      symm[i] += val;
    }
  }

  for (size_t i = 0; i < symm.size(); ++i) {
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}
    if (!(allowed_atoms[i])) {continue;}
    // On the off chance that no neighbors have been assigned to this atom,
    // we estimate the symmetry parameter to be the previous atom's value.
    // This assumes that atoms with ID's close to each other are close to
    // each other in the simulation, which may not always be the case.
    else if (iatom[0][i] == 0 && i != 0) {
      if (show_warnings) {
        cout << "\nWarning: no neighbors detected for atom " << atoms[i].getId()
        << ". Using the symmetry parameter of the previous atom = " << symm[i - 1] << "\n";
      }
      symm[i] = symm[i - 1];
    }
    else {symm[i] /= iatom[0][i];}

    if (symm[i] <= input.fi_cut) {atoms[i].setMark(1);}
    else {atoms[i].setMark(2);}
  }

  return symm;
}

vector <double> janssensSymmetryParameter(vector <Atom>& atoms,
                                          const vector <vector <int> >& iatom,
                                          const vector <bool> allowed_atoms,
                                          const inputData& input,
                                          const boxData& box) {
  vector <double> symm (atoms.size(), 0.0);
  vector <Position> neighbor_positions;
  double x, y, z, rxij, ryij, rzij, drij_sq, rcut, val;
  vector <double> distances;

  if (input.fi_cut <= 0 || input.fi_cut >= 1.0) {
    cerr << "Error: the Janssens method requires an orientation parameter cutoff in the range 0 < fi_cut < 1\n";
    exit(BOUNDS_ERROR);
  }

  rcut = input.a0 * input.r_cut;

  for (size_t i = 0; i < atoms.size(); ++i) {
    // move to the next atom if ignoring this atom type
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}

    // move to the next atom if not calculating for this atom
    if (!allowed_atoms[i]) {continue;}

    // Otherwise, we need to identify the positions of the nearest neighbors
    x = atoms[i].getWrapped()[0];
    y = atoms[i].getWrapped()[1];
    z = atoms[i].getWrapped()[2];
    neighbor_positions.clear();

    for (int l = 1; l <= iatom[0][i]; ++l) {
      unsigned int id = iatom[l][i];
      // calculate the distance between the central atom and its neighbors
      rxij = atoms[id].getWrapped()[0] - x;
      ryij = atoms[id].getWrapped()[1] - y;
      rzij = atoms[id].getWrapped()[2] - z;

      // Apply PBCs
      rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
      ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
      rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

      drij_sq = rxij * rxij + ryij * ryij + rzij * rzij;

      if (drij_sq < rcut * rcut) {neighbor_positions.push_back(Position(rxij, ryij, rzij));}
    }

    // Now that we've found the neighbor positions that are within the cutoff,
    // we need to calculate the orientation parameter
    // Only do this for the N < first_nn_list[input.crystal_structure].size() atoms!
    size_t num_atoms = (neighbor_positions.size() < first_nn_list.at(input.crystal_structure).size()) ? neighbor_positions.size() : first_nn_list.at(input.crystal_structure).size();
    for (size_t j = 0; j < num_atoms; ++j) {
      x = neighbor_positions[j][0];
      y = neighbor_positions[j][1];
      z = neighbor_positions[j][2];
      distances.clear();
      for (size_t k = 0; k < input.perfect_matrix_neighbors.size(); ++k) {
        rxij = x - input.perfect_matrix_neighbors[k].getX();
        ryij = y - input.perfect_matrix_neighbors[k].getY();
        rzij = z - input.perfect_matrix_neighbors[k].getZ();
        distances.push_back(rxij * rxij + ryij * ryij + rzij * rzij);
      }

      // Add the minimum distance value
      val = sqrt(*min_element(distances.begin(), distances.end()));
      symm[i] += val;
    }
    symm[i] /= num_atoms;
  }

  for (size_t i = 0; i < symm.size(); ++i) {
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}
    if (!(allowed_atoms[i])) {continue;}
    else if (iatom[0][i] == 0 && i != 0) {
      if (show_warnings) {
        cout << "\nWarning: no neighbors detected for atom " << atoms[i].getId()
             << ". Using the symmetry parameter of the previous atom = " << symm[i - 1] << "\n";
      }
      symm[i] = symm[i - 1];
    }

    if (symm[i] > input.janssens_symm * input.fi_cut) {atoms[i].setMark(1);}
    else {atoms[i].setMark(2);}
  }
  return symm;
}

vector <double> ulomekSymmetryParameter(vector <Atom>& atoms,
                                        const vector <vector <int> >& iatom,
                                        const vector <bool> allowed_atoms,
                                        const inputData& input,
                                        const boxData& box) {
  vector <double> symm (atoms.size(), 0.0);
  vector <vector <double> > Q; // combined reciprocal lattice vectors
  vector <double> kappa = {1.0, 1.0, 1.0, -1.0, -1.0, -1.0};
  double x, y, z, dijx, dijy, dijz, dikx, diky, dikz, rsqij, rsqik, wij, wik, scalarprod;
  unsigned int id;

  double rcut = input.a0 * input.r_cut;
  double cutsq = rcut * rcut;
  // Store the reciprocal lattices
  for (size_t i = 0; i < input.reciprocal_1.size(); ++i) {
    Q.push_back(input.reciprocal_1[i]);
  }
  for (size_t i = 0; i < input.reciprocal_2.size(); ++i) {
    Q.push_back(input.reciprocal_2[i]);
  }

  for (size_t i = 0; i < atoms.size(); ++i) {
    // move to the next atom if ignoring this atom type
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}

    // move to the next atom if not calculating for this atom
    if (!allowed_atoms[i]) {continue;}

    // Store the position of this atom
    x = atoms[i].getWrapped()[0];
    y = atoms[i].getWrapped()[1];
    z = atoms[i].getWrapped()[2];

    for (int l = 1; l <= iatom[0][i]; ++l) {
      id = iatom[l][i];
      // calculate the distance between these atoms
      dijx = x - atoms[id].getWrapped()[0];
      dijy = y - atoms[id].getWrapped()[1];
      dijz = z - atoms[id].getWrapped()[2];

      // Apply PBCs
      dijx = dijx - anInt(dijx / box.Lx) * box.Lx;
      dijy = dijy - anInt(dijy / box.Ly) * box.Ly;
      dijz = dijz - anInt(dijz / box.Lz) * box.Lz;

      rsqij = (dijx * dijx + dijy * dijy + dijz * dijz) / cutsq;
      wij = rsqij * (rsqij - 2) + 1;

      for (int ll = 1; ll <= iatom[0][i]; ++ll) {
        id = iatom[ll][i];

        dikx = x - atoms[id].getWrapped()[0];
        diky = y - atoms[id].getWrapped()[1];
        dikz = z - atoms[id].getWrapped()[2];

        // Apply PBCs
        dikx = dikx - anInt(dikx / box.Lx) * box.Lx;
        diky = diky - anInt(diky / box.Ly) * box.Ly;
        dikz = dikz - anInt(dikz / box.Lz) * box.Lz;

        rsqik = (dikx * dikx + diky * diky + dikz * dikz) / cutsq;
        wik = rsqik * (rsqik - 2) + 1;
        for (int alpha = 0; alpha < 6; ++alpha) {
          scalarprod = Q[alpha][0] * (dikx - dijx) + Q[alpha][1] * (diky - dijy) + Q[alpha][2] * (dikz - dijz);
          symm[i] += kappa[alpha] * wij * wik * cos(scalarprod);
        } // alpha
      } // atom ll
    } // atom l
    symm[i] /= input.scalenorm;
  }

  // Note that an adjustment parameter can assign atoms to a grain boundary, rather than to a grain
  for (size_t i = 0; i < symm.size(); ++i) {
    if (symm[i] > input.fi_cut) {atoms[i].setMark(1);}
    else {atoms[i].setMark(2);}
  }
  return symm;
}

vector <double> trauttSymmetryParameter(vector <Atom>& atoms,
                                        const vector <vector <int> >& iatom,
                                        const vector <bool> allowed_atoms,
                                        const inputData& input,
                                        const boxData& box) {
  vector <double> symm (atoms.size(), 0.0);
  Position pref, pdiff, numerator; // position of reference atom, differences in position, top part of exponential in eq. 33
  unsigned int id;
  vector <double> symm_test(rotated_1nn.size(), 0.0);
  double a0_sq = input.a0 * input.a0;
  map <int, int> histogram; // key is the misorientation, value is the count
  VerifyNewFile verify("trautt_histogram.txt");
  // see Trautt & Mishin, Acta Mat 60 (2012) 2407-2424 equation (33)

  for (size_t i = 0; i < atoms.size(); ++i) {
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}
    if (!allowed_atoms[i]) {continue;}

    fill(symm_test.begin(), symm_test.end(), 0.0); // reset the values of symm_test
    pref = atoms[i].getWrapped();

    for (int l = 1; l <= iatom[0][i]; ++l) { // for each neighbor atom
      id = iatom[l][i];
      pdiff = pref - atoms[id].getWrapped(); // this is r_{ijk} in equation (33)
      pdiff[0] = pdiff[0] - anInt(pdiff[0] / box.Lx) * box.Lx; // PBCs
      pdiff[1] = pdiff[1] - anInt(pdiff[1] / box.Ly) * box.Ly;
      pdiff[2] = pdiff[2] - anInt(pdiff[2] / box.Lz) * box.Lz;

      for (size_t j = 0; j < rotated_1nn.size(); ++j) { // for each angle
        for (size_t k = 0; k < rotated_1nn[j].size(); ++k) { // for each perfect neighbor position
          numerator = pdiff - rotated_1nn[j][k];
          for (size_t m = 0; m < 3; ++m) {
            symm_test[j] += exp(-(numerator[m] * numerator[m]) / a0_sq);
          }
        }
      }
    }
    auto max_it = max_element(symm_test.begin(), symm_test.end());
    // cout << *max_it << " at index " << max_it - symm_test.begin() << "\n";
    ++histogram[max_it - symm_test.begin()];
    atoms[i].setMark(max_it - symm_test.begin());
    symm[i] = *max_it;
  }

  ofstream fout(verify.validNewFile().c_str());
  checkFileStream(fout, verify.validNewFile());

  for (auto& item: histogram) {
    if (item.second == 0) {continue;}
    fout << item.first << " " << item.second << "\n";
  }

  return symm;
}

void writeAtomsToFile(const string& filename, const vector <Atom>& atoms,
                      const vector <bool> allowed_atoms, const vector <double>& symm,
                      const inputData& input, const boxData& box) {
  int n_atoms_written = 0;

  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);

  // This writes things in a Tecplot-readable file format
  fout << "VARIABLES = \"Atom ID\", \"Atom Type\", ";
  if (has_charge) {fout << "\"Atom Charge\", ";}
  fout << "\"X\", \"Y\", \"Z\", \"Grain Number\", \"Orientation Parameter\"\n";

  for (size_t i = 0; i < atoms.size(); ++i) {
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}

    if (!(allowed_atoms[i])) {continue;}

    fout << atoms[i].getId() << " "
         << atoms[i].getType() << " ";
    if (has_charge) {fout << atoms[i].getCharge() << " ";}
    fout << (atoms[i].getWrapped()[0] + box.xlow) << " "
         << (atoms[i].getWrapped()[1] + box.ylow) << " "
         << (atoms[i].getWrapped()[2] + box.zlow)<< " "
         << atoms[i].getMark() << " "
         << symm[i] << "\n";
    ++n_atoms_written;
  }
  fout.close();

  int num_allowed_atoms = accumulate(allowed_atoms.begin(), allowed_atoms.end(), 0);
  if (num_allowed_atoms != n_atoms_written) {
    cerr << "Error: number of atoms written does not match number of atoms specified to be written\n"
         << "num_allowed_atoms = " << num_allowed_atoms << " != n_atoms_written = " << n_atoms_written << "\n";
    exit(ATOM_COUNT_ERROR);
  }
}

void processData(const inputData& input) {
  boxData box;
  int n_grain_1, n_grain_2; // number of atoms assigned to each grain
  int N; // number of atoms
  string str, filename2; // junk string variable, filename of the data with the grain assignment

  vector <Atom> atoms; // the atoms from the data file
  vector <bool> allowed_atoms; // the atoms that will be printed
  vector <vector <int> > iatom; // cell-linked list
  vector <double> symm; // the orientation factors

  VerifyNewFile verify1(input.outfile);
  ofstream fout_data(verify1.validNewFile().c_str());
  checkFileStream(fout_data, verify1.validNewFile());

  fout_data << "# Data consists of: [timestep/file, atoms in grain 1, atoms in grain 2]\n";

  // VerifyNewFile verify2("misorientation_data.txt");
  // ofstream fout_misorientation(verify2.validNewFile().c_str());
  // checkFileStream(fout_misorientation, verify2.validNewFile());

  int aa = 1;
  int x_index = -1, y_index = -1, z_index = -1; // indices for the x, y, and z positions of the atoms
  int id_index = -1, type_index = -1, charge_index = -1; // indices for the atom id, type, and charge
  for (vector <string>::const_iterator it = input.files.begin(); it != input.files.end(); ++it) {
    box.reset();
    n_grain_1 = 0;
    n_grain_2 = 0;
    ifstream fin((*it).c_str());
    checkFileStream(fin, *it);

    if ((*it).find("dump") != string::npos) {
      int timestep;
      getline(fin, str); // Gets ITEM: TIMESTEP
      fin >> timestep; // Gets the timestep number
      fin.ignore();
      if (aa == 1) {
        if (show_warnings && timestep != 0) {
          cout << "Warning: first data file is not at timestep 0\n";
        }
      }

      getline(fin, str); // Gets ITEM: NUMBER OF ATOMS
      fin >> N; // number of atoms
      fin.ignore();

      getline(fin, str); // Gets ITEM: BOX BOUNDS
      if (str.find("xy") == string::npos) {
        fin >> box.xlow >> box.xhigh
            >> box.ylow >> box.yhigh
            >> box.zlow >> box.zhigh;
      }
      else { // We also get the triclinic tilt factors for triclinic systems
        fin >> box.xlow >> box.xhigh >> box.xy
            >> box.ylow >> box.yhigh >> box.xz
            >> box.zlow >> box.zhigh >> box.yz;
      }
      fin.ignore();

      getline(fin, str); // Gets ITEM: ATOMS <data types>
      string data_types = str.substr(str.find("ATOMS") + 6);
      stringstream ss(data_types);
      string tmp;
      int idx = 0, xu_index = -1, yu_index = -1, zu_index = -1;
      while (ss >> tmp) {
        if (tmp.compare("id") == 0) {id_index = idx;}
        else if (tmp.compare("type") == 0) {type_index = idx;}
        else if (tmp.compare("q") == 0) {charge_index = idx;}
        else if (tmp.compare("x") == 0) {x_index = idx;}
        else if (tmp.compare("y") == 0) {y_index = idx;}
        else if (tmp.compare("z") == 0) {z_index = idx;}
        else if (tmp.compare("xu") == 0) {xu_index = idx;}
        else if (tmp.compare("yu") == 0) {yu_index = idx;}
        else if (tmp.compare("zu") == 0) {zu_index = idx;}
        ++idx;
      }
      // Use unwrapped positions if regular positions are not available
      if (x_index == -1) {
        if (xu_index == -1) {cout << "Error finding x position index"; exit(INPUT_FORMAT_ERROR);}
        else {x_index = xu_index;}
      }
      if (y_index == -1) {
        if (yu_index == -1) {cout << "Error finding y position index"; exit(INPUT_FORMAT_ERROR);}
        else {y_index = yu_index;}
      }
      if (z_index == -1) {
        if (zu_index == -1) {cout << "Error finding z position index"; exit(INPUT_FORMAT_ERROR);}
        else {y_index = yu_index;}
      }
      filename2 = (*it).substr(0,(*it).find(".dump")) + "_interface.dat";

      fout_data << timestep << " ";
    }
    else if ((*it).find(".dat") != string::npos) {
      getline(fin, str); // gets comment line TODO: This is where I need to parse what data values are included
      string data_types = str.substr(str.find("[") + 1, str.find("]") - str.find("[") - 1);
      stringstream ss(data_types);
      string tmp;
      int idx = 0;
      while (ss >> tmp) {
        if (tmp.compare("ID") == 0 || tmp.compare("id") == 0) {id_index = idx;}
        else if (tmp.compare("type") == 0) {type_index = idx;}
        else if (tmp.compare("charge") == 0 || tmp.compare("q") == 0) {charge_index = idx;}
        else if (tmp.compare("x") == 0) {x_index = idx;}
        else if (tmp.compare("y") == 0) {y_index = idx;}
        else if (tmp.compare("z") == 0) {z_index = idx;}
        ++idx;
      }
      if (x_index == -1) {cout << "Error finding x position index"; exit(INPUT_FORMAT_ERROR);}
      if (y_index == -1) {cout << "Error finding y position index"; exit(INPUT_FORMAT_ERROR);}
      if (z_index == -1) {cout << "Error finding z position index"; exit(INPUT_FORMAT_ERROR);}
      fin >> N >> str; // number of atoms
      fin >> str >> str >> str; // atom types n_types - specified in the input file
      fin >> box.xlow >> box.xhigh >> str >> str
          >> box.ylow >> box.yhigh >> str >> str
          >> box.zlow >> box.zhigh >> str >> str;
      fin.ignore();

      getline(fin, str);
      if (str.find("xy") != string::npos) { // Triclinic tilt factors
        fin >> box.xy >> box.xz >> box.yz >> str >> str >> str;
        fin.ignore();
      }

      getline(fin, str);
      getline(fin, str);
      filename2 = (*it).substr(0,(*it).find(".dat")) + "_interface.dat";

      fout_data << (*it).substr(0, (*it).find(".dat")) << " ";
    }
    else {
      cerr << "Error: Unknown file type\n";
      exit(FILE_FORMAT_ERROR);
    }
    data_indices.id = id_index;
    data_indices.type = type_index;
    data_indices.q = charge_index;
    data_indices.x = x_index;
    data_indices.y = y_index;
    data_indices.z = z_index;
    if (data_indices.q != -1) {has_charge = true;}
    box.calculateBoxLengths();

    allowed_atoms.resize(N,false);
    atoms = getAtomData(fin, input, box);

    iatom = generateCellLinkedList(atoms, input, box);

    if (input.atom_ids_list.size() > 0) {
      allowed_atoms.assign(N, false);
      for (size_t i = 0; i < input.atom_ids_list.size(); ++i) {
        allowed_atoms[input.atom_ids_list[i] - 1] = true;
      }
    }
    else {
      allowed_atoms.assign(N, true);
      if (input.ignored_atoms.size() > 0) {
        for (size_t i = 0; i < atoms.size(); ++i)
        {
          if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {
            allowed_atoms[i] = false;
          }
        }
      }
    }

    if (print_nearest_neighbors) {
      string neighbor_filename = input.generateNNFilename(aa);
      ofstream fout_neighbors(neighbor_filename.c_str());
      checkFileStream(fout_neighbors, neighbor_filename);

      for (size_t i = 0; i < atoms.size(); ++i) {
        if (!allowed_atoms[i]) {continue;} // only print the specified atoms
        else {
          vector <int> neighs; // neighbors of atom i
          fout_neighbors << "Atom " << i << " has " << iatom[0][i] << "neighbors:\n";
          for (int l = 1; l <= iatom[0][i]; ++l) {neighs.push_back(iatom[l][i]);}
          sort (neighs.begin(), neighs.end());
          for (size_t j = 0; j < neighs.size(); ++j) {fout_neighbors << "  Atom " << neighs[j] << "\n";}
        }
      }
    }

    // Now that we have the atoms safely stored, we can process them
    if (input.algorithm == 'j') {symm = janssensSymmetryParameter(atoms, iatom, allowed_atoms, input, box);}
    else if (input.algorithm == 'z') {symm = zhangSymmetryParameter(atoms, iatom, allowed_atoms, input, box);}
    else if (input.algorithm == 'u') {symm = ulomekSymmetryParameter(atoms, iatom, allowed_atoms, input, box);}
    else if (input.algorithm == 't') {
      initializeRotatedPerfectNeighborList(input.crystal_structure, input.perfect_matrix_neighbors, input.r_axis, input.a0);
      symm = trauttSymmetryParameter(atoms, iatom, allowed_atoms, input, box);
    }
    else {
      cerr << "Unknown orientation parameter algorithm. Exiting...\n";
      exit(ERROR_CODE_NOT_DEFINED);
    }

    for (size_t i = 0; i < atoms.size(); ++i) {
      if (input.algorithm == 't') {break;}
      if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}
      if (!allowed_atoms[i]) {continue;}
      if (atoms[i].getMark() == 1) {++n_grain_1;}
      else if (atoms[i].getMark() == 2) {++n_grain_2;}
      else {
        cerr << "Error: Unrecognized grain assignment.\n";
        exit(BOUNDS_ERROR);
      }
    }

    fout_data << n_grain_1 << " " << n_grain_2 << "\n";

    if (input.print_files_every != 0 && ((aa - 1) % input.print_files_every == 0 || aa == input.files.size())) {
      writeAtomsToFile(filename2, atoms, allowed_atoms, symm, input, box);
    }

    fin.close();

    cout << "\r";
    cout << "Processing of file \"" << (*it) << "\" completed.";
    cout.flush();
    ++aa;
  }
  fout_data.close();
}

vector <int> getAtomIdsList(const string& file) {
  string str;
  vector <int> id_list;
  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  while (getline(fin, str)) {
    stringstream ss(str);
    int tmp;
    if (!(ss >> tmp)) {
      cerr << "Error reading atom IDs.\n";
      exit(FILE_FORMAT_ERROR);
    }
    else {id_list.push_back(tmp);}
  }

  sort(id_list.begin(), id_list.end());
  return id_list;
}

void printAlgorithmCitations() {
  cout << "Trautt algorithm:\n"
       << "  Trautt Z.T. and Mishin, Y. Acta Materialia 60(5) 2407-2424 (2012)\n"
       << "Zhang algorithm:\n"
       << "  Zhang H. et al., Acta Materialia 53(1) 79-86 (2005)\n"
       << "Janssens algorithm:\n"
       << "  Janssens K. et al., Nature Materials 5(2) 124-127 (2006)\n"
       << "Ulomek algorithm:\n"
       << "  Ulomek F. et al., Modelling and Simulation in Materials Science and Engineering 23(2) 025007 (2015)\n\n";
}

int main(int argc, char** argv)
{
  inputData input;
  string input_file;

  try {
    cxxopts::Options options(argv[0], "Determines which grain atoms belong to.");
    options
      .positional_help("File")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Input file", cxxopts::value<string>(input_file), "file")
        ("o,output", "Output file. Uniqueness is assured.", cxxopts::value<string>(input.outfile)->default_value("data.txt"), "file")
        ("e,every", "Option to specify how frequently interface files should be written. 1 specifies every file, 0 specifies none, any other integer specifies than every n files should be written. First and last interface files are always written if n != 0.", cxxopts::value<unsigned int>(input.print_files_every)->default_value("1"), "n")
        ("i,ignore", "Ignore atoms of a specific type (number) during neighbor list creation - each type must be specified separately", cxxopts::value<vector <int> >(), "atom_type")
        ("g,algorithm", "The algorithm for assigning atoms to grains. u - Ulomek, z - Zhang, j - Janssens, t - Trautt", cxxopts::value<char>(input.algorithm)->default_value("u"))
        ("c,citations", "Show the citations for the different algorithms and exit", cxxopts::value<bool>()->default_value("false"))
        ("q,quiet", "Suppress warnings", cxxopts::value<bool>(show_warnings)->implicit_value("false"))
        ("atom-ids", "Only output the results for the specified IDs from the file (containing a list of atom ids)", cxxopts::value<string>(), "file")
        ("print-nearest-neighbors", "Print the nearest neighbor list to a file", cxxopts::value<string>(input.nn_filebase)->default_value("nearest_neighbors_*.txt"), "file")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("citations")) {
      printAlgorithmCitations();
      return EXIT_SUCCESS;
    }

    if (result.count("help") || result.count("file") == 0) {
      cout << options.help();
      printInputFileHelp();
      return EXIT_SUCCESS;
    }

    if (result.count("ignore")) {
      input.ignored_atoms = result["ignore"].as<vector <int> >();
      cout << "Ignoring atoms of type: ";
      for (vector <int>::iterator it = input.ignored_atoms.begin(); it != input.ignored_atoms.end();) {
        cout << *it;
        if (++it != input.ignored_atoms.end()) {cout << ", ";}
      }
      cout << "\n";
    }
    else {cout << "Processing all atom types.\n";}

    if (result.count("print-nearest-neighbors")) {print_nearest_neighbors = true;}

    if (result.count("atom-ids")) {input.atom_ids_list = getAtomIdsList(result["atom-ids"].as<string>());}

    if (result.count("file")) {
      parseInputFile(input, input_file);
      processData(input);
      cout << "\n";
    }
  }
  catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << "\n";
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath> // for sin, cos
#include <numeric> // for iota
#include <omp.h>

#include <cxxopts.hpp>

#include "atom.h"
#include "position.h"
#include "error_code_defines.h"
#include "verifyNewFile.h"

using namespace std;

#define PI 3.141592653589793
#define THREAD_NUM 4

// TODO: Fix race condition(s?) for parallelized loops

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

double anInt(double x) {
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) {x += 0.5;}
  if (x < 0.0) {x -= 0.5;}
  temp = (int)(x);
  return (double)(temp);
}


vector <vector <double> > dotp(const vector <vector <double> > &m1, const vector <vector <double> > &m2) {
  if (m1.size() == 0 || m2.size() == 0) {
    cerr << "Cannot multiply non-existent matrices\n";
    exit(MATH_ERROR);
  }

  if (m1[0].size() != m2.size()) {
    cerr << "Inner sizes much match: m1[0].size (" << m1[0].size() << ") != m2.size (" << m2.size() << ")\n";
    exit(MATH_ERROR);
  }

  vector <vector <double> > result (m1.size(), vector <double> (m2[0].size(), 0.0));

  for (size_t row = 0; row < m1.size(); ++row) {
    for (size_t col = 0; col < m2[0].size(); ++col) {
      double sum = 0;
      for (size_t n = 0; n < m1[0].size(); ++n) sum += m1[row][n] * m2[n][col];
      result[row][col] = sum;
    }
  }

  return result;
}


vector <vector <double> > dotp (const vector <vector <double> > &m1, const vector <double> &m2) {
  vector <vector <double> > col_vec (1, vector <double> (m2.size(), 0.0));
  for (size_t i = 0; i < m2.size(); ++i) {
    col_vec[0][i] = m2[i];
  }

  return dotp(m1,col_vec);
}


double dotp (const vector <double> &v1, const vector <double> &v2) {
  if (v1.size() != v2.size()) {
    cerr << "Vectors must be the same size! v1.size (" << v1.size() << " != v2.size(" << v2.size() << ")\n";
    exit(MATH_ERROR);
  }
  double result = 0.0;

  for (size_t i = 0; i < v1.size(); ++i) {
    result += v1[i] * v2[i];
  }

  return result;
}


vector <vector <double> > scale (const vector <vector <double> > &m, const double &scalar) {
  vector <vector <double> > result (m.size(), vector <double> (m[0].size()));

  for (size_t i = 0; i < m.size(); ++i) {
    for (size_t j = 0; j < m[0].size(); ++j) {
      result[i][j] = scalar * m[i][j];
    }
  }

  return result;
}


vector <vector <double> > transpose (const vector <vector <double> > &m) {
  vector <vector <double> > result (m[0].size(), vector <double> (m.size(), 0.0));

  for (size_t i = 0; i < m.size(); ++i) {
    for (size_t j = 0; j < m[0].size(); ++j) {
      result[j][i] = m[i][j];
    }
  }

  return result;
}


template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}


struct inputData {
  // directly from the input file
  int n_types;
  double theta, r_cut, a0, eta;
  string crystal_structure;
  vector <double> new_x_axis {0.0,0.0,0.0}, new_y_axis{0.0,0.0,0.0}, new_z_axis{0.0,0.0,0.0};
  vector <double> old_x_axis {0.0,0.0,0.0}, old_y_axis{0.0,0.0,0.0}, old_z_axis{0.0,0.0,0.0};
  vector <string> files; // the data files to examine
  string r_axis; // the string representation of the orientation matrix, given as zzz

  // from the command line
  vector <int> ignored_atoms, atom_ids_list;
  string outfile, nn_filebase;
  unsigned int print_files_every = 0;

  // calculated from input file
  double sinin, cosin; // sin and cos of the angle
  vector <vector <double> > dir_vec = vector <vector <double> > (6, vector <double>(3, 0.0));
  vector <vector <vector <double> > > reciprocal_vectors = vector <vector <vector <double> > > (2, vector <vector <double> > (3, vector <double> (3, 0.0)));
  double norm_fac = 0;

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

  // for Ulomek orientation parameter
  void calculateReciprocalLattices() {
    // hard coded numbers because these sizes do not change - they form the 3x3
    // orientation matrix of the two grains
    vector <vector <double> > rot_mat = {
      {cosin, -sinin, 0.0},
      {sinin, cosin, 0.0},
      {0.0, 0.0, 1.0}
    };

    vector <vector <double> > orient_mat (6,vector <double> (3, 0.0));
    orient_mat[0] = normalize(new_x_axis);
    orient_mat[1] = normalize(new_y_axis);
    orient_mat[2] = normalize(new_z_axis);
    vector <vector <double> > tmp = dotp(rot_mat, vector <vector <double> >(orient_mat.begin(), orient_mat.begin() + 3));
    orient_mat[3] = tmp[0];
    orient_mat[4] = tmp[1];
    orient_mat[5] = tmp[2];

    // The basis vectors need to have their transposes applied independently
    vector <vector <double> > dir_vec = dotp(orient_mat, scale(basis_vectors.at(crystal_structure), a0));
    tmp = vector <vector <double> > (dir_vec.begin(), dir_vec.begin() + 3);
    tmp = transpose(tmp);
    for (size_t i = 0; i < tmp.size(); ++i) dir_vec[i] = tmp[i];
    tmp = vector <vector <double> > (dir_vec.begin() + 3, dir_vec.end());
    tmp = transpose(tmp);
    for (size_t i = 0; i < tmp.size(); ++i) dir_vec[i+3] = tmp[i];


    // This code is adapted from LAMMPS fix_eco_force.cpp.
    double vol = 0.5 / PI *
              (dir_vec[0][0] * (dir_vec[1][1] * dir_vec[2][2] - dir_vec[2][1] * dir_vec[1][2]) + 
               dir_vec[1][0] * (dir_vec[2][1] * dir_vec[0][2] - dir_vec[0][1] * dir_vec[2][2]) + 
               dir_vec[2][0] * (dir_vec[0][1] * dir_vec[1][2] - dir_vec[1][1] * dir_vec[0][2]));
    double i_vol = 1.0 / vol;


    reciprocal_vectors[0][0][0] = (dir_vec[1][1] * dir_vec[2][2] - dir_vec[2][1] * dir_vec[1][2]) * i_vol;
    reciprocal_vectors[0][0][1] = (dir_vec[1][2] * dir_vec[2][0] - dir_vec[2][2] * dir_vec[1][0]) * i_vol;
    reciprocal_vectors[0][0][2] = (dir_vec[1][0] * dir_vec[2][1] - dir_vec[2][0] * dir_vec[1][1]) * i_vol;
    reciprocal_vectors[0][1][0] = (dir_vec[2][1] * dir_vec[0][2] - dir_vec[0][1] * dir_vec[2][2]) * i_vol;
    reciprocal_vectors[0][1][1] = (dir_vec[2][2] * dir_vec[0][0] - dir_vec[0][2] * dir_vec[2][0]) * i_vol;
    reciprocal_vectors[0][1][2] = (dir_vec[2][0] * dir_vec[0][1] - dir_vec[0][0] * dir_vec[2][1]) * i_vol;
    reciprocal_vectors[0][2][0] = (dir_vec[0][1] * dir_vec[1][2] - dir_vec[1][1] * dir_vec[0][2]) * i_vol;
    reciprocal_vectors[0][2][1] = (dir_vec[0][2] * dir_vec[1][0] - dir_vec[1][2] * dir_vec[0][0]) * i_vol;
    reciprocal_vectors[0][2][2] = (dir_vec[0][0] * dir_vec[1][1] - dir_vec[1][0] * dir_vec[0][1]) * i_vol;

    vol = 0.5  / PI * 
        (dir_vec[3][0] * (dir_vec[4][1] * dir_vec[5][2] - dir_vec[5][1] * dir_vec[4][2]) + 
         dir_vec[4][0] * (dir_vec[5][1] * dir_vec[3][2] - dir_vec[3][1] * dir_vec[5][2]) + 
         dir_vec[5][0] * (dir_vec[3][1] * dir_vec[4][2] - dir_vec[4][1] * dir_vec[3][2]));
    i_vol = 1.0 / vol;
    reciprocal_vectors[1][0][0] = (dir_vec[4][1] * dir_vec[5][2] - dir_vec[5][1] * dir_vec[4][2]) * i_vol;
    reciprocal_vectors[1][0][1] = (dir_vec[4][2] * dir_vec[5][0] - dir_vec[5][2] * dir_vec[4][0]) * i_vol;
    reciprocal_vectors[1][0][2] = (dir_vec[4][0] * dir_vec[5][1] - dir_vec[5][0] * dir_vec[4][1]) * i_vol;
    reciprocal_vectors[1][1][0] = (dir_vec[5][1] * dir_vec[3][2] - dir_vec[3][1] * dir_vec[5][2]) * i_vol;
    reciprocal_vectors[1][1][1] = (dir_vec[5][2] * dir_vec[3][0] - dir_vec[3][2] * dir_vec[5][0]) * i_vol;
    reciprocal_vectors[1][1][2] = (dir_vec[5][0] * dir_vec[3][1] - dir_vec[3][0] * dir_vec[5][1]) * i_vol;
    reciprocal_vectors[1][2][0] = (dir_vec[3][1] * dir_vec[4][2] - dir_vec[4][1] * dir_vec[3][2]) * i_vol;
    reciprocal_vectors[1][2][1] = (dir_vec[3][2] * dir_vec[4][0] - dir_vec[4][2] * dir_vec[3][0]) * i_vol;
    reciprocal_vectors[1][2][2] = (dir_vec[3][0] * dir_vec[4][1] - dir_vec[4][0] * dir_vec[3][1]) * i_vol;

    double rsq, weight, wsum, scalar_product;
    vector <double> dr = vector <double> (3, 0.0); // displacement vector
    vector <double> reesum = vector <double> (3, 0.0); // sum of real part
    vector <double> imesum = vector <double> (3, 0.0); // sum of imaginary part
    int max_co = 4; // will produce wrong results for rcut > 3 * lattice constant
    int neigh = 0; // number of neighbors
    double cut_sq = r_cut * r_cut * a0 * a0;
    double inv_cut_sq = 1.0 / cut_sq;

    int i, k, idx[3];
    for (idx[0] = -max_co; idx[0] <= max_co; ++idx[0]) {
      for (idx[1] = -max_co; idx[1] <= max_co; ++idx[1]) {
        for (idx[2] = -max_co; idx[2] <= max_co; ++idx[2]) {
          // distance of atoms
          for (i = 0; i < 3; ++i) {
            dr[i] = dir_vec[0][i] * idx[0] + dir_vec[1][i] * idx[1] + dir_vec[2][i] * idx[2];
          }
          rsq = dotp(dr, dr);

          if ((rsq != 0.0) && (rsq < cut_sq)) {
            ++neigh;
            rsq *= inv_cut_sq;
          
            weight =  rsq * (rsq - 2.0) + 1.0;
            wsum += weight;

            for (k = 0; k < 3; ++k) {
              scalar_product = dotp(reciprocal_vectors[1][k], dr);
              reesum[k] += weight * cos(scalar_product);
              imesum[k] -= weight * sin(scalar_product);
            } // k
          } // within rsq
        } // idx[2]
      } // idx[1]
    } // idx[0]
    
    norm_fac = 3.0 * wsum * wsum;
    for (k = 0; k < 3; ++k) {
      norm_fac -= reesum[k] * reesum[k] + imesum[k] * imesum[k];
    }

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


// function prototypes
vector <double> getAxisData(istream&);
void parseInputFile(inputData&, const string&);
vector <int> getAtomIdsList(const string&);
void printInputFileHelp();
void processData(const inputData&);
vector <Atom> getAtomData(ifstream&, const inputData&, const boxData&);
vector <vector <int> > generateCellLinkedList(const vector <Atom>&,
                                              const inputData&,
                                              const boxData&);
vector <double> calculateSymmetryParameter(vector <Atom>&, 
                                        const vector <vector <int> >&,
                                        const vector <bool>,
                                        const inputData&,
                                        const boxData&);
void writeAtomsToFile(const string&, const vector <Atom>&,
                      const vector <bool>, const vector <double>&,
                      const inputData&, const boxData&);


int main(int argc, char** argv)
{
  omp_set_num_threads(THREAD_NUM);
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
        ("q,quiet", "Suppress warnings", cxxopts::value<bool>(show_warnings)->implicit_value("false"))
        ("atom-ids", "Only output the results for the specified IDs from the file (containing a list of atom ids)", cxxopts::value<string>(), "file")
        ("print-nearest-neighbors", "Print the nearest neighbor list to a file", cxxopts::value<string>(input.nn_filebase)->default_value("nearest_neighbors_*.txt"), "file")
        ("n,eta", "Cutoff parameter defined in Ulomek et al.", cxxopts::value<double>(input.eta)->default_value("0.25"), "eta")
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

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


void processData(const inputData& input) {
  boxData box;
  int n_grain_1, n_grain_2, n_gb; // number of atoms assigned to each grain
  int N; // number of atoms
  string str, filename2; // junk string variable, filename of the data with the grain assignment

  vector <Atom> atoms; // the atoms from the data file
  vector <bool> allowed_atoms; // the atoms that will be printed
  vector <vector <int> > iatom; // cell-linked list
  vector <double> symm; // the orientation factors

  VerifyNewFile verify1(input.outfile);
  ofstream fout_data(verify1.validNewFile().c_str());
  checkFileStream(fout_data, verify1.validNewFile());

  fout_data << "# Data consists of: [timestep/file, atoms in grain 1, atoms in grain 2"
            << ", atom in GB]\n";

  // VerifyNewFile verify2("misorientation_data.txt");
  // ofstream fout_misorientation(verify2.validNewFile().c_str());
  // checkFileStream(fout_misorientation, verify2.validNewFile());

  int aa = 1;
  int x_index = -1, y_index = -1, z_index = -1; // indices for the x, y, and z positions of the atoms
  int id_index = -1, type_index = -1, charge_index = -1; // indices for the atom id, type, and charge
  for (vector <string>::const_iterator it = input.files.begin(); it != input.files.end(); ++it) { // TODO: parallelize this
    box.reset();
    n_grain_1 = 0;
    n_grain_2 = 0;
    n_gb = 0;
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
    symm = calculateSymmetryParameter(atoms, iatom, allowed_atoms, input, box);

    for (size_t i = 0; i < atoms.size(); ++i) {
      if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}
      if (!allowed_atoms[i]) {continue;}
      if (atoms[i].getMark() == 1) {++n_grain_1;}
      else if (atoms[i].getMark() == 2) {++n_grain_2;}
      else if (atoms[i].getMark() == 3) {++n_gb;}
      else {
        cerr << "Error: Unrecognized grain assignment.\n";
        exit(BOUNDS_ERROR);
      }
    }

    fout_data << n_grain_1 << " " << n_grain_2 << " " << n_gb << "\n";

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
         << "N = " << atoms.size() << " != n_atoms_read = " << n_atoms_read << "\n";
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

  #pragma omp parallel for collapse(3)
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


vector <double> calculateSymmetryParameter(vector <Atom>& atoms,
                                        const vector <vector <int> >& iatom,
                                        const vector <bool> allowed_atoms,
                                        const inputData& input,
                                        const boxData& box) {
  double x, y, z; // position of the central atom
  int id; // neighbor atom id;
  vector <double> symm (atoms.size(), 0.0);
  int k; // variable to loop over 3 reciprocal directions
  int lambda; // variable to loop over 2 crystals
  int dim; // variable to loop over 3 spatial components
  double dx, dy, dz; // stores current interatomic vector
  double squared_distance; // stores current squared distance
  // double chi; // stores current order parameter
  double weight; // stores current weight function
  double scalar_product; // stores current scalar product
  // double omega; // phase of sine transition
  // double omega_pre = (PI / 2.0) * (1.0 / input.eta); // prefactor for omega
  vector <double> dr (3, 0.0);

  double rcut = input.a0 * input.r_cut;
  double cutsq = rcut * rcut;
  double inv_cutsq = 1.0 / cutsq;

  #pragma omp parallel for
  for (size_t i = 0; i < atoms.size(); ++i) {
    // move to the next atom if ignoring this atom type
    if (find(input.ignored_atoms.begin(), input.ignored_atoms.end(), atoms[i].getType()) != input.ignored_atoms.end()) {continue;}

    // move to the next atom if not calculating for this atom
    if (!allowed_atoms[i]) {continue;}

    vector <vector <double> > real_phi (2, vector <double> (3, 0.0));
    vector <vector <double> > imag_phi (2, vector <double> (3, 0.0));;

    // Store the position of this atom
    x = atoms[i].getWrapped()[0];
    y = atoms[i].getWrapped()[1];
    z = atoms[i].getWrapped()[2];

    for (int l = 1; l <= iatom[0][i]; ++l) {
      id = iatom[l][i];
      // calculate the distance between these atoms
      dx = x - atoms[id].getWrapped()[0];
      dy = y - atoms[id].getWrapped()[1];
      dz = z - atoms[id].getWrapped()[2];

      // Apply PBCs
      dx = dx - anInt(dx / box.Lx) * box.Lx;
      dy = dy - anInt(dy / box.Ly) * box.Ly;
      dz = dz - anInt(dz / box.Lz) * box.Lz;

      dr = {dx, dy, dz};

      squared_distance = (dx * dx + dy * dy + dz * dz) / cutsq;
      if (squared_distance < cutsq) {
        squared_distance *= inv_cutsq;
        weight = squared_distance * (squared_distance - 2.0) + 1.0;
        
        for (lambda = 0; lambda < 2; ++lambda) {
          for (k = 0; k < 3; ++k) {
            scalar_product = dotp(input.reciprocal_vectors[lambda][k], dr);
            real_phi[lambda][k] += weight * cos(scalar_product);
            imag_phi[lambda][k] += weight * sin(scalar_product);
          } // k
        } // lambda
      } // below cutoff distance
    } // neighbor atoms

    for (k = 0; k < 3; ++k) {
      symm[i] += (real_phi[0][k] * real_phi[0][k] + imag_phi[0][k] * imag_phi[0][k] -
                  real_phi[1][k] * real_phi[1][k] - imag_phi[1][k] * imag_phi[1][k]);
    }

    symm[i] /= input.norm_fac;
    if (symm[i] > input.eta) {atoms[i].setMark(1);}
    else if (symm[i] < -input.eta) {atoms[i].setMark(2);}
    else {atoms[i].setMark(3);}
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


void parseInputFile(inputData& input, const string& input_file) {
  string str; // junk string variable

  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);

  getline(fin, str);
  stringstream ss(str);
  if (!(ss >> input.theta >> input.n_types >> input.r_cut >> input.a0 >> input.crystal_structure)) {
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
       << "\n  theta = " << input.theta
       << "\n  n_types = " << input.n_types
       << "\n  r_cut = " << input.r_cut
       << "\n  a0 = " << input.a0
       << "\n  eta = " << input.eta
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

  input.calculateReciprocalLattices();

  // get the files
  while (getline(fin, str)) {
    if (!str.empty()) {
      input.files.push_back(str);
    }
  }
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


void printInputFileHelp() {
  cout << "\n\nThe first line of the input file must contain the following items in order:\n"
       << "  1. The misorientation angle of the embedded grain\n"
       << "  2. The number of atom types in the simulation.\n"
       << "  3. The cutoff distance (r_cut) for generating the nearest neighbor list\n"
       << "     in terms of the lattice parameter.\n"
       << "  4. The lattice parameter in Angstroms.\n"
       << "  5. The crystal structure of the lattice (fcc, bcc, sc).\n\n"
       << "These lines are then followed by the orientation matrix of the system, i.e.\n"
       << "  1 0 0 (orientation of the x axis)\n"
       << "  0 1 0 (orientation of the y axis)\n"
       << "  0 0 1 (orientation of the z axis)\n\n"
       << "The list of files (one file per line) should be placed afterwards.\n\n";
}

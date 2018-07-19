#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath> // for cos, sin
#include <cstdlib>
#include <cstdio>
#include <cxxopts.hpp>
#include "atom.h"
#include "error_code_defines.h"

using namespace std;

#define PI 3.141592653589793

bool warnings = true;
bool print_nearest_neighbors = false;

struct inputVars
{
  int n_files, n_types;
  double theta, r_cut, fi_cut, a0;
  double r_cut_sq;
  vector <double> new_x_axis{0,0,0}, new_y_axis{0,0,0}, new_z_axis{0,0,0};
  vector <double> old_x_axis{0,0,0}, old_z_axis{0,0,0};

  void calculateOldX()
  {
    old_x_axis[0] = new_x_axis[0];
    old_x_axis[1] = new_y_axis[0];
    old_x_axis[2] = new_z_axis[0];
  }

  void calculateOldZ()
  {
    old_z_axis[0] = new_x_axis[2];
    old_z_axis[1] = new_y_axis[2];
    old_z_axis[2] = new_z_axis[2];
  }

  void calculateRCutSq()
  {
    r_cut_sq = r_cut * r_cut;
  }
} input;

struct boxData
{
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double Lx, Ly, Lz;

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
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

double calculateAxisMagnitude(vector <double>& axis)
{
  return sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
}

void normalize(vector <double>& axis)
{
  double sum = 0;
  double magnitude = calculateAxisMagnitude(axis);

  if (axis.size() > 3)
  {
    cout << "Vector size is too large (" << axis.size() << ").  Cannot normalize.\n";
    exit(VECTOR_SIZE_ERROR);
  }

  for (unsigned int i = 0; i < axis.size(); ++i)
  {
    axis[i] /= magnitude;
  }
}

vector <double> getAxisData(istream& fin)
{
  string str; // junk string variable
  vector <double> axis (3,0);
  getline(fin, str);
  stringstream ss(str);
  if (!(ss >> axis[0] >> axis[1] >> axis[2]))
  {
    cout << "Error reading axis.\n";
    exit(AXIS_READ_ERROR);
  }
  return axis;
}

vector <string> parseInputFile(const string& input_file)
{
  string str; // junk string variable
  vector <string> files; // filenames to be processed

  ifstream fin_input(input_file.c_str());
  checkFileStream(fin_input, input_file);

  getline(fin_input, str);
  stringstream ss_input(str);
  if (!(ss_input >> input.n_files >> input.theta >> input.n_types >> input.r_cut >> input.fi_cut >> input.a0))
  {
    cout << "Error reading the input file.  Did you forget a value?\n"
         << "The first line of the input file must contain the following six items:\n"
         << "\t1. The number of files to be processed.\n"
         << "\t2. The misorientation angle of the grains with respect to each other.\n"
         << "\t3. The number of atom types in the simulation.\n"
         << "\t4. The cutoff distance for determining grain assignment in terms of the lattice parameter.\n"
         << "\t5. The cutoff value for which orientation parameters are assigned to each grain.\n"
         << "\t6. The lattice parameter in Angstroms.\n";
    exit(INPUT_FORMAT_ERROR);
  }
  input.calculateRCutSq();
  // This stores the transformation matrix from the <100> orientation to the new orientation
  input.new_x_axis = getAxisData(fin_input);
  input.new_y_axis = getAxisData(fin_input);
  input.new_z_axis = getAxisData(fin_input);

  cout << "Input parameters:\n"
       << "\tn_files = " << input.n_files << endl
       << "\ttheta = " << input.theta << endl
       << "\tn_types = " << input.n_types << endl
       << "\tr_cut = " << input.r_cut << endl
       << "\tfi_cut = " << input.fi_cut << endl
       << "\ta0 = " << input.a0 << endl;
  cout << "\tRotated coordinate system:\n"
       << "\t  x = " << input.new_x_axis[0] << " " << input.new_x_axis[1] << " " << input.new_x_axis[2] << endl
       << "\t  y = " << input.new_y_axis[0] << " " << input.new_y_axis[1] << " " << input.new_y_axis[2] << endl
       << "\t  z = " << input.new_z_axis[0] << " " << input.new_z_axis[1] << " " << input.new_z_axis[2] << endl;

  normalize(input.new_x_axis);
  normalize(input.new_y_axis);
  normalize(input.new_z_axis);

  // This calculates the transformation from the new orientation to the <100> orientation
  // Used later in the program.
  input.calculateOldX();
  input.calculateOldZ();
  normalize(input.old_x_axis);
  normalize(input.old_z_axis);

  int aa = 1;
  while (getline(fin_input, str) && aa <= input.n_files) // Get the files desired
  {
    if (!str.empty()) // empty lines are ignored.
    {
      files.push_back(str);
      ++aa;
    }
  }

  // warn if n_files and files.size() are not equal.
  if (warnings && files.size() != input.n_files)
  {
    cout << "Warning: " << input.n_files << " files specified in input file, but only " << files.size() << " found.\n";
  }

  return files;
}

void getAtomData(ifstream& fin, vector <Atom>& atoms)
{
  int n_atoms_read = 0; // number of atoms read
  string str; // holds the info
  vector <double> data; // holds the number of elements
  int atom_id, type;
  double charge, x, y, z, xu, yu, zu;
  while (getline(fin, str))
  {
    stringstream ss(str);
    double dummy;
    while (ss >> dummy)
    {
      data.push_back(dummy); // a way to count the number of data elements
    }
    atom_id = (int)(data[0]);
    type = (int)(data[1]);
    switch (data.size())
    {
      // atom has charge, and unwrapped as well as wrapped coordinates
      case 9: charge = data[2];
              x = data[3]; y = data[4]; z = data[5];
              xu = data[6]; yu = data[7]; zu = data[8];
              break;

      // atom does not have charge, and unwrapped as well as wrapped coordinates
      case 8: charge = 0.0;
              x = data[2]; y = data[3]; z = data[4];
              xu = data[5]; yu = data[6]; zu = data[7];
              break;

      // atom has charge, and only wrapped coordinates
      case 6: charge = data[2];
              x = data[3]; y = data[4]; z = data[5];
              xu = 0.0; yu = 0.0; zu = 0.0;
              break;
      // atom does not have charge, and only wrapped coordinates
      case 5: charge = 0.0;
              x = data[2]; y = data[3]; z = data[4];
              xu = 0.0; yu = 0.0; zu = 0.0;
              break;
      default: cout << "Unrecognized file format.  Expected format: id type charge* x y z xu* yu* zu*.  Note, elements marked with * are optional.\n";
               exit(FILE_FORMAT_ERROR);
    }
    data.clear();

    if (type > input.n_types)
    {
      cout << "Error: unexpected atom type.\n"
           << "n_types = " << input.n_types << " < this atom's type = " << type << endl;
      exit(ATOM_TYPE_ERROR);
    }

    // We adjust the positions so that we start at the origin so we can easily
    // assign to cells
    x = x / input.a0 - box.xlow;
    y = y / input.a0 - box.ylow;
    z = z / input.a0 - box.zlow;

    // We make the atom id (almost) match thre index.  There is a difference of 1
    // because C++ indexes from 0
    atoms[atom_id - 1] = Atom(atom_id, type, charge, x, y, z);
    // The unwrapped coordinates we do not do anything with here, so we simply
    // transcribe them.
    atoms[atom_id - 1].setXu(xu);
    atoms[atom_id - 1].setYu(yu);
    atoms[atom_id - 1].setZu(zu);
    ++n_atoms_read;
  }
  fin.close();

  if (n_atoms_read != atoms.size())
  {
    cout << "Error: number of atoms read does not match number of atoms in the simulation.\n"
         << "N = " << atoms.size() << " != n_atoms_read = " << n_atoms_read << endl;
    exit(ATOM_COUNT_ERROR);
  }
}

void generateCellLinkedList(const vector <Atom>& atoms, vector <vector <int> >& iatom)
{
  int ncellx, ncelly, ncellz; // number of cells in each direction
  int idx, idy, idz; // cell number in each direction
  double lcellx, lcelly, lcellz; // length of cells in each direction
  int n_atoms_per_cell; // number of atoms allowed per cell
  double drij_sq, rxij, ryij, rzij; // square of distance, x, y, and z separation.
  vector <vector <vector <int> > > icell; // cell index
  vector <vector <vector <vector <int> > > > pcell; // atom index in each cell

  // First we generate the number of cells in each direction
  ncellx = (int)(box.Lx / input.r_cut) + 1;
  ncelly = (int)(box.Ly / input.r_cut) + 1;
  ncellz = (int)(box.Lz / input.r_cut) + 1;

  // Length of cells in each direction
  lcellx = box.Lx / ncellx;
  lcelly = box.Ly / ncelly;
  lcellz = box.Lz / ncellz;

  // Minimum number of atoms allowed of 100
  n_atoms_per_cell = max((int)(atoms.size() / (double)(3 * ncellx * ncelly * ncellz)), 100);

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
  for (unsigned int i = 0; i < atoms.size(); ++i) // Look at each atom
  {
    //if (atoms[i].getType() == 2) continue; // ignoreing O atoms
    if (atoms[i].getType() != 1) {continue;} // Only want U atoms.  NOTE: This needs to change for substitutional defects!
    // Assign this atom to a cell
    // Rounds towards 0 with a type cast
    idx = (int)(atoms[i].getX() / lcellx); // assign the x cell
    idy = (int)(atoms[i].getY() / lcelly); // assign the y cell
    idz = (int)(atoms[i].getZ() / lcellz); // assign the z cell
    // Check if we went out of bounds
    // C++ indexes from 0, so we have to subtract 1 from the maximum value to
    // stay within our memory bounds
    if (idx >= ncellx) idx = ncellx - 1;
    if (idy >= ncelly) idy = ncelly - 1;
    if (idz >= ncellz) idz = ncellz - 1;

    // This is for unwrapped coordinates
    while (idx < 0.0) idx = ncellx + idx; // Note that this keeps things within the bounds set by lcellx
    while (idy < 0.0) idy = ncelly + idy; // Note that this keeps things within the bounds set by lcelly
    while (idz < 0.0) idz = ncellz + idz; // Note that this keeps things within the bounds set by lcellz

    ++icell[idx][idy][idz]; // increase the number of atoms in this cell
    // assign the atom number to this index.
    pcell[idx][idy][idz][icell[idx][idy][idz] - 1] = i;
  }

  for (int i = 0; i < ncellx; ++i) // For each x cell
  {
    for (int j = 0; j < ncelly; ++j) // For each y cell
    {
      for (int k = 0; k < ncellz; ++k) // For each z cell
      {
        for (int l = 0; l < icell[i][j][k]; ++l) // For each atom in this cell
        {
          int id = pcell[i][j][k][l]; // store this atom id
          // Now we check each sub cell around the current one
          for (int ii = -1; ii < 2; ++ii) // allowed values: -1, 0, and 1
          {
            for (int jj = -1; jj < 2; ++jj)
            {
              for (int kk = -1; kk < 2; ++kk)
              {
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
                for (int m = 0; m < icell[ia][ja][ka]; ++m)
                {
                  int jd = pcell[ia][ja][ka][m];
                  // If jd <= id, we've already dealt with this interaction
                  if (jd <= id)
                  {
                    continue;
                  }

                  // Now the actual calculations!
                  rxij = atoms[id].getX() - atoms[jd].getX();
                  ryij = atoms[id].getY() - atoms[jd].getY();
                  rzij = atoms[id].getZ() - atoms[jd].getZ();

                  // Apply PBCs
                  rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
                  ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
                  rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

                  // Now calculate the distance
                  drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                  if (drij_sq > input.r_cut_sq)
                  {
                    continue; // move to the next atom if we're too far away
                  }

                  // Create the neighbor list
                  iatom[0][id] += 1; //for atom id - number of neighbors
                  iatom[(iatom[0][id])][id] = jd; // point to the next atom
                  iatom[0][jd] += 1; // for atom jd
                  iatom[(iatom[0][jd])][jd] = id;
                } // m
              } //kk
            } //jj
          } //ii
        } // l
      } // k
    } // j
  } // i
}

void writeAtomsToFile(const string& filename, const vector <Atom>& atoms, const vector <double>& symm)
{
  int n_atoms_read = 0;

  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);

  fout << "VARIABLES = \"Atom ID\", \"Atom Type\", ";
  if (input.n_types == 2) {fout << "\"Atom Charge\", ";}
  fout << "\"X\", \"Y\", \"Z\", \"Grain Number\", \"Orientation Parameter\",\"Xu\", \"Yu\", \"Zu\"\n";

  // This writes things in a Tecplot-readable FILE_FORMAT_ERROR
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getType() != 1) // only type 1 atoms... TODO: Change this later.
    {
      continue;
    }
    fout << atoms[i].getId() << " "
         << atoms[i].getType() << " ";
    if (input.n_types == 2) {fout << atoms[i].getCharge() << " ";}
    fout << (atoms[i].getX() + box.xlow) * input.a0 << " "
         << (atoms[i].getY() + box.ylow) * input.a0 << " "
         << (atoms[i].getZ() + box.zlow) * input.a0 << " "
         << atoms[i].getMark() << " "
         << symm[i] << " " << atoms[i].getXu() << " " << atoms[i].getYu()
         << " " << atoms[i].getZu() << endl;
    ++n_atoms_read;
  }
  fout.close();
  if (input.n_types == 2) {n_atoms_read *= 3;}
  if (n_atoms_read != atoms.size())
  {
    cout << "Error: number of atoms written does not match number of atoms in the simulation.\n"
         << "N = " << atoms.size() << " != n_atoms_read = " << n_atoms_read << endl;
    exit(ATOM_COUNT_ERROR);
  }
}

void processData(vector <string>& files, const cxxopts::ParseResult& result)
{
  int N; // Number of atoms
  double xlow, xhigh, ylow, yhigh, zlow, zhigh; // Box bounds from the data files
  double xy, xz, yz; // tilt factors
  string filename2, str; // file to write atom data to, junk string variable
  vector <Atom> atoms;
  vector <vector <int> > iatom; // cell-linked list
  vector <double> symm;
  double x, y, z, rxij, ryij, rzij, drij_sq, xtemp, ytemp, ztemp;
  double costheta_sq;
  double coeffs [2] = {3.0, 2.0};
  int n_grain_1, n_grain_2;

  ofstream fout_data((result["output"].as<string>()).c_str());
  checkFileStream(fout_data, result["output"].as<string>());

  fout_data << "# Data consists of: [timestep, atoms in grain 1, atoms in grain 2]\n";

  int aa = 1;
  for (vector <string>::iterator it = files.begin(); it != files.end(); ++it)
  {

    ifstream fin((*it).c_str());
    checkFileStream(fin, *it);

    // Parse the file based on it's contents.  Note that if "dump" appears anywhere
    // in the filename, this script will assume it's a dump file.  If dump is not
    // found, it will look for "dat", and assume it is a LAMMPS input file.
    if ((result["type"].as<string>()).compare("dump") == 0)
    {
      int timestep;
      getline(fin, str); // Gets ITEM: TIMESTEP
      fin >> timestep; // Gets the timestep number
      fin.ignore();
      if (aa == 1)
      {
        if (warnings && timestep != 0)
        {
          cout << "Warning: first data file is not at timestep 0! "
               << "Ignore this warning if this is intentional.\n";
        }
      }
      getline(fin, str); // Gets ITEM: NUMBER OF ATOMS
      fin >> N;
      fin.ignore();
      getline(fin, str); //get ITEM: BOX BOUNDS
      if (str.find("xy") == string::npos)
      {
        fin >> xlow >> xhigh;
        fin >> ylow >> yhigh;
        fin >> zlow >> zhigh;
      }
      else  // We also get the triclinic tilt factors for triclinic systems
      {
        fin >> xlow >> xhigh >> xy;
        fin >> ylow >> yhigh >> xz;
        fin >> zlow >> zhigh >> yz;
      }
      fin.ignore();
      getline(fin, str); // Gets ITEM: ATOMS <data types>
      if (!result.count("append"))
      {
        filename2 = (*it).substr(0,(*it).find(".dump")) + "_interface.dat";
      }
      else
      {
        filename2 = "temp.dat";
      }

      fout_data << timestep << " ";
    }
    else if ((result["type"].as<string>()).compare("dat") == 0)
    {
      getline(fin, str); // gets comment line
      fin >> N >> str; // number of atoms
      fin >> str >> str >> str; // Number of types is specified in the input file
      fin >> xlow >> xhigh >> str >> str;
      fin >> ylow >> yhigh >> str >> str;
      fin >> zlow >> zhigh >> str >> str;
      fin.ignore();
      getline(fin,str);
      if (str.find("xy") != string::npos) // Triclinic tilt factors
      {
        fin >> xy >> xz >> yz >> str >> str >> str;
        fin.ignore();
      }
      getline(fin, str);
      if (!result.count("append"))
      {
        filename2 = (*it).substr(0,(*it).find(".dat")) + "_interface.dat";
      }
      else
      {
        cout << "This will be implemented later.\n";
        filename2 = (*it).substr(0,(*it).find(".dat")) + "_interface.dat";
        //filename2 = "temp.dat";
      }

      fout_data << (*it).substr(0, (*it).find(".dat"));
    }
    else
    {
      cout << "File type error";
      exit(FILE_FORMAT_ERROR);
    }

    // Convert the bounds in terms of a0
    xlow /= input.a0;
    xhigh /= input.a0;
    ylow /= input.a0;
    yhigh /= input.a0;
    zlow /= input.a0;
    zhigh /= input.a0;

    box.xlow = xlow; box.xhigh = xhigh;
    box.ylow = ylow; box.yhigh = yhigh;
    box.zlow = zlow; box.zhigh = zhigh;
    box.calculateBoxLengths();

    atoms.resize(N, Atom());
    iatom.clear();
    getAtomData(fin, atoms);
    generateCellLinkedList(atoms, iatom);
    if (print_nearest_neighbors)
    {
      stringstream neighbor;
      string neighbor_filename;
      neighbor << "nearest_neighbor_list_" << aa << ".txt";
      neighbor >> neighbor_filename;
      ofstream fout_neighbors(neighbor_filename.c_str());
      checkFileStream(fout_neighbors, neighbor_filename);
      for (unsigned int i = 0; i < atoms.size(); ++i)
      {
        vector <int> neighs;
        fout_neighbors << "Atom " << i << " has " << iatom[0][i] << " neighbors:\n";
        for (int l = 1; l <= iatom[0][i]; ++l)
        {
          neighs.push_back(iatom[l][i]);
        }
        sort(neighs.begin(), neighs.end());
        for (unsigned int i = 0; i < neighs.size(); ++i)
        {
          fout_neighbors << "\tAtom " << neighs[i] << endl;
        }
      }
    }

    // Now that we have the atoms safely stored, we can process them.
    symm.resize(atoms.size(), 0); // Assign each atom a symmetry parameter
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      // if (atoms[i].getType() == 2) // ignore O atoms.
      if (atoms[i].getType() != 1)
      {
        continue; // Only focus on U (or the single atom type)
      }

      x = atoms[i].getX();
      y = atoms[i].getY();
      z = atoms[i].getZ();

      // We start at l = 1 because if we start at l = 0, we just re-use the same
      // atom over and over.
      for (int l = 1; l <= iatom[0][i]; ++l)
      {
        unsigned int id = iatom[l][i];

        rxij = atoms[id].getX() - x;
        ryij = atoms[id].getY() - y;
        rzij = atoms[id].getZ() - z;

        // Apply PBCs.  Note that applying PBCs with the positions projected
        // in the <100> reference frame messes up the calculations!
        rxij = rxij - anInt(rxij / box.Lx) * box.Lx;
        ryij = ryij - anInt(ryij / box.Ly) * box.Ly;
        rzij = rzij - anInt(rzij / box.Lz) * box.Lz;

        // Project this vector onto the (001) plane
        // NOTE: In order for this method to work well, the correct cutoff distance needs to be used!
        xtemp = (rxij * input.old_z_axis[0] + ryij * input.old_z_axis[1] + rzij * input.old_z_axis[2]) * input.old_z_axis[0];
        ytemp = (rxij * input.old_z_axis[0] + ryij * input.old_z_axis[1] + rzij * input.old_z_axis[2]) * input.old_z_axis[1];
        ztemp = (rxij * input.old_z_axis[0] + ryij * input.old_z_axis[1] + rzij * input.old_z_axis[2]) * input.old_z_axis[2];

        rxij -= xtemp;
        ryij -= ytemp;
        rzij -= ztemp;


        // Calculate the magnitude of the distance
        drij_sq = rxij * rxij + ryij * ryij + rzij * rzij;

        if (drij_sq < 1.0E-8) // Handles the case where the projected position of the atom is right on top of the current atom.
        {
          symm[i] += 1;
          //cout << "Note: drij_sq = 0.0\n";
          continue;
        }
        // cos = dot(A,B) / (|A|*|B|)
        double dotp = (rxij * input.old_x_axis[0] + ryij * input.old_x_axis[1] + rzij * input.old_x_axis[2]);
        costheta_sq = (dotp * dotp) / (drij_sq);
        double val = (coeffs[0] - coeffs[1] * costheta_sq) * (coeffs[0] - coeffs[1] * costheta_sq) * costheta_sq;
        symm[i] += val;
      }
    }
    for (unsigned int i = 0; i < symm.size(); ++i)
    {
      if (atoms[i].getType() != 1) // We are making a major assumption here that the atom assigned as type 1 has the important crystallographic structure
      {
        continue;
      }
      // On the off chance that no neighbors have been assigned to this atom,
      // we estimate the symmetry parameter to be the previous atom's value.
      // This assumes that atoms with ID's close to each other are close to
      // each other in the simulation, which may not always be the case.
      else if (iatom[0][i] == 0 && i != 0)
      {
        if (warnings)
        {
          cout << "\nWarning: no neighbors detected for atom " << atoms[i].getId() << ".  Using the symmetry parameter of the previous atom = " << symm[i - 1] << endl;
        }
        symm[i] = symm[i - 1];
      }
      else
      {
        symm[i] /= iatom[0][i];
      }

      if (atoms[i].getType() != 1)
      {
        continue;
      }
      if (symm[i] <= input.fi_cut)
      {
        atoms[i].setMark(1);
        ++n_grain_1;
      }
      else // greater than the cutoff
      {
        atoms[i].setMark(2);
        ++n_grain_2;
      }
    }

    fout_data << n_grain_1 << " " << n_grain_2 << endl;

    writeAtomsToFile(filename2, atoms, symm);

    fin.close();
    if (result.count("append"))
    {
      cout << "This will be implemented later.\n";
//      rename((*it).c_str(), (*it + ".bak").c_str()); // make a backup of the original file
//      rename(filename2.c_str(), (*it + "_v2").c_str()); // change the temp file to the original filename
    }
    cout << "\r";
    cout << "Processing of file\"" << (*it) << "\" completed.";
    cout.flush();
    ++aa;
  }
  fout_data.close();
}

int main(int argc, char** argv)
{
  bool append, quiet, calculate_microrotation;
  string input_file, output_file, filetype;
  vector <string> filenames;
  try
  {
    cxxopts::Options options(argv[0], "Script that can determine various properties per atom.");
    options
      .positional_help("File")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "Input file", cxxopts::value<string>(input_file), "file")
        ("a,append", "Append to the processed data file.", cxxopts::value<bool>(append)->default_value("false"))
        ("t,type", "Flag specifying the input file type - LAMMPS input file (dat), or LAMMPS dump file (dump)", cxxopts::value<string>(filetype)->default_value("dump"))
        ("m,microrotation", "Flag specifying that the user wants the microrotation parameter calculated. Not yet implemented", cxxopts::value<bool>(calculate_microrotation)->default_value("false")->implicit_value("true"))
        ("o,output", "Output file for calculated data", cxxopts::value<string>(output_file)->default_value("data.txt"), "file")
        ("print-nearest-neighbors", "Print the nearest neighbor list to a file", cxxopts::value<string>()->implicit_value("nearest_neighbors_*.txt"), "file")
        ("q,quiet", "Suppress warnings from the code", cxxopts::value<bool>(quiet)->implicit_value("true")->default_value("false"))
        ("h,help", "Show the help");

    options.parse_positional({"file"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("file") == 0)
    {
      cout << options.help() << endl << endl
           << "The first line of the input file must contain the following six items:\n"
           << "\t1. The number of files to be processed.\n"
           << "\t2. The misorientation angle of the grains with respect to each other.\n"
           << "\t3. The number of atom types in the simulation.\n"
           << "\t4. The cutoff distance (r_cut) for determining grain assignment in \n\t   terms of the lattice parameter.\n"
           << "\t5. The cutoff value (fi_cut) for which orientation parameters are \n\t   assigned to each grain.\n"
           << "\t6. The lattice parameter in Angstroms.\n";
      return EXIT_SUCCESS;
    }

    if (result.count("quiet"))
    {
      warnings = false;
    }

    if (result.count("print-nearest-neighbors"))
    {
      print_nearest_neighbors = true;
    }

    if (result.count("file"))
    {
      if (result.count("type"))
      {
        string val = result["type"].as<string>();
        if ((val.compare("dump") != 0 || val.compare("dat") != 0))
        {
          cout << "Data file must be either a LAMMPS dump file (\"*.dump\") or a LAMMPS input data file (\"*.dat\")\n";
          exit(INPUT_FORMAT_ERROR);
        }
      }
      if (result.count("microrotation"))
      {
        cout << "Microrotation not implemented yet.\n";
      }
      filenames = parseInputFile(input_file);
      processData(filenames, result);
      cout << endl;
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath> // for acos, sin
#include "atom.h"

using namespace std;

#define PI 3.141592653589793

// Cutoff distances
#define UU_CUT 4.8496 // from Bai et al. Acta Materialia 85 (2015) 95-106: 0.866 * a0
#define UO_CUT 2.5 // actual value is 0.43301
#define OO_CUT 3.7 // actual value is 2.7625

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

// Conversion functions
double deg2rad(double x)
{
  return x * PI / 180.0;
}

double rad2deg(double x)
{
  return x * 180.0 / PI;
}

bool pairCmp(pair<int, double> &a, pair<int, double> &b)
{
  // round to the nearest 1000th
  double a1 = anInt(a.second * 1000) / 1000.0;
  double b1 = anInt(b.second * 1000) / 1000.0;
  return (a1 == b1);
}

bool pairSort(pair<int, double> &a, pair <int, double> &b)
{
  return (a.second < b.second);
}

int main(int argc, char** argv)
{
  string filename1, filename2, str; // filenames read from and written to, junk variable
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounds variables
  int N, n_type, n_atoms_read = 0; // number of atoms, atom types, number of atoms read
  vector <Atom> atoms; // all of the atoms from the file
  vector <int> temp_atoms, counts; // Just the indices of the nn_atoms to atom i, counts for the symmetry parameter
  vector <pair<int, vector<int> > > nn_atoms; // nearest neighbors
  vector <double> angles, symm; // a vector to hold the calculated angles, calculated symmetry parameters
  vector <pair <int, double> > symm_param, temp, unique_param; // vector holding the symmetry parameter for atom i
  int id, type, counter = 0, max1_index = 0, max2_index = 0; // id and type number of atom
  double charge, x, y, z, x1, y1, z1; // charge and position of atom
  double rxij, ryij, rzij, drij_sq, magnitude, theta;
  double uu_cut_sq = UU_CUT * UU_CUT, oo_cut_sq = OO_CUT * OO_CUT;
  double sintheta, sintheta_sq, total, max1 = 0, max2 = 0;

  if (argc == 2)
  {
    filename1 = argv[1];
    filename2 = filename1.substr(0,filename1.find(".dat")) + "_interface.dat";
  }
  else if (argc == 3)
  {
    filename1 = argv[1];
    filename2 = argv[2];
  }
  else
  {
    cout << "Please enter the data file to be read: ";
    cin  >> filename1;
    filename2 = filename1.substr(0,filename1.find(".dat")) + "_interface.dat";
  }

  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error opening file " << filename1 << endl;
    return 1;
  }

  ofstream fout(filename2.c_str());
  if (fout.fail())
  {
    cout << "Error opening file " << filename2 << endl;
    return 1;
  }

  getline(fin, str); // Gets the comment line;
  fin >> N >> str; // Gets the number of atoms
  fin >> n_type >> str >> str; // gets the number of atom types
  fin >> xlow >> xhigh >> str >> str;
  fin >> ylow >> yhigh >> str >> str;
  fin >> zlow >> zhigh >> str >> str;
  fin >> str; // Gets the Atoms line
  Lx = xhigh - xlow;
  Ly = yhigh - ylow;
  Lz = zhigh - zlow;

  cout << "Reading atoms...\n";
  while (fin >> id >> type >> charge >> x >> y >> z)
  {
    Atom temp(id, type, charge, x, y, z);
    atoms.push_back(temp);
    ++n_atoms_read;
  }

  if (n_atoms_read != N)
  {
    cout << "Error: number of atoms read does not match number of atoms in the simulation.\n"
         << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
    return 2;
  }

  // Now that we have the atoms safely stored, we can process them.  Restarting
  // the counter to make sure we calculate the correct number of angles!
  n_atoms_read = 0;
  cout << "Calculating nearest neighbors...\n";
  for (int i = 0; i < atoms.size(); ++i)
  {
    temp_atoms.clear();

    x1 = atoms[i].getX();
    y1 = atoms[i].getY();
    z1 = atoms[i].getZ();
    if (atoms[i].getType() == 1) // looking at U atoms
    {
      for (int j = 0; j < atoms.size(); ++j)
      {
        // Look at all of the U atoms that are within UU_CUT
        if (i == j)
          continue;
        if (atoms[j].getType() == 1)
        {
          // calculate the distances
          rxij = x1 - atoms[j].getX();
          ryij = y1 - atoms[j].getY();
          rzij = z1 - atoms[j].getZ();

          // Apply PBCs
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;

          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

          if (drij_sq < uu_cut_sq)
          {
            // Place them in a vector for later analysis
            // Note that this stores the atom ID, NOT the atom index!
            temp_atoms.push_back(atoms[j].getId());
          }
        }
      }
    nn_atoms.push_back(make_pair(atoms[i].getId(), temp_atoms));
    }
    ++n_atoms_read;
  }

  if (n_atoms_read != N)
  {
    cout << "Error: number of nearest neigbor sets does not match number of atoms in the simulation.\n"
         << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
    return 3;
  }

  // Now calculate the angles.  See the orientation parameter calculated as in
  // Bai et al. above, and in Zhang et al. Acta Materialia 53 (2005) 79-86
  n_atoms_read = 0;
  cout << "Calculating symmetry parameters...\n";

  for (int i = 0; i < nn_atoms.size(); ++i)
  {

    // for each set, calculate the vector distance between the atom pairs
    // The -1 is here because the IDs are (assumed to be) indexed from 1, not 0
    // We normalize everything to a single x-y plane, so the z distances are not
    // calculated.
    x = atoms[(nn_atoms[i].first) - 1].getX();
    y = atoms[(nn_atoms[i].first) - 1].getY();
    // The distances are projected onto the x-y plane, so the z component is always 0

    angles.clear(); // clear our angles vector for the next atom
    for (int j = 0; j < nn_atoms[i].second.size(); ++j)
    {
      // Calculate the angle for this pair
      rxij = atoms[(nn_atoms[i].second)[j] - 1].getX() - x;
      ryij = atoms[(nn_atoms[i].second)[j] - 1].getY() - y;

      // Apply PBCs
      rxij = rxij - anInt(rxij / Lx) * Lx;
      ryij = ryij - anInt(ryij / Ly) * Ly;

      magnitude = sqrt(rxij * rxij + ryij * ryij);

      theta = acos(rxij / magnitude); // gives the result in radians

      angles.push_back(theta);
      // Double check to see if we have already added in an angle for this atom
      // Add it into the map if we don't already have it
    }

    if (angles.size() == 0)
    {
      cout << "Error!  No angles calculated for atom " << nn_atoms[i].first << endl;
      return 4;
    }

    // Now we figure out the symmetry parameter of our atom
    total = 0.0;
    for (int j = 0; j < angles.size(); ++j)
    {
      sintheta = angles[j];
      sintheta_sq = sintheta * sintheta;
      total += (3 - 2 * sintheta_sq) * (3 - 2 * sintheta_sq) * sintheta_sq;
    }
    total /= angles.size();
    symm_param.push_back(make_pair(nn_atoms[i].first, total));
  }

  temp = symm_param;
  sort(symm_param.begin(), symm_param.end(), pairSort);
  vector <pair <int, double> >::iterator it;
  it = unique(symm_param.begin(), symm_param.end(), pairCmp);
  symm_param.resize(distance(symm_param.begin(), it));

  for (int i = 0; i < symm_param.size(); ++i)
  {
    symm_param[i].second = anInt(symm_param[i].second * 1000) / 1000.0;
  }
  unique_param = symm_param;
  symm_param = temp;

  for (int i = 0; i < symm_param.size(); ++i)
  {
    symm.push_back(anInt(symm_param[i].second * 1000) / 1000.0);
  }

  for (int i = 0; i < unique_param.size(); ++i)
  {
    counts.push_back(count(symm.begin(), symm.end(), unique_param[i].second));
    counter += counts[i];
    if (counts[i] > max1)
    {
      max2 = max1;
      max1 = unique_param[i].second;

      max2_index = max1_index;
      max1_index = i;
    }
    else if (counts[i] > max2)
    {
      max2 = unique_param[i].second;
      max2_index = i;
    }
  }

  cout << "The max number of counts occurs at " << max1_index << " with a value of " << max1 << " occurring " << counts[max1_index] << " times\n"
       << "The second max number of counts occurs at " << max2_index << " with a value of " << max2 << " occurring " << counts[max2_index] << " times\n";

  if (counter != nn_atoms.size())
  {
    cout << "Error: total number of counts is not equal to number of atoms calculated!\n"
         << "counter = " << counter << " != nn_atoms.size() = " << nn_atoms.size() << endl;
  }

  /*// Make sure we write the entire set of atoms
  n_atoms_read = 0;
  cout << "Writing results...\n";
  for (int i = 0; i < atoms.size(); ++i)
  {
    fout << atoms[i].getId() << " " << atoms[i].getType() << " "
         << atoms[i].getCharge() << " " << atoms[i].getX() << " "
         << atoms[i].getY() << " " << atoms[i].getZ() << " " << atoms[i].getMark() << endl;
    ++n_atoms_read;
  }

  if (n_atoms_read != N)
  {
    cout << "Error: number of atoms written does not match number of atoms in the simulation.\n"
         << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
    return 6;
  }*/


  fin.close();
  fout.close();
  return 0;
}

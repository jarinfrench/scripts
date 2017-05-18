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

bool pairGreater(pair<int, double> &a, pair <int, double> &b)
{
  return (a.second > b.second);
}

int main(int argc, char** argv)
{
  string filename1, filename2, str; // filenames read from and written to, junk variable
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounds variables
  int N, n_type, n_atoms_read = 0; // number of atoms, atom types, number of atoms read
  vector <Atom> atoms; // all of the atoms from the file
  vector <int> counts; // Just the indices of the nn_atoms to atom i, counts for the symmetry parameter
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

  // Open up the files for reading and writing.
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

  // Pull out the relevant information from the heading
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

  while (fin >> id >> type >> charge >> x >> y >> z)
  {
    if (type > n_type)
    {
      cout << "Error: unexpected atom type.\n"
           << "n_types = " << n_type << " < this atom's type = " << type << endl;
      return 2;
    }
    // Read the atoms line by line, and put them into a vector for analysis.
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

  for (int i = 0; i < atoms.size(); ++i)
  {
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

          // Calculate the distance
          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

          // If we are within the cutoff distance, store it atoms for later.
          if (drij_sq < uu_cut_sq)
          {
            /* Now calculate the angle.  See the orientation parameter
            ** calculated as in Bai et al. above, and in Zhang et al., Acta
            ** Materialia, 53 (2005), 79-86
            */
            // We only use the x and y for the magnitude because we project
            // the vector onto the x-y plane, effectively setting z = 0.
            magnitude = sqrt(rxij * rxij + ryij * ryij);
            theta = acos(rxij / magnitude); // gives the result in radians

            // Store the angles for analysis after this loop
            angles.push_back(theta);
          }
        }
      }

      if (angles.size() == 0)
      {
        cout << "Error!  No angles calculated for atom " << atoms[i].getId() << endl;
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
      symm_param.push_back(make_pair(atoms[i].getId(), total));
      symm.push_back(anInt(total * 1000) / 1000.0);
    }
  }

  // Find the unique values of the symmetry parameter
  temp = symm_param; // save the original data set
  sort(symm_param.begin(), symm_param.end(), pairSort); // sort it
  vector <pair <int, double> >::iterator it; // initialize an iterator
  it = unique(symm_param.begin(), symm_param.end(), pairCmp); // determine where the unique values stop
  symm_param.resize(distance(symm_param.begin(), it)); // resize the sorted vector

  for (int i = 0; i < symm_param.size(); ++i)
  {
    symm_param[i].second = anInt(symm_param[i].second * 1000) / 1000.0;
  }
  unique_param = symm_param;
  symm_param = temp;
  // FIXME: There shouldn't be this many unique symmetry parameters
  cout << "There are " << unique_param.size() << " unique symmetry parameters.\n";

  for (int i = 0; i < unique_param.size(); ++i)
  {
    counts.push_back(count(symm.begin(), symm.end(), unique_param[i].second));
    counter += counts[i];
    //cout << "unique_param[" << i << "].second = " << unique_param[i].second << " occurs " << counts[i] << " times.\n";
  }

  if (counter != symm_param.size())
  {
    cout << "Error in counting the elements\n";
    return 6;
  }

  max1 = *max_element(counts.begin(), counts.end());
  max1_index = distance(counts.begin(), max_element(counts.begin(), counts.end()));

  vector <int> count_temp = counts;
  nth_element(counts.begin(), counts.begin()+1, counts.end(), greater<int>());
  vector <int> sorted_counts = counts;
  counts = count_temp;
  max2 = *find(counts.begin(), counts.end(), sorted_counts[1]);
  max2_index = distance(counts.begin(), find(counts.begin(), counts.end(), sorted_counts[1]));

  cout << "The max number of counts occurs at " << max1_index << " with a value of " << unique_param[max1_index].second << " occurring " << counts[max1_index] << " times\n"
       << "The second max number of counts occurs at " << max2_index << " with a value of " << unique_param[max2_index].second << " occurring " << counts[max2_index] << " times\n";

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

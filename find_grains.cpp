#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath> // for acos, sin
#include "atom.h"

using namespace std;

#define PI 3.141592653589793

// Cutoff distances
#define UU_CUT 4.8496 // from Bai et al. Acta Materialia 85 (2015) 95-106: 0.866 * a0
//#define UO_CUT 2.5 // actual value is 0.43301
//#define OO_CUT 3.7 // actual value is 2.7625

// Calculate the rounded value of x
double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

// Comparison functions for sorting a vector of pairs
bool pairCmp(pair<int, double> &a, pair<int, double> &b)
{
  return (a.second == b.second);
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
  vector <int> counts; // Counts for the symmetry parameter
  vector <double> symm; // a vector to hold the calculated symmetry parameters
  vector <pair <int, double> > symm_param, temp, unique_param; // vector holding the symmetry parameter for atom i
  int id, type, max1_index, max2_index; // id and type number of atom
  unsigned int counter = 0;
  double charge, x, y, z; // charge and position of atom
  double rxij, ryij, rzij, drij_sq, magnitude, theta;
  double uu_cut_sq = UU_CUT * UU_CUT;//, oo_cut_sq = OO_CUT * OO_CUT;
  double sintheta, sintheta_sq, total;

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
  // This is for a LAMMPS input file.
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

  // Now that we have the atoms safely stored, we can process them.
  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    x = atoms[i].getX();
    y = atoms[i].getY();
    z = atoms[i].getZ();
    if (atoms[i].getType() == 1) // looking at U atoms
    {
      counter = 0; // reset the counter
      total = 0.0; // reset the total
      for (unsigned int j = 0; j < atoms.size(); ++j)
      {
        // Look at all of the U atoms that are within UU_CUT
        // Don't compare the atom to itself
        if (i == j)
          continue;
        if (atoms[j].getType() == 1)
        {
          // calculate the distances
          rxij = x - atoms[j].getX();
          ryij = y - atoms[j].getY();
          rzij = z - atoms[j].getZ();

          // Apply PBCs
          rxij = rxij - anInt(rxij / Lx) * Lx;
          ryij = ryij - anInt(ryij / Ly) * Ly;
          rzij = rzij - anInt(rzij / Lz) * Lz;

          // Calculate the magnitude of the distance
          drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

          // If we are within the cutoff distance, do the work
          if (drij_sq < uu_cut_sq)
          {
            /* Now calculate the angle.  See the orientation parameter
            ** calculated as in Bai et al. above, and in Zhang et al., Acta
            ** Materialia, 53 (2005), 79-86
            */
            magnitude = sqrt(drij_sq);

            theta = acos(rxij / magnitude); // gives the result in radians
            ++counter;

            // Calculate the symmetry parameter
            sintheta = sin(theta);
            sintheta_sq = sintheta * sintheta;
            total += (3 - 4 * sintheta_sq) * (3 - 4 * sintheta_sq) * sintheta_sq - sintheta_sq;
          }
        }
      }
      total /= counter;
      total = anInt(total * 100.0) / 100.0; // Round
      symm_param.push_back(make_pair(i, total)); // Store them for analysis
      symm.push_back(total);
    }
  }

  // Find the unique values of the symmetry parameter
  temp = symm_param; // save the original data set

  sort(symm_param.begin(), symm_param.end(), pairSort); // sort it

  vector <pair <int, double> >::iterator it; // initialize an iterator

  it = unique(symm_param.begin(), symm_param.end(), pairCmp); // determine where the unique values stop

  symm_param.resize(distance(symm_param.begin(), it)); // resize the sorted unique vector

  unique_param = symm_param; // Store the unique values

  symm_param = temp; // put the original values back

  // Count the occurrence of each unique value
  counter = 0; // Resets our counter
  for (unsigned int i = 0; i < unique_param.size(); ++i)
  {
    counts.push_back(count(symm.begin(), symm.end(), unique_param[i].second));
    counter += counts[i];
  }

  if (counter != symm_param.size())
  {
    cout << "Error in counting the elements\n"
         << "counter = " << counter << " != symm_param.size() = " << symm_param.size() << endl;
    return 6;
  }

  //max1 = *max_element(counts.begin(), counts.end());
  max1_index = distance(counts.begin(), max_element(counts.begin(), counts.end()));

  vector <int> count_temp = counts;
  nth_element(counts.begin(), counts.begin()+1, counts.end(), greater<int>());
  vector <int> sorted_counts = counts;
  counts = count_temp;
  //max2 = *find(counts.begin(), counts.end(), sorted_counts[1]);
  max2_index = distance(counts.begin(), find(counts.begin(), counts.end(), sorted_counts[1]));
  double max1 = unique_param[max1_index].second;
  double max2 = unique_param[max2_index].second;

  cout << "The highest occurring value is " << max1 << " and the second highest occurring value is " << max2 << endl;

  // Create a histogram of the number of counts for each value in unique_param
  vector <pair <double, int> > hist;
  ofstream fout2((filename1.substr(0,filename1.find(".dat")) + "_histogram.csv").c_str());
  if (fout2.fail())
  {
    cout << "Error opening file " << filename1.substr(0,filename1.find(".dat")) + "_histogram.csv" << endl;
    cout << "Line 280\n";
    return 1;
  }
  for (unsigned int i = 0; i < unique_param.size(); ++i)
  {
    fout2 << unique_param[i].second << "," << counts[i] << endl;
  }
  fout2.close();

  // Calculate the cutoff distance - we are specifying that halfway between the
  // two most occurring values is the cutoff.
  double cutoff = (unique_param[max1_index].second + unique_param[max2_index].second) / 2.0;

  // Now we mark the atoms based its counts
  for (unsigned int i = 0; i < symm_param.size(); ++i)
  {
    if (symm_param[i].second <= cutoff)
    {
      atoms[symm_param[i].first].setMark(1);
    }
    else
    {
      atoms[symm_param[i].first].setMark(2);
    }
  }

  // Make sure we write the entire set of atoms
  // This writes things in a tecplot-readable format.
  n_atoms_read = 0;
  for (unsigned int i = 0; i < symm_param.size(); ++i)
  {
    fout << atoms[symm_param[i].first].getId() << " "
         << atoms[symm_param[i].first].getType() << " "
         << atoms[symm_param[i].first].getCharge() << " "
         << atoms[symm_param[i].first].getX() << " "
         << atoms[symm_param[i].first].getY() << " "
         << atoms[symm_param[i].first].getZ() << " "
         << atoms[symm_param[i].first].getMark() << " "
         << symm_param[i].second << endl;
    ++n_atoms_read;
  }

  for (unsigned int i = 0; i < atoms.size(); ++i)
  {
    if (atoms[i].getType() == 1)
      continue;
    else // O atoms get assigned a value of 1 by default, which the U atoms will never get.
    {
      fout << atoms[i].getId() << " "
           << atoms[i].getType() << " "
           << atoms[i].getCharge() << " "
           << atoms[i].getX() << " "
           << atoms[i].getY() << " "
           << atoms[i].getZ() << " "
           << atoms[i].getMark() << " "
           << 1 << endl;
      ++n_atoms_read;
    }
  }

  if (n_atoms_read != N)
  {
    cout << "Error: number of atoms written does not match number of atoms in the simulation.\n"
         << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
    return 6;
  }


  fin.close();
  fout.close();
  return 0;
}

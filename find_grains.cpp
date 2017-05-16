#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <cmath> // for acos
#include "atom.h"

using namespace std;

#define PI 3.141592653589793

// Cutoff distances
#define UU_CUT 5.5 // actual value is 3.856
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

// Counting struct for couting number of occurences of a value in a map.
struct Compare
{
  double num;
  Compare(const double& num) : num(num) {}
};

bool operator==(const std::pair<int, double>&p, const Compare& c)
{
  return c.num == p.second;
}

bool operator==(const Compare& c, const std::pair<int, double>&p)
{
  return c.num == p.second;
}

int main(int argc, char** argv)
{
  string filename1, filename2, str; // filenames read from and written to, junk variable
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, Lx, Ly, Lz; // bounds variables
  int N, n_type, n_atoms_read = 0; // number of atoms, atom types, number of atoms read
  vector <Atom> atoms; // all of the atoms from the file
  vector <int> temp_atoms; // Just the indices of the nn_atoms to atom i
  vector <pair<int, vector<int> > > nn_atoms; // nearest neighbors
  vector <int> angles_count; // holds the number of occurences of each angle (given in order of the set unique_angles)
  vector <int> sorted_angles; // sorted angles
  map <int, double> angles; // map of the angles for the atom ID
  set <double> unique_angles;
  int id, type, max1, max2, max1_index, max2_index, counter; // id and type number of atom, counter variable
  double charge, x, y, z, x1, y1, z1; // charge and position of atom
  double rxij, ryij, rzij, drij_sq, magnitude, theta;
  double uu_cut_sq = UU_CUT * UU_CUT, oo_cut_sq = OO_CUT * OO_CUT;
  double max1_angle, max2_angle;

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
  cout << "[" << flush; // Progress bar
  for (int i = 0; i < atoms.size(); ++i)
  {
    if (i % (N / 100) == 0 )
    {
      cout << "." << flush;
    }
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
    }
    else if (atoms[i].getType() == 2) // looking at O atoms
    {
      for (int j = 0; j < atoms.size(); ++j)
      {
        // Look at all of the O atoms that are within OO_CUT
        if (i == j)
          continue;
        if (atoms[j].getType() == 2)
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

          if (drij_sq < oo_cut_sq)
          {
            // Place them in a vector for later analysis
            // Note that this stores the atom ID, NOT the atom index!
            temp_atoms.push_back(atoms[j].getId());
          }
        }
      }
    }
    if (temp_atoms.size() < 1)
    {
      cout << "Error determining nearest neighbors for atom " << atoms[i].getId() << endl;
    }
    nn_atoms.push_back(make_pair(atoms[i].getId(), temp_atoms));
    ++n_atoms_read;
  }

  if (n_atoms_read != N)
  {
    cout << "Error: number of nearest neigbor sets does not match number of atoms in the simulation.\n"
         << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
    return 3;
  }

  cout << ".] COMPLETE\n";

  // Now calculate the angles.  We only use the particle that is the farthest
  // away in the x direction.  In the ideal case, this means that we are just shy
  // of the second nearest neighbor for each atom (in a perfect lattice).
  // Reset the counter again to make sure we get the right number of angles...
  n_atoms_read = 0;
  cout << "Calculating angles...\n";
  cout << "[" << flush;
  for (int i = 0; i < nn_atoms.size(); ++i)
  {
    cout << "Index " << i << " is for atom id " << nn_atoms[i].first << endl;
  }
  for (int i = 0; i < nn_atoms.size(); ++i)
  {
    if (i % (N / 100) == 0 )
    {
      cout << "." << flush;
    }
    // for each set, calculate the vector distance between the atom pairs
    // The -1 is here because the IDs are (assumed to be) indexed from 1, not 0
    // We normalize everything to a single x-y plane, so the z distances are not
    // calculated.
    x = atoms[(nn_atoms[i].first) - 1].getX();
    y = atoms[(nn_atoms[i].first) - 1].getY();

    rxij = 0;
    for (int j = 0; j < nn_atoms[i].second.size(); ++j)
    {
      // Calculate the angle for this pair
      x1 = atoms[(nn_atoms[i].second)[j] - 1].getX();
      y1 = atoms[(nn_atoms[i].second)[j] - 1].getY();

      // if our current distance is greater than (or equal to, just to double
      // check) the last saved one, (re)do the calculation
      if (((x1 - x) - anInt((x1 - x) / Lx) * Lx) >= rxij)
      {
        rxij = x1 - x;
        ryij = y1 - y;

        // Apply PBCs
        rxij = rxij - anInt(rxij / Lx) * Lx;
        ryij = ryij - anInt(ryij / Ly) * Ly;

        magnitude = sqrt(rxij * rxij + ryij * ryij);
        theta = acos(rxij / magnitude); // gives the result in radians

        // Double check to see if we have already added in an angle for this atom
        // Add it into the map if we don't already have it
        if (angles.count(nn_atoms[i].first) == 0)
        {
          ++n_atoms_read;
          cout << "Initial angle written for atom " << nn_atoms[i].first << endl;
        }
        angles[nn_atoms[i].first] = anInt(abs(rad2deg(theta)));
      }
    }
  }
  cout << ".] COMPLETE\n";

  if (n_atoms_read != N)
  {
    cout << "Error: number of calculated angles does not match number of atoms in the simulation.\n"
         << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
    return 4;
  }

  // Determine the number of unique angles
  for (map <int, double>::iterator it = angles.begin();
       it != angles.end();
       ++it)
  {
    unique_angles.insert(it->second);
  }

  // count the occurence of each angle, make sure we have the same number of
  // angles that we have atoms!
  n_atoms_read = 0;
  set <double>::iterator set_it = unique_angles.begin();
  for (int i = 0; i < unique_angles.size(); ++i)
  {
    counter = count(angles.begin(), angles.end(), Compare(*set_it));
    n_atoms_read += counter;
    angles_count.push_back(counter);
    cout << count(angles.begin(), angles.end(), Compare(*set_it)) << " have angle " << setprecision(10) << *set_it << endl;
    ++set_it;
  }

  if (n_atoms_read != N)
  {
    cout << "Error: number of angles does not match number of atoms in the simulation.\n"
         << "N = " << N << " != n_atoms_read = " << n_atoms_read << endl;
    return 5;
  }
  //sort(angles_count.rbegin(), angles_count.rend());
  //max1 = angles_count[0];
  //max2 = angles_count[1];

  cout << "The index with the highest angle is " << distance(angles_count.begin(), max_element(angles_count.begin(), angles_count.end()))
       << " which occurs " << *max_element(angles_count.begin(), angles_count.end()) << " times.\n";
  max1 = *max_element(angles_count.begin(), angles_count.end());
  max1_index = distance(angles_count.begin(), max_element(angles_count.begin(), angles_count.end()));
  set_it = unique_angles.begin();
  advance(set_it, max1_index);
  max1_angle = *set_it;

  vector <int> temp = angles_count;

  // This "sorts" the array such that the second highest value is in the second position.
  nth_element(angles_count.begin(), angles_count.begin()+1, angles_count.end(), greater<int>());
  sorted_angles = angles_count;
  angles_count = temp;
  cout << "The index with the second highest angle is " << distance(angles_count.begin(), find(angles_count.begin(), angles_count.end(), sorted_angles[1]))
       << " which occurs " << *find(angles_count.begin(), angles_count.end(), sorted_angles[1]) << " times.\n";
  max2 = *find(angles_count.begin(), angles_count.end(), sorted_angles[1]);
  max2_index = distance(angles_count.begin(), find(angles_count.begin(), angles_count.end(), sorted_angles[1]));
  set_it = unique_angles.begin();
  advance(set_it, max2_index);
  max2_angle = *set_it;

  //cout << max1 << "\t" << max2 << endl;
  //cout << max1_angle << "\t" << max2_angle << endl;
  for (int i = 0; i < nn_atoms.size(); ++i)
  {
    if (angles[nn_atoms[i].first] == max1_angle)
    {
      atoms[nn_atoms[i].first - 1].setMark(1);
    }
    else if (angles[nn_atoms[i].first] == max2_angle)
    {
      atoms[nn_atoms[i].first - 1].setMark(2);
    }
    else
    {
      atoms[nn_atoms[i].first - 1].setMark(3);
    }
  }

  // Make sure we write the entire set of atoms
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
  }


  fin.close();
  fout.close();
  return 0;
}

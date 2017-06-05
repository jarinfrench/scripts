/******************************************************************************
* This script converts a file in the format .xyz to the .dat format.  Note that
* if copying data from GBStudio, the output format MUST be xmol, otherwise there
* will be extra information that this script will not handle properly.
******************************************************************************/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
  string filename1, filename2, filename3, str;
  int n_atoms, N, type;
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double x, y, z, charge, temp = 1.0;
  char atom_type;

  if (argc == 1)
  {
    cout << "Please enter the file to be converted: ";
    cin >> filename1;
  }
  else
  {
    filename1 = argv[1];
  }

  filename2 = "temp.dat"; // temporary filename until we determine what it needs to be
  filename3 = filename1.substr(0, filename1.find(".")) + ".dat";

  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Unable to open file " << filename1 << endl;
    return 1;
  }

  ofstream fout(filename2.c_str());
  if (fout.fail())
  {
    cout << "Unable to open file " << filename2 << endl;
    return 1;
  }

  fout << "This bulk UO2 coordinates format: [ID type charge x y z]\n\n";

  fin >> N; // Number of atoms
  fout << N << "  atoms\n";

  fout << "2   atom types\n"; // This is assuming UO2

  fin >> xhigh >> yhigh >> zhigh >> str >> str >> str; // Boundaries
  xlow = 0.0;
  ylow = 0.0;
  zlow = 0.0;
  fout << fixed;
  fout << setprecision(6) << xlow << " " << xhigh << " xlo xhi\n";
  fout << ylow << " " << yhigh << " ylo yhi\n";
  fout << zlow << " " <<  zhigh << " zlo zhi\n";

  fout << "\nAtoms\n\n";

  n_atoms = 0;
  while (fin >> atom_type >> x >> y >> z) // Read in the data
  {
    if (atom_type == 'U') // Given a U atom, write the correct type and charge
    {
      type = 1;
      charge = 2.4;
    }
    else if (atom_type == 'O') // Same with the O atom
    {
      type = 2;
      charge = -1.2;
    }
    else // We are only dealing with U and O, so if we have something else, there's a problem
    {
      cout << "Unrecognized atom type.\n";
      return 2;
    }

    // This just double checks that our bounds are set correctly.  If not, they
    // will be correct later.
    if (x < xlow)
    {
      temp = -xhigh;
      xhigh = 0.0;
      xlow = temp;
    }
    if (y < ylow)
    {
      temp = -yhigh;
      yhigh = 0.0;
      ylow = temp;
    }
    if (z < zlow)
    {
      temp = -zhigh;
      zhigh = 0.0;
      zlow = temp;
    }
    fout << setprecision(0) << ++n_atoms << " " << type << " "
         << setprecision(1) << charge << " "
         << setprecision(6) << x << " " << y << " " << z << endl;

  }

  // Close the file streams
  fin.close();
  fout.close();
  if (n_atoms != N)
  {
    cout << "Error counting atoms: n_atoms = " << n_atoms << " != N = " << N << endl;
    return 3;
  }

  // Now fix the bounds if necessary
  if (temp < 0) // This works because we initialize temp as 1.0, and throughout this code it is ONLY negative
  {
    ifstream fin(filename2.c_str());
    if (fin.fail())
    {
      cout << "Unable to open file " << filename2 << endl;
      return 1;
    }
    ofstream fout(filename3.c_str());
    if (fout.fail())
    {
      cout << "Unable to open file " << filename3 << endl;
      return 1;
    }
    fout << fixed;

    // We need to re-write everything to the new file.
    getline(fin, str); // Comment line
    fout << str << "\n\n";
    getline(fin,str); // blank line
    getline(fin, str); // N atoms line
    fout << str << "\n";
    getline(fin, str); // n atom types line
    fout << str << "\n";
    fin >> str >> str >> str >> str; // x bounds
    fout << setprecision(6) << xlow << " " << xhigh << " xlo xhi\n";
    fin >> str >> str >> str >> str; // y bounds
    fout << ylow << " " << yhigh << " ylo yhi\n";
    fin >> str >> str >> str >> str; // z bounds
    fout << zlow << " " << zhigh << " zlo zhi";
    fin.clear();
    while (getline(fin,str))
    {
      fout << str << "\n";
    }
    if (remove(filename2.c_str()) != 0)
    {
      cout << "Error deleting temp file \"temp.dat\"\n";
    }
    // Close the file streams
    fin.close();
    fout.close();
  }
  // If the bounds haven't changed, just rename the file (a LOT easier!)
  else // rename the temp file.
  {
    rename(filename2.c_str(), filename3.c_str());
  }

  return 0;
}

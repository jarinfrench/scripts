#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
  string filename1, filename2, str;
  int n_atoms, n_atom_types, N, type;
  double xlow, xhigh, ylow, yhigh, zlow, zhigh;
  double x, y, z, charge;
  char atom_type;

  if (argc == 1)
  {
    cout << "Please enter the filename to be converted: ";
    cin >> filename1;
  }
  else
  {
    filename1 = argv[1];
  }

  filename2 = filename1.substr(0, filename1.find(".")) + ".dat";

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

  getline(fin, str); // initial comment line
  fout << str << "\n\n";

  fin >> N >> str; // Number of atoms
  fout << N << "  atoms\n";

  fin >> n_atom_types >> str >> str; // Number of atom types
  fout << n_atom_types << "  atom types\n";

  fin >> xlow >> xhigh >> str >> str; // Boundaries
  fout << fixed;
  fout << setprecision(6) << xlow << " " << xhigh << " xlo xhi\n";

  fin >> ylow >> yhigh >> str >> str;
  fout << ylow << " " << yhigh << " ylo yhi\n";

  fin >> zlow >> zhigh >> str >> str;
  fout << zlow << " " <<  zhigh << " zlo zhi\n";

  fin >> str; // Atoms line
  fout << "\nAtoms\n\n";

  n_atoms = 0;
  while (fin >> atom_type >> x >> y >> z)
  {
    if (atom_type == 'U')
    {
      type = 1;
      charge = 2.4;
    }
    else if (atom_type = 'O')
    {
      type = 2;
      charge = -1.2;
    }
    else
    {
      cout << "Unrecognized atom type.\n";
      return 2;
    }

    fout << setprecision(0) << ++n_atoms << " " << type << " "
         << setprecision(1) << charge << " "
         << setprecision(6) << x << " " << y << " " << z << endl;

  }

  fin.close();
  fout.close();
  if (n_atoms != N)
  {
    cout << "Error counting atoms: n_atoms = " << n_atoms << " != N = " << N << endl;
    return 3;
  }

  return 0;

}

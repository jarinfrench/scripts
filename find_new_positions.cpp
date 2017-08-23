#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include "atom.h"

using namespace std;

int main(int argc, char **argv)
{
  // External values
  string file1, file2, file3, str;
  vector <Atom> reference, compared; // List of atoms to start, and compare to
  int id, type, N;
  double charge, x, y, z, xlo, xhi, ylo, yhi, zlo, zhi, Lx, Ly, Lz;

  if (argc != 3)
  {
    cout << "Please enter the reference data: ";
    cin  >> file1;
    cout << "Please enter the filename of the data to process: ";
    cin  >> file2;
  }
  else
  {
    file1 = argv[1];
    file2 = argv[2];
  }

  ifstream fin(file1.c_str());
  if (fin.fail())
  {
    cout << "Error opening file: " << file1 << endl;
    return 1;
  }

  ifstream fin2(file2.c_str());
  if (fin2.fail())
  {
    cout << "Error opening file: " << file2 << endl;
    return 1;
  }

  file3 = file1.substr(0,file1.find(".dump"))+"_tracked.dat";
  ofstream fout(file3.c_str());
  if (fout.fail())
  {
    cout << "Error opening file: " << file3 << endl;
    return 1;
  }

  fin >> str >> str; // ITEM: TIMESTEP
  fin >> str; // <timestep>
  fin.ignore();
  getline(fin, str); // ITEM: NUMBER OF ATOMS
  fin >> N;
  fin.ignore();
  getline(fin, str); // ITEM: BOX BOUNDS
  fin >> xlo >> xhi;
  fin >> ylo >> yhi;
  fin >> zlo >> zhi;
  fin.ignore();
  getline(fin, str); // ITEM: ATOMS
  Lx = xhi - xlo;
  Ly = yhi - ylo;
  Lz = zhi - zlo;

  for (int i = 0; i < 9; ++i)
  {
    getline(fin2, str);
    // Ignore the first 9 lines of the second file
  }

  reference.resize(N, Atom());
  compared.resize(N, Atom());

  while (fin >> id >> type >> charge >> x >> y >> z)
  {
    reference[id - 1] = Atom(id, type, charge, x, y, z);
    if (x < 11 * Lx / 20 && x > 9 * Lx / 20)
    {
      reference[id - 1].setMark(1);
    }
    fin2 >> id >> type >> charge >> x >> y >> z;
    compared[id - 1] = Atom(id, type, charge, x, y, z);
  }

  for (unsigned int i = 0; i < reference.size(); ++i)
  {
    if (reference[i].getMark() == 1)
    {
      compared[i].setMark(1);
    }

    fout << compared[i].getId() << " " << compared[i].getType() << " "
         << compared[i].getCharge() << " " << compared[i].getX() << " "
         << compared[i].getY() << " " << compared[i].getZ() << " "
         << compared[i].getMark() << endl;
  }


  return 0;
}

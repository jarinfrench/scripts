#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char** argv)
{
  string filename1, filename2, str; // filenames to be converted, written to, string variables

  if (argc == 1)
  {
    // prompt for the filename to be converted.
    cout << "Please enter the filename of the CSV file to be converted: ";
    cin  >> filename1;
  }
  else
  {
    filename1 = argv[1];
  }

  if (filename1.substr(filename1.length() - 3) != "csv")
  {
    cout << "Error: please enter a file of type CSV (ending in .csv).\n";
    return 1;
  }

  // The new file will have the same name, just a different file extension.
  filename2 = filename1.substr(0, filename1.find(".")) + ".tec";

  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error: unable to open file " << filename1 << endl;
    return 2;
  }

  ofstream fout(filename2.c_str());
  if (fout.fail())
  {
    cout << "Error: unable to open file " << filename2 << endl;
    return 2;
  }

  // The first four lines are helpful for python processing, but not here, so
  // skip them

  getline(fin, str); // LAMMPS version
  getline(fin, str); // Box size
  getline(fin, str); // Number of atoms
  getline(fin, str); // Unit style

  // Now we read the rest of the file
  while (getline(fin, str, ','))
  {
    fout << str << " ";
  }

  return 0;
}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char **argv)
{
  string filename1, filename2, str, str2, str3; // filenames and line variablea
  int N; // number of atoms
  vector<long double> energies; // all of the energies we've read
  long double en1, en2, en3; // current energies
  double theta;

  // Get the filename
  if (argc != 2)
  {
    cout << "Please enter the filename that contains the energy from the simulation: ";
    cin  >> filename1;
  }
  else
    filename1 = argv[1];

  filename2 = "100Twist_total_energy.txt"; // the file we will write to.
  str2 = "  Energy initial, next-to-last, final =";
  str3 = "ERROR: Lost atoms: ";

  stringstream th; // this defaults to theta being 0 if it can't find a value
  th << filename1.substr(filename1.find("_") + 1, filename1.find("degree") - filename1.find("_") - 1);
  th >> theta;

  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error! Cannot open file " << filename1 << endl;
    return -1;
  }
  getline(fin, str); // LAMMPS (17 Nov 2016)
  getline(fin, str); // Reading data file ...
  getline(fin, str); //   orthogonal box = (0 0 0) to (272.651 272.651 27.2651)
  getline(fin, str); //   8 by 4 by 1 MPI processor grid
  getline(fin, str); //   reading atoms ...
  fin >> N >> str;
  fin.ignore(); // clears the newline character

  while (getline(fin,str))
  {
    if (fin.fail())
    {
      cout << "Read error!\n";
    }
    if (str == str2)
    {
      fin >> en1 >> en2 >> en3;
      energies.push_back(en1);
      energies.push_back(en2);
      energies.push_back(en3);
      fin.ignore(); // clear the newline character
    }
    if (str.substr(0,19) == str3) // Simulation ended due to lost atoms
    {
      energies.clear(); // clear all of the energies
      cout << "Simulation error: lost atoms.\n";
      break;
    }
  }

  if (energies.size() == 0)
  {
    cout << "Did not find any energies to return.\n\n";
    fin.close();
    return 0;
  }

  ofstream fout(filename2.c_str(), ofstream::app);
  if (fout.fail())
  {
    cout << "Error! Cannot open file " << filename2 << endl;
    return -1;
  }

  fout << theta << " " << setprecision(15)
       << *min_element(energies.begin(), energies.end()) << " "
       << N << endl;


  fin.close();
  fout.close();

  return 0;
}

/******************************************************************************
* This script takes as input the LAMMPS output file that contains the energies
* from the simulations, and the destination filename.  If the file already exists,
* the data is simply appended to it.
******************************************************************************/
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
  string filename1, filename2, str, str2, str3; // filenames and line variables
  int N; // number of atoms
  vector<long double> energies; // all of the energies we've read
  long double en1, en2, en3; // current energies
  double theta;

  // Get the filename
  if (argc != 3)
  {
    // Prompt for the filenames if we don't already have them
    cout << "Please enter the filename that contains the energy from the simulation: ";
    cin  >> filename1;

    cout << "Please enter the filename to write to: ";
    cin  >> filename2;
  }
  else
  {
    filename1 = argv[1];
    filename2 = argv[2];
  }

  // These are the strings that we are looking for.
  str2 = "  Energy initial, next-to-last, final = ";
  str3 = "ERROR: Lost atoms: ";

  stringstream th; // this defaults to theta being 0 if it can't find a value

  // This line is rather arbitrary, but looks for the first occurrence of the "_"
  // character after the 9th character.  Given a file of the format
  // minimize_100_*degree.dat means that the first occurrence of the "_"
  // character is the one right before the * - specifying the angle.  Note that
  // if the phrase "no_GB" is found instead of "***degree", the theta value will
  // automatically set the angle to 0, which is what we want.
  th << filename1.substr(filename1.find("_", 9) + 1, filename1.find("degree") - filename1.find("_") - 1);
  th >> theta;

  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error! Cannot open file \"" << filename1 << "\"" << endl;
    return 1;
  }
  // The next few lines is just to ignore the first few lines in the output file.
  getline(fin, str); // LAMMPS (17 Nov 2016)
  getline(fin, str); // Reading data file ...
  getline(fin, str); //   orthogonal box = (0 0 0) to (272.651 272.651 27.2651)
  getline(fin, str); //   8 by 4 by 1 MPI processor grid
  getline(fin, str); //   reading atoms ...
  fin >> N >> str; // Having N is useful!  We store it here.
  fin.ignore(); // clears the newline character

  while (getline(fin,str)) // Check every line
  {
    if (fin.fail())
    {
      cout << "Read error!\n";
    }
    if (str == str2) // If we found the line that matches our target:
    {
      // Extract the energies, store them in a vector, and move ahead.
      fin >> en1 >> en2 >> en3;
      energies.push_back(en1);
      energies.push_back(en2);
      energies.push_back(en3);
      fin.ignore(); // clear the newline character
    }
    // This checks to see if the simulation ended due to lost atoms.  If so, we
    // aren't going to find any energies, so we just end the loop here.
    if (str.substr(0,19) == str3) // Simulation ended due to lost atoms
    {
      energies.clear(); // clear all of the energies
      cout << "Simulation error: lost atoms.\n";
      break;
    }
  }

  // If we didn't find any energies, no point in doing anything else.
  if (energies.size() == 0)
  {
    cout << "Did not find any energies to return in file " << filename1 << "\n\n";
    fin.close();
    return 0;
  }

  // Write the angle, the SMALLEST value in the vector, and the number of atoms.
  // Generally speaking, the minimum value should be the last value written to
  // the vector, but we use *min_element() just in case.  Note that if the
  // simulation annealed the structure, there is additional, manual analysis
  // that must be done.
  ofstream fout(filename2.c_str(), ofstream::app);
  if (fout.fail())
  {
    cout << "Error! Cannot open file \"" << filename2 << "\"" << endl;
    return 1;
  }

  fout << theta << " " << setprecision(15)
       << *min_element(energies.begin(), energies.end()) << " "  // min_element is a pointer, so we dereference with *
       << N << endl;

  // Close the file streams.
  fin.close();
  fout.close();

  return 0;
}

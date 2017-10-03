/******************************************************************************
* This program parses output from LAMMPS minimized structures and formats it
* so LAMMPS can run a simulation with the results.
*******************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>

using namespace std;

int main(int argc, char **argv)
{
  string filename1, filename2, str, element; // Filename of read and written file.  Also a junk variable to hold stuff
  int N,  n_total = 0; // Number of atoms, number read
  double xlow, xhigh, ylow, yhigh, zlow, zhigh; // Box bounds
  int id, type; // Atom id number and type number
  double charge, x, y, z; // atom charge, and x, y and z positions
  double r_grain;

  if (argc != 4)
  {
    cout << "Please enter the filename to be parsed: ";
    cin  >> filename1;

    cout << "Please enter the radius of the rotated grain: ";
    cin  >> r_grain;

    cout << "Please enter the element: ";
    cin  >> element;
  }
  else
  {
    filename1 = argv[1]; // First value of argv is the script name
    r_grain = strtod(argv[2], NULL); // Second parameter specifies the endptr - ignored if NULL
    element = argv[3];
  }

  transform(element.begin(), element.end(), element.begin(), ::toupper);
  if (element.compare("UO2") != 0)
  {
    transform(element.begin(), element.end(), element.begin(), ::tolower);
    element[0] = toupper(element[0]);
  }
  ostringstream fn2;
  // The +8 moves the starting point to where the filename indicates the degree
  // The format of the filename SHOULD be: dump3.pos.<axis>.###degree.<step #>.dat
  // The -14 comes from calculating the number of characters after the first
  // number indicating the angle.  A hard coded value will not always get the
  // right amount of chars.
  fn2 << element << "_minimized_" << filename1.substr(filename1.find("pos") + 8,
                         filename1.find(".dat") - 14 - filename1.find("pos"))
      << "_r" << r_grain << "A.dat";
  filename2 = fn2.str();

  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error opening file " << filename1 << endl;
    return 1;
  }

  ofstream fout(filename2.c_str());
  fout << fixed; // Always use the decimal!
  if (fout.fail())
  {
    cout << "Error opening file " << filename2 << endl;
    return 1;
  }

  // Parse the data so that we can write the relevant data to filename2
  getline(fin, str); //Gets the first line (TIMESTEP)
  fout << "These " << element << " coordinates have been annealed: ";
  if (element.compare("UO2") == 0)
  {
    fout << "[ID type charge x y z]\n\n";
  }
  else
  {
    fout << "[ID type x y z]\n\n";
  }

  fin >> str; // This line should be the timestep ended on
  fin.ignore(); // This is important!  There is an extra newline character in the stream that will be the first thing read by getline!
  getline(fin, str); // This should get the "ITEM: NUMBER OF ATOMS" line
  fin >> N;
  fout << N << "  atoms\n"; // Store this value in the new file
  if (element.compare("UO2") == 0)
  {
    fout << "2   atom types\n"; // WARNING: This is assuming we are working with UO2 only!
  }
  else
  {
    fout << "1   atom types\n";
  }

  fin.ignore();
  getline(fin, str); // Gets the line "ITEM: BOX BOUNDS pp pp pp"
  fin >> xlow >> xhigh;
  fin >> ylow >> yhigh;
  fin >> zlow >> zhigh;

  fout.precision(6);
  fout << xlow << "\t" << xhigh << "\txlo xhi\n";
  fout << ylow << "\t" << yhigh << "\tylo yhi\n";
  fout << zlow << "\t" << zhigh << "\tzlo zhi\n";

  fin.ignore();
  getline(fin, str); // Gets the line "ITEM: ATOMS id type q x y z"
  fout << "\nAtoms\n\n"; // Now we can start reading the atoms

  while (getline(fin,str))
  {
    if (fin.fail())
    {
      cout << "Read error.\n";
      break;
    }

    stringstream ss(str);
    if (element.compare("UO2") == 0)
    {
      ss >> id >> type >> charge >> x >> y >> z;
    }
    else
    {
      ss >> id >> type >> x >> y >> z;
    }

    ++n_total;

    fout.precision(0);
    fout << id << " " << type << " ";
    if (element.compare("UO2") == 0)
    {
      fout.precision(1);
      fout << charge << " ";
    }
    fout.precision(6);
    fout << x << " " << y  << " " << z << endl;
  }

  if (n_total != N)
  {
    cout << "Error: " << n_total << " atoms were read out of " << N << " atoms in the file.\n";
    return 2;
  }

  return 0;
}

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>

#include "atom.h"
#include "position.h"

using namespace std;

vector <double> compareFiles(vector <Atom> & reference_atoms, string current_file)
{
  bool dump;
  string str;
  unsigned int N2, n_atoms_read2 = 0;
  vector <Atom> current_atoms;
  vector <double> msd (4, 0.0);
  vector <double> ref_CoM(3,0.0), curr_CoM(3,0.0), diff_CoM(3,0.0);
  ifstream fin2(current_file.c_str());

  for (unsigned int i = 0; i < reference_atoms.size(); ++i)
  {
    ref_CoM[0] += reference_atoms[i].getWrapped()[0];
    ref_CoM[1] += reference_atoms[i].getWrapped()[1];
    ref_CoM[2] += reference_atoms[i].getWrapped()[2];
  }
  ref_CoM[0] /= reference_atoms.size();
  ref_CoM[1] /= reference_atoms.size();
  ref_CoM[2] /= reference_atoms.size();

  if (fin2.fail())
  {
    cerr << "Error opening file: \"" << current_file << "\"\n";
    exit(1);
  }



  if (current_file.find("dump") == string::npos) // if "dump" is not in the filename
  {
    dump = false; // it isn't a dump file, assume the file is an input file.
  }
  else
  {
    dump = true; // if it is in the filename, assume that it's a LAMMPS dump file
  }

  if (dump)
  {
    // This is for a LAMMPS dump file
    getline(fin2, str); // Gets ITEM: TIMESTEP
    getline(fin2, str); // Gets the timestep number

    getline(fin2, str); // Gets ITEM: NUMBER OF ATOMS
    fin2 >> N2;
    fin2.ignore();
    getline(fin2, str); //get ITEM: BOX BOUNDS
    getline(fin2, str); // get the x bounds
    getline(fin2, str); // get the y bounds
    getline(fin2, str); // get the z bounds
    getline(fin2, str); // Gets ITEM: ATOMS <data types>
  }
  else
  {
    // This is for a LAMMPS input file
    getline(fin2, str); // gets comment line
    fin2 >> N2 >> str; // number of atoms
    fin2.ignore();
    getline(fin2, str); // n atom types
    getline(fin2, str); // x bounds
    getline(fin2, str); // y bounds
    getline(fin2, str); // z bounds
    fin2 >> str;
    fin2.ignore();
  }

  if (reference_atoms.size() != N2)
  {
    cerr << "Error: number of atoms between data files do not match. N1 = " << reference_atoms.size() << " != N2 = " << N2 << "\n";
    exit(2);
  }

  current_atoms.resize(N2, Atom());

  while (getline(fin2, str))
  {
    int atom_id, atom_type;
    double charge, x, y, z;
    stringstream ss(str);
    if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z))
    {
      // If the line does NOT have 6 entries, it must have 5, so
      ss >> atom_id >> atom_type >> x >> y >> z;
      charge = 0.0;
    }
    Position p (x,y,z);
    current_atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, p);
    ++n_atoms_read2;
  }

  for (unsigned int i = 0; i < current_atoms.size(); ++i)
  {
    curr_CoM[0] += current_atoms[i].getWrapped()[0];
    curr_CoM[1] += current_atoms[i].getWrapped()[1];
    curr_CoM[2] += current_atoms[i].getWrapped()[2];
  }
  curr_CoM[0] /= current_atoms.size();
  curr_CoM[1] /= current_atoms.size();
  curr_CoM[2] /= current_atoms.size();
  diff_CoM[0] = curr_CoM[0] - ref_CoM[0];
  diff_CoM[1] = curr_CoM[1] - ref_CoM[1];
  diff_CoM[2] = curr_CoM[2] - ref_CoM[2];

  fin2.close();

  if (n_atoms_read2 != N2)
  {
    cerr << "Error reading current atoms. n_atoms_read = " << n_atoms_read2 << " != N = " << N2 << endl;
    exit(3);
  }

  // Now calculate the MSD values
  for (unsigned int i = 0; i < reference_atoms.size(); ++i)
  {
    msd[0] += ((current_atoms[i].getWrapped()[0] - reference_atoms[i].getWrapped()[0]) - diff_CoM[0]) * ((current_atoms[i].getWrapped()[0] - reference_atoms[i].getWrapped()[0]) - diff_CoM[0]);
    msd[1] += ((current_atoms[i].getWrapped()[1] - reference_atoms[i].getWrapped()[1]) - diff_CoM[1]) * ((current_atoms[i].getWrapped()[1] - reference_atoms[i].getWrapped()[1]) - diff_CoM[1]);
    msd[2] += ((current_atoms[i].getWrapped()[2] - reference_atoms[i].getWrapped()[2]) - diff_CoM[2]) * ((current_atoms[i].getWrapped()[2] - reference_atoms[i].getWrapped()[2]) - diff_CoM[2]);
  }

  msd[0] /= reference_atoms.size();
  msd[1] /= reference_atoms.size();
  msd[2] /= reference_atoms.size();
  msd[3] = msd[0] + msd[1] + msd[2];

  return msd;
}

int main(int argc, char **argv)
{
  string reference, current, input_filename, output_file, str; // reference and current position files, input filename, junk variable
  vector <double> msd; // MSD values for x, y, z, and combined
  vector <Atom> reference_atoms, current_atoms;
  bool dump, input_file; // Is or is not a dump file, is or is not an input file
  int N1, n_atoms_read1 = 0; // Number of atoms in the simulation

  // Atom parameters
  int atom_id, atom_type;
  double charge, x, y, z;

  // Check for input, prompt for values we need
  if (argc == 2)
  {
    // assume that we have been given an input file with a list of files to process.
    input_file = true;
    input_filename = argv[1];
  }
  else
  {
    input_file = false;
    if (argc != 3)
    {
      cout << "Please enter the filename of the reference system: ";
      cin  >> reference;

      cout << "Please enter the filename of the current system: ";
      cin  >> current;
    }
    else
    {
      reference = argv[1];
      current = argv[2];
    }
  }

  // Open up the file streams - check for errors
  ifstream fin_input;
  ofstream fout;
  if (input_file)
  {
    fin_input.open(input_filename.c_str());
    if (fin_input.fail())
    {
      cerr << "Error opening file: \"" << input_filename << "\"\n";
      return 1;
    }

    fout.open("msd_data.csv");
    if (fout.fail())
    {
      cerr << "Error opening up output file \"msd_data.csv\"\n";
      return 1;
    }

    // get the reference file
    fin_input >> reference;
  }

  ifstream fin1(reference.c_str());
  if (fin1.fail())
  {
    cerr << "Error opening file: \"" << reference << "\"\n";
    return 1;
  }

  // Pull out the relevant information from the heading
  // Determine if the file is a dump file or an input file
  if (reference.find("dump") == string::npos) // if "dump" is not in the filename
  {
    dump = false; // it isn't a dump file, assume the file is an input file.
  }
  else
  {
    dump = true; // if it is in the filename, assume that it's a LAMMPS dump file
  }

  if (dump)
  {
    // This is for a LAMMPS dump file
    getline(fin1, str); // Gets ITEM: TIMESTEP
    getline(fin1, str); // Gets the timestep number

    getline(fin1, str); // Gets ITEM: NUMBER OF ATOMS
    fin1 >> N1;
    fin1.ignore();
    getline(fin1, str); //get ITEM: BOX BOUNDS
    getline(fin1, str); // get the x bounds
    getline(fin1, str); // get the y bounds
    getline(fin1, str); // get the z bounds
    getline(fin1, str); // Gets ITEM: ATOMS <data types>
  }
  else
  {
    // This is for a LAMMPS input file
    getline(fin1, str); // gets comment line
    fin1 >> N1 >> str; // number of atoms
    fin1.ignore();
    getline(fin1, str); // n atom types
    getline(fin1, str); // x bounds
    getline(fin1, str); // y bounds
    getline(fin1, str); // z bounds
    fin1 >> str;
    fin1.ignore();
  }

  reference_atoms.resize(N1, Atom());

  while (getline(fin1, str))
  {
    stringstream ss(str);
    if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z))
    {
      // If the line does NOT have 6 entries, it must have 5, so
      ss >> atom_id >> atom_type >> x >> y >> z;
      charge = 0.0;
    }

    Position p(x,y,z);
    reference_atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, p);
    ++n_atoms_read1;
  }

  fin1.close();

  if (n_atoms_read1 != N1)
  {
    cerr << "Error reading reference atoms. n_atoms_read = " << n_atoms_read1 << " != N = " << N1 << endl;
    return 3;
  }

  if (input_file)
  {
    while (fin_input >> current)
    {
      msd = compareFiles(reference_atoms, current);
      fout << msd[0] << "," << msd[1] << "," << msd[2] << "," << msd[3] << "\n";
    }
    fin_input.close();
    fout.close();
  }
  else
  {
    msd = compareFiles(reference_atoms, current);
    cout << "The MSD values are:"
         << "\n\t MSD_x = " << msd[0]
         << "\n\t MSD_y = " << msd[1]
         << "\n\t MSD_z = " << msd[2]
         << "\n\t MSD_total = " << msd[3] << endl;
  }

  return 0;
}

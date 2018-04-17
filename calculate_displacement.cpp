#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "atom.h"

using namespace std;

bool checkUnwrapped(const vector <Atom> & a)
{
  for (unsigned int i = 0; i < a.size(); ++i)
  {
    if (a[i].getXu() != 0 || a[i].getYu() != 0 || a[i].getZu() != 0)
    {
      return false;
    }
  }
  return true;
}

void processFile(ifstream & fin, vector <Atom> & atoms)
{
  string str;
  unsigned int atom_id, n_read = 0;
  int atom_type, grain_num;
  bool warn = false;
  double charge, x, y, z, f_i, xu, yu, zu;
  getline(fin,str); // we ignore the first line
  while (getline(fin, str))
  {
    stringstream ss(str);
    stringstream::pos_type pos = ss.tellg();
    if (!(ss >> atom_id >> atom_type >> charge >> x >> y >> z >> grain_num >> f_i >> xu >> yu >> zu))
    {
      ss.clear();
      ss.seekg(pos, ss.beg);
      charge = 0.0;
      if (!(ss >> atom_id >> atom_type >> x >> y >> z >> grain_num >> f_i >> xu >> yu >> zu))
      {
        if (!warn)
        {
          cout << "WARNING: Unable to determine the unwrapped coordinates.  Proceeding using assumed wrapped coordinates.\n";
          warn = true;
        }
        xu = x;
        yu = y;
        zu = z;
        // unwrapped coordinates not included.  Cannot calculate displacement vectors
        /*cout << "Error reading data.  Please make sure the data file has the following information (in order):\n"
             << "\t<atom_id> <atom_type> <atom_charge>(optional) x y z <grain_number> <orientation_parameter> xu yu zu\n";
        exit(4);*/
      }
    }
    if (atom_id > atoms.size())
    {
      atoms.resize(atom_id, Atom());
    }
    atoms[atom_id - 1] = Atom(atom_id, atom_type, charge, x, y, z);
    atoms[atom_id - 1].setXu(xu);
    atoms[atom_id - 1].setYu(yu);
    atoms[atom_id - 1].setZu(zu);
    atoms[atom_id - 1].setMark(grain_num);
    ++n_read;
  }

  if (checkUnwrapped(atoms))
  {
    cout << "WARNING: All unwrapped coordinate values are zero.  Proceeding using wrapped coordinates.\n";
    for (unsigned int i = 0; i < atoms.size(); ++i)
    {
      atoms[i].setXu(atoms[i].getX());
      atoms[i].setYu(atoms[i].getY());
      atoms[i].setZu(atoms[i].getZ());
    }
  }

  if (n_read != atoms.size())
  {
    cout << "Error reading file.  n_read = " << n_read << " != atoms.size() = " << atoms.size() << endl;
    exit(3);
  }
}

void writeData(ofstream & fout, vector <Atom> const & atoms1, vector <Atom> const & atoms2)
{
  double disp_x, disp_y, disp_z, disp_mag;
  fout << "VARIABLES = \"Atom ID\", \"Atom Type\", \"Atom Charge\", \"Xu\", \"Yu\", \"Zu\", \"Changes Grain\", \"X(K)\", \"Y(K)\", \"Z(K)\", \"Magnitude\"\n";
  for (unsigned int i = 0; i < atoms1.size(); ++i)
  {
    disp_x = atoms1[i].getXu() - atoms2[i].getXu();
    disp_y = atoms1[i].getYu() - atoms2[i].getYu();
    disp_z = atoms1[i].getZu() - atoms2[i].getZu();
    int grain_change = 0;
    if (atoms1[i].getMark() != atoms2[i].getMark())
    {
      grain_change = 1;
    }
    disp_mag = sqrt(disp_x * disp_x + disp_y * disp_y + disp_z * disp_z);
    fout << atoms2[i].getId() << " " << atoms2[i].getType() << " "
         << atoms2[i].getCharge() << " " << atoms2[i].getXu() << " "
         << atoms2[i].getYu() << " " << atoms2[i].getZu() << " " << grain_change << " "
         << disp_x << " " << disp_y << " " << disp_z << " " << disp_mag << endl;
  }
}

int main(int argc, char **argv)
{
  string filename1, filename2, input_file, output_file, str;
  vector <Atom> atoms_1, atoms_2;
  bool input, compare_to_one_file = false;

  if (argc == 2)
  {
    input_file = argv[1];
    input = true;
    char comp;
    //cout << "Would you like to compare all files to the first file? ";
    //cin  >> comp;

    //if (comp == 'y' || comp == 'Y')
    //  compare_to_one_file = true;
  }
  else if (argc == 3)
  {
    filename1 = argv[1];
    filename2 = argv[2];
    input = false;
    output_file = filename1.substr(0,filename1.find("_")) + "to"+filename2.substr(0,filename2.find("_")) + "_displacement_data.dat";
  }
  else
  {
    cout << "Please enter the reference system: ";
    cin  >> filename1;

    cout << "Please enter the compared system: ";
    cin  >> filename2;

    output_file = filename1.substr(0,filename1.find("_")) + "to"+filename2.substr(0,filename2.find("_")) + "_displacement_data.dat";

    input = false;
  }

  if (input)
  {
    ifstream fin_input(input_file.c_str());
    fin_input >> filename1;

    ifstream fin(filename1.c_str());
    if (fin.fail())
    {
      cout << "Error opening file \"" << filename1 << "\"\n";
      return 1;
    }

    processFile(fin, atoms_1);
    fin.close();

    cout << "Processing of file \"" << filename1 << "\" complete.\n";

    while (fin_input >> filename2)
    {
      output_file = filename1.substr(0,filename1.find("_")) + "to"+filename2.substr(0,filename2.find("_")) + "_displacement_data.dat";
      ifstream fin2(filename2.c_str());
      //output_file = filename2.substr(0,filename2.find("_")) + "_displacement_data.dat";

      if (fin2.fail())
      {
        cout << "Error opening file \"" << filename2 << "\"\n";
        return 1;
      }

      processFile(fin2,atoms_2);
      fin2.close();
      cout << "Processing of file \"" << filename2 << "\" complete.\n";

      if (atoms_1.size() != atoms_2.size())
      {
        cout << "Error reading atom data.  atoms_1.size() = " << atoms_1.size() << " != atoms_2.size() = " << atoms_2.size() << endl;
        return 3;
      }

      ofstream fout(output_file.c_str());
      if (fout.fail())
      {
        cout << "Error opening file \"" << output_file << "\"\n";
        return 1;
      }

      writeData(fout, atoms_1, atoms_2);
      // Here we implement an option to compare with only one data file
      if (!compare_to_one_file)
      {
        atoms_1 = atoms_2;
        filename1 = filename2;
      }

      fout.close();
    }
  }
  else
  {
    ifstream fin(filename1.c_str());
    if (fin.fail())
    {
      cout << "Error opening file \"" << filename1 << "\"\n";
      return 1;
    }

    processFile(fin, atoms_1);
    fin.close();

    ifstream fin2(filename2.c_str());
    if (fin2.fail())
    {
      cout << "Error opening file \"" << filename2 << "\"\n";
      return 1;
    }

    processFile(fin2, atoms_2);
    fin2.close();

    ofstream fout(output_file.c_str());
    if (fout.fail())
    {
      cout << "Error opening file \"" << output_file << "\"\n";
      return 1;
    }
    writeData(fout, atoms_1, atoms_2);
    fout.close();
  }

  return 0;
}

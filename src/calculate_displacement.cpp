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

void writeData(ofstream & fout, vector <Atom> const & reference, vector <Atom> const & compared)
{
  double disp_x, disp_y, disp_z, disp_mag;
  fout << "VARIABLES = \"Atom ID\", \"Atom Type\", \"Atom Charge\", \"Xu\", \"Yu\", \"Zu\", \"Changes Grain\", \"X(K)\", \"Y(K)\", \"Z(K)\", \"Magnitude\"\n";
  for (unsigned int i = 0; i < reference.size(); ++i)
  {
    disp_x = compared[i].getXu() - reference[i].getXu();
    disp_y = compared[i].getYu() - reference[i].getYu();
    disp_z = compared[i].getZu() - reference[i].getZu();
    int grain_change = 0;
    if (reference[i].getMark() != compared[i].getMark())
    {
      grain_change = 1;
    }
    disp_mag = sqrt(disp_x * disp_x + disp_y * disp_y + disp_z * disp_z);
    fout << reference[i].getId() << " " << reference[i].getType() << " "
         << reference[i].getCharge() << " " << reference[i].getXu() << " "
         << reference[i].getYu() << " " << reference[i].getZu() << " " << grain_change << " "
         << disp_x << " " << disp_y << " " << disp_z << " " << disp_mag << endl;
  }
}

int main(int argc, char **argv)
{
  string filename1, filename2, input_file, output_file, str;
  vector <Atom> reference, compared;
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

    processFile(fin, reference);
    fin.close();

    cout << "Processing of file \"" << filename1 << "\" complete.";
    cout.flush();

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

      processFile(fin2,compared);
      fin2.close();
      cout << "\r";
      cout << "Processing of file \"" << filename2 << "\" complete.";

      if (compared.size() != reference.size())
      {
        cout << "Error reading atom data.  reference.size() = " << reference.size() << " != compared.size() = " << compared.size() << endl;
        return 3;
      }

      ofstream fout(output_file.c_str());
      if (fout.fail())
      {
        cout << "Error opening file \"" << output_file << "\"\n";
        return 1;
      }

      writeData(fout, reference, compared);
      // Here we implement an option to compare with only one data file
      if (!compare_to_one_file)
      {
        reference = compared;
        filename1 = filename2;
      }

      fout.close();
    }
    cout << endl;
  }
  else
  {
    ifstream fin(filename1.c_str());
    if (fin.fail())
    {
      cout << "Error opening file \"" << filename1 << "\"\n";
      return 1;
    }

    processFile(fin, reference);
    fin.close();

    ifstream fin2(filename2.c_str());
    if (fin2.fail())
    {
      cout << "Error opening file \"" << filename2 << "\"\n";
      return 1;
    }

    processFile(fin2, compared);
    fin2.close();

    ofstream fout(output_file.c_str());
    if (fout.fail())
    {
      cout << "Error opening file \"" << output_file << "\"\n";
      return 1;
    }
    writeData(fout, reference, compared);
    fout.close();
  }

  return 0;
}

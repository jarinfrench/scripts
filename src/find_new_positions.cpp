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
  string file1, file2, file3, input_file, str, str2;
  vector <Atom> reference, compared; // List of atoms to start, and compare to
  vector <string> vars;
  int id, id2, type, type2, N;
  double charge = 0, x, y, z, xu, yu, zu, num;
  double charge2 = 0, x2, y2, z2, xu2, yu2, zu2, num2;
  double xlo, xhi, ylo, yhi, zlo, zhi, Lx, Ly, Lz;
  bool input;

  if (argc == 2)
  {
    input = true;
    input_file = argv[1];
  }
  else if (argc == 3)
  {
    input = false;
    file1 = argv[1];
    file2 = argv[2];
    file3 = file2.substr(0,file2.find(".dump"))+"_tracked.dat";
  }
  else if (argc == 4)
  {
    file1 = argv[1];
    file2 = argv[2];
    file3 = argv[3];
  }
  else
  {
    cout << "Usage: " << argv[0] << " input_file OR " << argv[0] << " file1 file2 [output file]\n";
  }

  if (input)
  {
    ifstream fin_input(input_file.c_str());
    if (fin_input.fail())
    {
      cout << "Error opening file: " << input_file << endl;
      return 1;
    }
    getline(fin_input, file1);
    ifstream fin(file1.c_str());
    if (fin.fail())
    {
      cout << "Error opening file: " << file1 << endl;
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
    stringstream tmp(str);
    tmp >> str >> str;
    while (tmp >> str)
    {
      vars.push_back(str);
    }

    Lx = xhi - xlo;
    Ly = yhi - ylo;
    Lz = zhi - zlo;

    reference.resize(N, Atom());
    compared.resize(N, Atom());

    while (getline(fin, str))
    {
      stringstream ss(str);
      for (unsigned int i = 0; i < vars.size(); ++i)
      {
        ss >> num;
        if (vars[i].compare("id") == 0)
        {
          id = (int)(num);
        }
        else if (vars[i].compare("type") == 0)
        {
          type = (int)(num);
        }
        else if (vars[i].compare("q") == 0)
        {
          charge = num;
        }
        else if (vars[i].compare("x") == 0)
        {
          x = num;
        }
        else if (vars[i].compare("y") == 0)
        {
          y = num;
        }
        else if (vars[i].compare("z") == 0)
        {
          z = num;
        }
        else
        {
          continue;
        }
      }

      reference[id - 1] = Atom(id, type, charge, x, y, z);
      if (x < 11 * Lx / 20 && x > 9 * Lx / 20)
      {
        reference[id - 1].setMark(1);
      }
    }
    fin.close();

    while (getline(fin_input, file2))
    {
      ifstream fin2(file2.c_str());
      if (fin2.fail())
      {
        cout << "Error opening file: " << file2 << endl;
        return 1;
      }

      file3 = file2.substr(0,file2.find(".dump"))+"_tracked.dat";
      ofstream fout(file3.c_str());
      if (fout.fail())
      {
        cout << "Error opening file: " << file3 << endl;
        return 1;
      }

      for (int i = 0; i < 9; ++i)
      {
        getline(fin2, str);
        // Ignore the first 9 lines of the second file
      }

      while (getline(fin2, str))
      {
        stringstream ss(str);
        for (unsigned int i = 0; i < vars.size(); ++i)
        {
          ss >> num;
          if (vars[i].compare("id") == 0)
          {
            id = (int)(num);
          }
          else if (vars[i].compare("type") == 0)
          {
            type = (int)(num);
          }
          else if (vars[i].compare("q") == 0)
          {
            charge = num;
          }
          else if (vars[i].compare("x") == 0)
          {
            x = num;
          }
          else if (vars[i].compare("y") == 0)
          {
            y = num;
          }
          else if (vars[i].compare("z") == 0)
          {
            z = num;
          }
          else
          {
            continue;
          }
        }

        compared[id - 1] = Atom(id, type, charge, x, y, z);
      }
      for (unsigned int i = 0; i < reference.size(); ++i)
      {
        if (reference[i].getMark() == 1)
        {
          compared[i].setMark(1);
          fout << compared[i].getId() << " " << compared[i].getType() << " "
          << compared[i].getCharge() << " " << compared[i].getX() << " "
          << compared[i].getY() << " " << compared[i].getZ() << " "
          << compared[i].getMark() << endl;
        }
      }

      cout << "\rFile " << file2 << " processed.";
      cout.flush();
      fin2.close();
      fout.close();
    }
    fin_input.close();
    cout << endl;
  } // end if is input
  else
  {
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
    stringstream tmp(str);
    tmp >> str >> str;
    while (tmp >> str)
    {
      vars.push_back(str);
    }

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

    while (getline(fin,str) && getline(fin2, str2))
    {
      stringstream ss(str);
      stringstream ss2(str2);
      for (unsigned int i = 0; i < vars.size(); ++i)
      {
        ss >> num;
        ss2 >> num2;
        if (vars[i].compare("id") == 0)
        {
          id = (int)(num);
          id2 = (int)(num2);
        }
        else if (vars[i].compare("type") == 0)
        {
          type = (int)(num);
          type2 = (int)(num2);
        }
        else if (vars[i].compare("q") == 0)
        {
          charge = num;
          charge2 = num2;
        }
        else if (vars[i].compare("x") == 0)
        {
          x = num;
          x2 = num2;
        }
        else if (vars[i].compare("y") == 0)
        {
          y = num;
          y2 = num2;
        }
        else if (vars[i].compare("z") == 0)
        {
          z = num;
          z2 = num2;
        }
        else
        {
          continue;
        }
      }

      reference[id - 1] = Atom(id, type, charge, x, y, z);
      if (x < 11 * Lx / 20 && x > 9 * Lx / 20)
      {
        reference[id - 1].setMark(1);
      }

      compared[id2 - 1] = Atom(id2, type2, charge2, x2, y2, z2);
    }

    for (unsigned int i = 0; i < reference.size(); ++i)
    {
      if (reference[i].getMark() == 1)
      {
        compared[i].setMark(1);
        fout << compared[i].getId() << " " << compared[i].getType() << " "
        << compared[i].getCharge() << " " << compared[i].getX() << " "
        << compared[i].getY() << " " << compared[i].getZ() << " "
        << compared[i].getMark() << endl;
      }
    }
    fin.close();
    fin2.close();
    fout.close();
  }

  return 0;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

// Simple method to switch the key-value pair in a map.
template <typename A, typename B>
pair <B,A> flip_pair(const pair <A,B> &p)
{
  return pair<B,A>(p.second, p.first);
}

template <typename A, typename B>
map <B,A> flip_map(const map <A,B> &src)
{
  map <B, A> dst;
  transform(src.begin(), src.end(), inserter(dst, dst.begin()), flip_pair<A,B>);

  return dst;
}

int main(int argc, char** argv)
{
  string filename1, filename2, str, str2; // file to be parsed, written to, junk var
  string date1, date2, date3, date; // Different parts of the date
  string indicator; // line that specifies that we have important information coming up next.
  string indicator_end; // line that specifies the end of the important info.
  string unit_indicator, unit_style; // line that specifies what unit style we are using
  double Lx = 0.0, Ly = 0.0, Lz = 0.0; // box size in each direction.
  double a1, a2, a3, a4, a5, a6; // dummy variables
  int N = -1, n_labels = 0;; // Number of atoms, number of labels
  map <string, int> labels;
  map <int, string> flipped_labels;
  bool written_labels = false; // boolean determining if we've written the labels
  if (argc == 1)
  {
    cout << "Please enter the name of the LAMMPS txt output file to be parsed: ";
    cin  >> filename1;

    filename2 = "parsed_output.csv";
  }
  else if (argc == 2)
  {
    filename1 = argv[1];
    filename2 = "parsed_output.csv";
  }
  else if (argc == 3)
  {
    filename1 = argv[1];
    filename2 = argv[2];
  }

  indicator = "Per MPI rank memory allocation";
  indicator_end = "Loop time of";
  unit_indicator = "Unit style";

  // open up the filestreams
  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error: unable to open file " << filename1 << endl;
    return 1;
  }

  ofstream fout(filename2.c_str());
  if (fout.fail())
  {
    cout << "Error: unable to open file " << filename2 << endl;
    return 1;
  }

  // Read the first line of the file.  If it doesn't say LAMMPS (date), return
  // an error message saying that this doesn't look like a LAMMPS output.txt file
  fin >> str >> date1 >> date2 >> date3;
  if (str != "LAMMPS")
  {
    cout << "Error: This file does not look like a LAMMPS output.txt file.\n"
         << "If you're sure that this is a LAMMPS output file, check that the "
         << "beginning of your file has the format \"LAMMPS ([date])\"\n";
    return 2;
  }
  else
  {
    date = date1.substr(1) + " " + date2 + " " + date3.substr(0,date3.find(")"));
    fout << "\"LAMMPS version:\",\"" << date << "\"" << endl;
  }

  // Read lines until we can find the box length
  // Note that there is a leftover "\n" character that is in fin, so we use
  // the ignore function to clear out the newline character.
  fin.ignore();

  while (getline(fin, str))
  {
    if (str.find("orthogonal box") != string::npos)
    {
      // This is a bit convoluted... here I am extracting the box sizes using
      // stringstreams.  The difficulty comes in that the expected line looks
      // like (a1 a2 a3) to (a4 a5 a6), so I have to deal with the parentheses.
      // There is definitely a better way to do it, I'm just not focused on
      // what that way is right now.
      stringstream ss(str.substr(str.find("(")+1));
      ss >> a1 >> a2 >> a3 >> str2 >> str2 >> str2 >> a5 >> a6;
      stringstream ss2(str2.substr(1));
      ss2 >> a4;

      Lx = a4 - a1;
      Ly = a5 - a2;
      Lz = a6 - a3;
      fout << "\"Lx Ly Lz =\"," << Lx << "," << Ly << "," << Lz << endl;
      break;
    }
  }

  if ((Lx == 0.0) && (Ly == 0.0) & (Lz == 0.0))
  {
    cout << "Error: unable to determine box size.\n";
    return 3;
  }

  // Now get the number of atoms.  The LAMMPS output seems to specify the number
  // of atoms either by saying "Created N atoms" or by saying
  // "reading atoms ...\nN atoms".  So here I look for the line that contains
  // the word "atoms", but does not contain the character sequence "..."
  while (getline(fin, str))
  {
    if (str.find("atoms") != string::npos && str.find("...") == string::npos)
    {
      stringstream ss(str);
      if (str.find("Created") == string::npos)
      {
        ss >> N >> str2;
      }
      else
      {
        ss >> str2 >> N >> str2;
      }

      fout << "\"N =\"," << N << endl;
      break;
    }
  }

  if (N < 0)
  {
    cout << "Error: unable to determine number of atoms.\n";
    return 4;
  }

  while (getline(fin, str))
  {
    if (str.find(unit_indicator) != string::npos)
    {
      stringstream ss(str);
      ss >> str2 >> str2 >> str2 >> unit_style;
      fout << "\"Unit style:\"," << unit_style << endl;
      break; // once we've found the unit style, we can move on
    }
  }

  // Now we read through the file line by line until we reach an indicator,
  // which tells us that we have important information coming in.
  while (getline(fin, str))
  {
    if (str.find(indicator) != string::npos) // found our indicator line
    {
      getline(fin, str); // get the line containing our data labels
      stringstream ss(str);
      while (ss >> str2) // get the labels individually, but only once!
      {
        if (labels.find(str2) == labels.end()) // we don't have this label yet
        {
          // if we have already written the labels, but we don't yet have this
          // label, we need to re-write the labels.
          if (written_labels)
          {
            written_labels = false;
          }
          labels.insert(pair <string, int> (str2, ++n_labels));
        }
      }
      if (!written_labels)
      {
        flipped_labels = flip_map(labels);
        for (map<int, string>::const_iterator it = flipped_labels.begin();
             it != flipped_labels.end();
             ++it)
        {
          if (it == --flipped_labels.end())
          {
            fout << it -> second;
          }
          else
          {
            fout << it -> second << ",";
          }
        }
        fout << endl;
        written_labels = true;
      }
      while (getline(fin, str2))
      {
        if (str2.find(indicator_end) != string::npos) // We reached the end of the data
        {
          break;
        }
        else
        {
          stringstream ss(str2);
          ss >> a1;
          fout << a1;
          while (ss >> a1)
          {
            fout << "," << a1;
          }
          fout << endl;
        }
      }
    }
  }

  return 0;
}

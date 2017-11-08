#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
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

bool check4Double(const vector <string> & standard, const string & user)
{
  if (find(standard.begin(), standard.end(), user) != standard.end())
  {
    return true; // we do have a double
  }
  return false; // we don't have a double
}

int main(int argc, char** argv)
{
  string filename1, filename2, str, str2; // file to be parsed, written to, junk var
  string date1, date2, date3, date; // Different parts of the date
  string unit_style; // line that specifies what unit style we are using
  vector <string> indicators; // the indicators we want to look for
  vector <string> default_indicators; // indicators we already look for
  double Lx = 0.0, Ly = 0.0, Lz = 0.0; // box size in each direction.
  double time_step; // time step value
  double a1, a2, a3, a4, a5, a6; // dummy variables
  int N = -1, n_labels = 0;; // Number of atoms, number of labels
  int n_indicators; // number of indicators
  map <string, int> labels;
  map <int, string> flipped_labels;
  bool written_labels = false; // boolean determining if we've written the labels
  bool set_unit_style = false; // have we set the unit style or no?

  // Default indicators:
  default_indicators.push_back("orthogonal box"); // 0 - indicates the bounds of the system (same line)
  default_indicators.push_back("atoms"); // 1 - Number of atoms indicator line
  default_indicators.push_back("..."); // 2 - same
  default_indicators.push_back("Created"); // 3 - same
  default_indicators.push_back("in group"); // 4 - avoid adding atoms that are put in groups
  default_indicators.push_back("Per MPI rank memory allocation"); // 5 - indicates the next line contains the names of the data
  default_indicators.push_back("Memory usage per processor"); // 6 - same thing as above, but for a slightly different process
  default_indicators.push_back("Loop time of"); // 7 - Indicates the end of data
  default_indicators.push_back("Unit style"); // 8 - Indicates the unit style being used
  default_indicators.push_back("  Time step     :"); // 9 - Indicates the size of a time step

 // Argument checking
  if (argc == 1)
  {
    cout << "Please enter the name of the LAMMPS txt output file to be parsed: ";
    cin  >> filename1;

    filename2 = "parsed_output.csv";

    cout << "Please enter the number of indicators you want to check for: ";
    cin  >> n_indicators;
    indicators.resize(n_indicators,"not assigned");
    for (int i = 0; i < n_indicators; ++i)
    {
      cout << "Please enter indicator " << i << ": ";
      cin  >> indicators[i];
      if (check4Double(default_indicators, indicators[i]))
      {
        cout << "Indicator \"" << indicators[i] << "\" is already being checked.\n";
      }
    }
  }
  else if (argc == 2)
  {
    filename1 = argv[1];
    filename2 = "parsed_output.csv";

    cout << "Please enter the number of indicators you want to check for: ";
    cin  >> n_indicators;
    indicators.resize(n_indicators,"not assigned");
    for (int i = 0; i < n_indicators; ++i)
    {
      cout << "Please enter indicator " << i << ": ";
      cin  >> indicators[i];
      if (check4Double(default_indicators, indicators[i]))
      {
        cout << "Indicator \"" << indicators[i] << "\" is already being checked.\n";
      }
    }
  }
  else if (argc == 3)
  {
    filename1 = argv[1];
    filename2 = argv[2];

    cout << "Please enter the number of indicators you want to check for: ";
    cin  >> n_indicators;
    indicators.resize(n_indicators,"not assigned");
    for (int i = 0; i < n_indicators; ++i)
    {
      cout << "Please enter indicator " << i << ": ";
      cin  >> indicators[i];
      if (check4Double(default_indicators, indicators[i]))
      {
        cout << "Indicator \"" << indicators[i] << "\" is already being checked.\n";
      }
    }
  }
  else if (argc > 3)
  {
    filename1 = argv[1];
    filename2 = argv[2];

    n_indicators = argc - 3;
    indicators.resize(n_indicators, "not assigned");

    for (int i = 0; i < n_indicators; ++i)
    {
      indicators[i] = argv[i+3];
      if (check4Double(default_indicators, indicators[i]))
      {
        cout << "Indicator \"" << indicators[i] << "\" is already being checked.\n";
      }
    }
  }

  if (n_indicators > 0)
  {
    cout << "Note: user-specified indicators are not currently implemented beyond prompting for them.\n";
  }

  // open up the filestreams
  ifstream fin(filename1.c_str());
  if (fin.fail())
  {
    cout << "Error: unable to open file \"" << filename1 << "\"" << endl;
    return 1;
  }

  fstream fout(filename2.c_str(), fstream::out);
  if (fout.fail())
  {
    cout << "Error: unable to open file \"" << filename2 << "\"" << endl;
    return 1;
  }
  fout.close();

  fout.open(filename2.c_str(), fstream::in | fstream::out);
  fout.precision(5);
  fout.setf(ios::fixed);
  fout.setf(ios::showpoint);

  // Read the first line of the file.  If it doesn't say LAMMPS (date), return
  // an error message saying that this doesn't look like a LAMMPS output.txt file
  fin >> str >> date1 >> date2 >> date3;
  if (str != "LAMMPS")
  {
    cout << "Error: This file does not look like a LAMMPS output.txt file.\n"
         << "If you are sure that this is a LAMMPS output file, check that the "
         << "beginning of your file has the format \"LAMMPS ([date])\"\n";
    return 2;
  }
  else
  {
    date = date1.substr(1) + " " + date2 + " " + date3.substr(0,date3.find(")"));
    fout << "\"LAMMPS version:\",\"" << date << "\"" << endl;
  }

  // Note that there is a leftover "\n" character that is in fin, so we use
  // the ignore function to clear out the newline character.
  fin.ignore();

  while (getline(fin, str))
  {
    // find the box length
    if (str.find(default_indicators[0]) != string::npos)
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
      continue;
    }

    // Get the number of atoms.  The LAMMPS output seems to specify the number
    // of atoms either by saying "Created N atoms" or by saying
    // "reading atoms ...\nN atoms".  So here I look for the line that contains
    // the word "atoms", but does not contain the character sequence "..."
    if (str.find(default_indicators[1]) != string::npos && str.find(default_indicators[2]) == string::npos)
    { // found the word "atoms", did not find the string "..."
      stringstream ss(str);
      if (str.find(default_indicators[3]) == string::npos)
      { // Did not find the word "Created"
        if (str.find(default_indicators[4]) == string::npos)
        { // did not find the words "in group"
          if (N != -1)
          {
            int temp;
            ss >> temp >> str2;
            N += temp;
            long last_pos = fout.tellp(); // Current position of the output stream
            fout.seekg(0); // Go back to the beginning of the file we are writing to
            getline(fout,str2); // Get the first line
            getline(fout,str2); // get the second line
            long pos = fout.tellg(); // get the position of the next line
            fout.seekp(pos); // move the output position to that position
            fout << "\"N =\"," << N << endl;
            fout.seekp(last_pos+1); // Go to the end of the file, and add one more position to it.
          }
          else
          {
            ss >> N >> str2;
            fout << "\"N =\"," << N << endl;
          }
        }
      }
      else
      { // found the word "Created"
        if (str.find(default_indicators[4]) == string::npos)
        { // did not find the words "in group"
          if (N != -1)
          {
            int temp;
            ss >> str2 >> temp >> str2;
            N += temp;
            long last_pos = fout.tellp();
            fout.seekg(0); // Go back to the beginning of the file we are writing to
            getline(fout,str2); // Get the first line
            getline(fout,str2);
            long pos = fout.tellg();
            fout.seekp(pos);
            fout << "\"N =\"," << N << endl;
            fout.seekp(last_pos+1);
          }
          else
          {
            ss >> str2 >> N >> str2;
            fout << "\"N =\"," << N << endl;
          }
        }
      }
      continue;
    }

    if (str.find(default_indicators[8]) != string::npos && !set_unit_style)
    {
      stringstream ss(str);
      ss >> str2 >> str2 >> str2 >> unit_style;
      fout << "\"Unit style:\"," << unit_style << endl;
      set_unit_style = true;
      continue; // once we've found the unit style, we can move on
    }

    if (str.find(default_indicators[9]) != string::npos)
    {
      stringstream ss(str);
      ss >> str >> str >> str >> time_step;
      fout << "\"Time step:\"," << time_step << endl;
      continue;
    }

    // Now we read through the file line by line until we reach an indicator,
    // which tells us that we have important information coming in.
    if (str.find(default_indicators[5]) != string::npos || str.find(default_indicators[6]) != string::npos) // found our indicator line
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
        if (str2.find(default_indicators[7]) != string::npos) // We reached the end of the data
        {
          break;
        }
        else
        {
          stringstream ss(str2);
          ss >> a1;
          fout << setprecision(5) << a1;
          while (ss >> a1)
          {
            fout << "," << setprecision(5) << a1;
          }
          fout << endl;
        }
      }
    }
  }

  fin.close();
  fout.close();

  // Error check for determining box size
  if ((Lx == 0.0) && (Ly == 0.0) & (Lz == 0.0))
  {
    cout << "Error: unable to determine box size.\n";
    return 3;
  }

  if (N < 0)
  {
    cout << "Error: unable to determine number of atoms.\n";
    return 4;
  }



  return 0;
}

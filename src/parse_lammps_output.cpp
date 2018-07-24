#include <iostream>
// #include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cxxopts.hpp>
#include "error_code_defines.h"

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

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

bool check4Double(const vector <string> & standard, const string & user)
{
  if (find(standard.begin(), standard.end(), user) != standard.end())
  {
    return true; // we do have a double
  }
  return false; // we don't have a double
}

string removeParentheses(string str)
{
  str.erase(remove_if(str.begin(), str.end(),
                      [](char ch){return ch=='(' || ch == ')'; }), str.end());
  return str;
}

void setDefaultIndicators(vector <string>& indicator)
{
  indicator.push_back("orthogonal box"); // 0 - indicators the bounds of the system (same line)
  indicator.push_back("Reading atoms..."); // 1 - number of atoms indicator line (number of atoms on next line)
  indicator.push_back("Created"); // 2 - atoms indicator line (same line)
  indicator.push_back("Per MPI rank memory allocation"); // 3 - indicates the next line contains data labels
  indicator.push_back("Memory usage per processor"); // 4 - same as 4 but for a slightly different process
  indicator.push_back("Loop time of"); // 5 - indicates the end of data
  indicator.push_back("Unit style"); // 6 - indicates the unit style being used
  indicator.push_back("  Time step     :"); // 7 - Indicates the size of a time step
}

void parseFile(const string& input_file, const string& output_file,
               const vector <string>& indicators, const string& separator)
{
  string str, date1, date2, date3, date, str2;
  double Lx, Ly, Lz;
  int N, time_step;
  string unit_style;
  bool set_unit_style = false;
  vector <string> default_indicators;
  stringstream ss;

  setDefaultIndicators(default_indicators);
  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);

  ofstream fout(output_file.c_str()); // creates the file if it doesn't exist
  checkFileStream(fout, output_file);
  fout.precision(5);
  fout << fixed;

  // Read the first line of the file.  If it doesn't say LAMMPS (date), return
  // an error message saying it doesn't look like a LAMMPS output.txt file
  fin >> str >> date1 >> date2 >> date3;
  if (str.compare("LAMMPS") != 0)
  {
    cout << "Error: this file does not look like a LAMMPS output file.\n"
         << "If you are sure that this is a LAMMPS output file, check that the "
         << "beginning of your file has the format \"LAMMPS ([date])\"\n";
    exit(FILE_FORMAT_ERROR);
  }
  else
  {
    date = date1.substr(1) + " " + date2 + " " + date3.substr(0,date3.find(")"));
  }
  fin.ignore();

  while (getline(fin, str))
  {

    // box length
    if (str.find(default_indicators[0]) != string::npos)
    {
      double a1, a2, a3, a4, a5, a6;
      str2 = str.substr(str.find("("));
      str2 = removeParentheses(str2);
      ss.str(str2); // This will now have the format "xmin ymin zmin to xmax ymax zmax"
      ss >> a1 >> a2 >> a3 >> str2 >> a4 >> a5 >> a6;
      Lx = a4 - a1;
      Ly = a5 - a2;
      Lz = a6 - a3;
      continue; // we're done with this line
    }

    // Number of atoms.
    if (str.find(default_indicators[1]) != string::npos)
    { // Found "Reading atoms..."
      getline(fin, str);
      ss.str(str);
      ss >> N >> str2;
      continue;
    }

    if (str.find(default_indicators[2]) != string::npos)
    { // found "Created"
      ss.str(str);
      ss >> str2 >> N >> str2;
      continue;
    }

    // Unit style
    if (str.find(default_indicators[6]) != string::npos && !set_unit_style)
    {
      ss.str(str);
      ss >> str2 >> str2 >> str2 >> unit_style;
      set_unit_style = true;
      continue;
    }

    // Time step
    if (str.find(default_indicators[7]) != string::npos)
    {
      ss.str(str);
      ss >> str2 >> str2 >> str2 >> time_step;
      continue;
    }

    // Data labels, and subsequent data
    if (str.find(default_indicators[3]) != string::npos ||
        str.find(default_indicators[4]) != string::npos)
    {
      getline(fin, str);
      ss.str(str);
      while (ss >> str2) // read each label
      {
        if (labels.find(str2) == labels.end()) // Don't have this label yet
        {
          // if we have already written the labels, but we don't yet have this
          // label, we need to rewrite the labels.
        }
      }
    }

  }
}

int main(int argc, char** argv)
{
  vector <string> indicators;
  try
  {
    cxxopts::Options options(argv[0], "LAMMPS ouput parsing program.");
    options
      .positional_help("Input")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("i,input", "Input file to be processed.", cxxopts::value<string>(input_file), "Input")
        ("o,output", "Output filename.", cxxopts::value<string>(output_file)->default("parsed_output.txt"), "Output")
        ("add-indicator", "Add an indicator to check for in the input file", cxxopts::value<vector <string> >(), "String")
        ("s,separator", "Specify the separator to use", cxxopts::value<string>(separator)->default_value(",")->implicit_value(" "))
        ("h,help", "Show the help");

    options.parse_positional({"input"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("file") == 0)
    {
      cout << options.help() << endl;
    }

    if (result.count("add-indicator"))
    {
      cout << "Note that user-defined indicators are not currently implemented.\n";
      indicators = result["add-indicator"].as<vector <string> >();
    }

    if (result.count("input"))
    {
      parseFile(input_file, output_file, indicators, separator);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }
  return EXIT_SUCCESS;
}

/*int main(int argc, char** argv)
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







  // This makes sure the file exists
  fstream fout(filename2.c_str(), fstream::out);
  if (fout.fail())
  {
    cout << "Error: unable to open file \"" << filename2 << "\"" << endl;
    return 1;
  }
  fout.close(); // We close this here...

  // so we can reopen the file for both reading and writing at the same time
  fout.open(filename2.c_str(), fstream::in | fstream::out);
  fout.precision(5);
  fout.setf(ios::fixed);
  fout.setf(ios::showpoint);

  fout << "\"LAMMPS version:\",\"" << date << "\"" << endl;


  // Note that there is a leftover "\n" character that is in fin, so we use
  // the ignore function to clear out the newline character.
  fin.ignore();

  while (getline(fin, str))
  {

    fout << "\"N = \"," << N << endl;


      fout << "\"Unit style:\"," << unit_style << endl;




      fout << "\"Time step:\"," << time_step << endl;


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
}*/

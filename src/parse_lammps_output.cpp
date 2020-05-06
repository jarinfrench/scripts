#include <iostream>
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
    cerr << "Error opening file \"" << file << "\"\n";
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
  indicator.push_back("reading atoms ..."); // 1 - number of atoms indicator line (number of atoms on next line)
  indicator.push_back("Created"); // 2 - atoms indicator line (same line)
  indicator.push_back("Per MPI rank memory allocation"); // 3 - indicates the next line contains data labels
  indicator.push_back("Memory usage per processor"); // 4 - same as 4 but for a slightly different process
  indicator.push_back("Loop time of"); // 5 - indicates the end of data
  indicator.push_back("Unit style"); // 6 - indicates the unit style being used
  indicator.push_back("  Time step     :"); // 7 - Indicates the size of a time step
}

void writeLabels(ofstream& fout, const map <string, int>& labels, const string& separator)
{
  map <int, string> flipped_labels = flip_map(labels);
  for (map <int, string>::const_iterator it = flipped_labels.begin();
       it != flipped_labels.end();
       ++it)
  {
    if (it == --flipped_labels.end()) {fout << it -> second << endl;}
    else {fout << it -> second << separator;}
  }
}

void parseFile(const string& input_file, const string& output_file,
               const vector <string>& indicators, const string& separator)
{
  string str, date1, date2, date3, date, str2;
  double Lx, Ly, Lz;
  double a1;
  double time_step;
  int N, n_labels = 0;
  string unit_style;
  vector <string> default_indicators;
  map <string, int> labels;
  map <int, string> flipped_labels;
  stringstream ss;

  setDefaultIndicators(default_indicators);
  ifstream fin(input_file.c_str());
  checkFileStream(fin, input_file);

  ofstream fout(output_file.c_str());
  checkFileStream(fout, output_file);
  fout.precision(5);
  fout << fixed;

  // Read the first line of the file.  If it doesn't say LAMMPS (date), return
  // an error message saying it doesn't look like a LAMMPS output.txt file
  fin >> str >> date1 >> date2 >> date3;
  if (str.compare("LAMMPS") != 0)
  {
    cerr << "Error: this file does not look like a LAMMPS output file.\n"
         << "If you are sure that this is a LAMMPS output file, check that the "
         << "beginning of your file has the format \"LAMMPS ([date])\"\n";
    exit(FILE_FORMAT_ERROR);
  }
  else
  {
    date = date1.substr(1) + " " + date2 + " " + date3.substr(0,date3.find(")"));
  }
  fin.ignore();

  // Write the LAMMPS version to the file
  fout << "\"LAMMPS version:\"" << separator << "\"" << date << "\"" << endl;

  bool box_size_found = false;
  bool n_atoms_found = false;
  bool time_step_found = false;
  bool unit_style_set = false;
  bool labels_written = false;

  while (getline(fin, str))
  {

    // box length
    if (str.find(default_indicators[0]) != string::npos && !box_size_found)
    {
      double a1, a2, a3, a4, a5, a6;
      str2 = str.substr(str.find("("));
      str2 = removeParentheses(str2);
      ss.str(str2); // This will now have the format "xmin ymin zmin to xmax ymax zmax"
      ss >> a1 >> a2 >> a3 >> str2 >> a4 >> a5 >> a6;
      Lx = a4 - a1;
      Ly = a5 - a2;
      Lz = a6 - a3;

      fout << "\"Lx Ly Lz =\"" << separator << Lx << separator << Ly << separator << Lz << endl;
      box_size_found = true;
      continue; // we're done with this line
    }

    // Number of atoms.
    if (str.find(default_indicators[1]) != string::npos && !n_atoms_found)
    { // Found "reading atoms ..."
      ss.clear();
      getline(fin, str);
      ss.str(str);
      ss >> N >> str2;
      n_atoms_found = true;

      fout << "\"N =\"" << separator << N << endl;
      continue;
    }
    else if (str.find(default_indicators[1]) != string::npos && n_atoms_found)
    {
      ss.clear();
      getline(fin, str);
      ss.str(str);
      int tmp;
      ss >> tmp >> str2;

      if (tmp != N)
      {
        fout << "\"N =\"" << separator << tmp << endl;
      }
      continue;
    }

    if (str.find(default_indicators[2]) != string::npos && !n_atoms_found)
    { // found "Created"
      ss.clear();
      ss.str(str);
      ss >> str2 >> N >> str2;
      n_atoms_found = true;

      fout << "\"N =\"" << separator << N << endl; // TODO: Make sure this only occurs once!
      continue;
    }
    else if (str.find(default_indicators[2]) != string::npos && n_atoms_found)
    {
      ss.clear();
      ss.str(str);
      int tmp;
      ss >> str2 >> tmp >> str2;

      if (tmp != N)
      {
        fout << "\"N =\"" << separator << N << endl;
      }
      continue;
    }

    // Unit style
    if (str.find(default_indicators[6]) != string::npos && !unit_style_set)
    {
      ss.clear();
      ss.str(str);
      ss >> str2 >> str2 >> str2 >> unit_style;
      unit_style_set = true;

      fout << "\"Unit style:\"" << separator << unit_style << endl;
      continue;
    }

    // Time step
    if (str.find(default_indicators[7]) != string::npos && !time_step_found)
    {
      ss.clear();
      ss.str(str);
      ss >> str2 >> str2 >> str2 >> time_step;
      time_step_found = true;

      fout << "\"Time step:\"" << separator << time_step << endl;
      continue;
    }

    // Data labels, and subsequent data
    if (str.find(default_indicators[3]) != string::npos ||
        str.find(default_indicators[4]) != string::npos)
    {
      getline(fin, str);
      ss.clear();
      ss.str(str);
      while (ss >> str2) // read each label
      {
        if (labels.find(str2) == labels.end()) // Don't have this label yet
        {
          // if we have already written the labels, but we don't yet have this
          // label, we need to rewrite the labels.
          if (labels_written)
          {
            labels_written = false;
          }
          labels.insert(pair <string, int> (str2, ++n_labels));
        }
      }
      if (!labels_written)
      {
        if (!time_step_found)
        {
          fout << "\"This data is (assumed) to be from a minimization\"\n";
        }
        writeLabels(fout, labels, separator);
        labels_written = true;
      }

      while (getline(fin, str2)) // Now get the data
      {
        if (str2.find(default_indicators[5]) != string::npos)
        {
          fout << endl; // break up the data
          break;
        }
        else
        {
          ss.clear();
          ss.str(str2);
          ss >> a1;
          fout.precision(5);
          fout << a1;
          while (ss >> a1)
          {
            fout << separator << a1;
          }
          fout << endl;
        }
      }
    }
  }

  fin.close();
  fout.close();

  if (!box_size_found)
  {
    cout << "Warning: unable to determine box size.\n";
  }

  if (!n_atoms_found)
  {
    cout << "Warning: unable to determine the number of atoms.\n";
  }
}

int main(int argc, char** argv)
{
  vector <string> indicators;
  string input_file, output_file, separator;
  try
  {
    cxxopts::Options options(argv[0], "LAMMPS output parsing program.");
    options
      .positional_help("Input")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("i,input", "Input file to be processed", cxxopts::value<string>(input_file), "Input")
        ("o,output", "Output filename", cxxopts::value<string>(output_file)->default_value("parsed_output.txt"), "Output")
        ("s,separator", "Specify the separator to use", cxxopts::value<string>(separator)->default_value(",")->implicit_value(" "))
        ("h,help", "Show the help");

    options.parse_positional({"input"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || result.count("input") == 0)
    {
      cout << options.help() << endl;
    }

    if (result.count("input"))
    {
      parseFile(input_file, output_file, indicators, separator);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }
  return EXIT_SUCCESS;
}

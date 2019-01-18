#include <iostream>
#include <fstream>
#include <string>
#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

void convertFile(const string& infile, const string& outfile)
{
  string str;

  ifstream fin(infile.c_str());
  checkFileStream(fin, infile);

  ofstream fout(outfile.c_str());
  checkFileStream(fout, outfile);

  for (int i = 0; i < 9; ++i)
  {
    getline(fin, str);
  }

  // read through the file one item at a time (stopping at ','), then output that item with a space.
  while (getline(fin, str)) {fout << str << "\n";}

  cout << "File " << infile << " processed.\r";
  cout.flush();
}

int main(int argc, char** argv)
{
  vector <string> infile;
  string outfile;

  try
  {
    cxxopts::Options options(argv[0], "Converts a LAMMPS dump file to a .tec file");
    options
    .positional_help("file file")
    .show_positional_help();

    options
    .allow_unrecognised_options()
    .add_options()
      ("i,input", "Input file", cxxopts::value<vector <string> >(infile), "file")
      ("o,output", "Name of converted file", cxxopts::value<string>(outfile)->default_value("*.tec"), "file")
      ("h,help", "Show the help");

    options.parse_positional({"input"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !(result.count("input")))
    {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (result.count("input"))
    {
      for (unsigned int i = 0; i < infile.size(); ++i)
      {
        if (!(result.count("output")))
        {
          outfile = infile[i].substr(0, infile[i].find(".dump")) + ".tec";
        }
        convertFile(infile[i], outfile);
      }
      cout << endl;
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

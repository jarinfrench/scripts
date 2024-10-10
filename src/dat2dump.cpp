#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file\"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

void convertFile(const string& infile, const string& outfile) {
  string str, str2; // tmp variables

  ifstream fin(infile.c_str());
  checkFileStream(fin, infile);

  ofstream fout(outfile.c_str());
  checkFileStream(fout, outfile);

  // TODO: Auto-populate the variables based on the comment after atoms or a format string in the first line
  while (getline(fin, str)) {
    stringstream ss(str);
    ss >> str2;
    if (str2.compare("Atoms") == 0) break;
  }
  getline(fin, str); // blank line after "Atoms" specification

  while (getline(fin,str)) {
    if (all_of(str.begin(), str.end(), ::isspace)) { // Check if we've reached the end of the atoms
      break;
    }
    fout << str << "\n";
  }
}

int main(int argc, char** argv) {
  vector <string> infiles;
  string outfile;

  try {
    cxxopts::Options options(argv[0], "Converts a LAMMPS input file to a Tecplot-readable file");
    options
    .positional_help("infile outfile")
    .show_positional_help();

    options
    .add_options()
      ("i,input", "Input file", cxxopts::value<vector <string> >(infiles), "file")
      ("o,output", "Output file", cxxopts::value<string>(outfile)->default_value("*.tec"), "file")
      ("h,help", "Show the help");
    
    options.parse_positional({"input", "output"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !(result.count("input"))) {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (infiles.size() > 1 && result.count("output")) {
      cout << "WARNING: output file will only be for the last file!\n";
    }

    if (result.count("input")) {
      for (unsigned int i = 0; i < infiles.size(); ++i) {
        if (!(result.count("output"))) {
          size_t end_pos = infiles[i].find(".dat");
          if (end_pos == string::npos) {
            end_pos = infiles[i].find(".lmp");
          }
          if (end_pos == string::npos) {
            cerr << "Unrecognized input file format. Expected '*.dat' or '*.lmp'\n";
            return FILE_NAME_ERROR;
          }
          outfile = infiles[i].substr(0, end_pos) + ".tec";
        }
        convertFile(infiles[i], outfile);
      }
    }
  } catch (const cxxopts::OptionException& e) {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}
/******************************************************************************
* This program parses output from LAMMPS minimized structures and formats it
* so LAMMPS can run a simulation with the results.
*******************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cxxopts.hpp>
#include "atom.h"
#include "error_code_defines.h"

using namespace std;

struct dataValues
{
  int N, n_types = 0; // number of atoms, number of atom types
  double xlow, xhigh, ylow, yhigh, zlow, zhigh, xy, xz, yz; // box bounds and tilt factors
  vector <string> data_labels;
  vector <Atom> atoms;
  bool using_charge = false;
  bool has_tilt = false;
};

template <typename T>
void checkFileStream(T& stream, const string& file)
{
  if (stream.fail())
  {
    cout << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

dataValues readData(const string& infile)
{
  string str; // junk string variable
  dataValues data;
  int id = 0, type = 0, n_total = 0;
  double charge = 0.0, x = 0.0, y = 0.0, z = 0.0;

  ifstream fin(infile.c_str());
  checkFileStream(fin, infile);

  for (int i = 0; i < 3; ++i) {getline(fin, str);} // gets some header info we don't need
  fin >> data.N;

  data.atoms.resize(data.N, Atom());

  fin.ignore();
  getline(fin, str); // gets the box bounds info line
  if (str.find("xy") == string::npos)
  {
    fin >> data.xlow >> data.xhigh >> data.ylow >> data.yhigh >> data.zlow >> data.zhigh;
  }
  else
  {
    data.has_tilt = true;
    fin >> data.xlow >> data.xhigh >> data.xy >> data.ylow >> data.yhigh >> data.xz >> data.zlow >> data.zhigh >> data.yz;
  }

  fin.ignore();
  getline(fin, str); // gets the line "ITEM: ATOMS <atom data labels>"
  stringstream ss(str);
  ss >> str >> str;
  while (ss >> str) {data.data_labels.push_back(str);} // get the labels

  // Now get the atom data
  while (getline(fin, str))
  {
    stringstream ss2(str);
    for (int i = 0; i < data.data_labels.size(); ++i)
    {
      if (data.data_labels[i] == "id") {ss >> id;}
      else if (data.data_labels[i] == "type")
      {
        ss >> type;
        if (type > data.n_types) {data.n_types = type;}
      }
      else if (data.data_labels[i] == "q") {ss >> charge; data.using_charge = true;}
      else if (data.data_labels[i] == "x") {ss >> x;}
      else if (data.data_labels[i] == "y") {ss >> y;}
      else if (data.data_labels[i] == "z") {ss >> z;}
      else {continue;}
    }

    data.atoms[id - 1] = Atom(id, type, charge, x, y, z);
    ++n_total;
  }

  fin.close();

  if (n_total != data.N)
  {
    cout << "Error: " << n_total << " atoms were read out of " << data.N << " atoms in the file.\n";
    exit(ATOM_COUNT_ERROR);
  }

  return data;
}

void writeData(const string& outfile, const string& chem_formula, const dataValues& data)
{
  ofstream fout(outfile.c_str());
  checkFileStream(fout, outfile);

  fout << "These " << chem_formula << " coordinates are from a LAMMPS dump file: [id type ";
  if (find(data.data_labels.begin(), data.data_labels.end(), "q") != data.data_labels.end())
  {
    fout << "charge ";
  }
  fout << "x y z]\n\n"
       << data.N << "      atoms\n"
       << data.n_types << " atom types";

  fout.precision(6);
  if (data.has_tilt)
  {
    fout << data.xlow << " " << data.xhigh << " xlo xhi\n"
         << data.ylow << " " << data.yhigh << " ylo yhi\n"
         << data.zlow << " " << data.zhigh << " zlo zhi\n"
         << data.xy << " " << data.xz << " " << data.yz << " xy xz yz\n\n";
  }
  else
  {
    fout << data.xlow << " " << data.xhigh << " xlo xhi\n"
         << data.ylow << " " << data.yhigh << " ylo yhi\n"
         << data.zlow << " " << data.zhigh << " zlo zhi\n\n";
  }

  fout << "Atoms\n\n";

  for (unsigned int i = 0; i < data.atoms.size(); ++i)
  {
    fout.precision(0);
    fout << data.atoms[i].getId() << " " << data.atoms[i].getType() << " ";
    if (data.using_charge)
    {
      fout.precision(1);
      fout << data.atoms[i].getCharge() << " ";
    }
    fout.precision(6);
    fout << data.atoms[i].getX() << " " << data.atoms[i].getY() << " "
         << data.atoms[i].getZ() << endl;
  }

  fout.close();
}

int main(int argc, char** argv)
{
  string infile, outfile, chem_formula; // input file, output file, chemical formula
  dataValues data;

  try
  {
    cxxopts::Options options(argv[0], "Parse a LAMMPS dump file for use as a LAMMPS input data file");
    options
      .positional_help("infile outfile")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("f,file", "File to convert", cxxopts::value<string>(infile), "file")
        ("o,output", "Name of output file", cxxopts::value<string>(outfile), "file")
        ("e,element", "The chemical formula of the system", cxxopts::value<string>(chem_formula)->default_value("<none specified>"), "element")
        ("h,help", "Show the help");

    options.parse_positional({"file", "output"});
    auto result = options.parse(argc, argv);

    if (result.count("help") || !(result.count("file")))
    {
      cout << options.help() << endl;
      return EXIT_SUCCESS;
    }

    if (!(result.count("output")))
    {
      outfile = infile.substr(0, infile.rfind(".")) + ".dat"; // strips the last period off the original filename, and adds ".dat"
    }

    if (result.count("file"))
    {
      data = readData(infile);
      writeData(outfile, chem_formula, data);
    }

  }
  catch (const cxxopts::OptionException& e)
  {
    cout << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include <CLI11.hpp>

#include "atom.h"
#include "error_code_defines.h"

using namespace std;

struct FileData {
  vector <Atom> atoms;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  unsigned int timestep, N;
  int x_idx = -1, y_idx = -1, z_idx = -1;
  int id_idx = -1, type_idx = -1, charge_idx = -1;
};

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(FILE_OPEN_ERROR);
  }
}

string join(const vector<string>& strs, const string& separator = " ") {
  string result = "";
  if (strs.size() == 0) return result;
  if (strs.size() == 1) return strs[0];
  result = strs[0];
  for (size_t i = 1; i < strs.size(); ++i) result += separator + strs[i];
  return result;
}

// Taken shamelessly from https://stackoverflow.blog/2019/10/11/c-creator-bjarne-stroustrup-answers-our-top-five-c-questions/
// pointed there by https://stackoverflow.com/q/14265581, comment by Wais Kamal Mar 6 at 7:40
// the split() function is also taken from there
template <typename Delim>
string get_word(istream& ss, Delim d) {
  string word;
  for (char ch; ss.get(ch); ) { // skip delimiters
    if (!d(ch)) {
      word.push_back(ch);
      break;
    }
  }
  for (char ch; ss.get(ch); ) { // collect word
    if (!d(ch))
      word.push_back(ch);
    else
      break;
  }
  return word;
}

vector <string> split(const string& str, const string& delim = " ") {
  vector <string> words;
  stringstream ss(str);
  auto del = [&](char ch) {for (auto x : delim) if (x == ch) return true; return false;};
  for (string w; (w = get_word(ss, del)) != ""; ) words.push_back(w);
  return words;
}

// Taken from https://stackoverflow.com/a/447307
bool isFloat(const string& str) {
  istringstream ss(str);
  double f;
  ss >> noskipws >> f; // leading whitespace is _not_ ignored
  return ss.eof() && !ss.fail();
}

FileData parseFile(const string& file) {
  FileData data;
  string str, delim = ": ";
  vector <string> strs;
  string err_msg = "File " + file + " is not a valid LAMMPS dump file\n";

  ifstream fin(file.c_str());
  checkFileStream(fin, file);

  getline(fin, str); // Expected: ITEM: TIMESTEP
  strs = split(str, delim);
  if (strs.size() != 2) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  if (strs[0].compare("ITEM") != 0) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  if (strs[1].compare("TIMESTEP") != 0) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }

  getline(fin, str); // Expected: some integer
  strs = split(str, delim);
  if (strs.size() != 1) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  if (!isFloat(strs[0])) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  stringstream ss(strs[0]);
  ss >> data.timestep;

  getline(fin, str); // Expected: ITEM: NUMBER OF ATOMS
  strs = split(str, delim);
  if (strs.size() != 4) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  if (strs[0].compare("ITEM") != 0) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  if (join({strs[1], strs[2], strs[3]}).compare("NUMBER OF ATOMS") != 0) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }

  getline(fin, str); // Expected: some integer
  strs = split(str, delim);
  if (strs.size() != 1) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  if (!isFloat(strs[0])) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  ss.str(string()); ss.clear(); ss << strs[0];
  ss >> data.N;
  data.atoms.resize(data.N, Atom());

  getline(fin, str); // Expected: ITEM: BOX BOUNDS xy* xz* yz* xx yy zz, marked strings are optional
  strs = split(str, delim);
  if (strs.size() != 6 && strs.size() != 9) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  if (join({strs[1], strs[2]}).compare("BOX BOUNDS") != 0) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }

  for (unsigned int i = 0; i < 3; ++i) {
    getline(fin, str); // Expected: 2 or 3 floats
    strs = split(str, delim);
    if (strs.size() != 2 && strs.size() != 3) {
      cout << err_msg;
      exit(FILE_FORMAT_ERROR);
    }
    if (!isFloat(strs[0]) || !isFloat(strs[1])) {
      cout << err_msg;
      exit(FILE_FORMAT_ERROR);
    }
    if (i == 0) {
      ss.str(string()); ss.clear(); ss << strs[0];
      ss >> data.xlo;
      ss.str(string()); ss.clear(); ss << strs[1];
      ss >> data.xhi;
    } else if (i == 1) {
      ss.str(string()); ss.clear(); ss << strs[0];
      ss >> data.ylo;
      ss.str(string()); ss.clear(); ss << strs[1];
      ss >> data.yhi;
    } else { // i == 2
      ss.str(string()); ss.clear(); ss << strs[0];
      ss >> data.zlo;
      ss.str(string()); ss.clear(); ss << strs[1];
      ss >> data.zhi;
    }
  }

  getline(fin, str); // Expected: ITEM: ATOMS <data_types>
  strs = split(str, delim);
  int n_data = strs.size() - 2; // The number of data types
  if (strs[0].compare("ITEM") != 0) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  if (strs[1].compare("ATOMS") != 0) {
    cout << err_msg;
    exit(FILE_FORMAT_ERROR);
  }
  for (size_t i = 2; i < strs.size(); ++i) {
    if (strs[i].compare("x") == 0) data.x_idx = i - 2;
    else if (strs[i].compare("y") == 0) data.y_idx = i - 2;
    else if (strs[i].compare("z") == 0) data.z_idx = i - 2;
    else if (strs[i].compare("xu") == 0 && data.x_idx == -1) data.x_idx = i - 2;
    else if (strs[i].compare("yu") == 0 && data.y_idx == -1) data.y_idx = i - 2;
    else if (strs[i].compare("zu") == 0 && data.z_idx == -1) data.z_idx = i - 2;
    else if (strs[i].compare("xs") == 0 && data.x_idx == -1) data.x_idx = i - 2;
    else if (strs[i].compare("ys") == 0 && data.y_idx == -1) data.y_idx = i - 2;
    else if (strs[i].compare("zs") == 0 && data.z_idx == -1) data.z_idx = i - 2;
    else if (strs[i].compare("xsu") == 0 && data.x_idx == -1) data.x_idx = i - 2;
    else if (strs[i].compare("ysu") == 0 && data.y_idx == -1) data.y_idx = i - 2;
    else if (strs[i].compare("zsu") == 0 && data.z_idx == -1) data.z_idx = i - 2;
    else if (strs[i].compare("id") == 0) data.id_idx = i - 2;
    else if (strs[i].compare("type") == 0) data.type_idx = i - 2;
    else if (strs[i].compare("q") == 0) data.charge_idx = i - 2;
  }

  if (data.x_idx == -1 || data.y_idx == -1 || data.z_idx == -1) {
    cout << "Dump file does not contain the 3D coordinates of the atoms\n";
    exit(FILE_FORMAT_ERROR);
  } else if (data.id_idx == -1 || data.type_idx == -1) {
    cout << "Dump file does not contain the id and type of the atoms\n";
    exit(FILE_FORMAT_ERROR);
  }

  while(getline(fin, str)) {
    strs = split(str, delim);
    if (strs.size() != n_data) {
      cout << err_msg;
      exit(FILE_FORMAT_ERROR);
    }
    int id = 0, type = 0;
    double x = 0.0, y = 0.0, z = 0.0, charge = 0.0;
    for (size_t i = 0; i < n_data; ++i) {
      if (!isFloat(strs[i])) {
        cout << err_msg;
        exit(FILE_FORMAT_ERROR);
      } else {
        ss.str(string()); ss.clear(); ss << strs[i];
        if (i == data.id_idx) ss >> id;
        else if (i == data.type_idx) ss >> type;
        else if (i == data.x_idx) ss >> x;
        else if (i == data.y_idx) ss >> y;
        else if (i == data.z_idx) ss >> z;
        else if (i == data.charge_idx) ss >> charge;
      }
    }
    data.atoms[id - 1] = Atom(id, type, charge, Position(x, y, z));
  }

  fin.close();
  return data;
}

map <int, pair<int,int> > calculateDistribution(const FileData& data, const char& dist_type,
  const vector <int>& types, const double& binsize, const Position& center, int nbins) {
  map <int, pair<int, int> > result;

  // Break up the domain into binsize chunks depending on dist_type
  if (nbins == 0) {
    stringstream ss;
    string equation;
    switch(dist_type) {
      case 'x': nbins = floor((data.xhi - data.xlo) / binsize);
                ss << "floor((" << round(data.xhi) << " - " << round(data.xlo) << ") / " << binsize;
                break;
      case 'y': nbins = floor((data.yhi - data.ylo) / binsize);
                ss << "floor((" << round(data.yhi) << " - " << round(data.ylo) << ") / " << binsize;
                break;
      case 'z': nbins = floor((data.zhi - data.zlo) / binsize);
                ss << "floor((" << round(data.zhi) << " - " << round(data.zlo) << ") / " << binsize;
                break;
      case 'r': {
        double dlolo = distance2D(center, Position(data.xlo, data.ylo, 0.0));
        double dhilo = distance2D(center, Position(data.xhi, data.ylo, 0.0));
        double dlohi = distance2D(center, Position(data.xlo, data.yhi, 0.0));
        double dhihi = distance2D(center, Position(data.xhi, data.yhi, 0.0));
        nbins = floor(max(max(max(dlolo,dhilo),dlohi),dhihi) / binsize);
        ss << "floor(max(" << round(dlolo) << ", " << round(dhilo) << ", "
           << round(dlohi) << ", " << round(dhihi) << ") / " << binsize << ")";
        break;
      }
      default: cout << "Error: unknown distribution type\n";
               exit(OPTION_PARSING_ERROR);
    }
    equation = ss.str();
    cout << "Using binsize of " << binsize << ", nbins = " << nbins << " (" << equation << ")\n";
  }

  for (unsigned int i = 0; i < nbins; ++i) {
    result[i].first = 0; // set each type count to 0
    result[i].second = 0; // set each total count to 0
  }

  // Now we assign each atom (if specified in types) to a bin
  int bin; // the bin we assign the atom to
  for (size_t i = 0; i < data.atoms.size(); ++i) {


    switch (dist_type) {
      case 'x' : bin = floor((data.atoms[i].getWrapped().getX() - data.xlo) / binsize);
                 break;
      case 'y' : bin = floor((data.atoms[i].getWrapped().getY() - data.ylo) / binsize);
                 break;
      case 'z' : bin = floor((data.atoms[i].getWrapped().getZ() - data.zlo) / binsize);
                 break;
      case 'r' : bin = floor(distance2D(center, data.atoms[i].getWrapped()) / binsize);
                 break;
      default : cout << "Error: unknown distribution type\n";
                exit(OPTION_PARSING_ERROR);
    }
    if (bin >= nbins) bin = nbins - 1; // Force all 'far out' values to be put in the last bin
    ++result[bin].second; // increment the total number of atoms in this bin
    if (find(types.begin(), types.end(), data.atoms[i].getType()) != types.end()) {
      ++result[bin].first; // if the type is included in the list of target types, increment that counter
    }

  }
  return result;
}

void writeData(const map<int, pair <int, int> >& histogram, const string& infile, const string& outfile, const double& binsize) {
  ofstream fout(outfile.c_str(), ios_base::app);
  checkFileStream(fout, outfile);

  fout << "# Histogram data for file " << infile << " (binsize = " << binsize << ")\n";
  fout << "# bin type_count total_count\n";
  for (auto& h : histogram) {
    fout << h.first << " " << h.second.first << " " << h.second.second << "\n";
  }
  fout << "\n";

  fout.close();
}

int main(int argc, char** argv) {
  vector <string> files; // the files to analyze
  string outfile = "distribution.txt";
  char dist_type = 'x';
  double binsize = 1;
  vector <double> center{0,0};
  vector <int> types;
  CLI::App app{"Calculate the concentration profile of the specified snapshot(s)"};
  app.option_defaults()->always_capture_default();

  app.add_option("file", files, "The snapshot(s) to calculate the concentration profile(s) of")->required()->check(CLI::ExistingFile);
  app.add_option("-t,--types", types, "The atom type(s) to examine")->required();
  app.add_option("-s,--binsize", binsize, "The binsize used in the same units as the dump file (default = 1)")->check(CLI::PositiveNumber);
  app.add_option("-o,--outfile", outfile, "The name of the output file (default: distribution.txt)");
  app.add_option("-d,--dist-type", dist_type, "The distribution to calculate: x (default), y, z, or (r)adial")->check(CLI::IsMember({"x","y","z","r"}, CLI::ignore_case));
  app.add_option("-c,--center", center, "The center point to use for determining the radial concentration (default: (0,0))")->expected(2)->check(CLI::Number); //TODO: make this work for my Position class

  CLI11_PARSE(app,argc, argv);

  FileData data;
  map <int, pair <int, int> > histogram; // first int is the bin number, the pair is (first) the target atom type count, (second) the total atom type count
  for (auto& f : files) {
    data = parseFile(f);
    histogram = calculateDistribution(data, dist_type, types, binsize, Position(center[0], center[1], 0.0), histogram.size());
    writeData(histogram, f, outfile, binsize);

    cout << "Processed file " << f << "\r";
    cout << flush;
  }

  cout << "\n";

  return EXIT_SUCCESS;
}

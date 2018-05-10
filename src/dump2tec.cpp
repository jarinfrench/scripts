#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char **argv)
{
  string filename_in, filename_out, str; // input file name, output file name

  if (argc != 2)
  {
    cout << "Please enter the filename of the dump file to convert: ";
    cin  >> filename_in;
  }
  else
  {
    filename_in = argv[1];
  }

  filename_out = filename_in.substr(0,filename_in.find(".dump")) + ".tec";

  ifstream fin(filename_in.c_str());
  if (fin.fail())
  {
    cout << "Error opening file \"" << filename_in << "\"\n";
    return 1;
  }

  ofstream fout(filename_out.c_str());
  if (fout.fail())
  {
    cout << "Error opening file \"" << filename_out << "\"\n";
    return 1;
  }

  // Dump files contain some header information in the first 9 lines
  for (int i = 0; i < 9; ++i)
  {
    getline(fin, str);
  }

  while (getline(fin,str))
  {
    fout << str << endl;
  }

  return 0;
}

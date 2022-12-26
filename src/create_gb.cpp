#include <iostream>
#include <vector>

#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

int main(int argc, char** argv)
{
  vector <int> miller_indices;
  try
  {
    cxxopts::Options options(argv[0], "Creates an input file to create GBs using Atomsk.");
    options
      .positional_help("miller_indices")
      .show_positional_help();

    options
      .allow_unrecognised_options()
      .add_options()
        ("i,index", "all 9(?) Miller indices of the grain boundary", cxxopts::value<vector <int> >(miller_indices), "miller_index")
        ("h,help", "Show the help");

    options.parse_positional({"index"});
    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      cout << options.help() << "\n";
      return EXIT_SUCCESS;
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cerr << "Error parsing options: " << e.what() << "\n";
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

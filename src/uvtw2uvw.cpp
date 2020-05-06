#include <iostream>
#include <vector>
#include <algorithm>
#include <cxxopts.hpp>
#include "error_code_defines.h"

using namespace std;

bool no_bar = false;

int gcd (const int& u, const int& v)
{
  return (v != 0) ? gcd(v, u % v) : u;
}

void uvtw2uvw(const int& u, const int& v, const int& t, const int& w)
{
  int u2, v2, w2;
  u2 = u - t;
  v2 = v - t;
  w2 = w;

  // Simplify by finding the GCF of u2, v2, and w2 - Euclid's algorithm
  int res = gcd(gcd(abs(u2),abs(v2)),abs(w));
  u2 /= res;
  v2 /= res;
  w2 /= res;

  // simplify the original uvtw
  res = gcd(gcd(gcd(abs(u),abs(v)),abs(t)),abs(w));
  int u_sim = u / res;
  int v_sim = v / res;
  int t_sim = t / res;
  int w_sim = w / res;

  if (u_sim != u)
  {
    cout << "[";
    if (no_bar)
    {
      cout << u << " " << v << " " << t << " " << w << "] --> ";
    }
    else
    {
      (u < 0) ? cout << abs(u) << "\u0305" : cout << u; cout << " ";
      (v < 0) ? cout << abs(v) << "\u0305" : cout << v; cout << " ";
      (t < 0) ? cout << abs(t) << "\u0305" : cout << t; cout << " ";
      (w < 0) ? cout << abs(w) << "\u0305" : cout << w; cout << "] --> ";
    }
  }

  cout << "[";
  if (no_bar)
  {
    cout << u_sim << " " << v_sim << " " << t_sim << " " << w_sim << "] --> ["
         << u2 << " " << v2 << " " << w2 << "]\n";
  }
  else
  {
    // overbars for negative numbers - only overbars the last number
    (u < 0) ? cout << abs(u_sim) << "\u0305" : cout << u_sim; cout << " ";
    (v < 0) ? cout << abs(v_sim) << "\u0305" : cout << v_sim; cout << " ";
    (t < 0) ? cout << abs(t_sim) << "\u0305" : cout << t_sim; cout << " ";
    (w < 0) ? cout << abs(w_sim) << "\u0305" : cout << w_sim; cout << "] --> [";
    (u2 < 0) ? cout << abs(u2) << "\u0305" : cout << u2; cout << " ";
    (v2 < 0) ? cout << abs(v2) << "\u0305" : cout << v2; cout << " ";
    (w2 < 0) ? cout << abs(w2) << "\u0305" : cout << w2; cout << "]\n";
  }
}

int main(int argc, char** argv)
{
  int u, v, t, w;
  try
  {
    cxxopts::Options options(argv[0], "Convert the four-index (Miller-Bravais) notation to three-index notation for crystal directions");
    options
    .positional_help("u v t w")
    .show_positional_help();

    options
    .allow_unrecognised_options()
    .add_options()
    ("u", "First index", cxxopts::value<int>(u), "value")
    ("v", "Second index", cxxopts::value<int>(v), "value")
    ("t", "Third index", cxxopts::value<int>(t), "value")
    ("w", "Fourth index", cxxopts::value<int>(w), "value")
    ("no-bar", "Change output to show the negative number instead of the overbar", cxxopts::value<bool>(no_bar)->default_value("false")->implicit_value("true"));

    options.parse_positional({"u","v","t","w"});
    auto result = options.parse(argc, argv);

    if (result.count("help") ||
        result.count("u") == 0 || result.count("v") == 0 || result.count("t") == 0 || result.count("w") == 0)
    {
      cout << options.help() << endl;
      cout << "If you are trying to pass in negative indices, use the following format:\n"
           << "uvtw2uvw -- <u> <v> <t> <w>\n\n";
      return EXIT_SUCCESS;
    }
    else
    {
      uvtw2uvw(u,v,t,w);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

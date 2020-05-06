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

void uvw2uvtw(const int& u, const int& v, const int& w)
{
  int u2_num, v2_num, w2, t_num, den = 3, u2, v2, t;
  u2_num = 2 * u - v;
  v2_num = 2 * v - u;
  t_num = -(u2_num + v2_num);
  w2 = w;

  if (u2_num % 3 != 0 && v2_num % 3 != 0 && t_num % 3 != 0)
  {
    u2 = u2_num;
    v2 = v2_num;
    t = t_num;
    w2 *= den;
  }
  else
  {
    u2 = u2_num / den;
    v2 = v2_num / den;
    t = u2 + v2;
  }

  // Simplify by finding the GCF of u2, v2, t, and w2 - Euclid's algorithm
  int res = gcd(gcd(gcd(abs(u2),abs(v2)),abs(t)),abs(w));
  u2 /= res;
  v2 /= res;
  t /= res;
  w2 /= res;

  // simplify the original uvw
  res = gcd(gcd(abs(u),abs(v)),abs(w));
  int u_sim = u / res;
  int v_sim = v / res;
  int w_sim = w / res;

  if (u_sim != u)
  {
    cout << "[";
    if (no_bar)
    {
      cout << u << " " << v << " " << w << "] --> ";
    }
    else
    {
      (u < 0) ? cout << abs(u) << "\u0305" : cout << u; cout << " ";
      (v < 0) ? cout << abs(v) << "\u0305" : cout << v; cout << " ";
      (w < 0) ? cout << abs(w) << "\u0305" : cout << w; cout << "] --> ";
    }
  }

  cout << "[";
  if (no_bar)
  {
    cout << u_sim << " " << v_sim << " " << w_sim << "] --> ["
         << u2 << " " << v2 << " " << t << " " << w2 << "]\n";
  }
  else
  {
    // overbars for negative numbers - only overbars the last number
    (u < 0) ? cout << abs(u_sim) << "\u0305" : cout << u_sim; cout << " ";
    (v < 0) ? cout << abs(v_sim) << "\u0305" : cout << v_sim; cout << " ";
    (w < 0) ? cout << abs(w_sim) << "\u0305" : cout << w_sim; cout << "] --> [";
    (u2 < 0) ? cout << abs(u2) << "\u0305" : cout << u2; cout << " ";
    (v2 < 0) ? cout << abs(v2) << "\u0305" : cout << v2; cout << " ";
    (t < 0) ? cout << abs(t) << "\u0305" : cout << t; cout << " ";
    (w2 < 0) ? cout << abs(w2) << "\u0305" : cout << w2; cout << "]\n";
  }
}

int main(int argc, char** argv)
{
  int u, v, w;
  try
  {
    cxxopts::Options options(argv[0], "Convert the three-index notation to Miller-Bravais notation (four index) for crystal directions");
    options
    .positional_help("u v w")
    .show_positional_help();

    options
    .allow_unrecognised_options()
    .add_options()
    ("u", "First index", cxxopts::value<int>(u), "value")
    ("v", "Second index", cxxopts::value<int>(v), "value")
    ("w", "Third index", cxxopts::value<int>(w), "value")
    ("no-bar", "Change output to show the negative number instead of the overbar", cxxopts::value<bool>(no_bar)->default_value("false")->implicit_value("true"));

    options.parse_positional({"u","v","w"});
    auto result = options.parse(argc, argv);

    if (result.count("help") ||
        result.count("u") == 0 || result.count("v") == 0 || result.count("w") == 0)
    {
      cout << options.help() << endl;
      cout << "If you are trying to pass in negative indices, use the following format:\n"
           << "uvw2uvtw -- <u> <v> <w>\n\n";
      return EXIT_SUCCESS;
    }
    else
    {
      uvw2uvtw(u,v,w);
    }
  }
  catch (const cxxopts::OptionException& e)
  {
    cerr << "Error parsing options: " << e.what() << endl;
    return OPTION_PARSING_ERROR;
  }

  return EXIT_SUCCESS;
}

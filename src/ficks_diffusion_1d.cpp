#include <iostream>
#include <iomanip> // for setprecision
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace std;

int ENDTIME = 3000;
double C0 = 0.0; // initial concentration across the domain
double CS = 1.0; // concentration (constant) at one edge of the domain
double DIFFUSION = 1; // diffusion constant
double DX = 1.0; // Size of the step between grid points
double DX_SQ = DX*DX;
double DT = .001; // Timestep size
int TIMESTEPS = ENDTIME / DT;
int GRID_X = 100; // Number of grid points. The smaller this value, the less space is simulated
int NCHECK = TIMESTEPS / 100;

typedef vector <double> Field;

template <typename T>
void checkFileStream(T& stream, const string& file) {
  if (stream.fail()) {
    cerr << "Error opening file \"" << file << "\"\n";
    exit(2);
  }
}

string createProgressString(const int& point_progress){
  stringstream ss;

  for (unsigned int i = 0; i < point_progress; ++i) ss << ".";
  for (unsigned int i = 0; i < 70 - point_progress; ++i) ss << " ";

  return ss.str();
}

// This uses a 5-point stencil of {-2, -1, 0, 1, 2}, non-periodic boundary conditions
// the edge cases were separately calculated  using stencils of:
//   {0,1,2} for result[0] --> This is actually ignored as a boundary condition
//   {-1,0,1,2} for result[1]
//   {-2,-1,0,1} for result[s-2] (where s represents the end of the field, and thus out of bounds)
//   {-2,-1,0} for result[s-1]
//   function calculations were provided via https://web.media.mit.edu/~crtaylor/calculator.html
//   The case of using a stencil of {-1,0,1} was independently verified.
Field SecondOrderDerivative(Field field, double multiplier = 1.0) {
  unsigned int s = field.size();
  double denom = multiplier / (12 * DX_SQ);
  Field result(s,0.0);

  result[0] = (field[0] - 2 * field[1] + field[2]) / DX_SQ * multiplier;
  result[1] = (field[0] - 2 * field[1] + field[2]) / DX_SQ * multiplier;
  result[s - 2] = (field[s - 3] - 2 * field[s - 2] + field[s - 1]) / DX_SQ * multiplier;
  result[s - 1] = (field[s - 3] - 2 * field[s - 2] + field[s - 1]) / DX_SQ * multiplier;
  for (size_t i = 2; i < s - 2; ++i) {
    result[i] = (-field[i - 2] + 16 * field[i - 1] - 30 * field[i] + 16 * field[i + 1] - field[i + 2]) * denom;
  }
  return result;
}

// Similarly to SecondOrderDerivative, this uses the same five-point stencil, non-PBC.
// Edge cases were calculated separately (same stencils as above)
Field FirstOrderDerivative(Field field, double multiplier = 1.0) {
  unsigned int s = field.size();
  double denom = multiplier / (12 * DX);
  Field result(s,0.0);

  result[0] = (-3 * field[0] + 4 * field[1] - field[2]) / (2 * DX) * multiplier;
  result[1] = (-2 * field[0] - 3 * field[1] + 6 * field[2] - field[3]) / (6 * DX) * multiplier;
  result[s - 2] = (field[s - 4] - 6 * field[s - 3] + 3 * field[s - 2] + 2 * field[s - 1]) / (6 * DX) * multiplier;
  result[s - 1] = (field[s - 3] - 4 * field[s - 2] + 3 * field[s - 1]) / (2 * DX) * multiplier;

  for (size_t i = 2; i < s - 2; ++i) {
    result[i] = (field[i - 2] - 8 * field[i - 1] + 8 * field[i + 1] - field[i + 2]) * denom;
  }

  return result;
}

void printConcentration(const Field& field, const string& filename) {
  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);

  for (size_t i = 0; i < field.size(); ++i) fout << field[i] << "\n";

  fout.close();
}

int main(int argc, char** argv) {
  // generate the initial condition field of the initial concentration at t = 0
  Field concentration_2ndOrder(GRID_X, C0); // u(x,t)
  Field concentration_1stOrder(GRID_X, C0); // u(x,t)
  concentration_2ndOrder[0] = CS; // initial conditions
  concentration_1stOrder[0] = CS;
  Field flux(GRID_X, 0.0); // u_x(x,t)
  Field old_2ndOrder(GRID_X, C0);
  Field old_1stOrder(GRID_X, C0);
  Field du_dt_2ndOrder(GRID_X, 0.0); // u_t(x,t)
  Field du_dt_1stOrder(GRID_X, 0.0); // u_t(x,t)

  for (unsigned int a = 1; a <= TIMESTEPS; ++a) {
    // For the case of the second order, we can just extract the results directly from the function
    du_dt_2ndOrder = SecondOrderDerivative(concentration_2ndOrder, DIFFUSION);

    // For the case of 2 first order derivatives, we need to finangle a little bit
    flux = FirstOrderDerivative(concentration_1stOrder, -DIFFUSION);
    du_dt_1stOrder = FirstOrderDerivative(flux, -1);

    // move forward one timestep
    for (size_t i = 1; i < GRID_X; ++i) {
      concentration_2ndOrder[i] = old_2ndOrder[i] + du_dt_2ndOrder[i] * DT;
      concentration_1stOrder[i] = old_1stOrder[i] + du_dt_1stOrder[i] * DT;
    }

    if (a % NCHECK == 0) {
      string progress = createProgressString((int)(a / (double)(TIMESTEPS) * 70.0));
      cout << "\r" << progress
           << " (" << setprecision(1) << fixed
           << (double)(a) / TIMESTEPS * 100 << "%)"
           << flush;
      stringstream ss;
      string next_filename;
      ss << "out_" << a;
      ss >> next_filename;
      printConcentration(concentration_2ndOrder, next_filename + "_2ndOrder.txt");
      printConcentration(concentration_1stOrder, next_filename + "_1stOrder.txt");
      if (a == TIMESTEPS) cout << "\n";
    }
    old_2ndOrder = concentration_2ndOrder;
    old_1stOrder = concentration_1stOrder;
  }
}

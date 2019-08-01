#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <random> // 
#include <unistd.h>
#include <numeric> // iota
#include <algorithm> // sort
#include <omp.h> // for omp_get_cancellation

#include "error_code_defines.h"
#include "position.h"

using namespace std;

#define PI 3.141592653589793

// global variables - read in (or determined) by an input file
int GRID_X; // number of grid points in x
int GRID_Y; // number of grid points in y
double DX; // size of the step in each direction
double MAX_X; // maximum x coordinate (DX * GRID_X)
double MAX_Y; // maximum y coordinate (DX * GRID_Y)
unsigned int SEED; // random number generator seed
unsigned int N_ETA; // number of order parameters (~36 gives good results)
unsigned int NUMSTEPS; // number of timesteps to run
unsigned int NSKIP; // number of timesteps to run before outputting a file
int NCHECK; // number of timesteps to run before giving the user a progress update
double DT; // timestep
int EQUILIBRIUM = 500; // The number of timesteps for the sharp interface to equilibrate
string mobility_file; // the name of the txt file containing the L matrix

struct Field
{
  vector <vector <double> > values;
  Position center; // original center point
  bool is_active; // whether or not the field is active
  
  // constructor
  Field() 
  {
    values.assign(GRID_X, vector <double> (GRID_Y,0.0));
    is_active = true;
  }
};

struct MultiJunction
{
  vector <unsigned int> active_etas; // list of active etas for this junction
  vector <Position> points; // list of points that have the same active etas within the same general location
  Position position; // The position of the multi junction, calculated as an average of points
};

bool operator == (const MultiJunction& lhs, const MultiJunction& rhs)
{
  if (lhs.active_etas.size() != rhs.active_etas.size()) {return false;}
  for (unsigned int i = 0; i < lhs.active_etas.size(); ++i)
  {
    if ((lhs.active_etas)[i] != rhs.active_etas[i]) {return false;}
  }
  
  if (lhs.points.size() != rhs.points.size()) {return false;}
  for (unsigned int i = 0; i < lhs.points.size(); ++i)
  {
    if ((lhs.points)[i] != rhs.points[i]) {return false;}
  }
  
  if (lhs.position != rhs.position) {return false;}
  
  return true;
}

bool operator != (const MultiJunction& lhs, const MultiJunction& rhs)
{
  return !(lhs == rhs);
}

struct JunctionNeighbors
{
  MultiJunction junction;
  vector <MultiJunction> neighbors;
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

double anInt(double x)
{
  int temp; // temporary variable to hold the integer value of x
  // round to the largest absolute distance from zero
  if (x > 0.0) x += 0.5;
  if (x < 0.0) x -= 0.5;
  temp = (int)(x);
  return (double)(temp);
}

Position PBCs(Position p)
{
  double x = p.getX();
  double y = p.getY();
  
  x = x - anInt(x / MAX_X) * MAX_X;
  y = y - anInt(y / MAX_Y) * MAX_Y;
  
  return Position (x,y);
}

bool clockwiseSort(Position a, Position b /*, Position center*/)
{
  // see https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
  // WARNING: points a and b should have _already_ had the center subtracted from it!
  if (a.getX() /*- center.getX()*/ >= 0)
  {  
    if (b.getX() /*- center.getX()*/ < 0) {return true;}
    else {return false;}
  }
  if (abs(a.getX() /*- center.getX()*/) < 1.0e-8 && abs(b.getX() /*-center.getX()*/) < 1.0e-8)
  {
    return a.getY() > b.getY();
  }
  
  double det = (a.getX() /*- center.getX()*/) * (b.getY() /*- center.getY()*/) - (b.getX() /*- center.getX()*/) * (a.getY() /*- center.getY()*/);
  if (det < -1.0e-8) {return true;}
  else if (det > 1.0e-8) {return false;}
  
  // Points are on the same line, so check which is closer to the center
  double d1 = (a.getX() /*- center.getX()*/) * (a.getX() /*- center.getX()*/) + (a.getY() /*- center.getY()*/) * (a.getY() /*- center.getY()*/);
  double d2 = (b.getX() /*- center.getX()*/) * (b.getX() /*- center.getX()*/) + (b.getY() /*- center.getY()*/) * (b.getY() /*- center.getY()*/);
  return d1 > d2;
}

vector<size_t> sort_indices (vector<Position>& v)
{
  // From https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
  // Initialize original index locations
  vector <size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0); // a vector of 0 to v.size() - 1 in increments of 1
  
  // sort indices based on comparing values in v
  sort(idx.begin(), idx.end(),
      [&v](size_t i1, size_t i2) {return clockwiseSort (v[i1], v[i2]);});
  
  sort(v.begin(), v.end(), clockwiseSort);
  
  return idx;
}

Position averagePositions(const vector <Position>& positions)
{
  double x_avg = 0.0;
  double y_avg = 0.0;
  
  for (size_t i = 0; i < positions.size(); ++i)
  {
    x_avg += positions[i].getX();
    y_avg += positions[i].getY();
  }
  
  x_avg /= positions.size();
  y_avg /= positions.size();
  
  // This applies PBCs
  x_avg = (x_avg >= GRID_X) ? x_avg - GRID_X : x_avg;
  y_avg = (y_avg >= GRID_Y) ? y_avg - GRID_Y : y_avg;
  x_avg = (x_avg <= -1) ? x_avg + GRID_X : x_avg;
  y_avg = (y_avg <= -1) ? y_avg + GRID_Y : y_avg;
  
  // WARNING: Viewing the field i.e. in python (using meshgrid for the x and y data) _switches_ the x and y axes
  // plotting these points on the same plot will _not_ match up unless you switch the x and y!
  return Position (x_avg * DX, y_avg * DX); 
}

vector <Field> getFieldFromFile(const string& file)
{
  string str;
  vector <Field> etas (N_ETA, Field());
  ifstream fin(file.c_str());
  checkFileStream(fin, file);
  
  fin.close();
  
  return etas;
}

vector <Position> getCentroidsFromFile(const string& file)
{
  string str;
  int num_grains = 0;
  vector <Position> centroids (N_ETA, Position());
  ifstream fin(file.c_str());
  checkFileStream(fin, file);
  
  getline(fin, str); // first two lines are comment lines
  getline(fin, str);
  
  while (num_grains < N_ETA)
  {
    double x, y;
    getline(fin, str);
    stringstream ss(str);
    ss >> x >> y; // get the x and y positions of the centroids
    
    Position center(MAX_X * x, MAX_Y * y);
    centroids[num_grains++] = center;
  }
  
  fin.close();
  return centroids;
}

void assignGridToEta(vector <Field>& etas)
{
  for (unsigned int i = 0; i < GRID_X; ++i)
  {
    for (unsigned int j = 0; j < GRID_Y; ++j)
    {
      int min_id = -1; // invalid id number
      double min = 1e12; // some arbitrarily large number
      
      for (unsigned int k = 0; k < N_ETA; ++k)
      {
        // find the distance to each point
        double x = i * DX - etas[k].center.getX();
        double y = j * DX - etas[k].center.getY();
        
        x = x - anInt(x / MAX_X) * MAX_X;
        y = y - anInt(y / MAX_Y) * MAX_Y;
        
        double distance = (x * x) + (y * y);
        
        if (distance < min)
        {
          min = distance;
          min_id = k;
        }
      }

      etas[min_id].values[i][j] = 1.0;
    }
  }
}

// Voronoi tesselation initial condition.  Centers of the voronoi points are
// specified randomly.  Each grid point is compared to the voronoi centers,
// and assigned to the one nearest to them.  Current implementation is not very
// efficient - O(n^2)*O(m), where n is the number of grid points in one direction,
// and m is the number of order parameters.
void generateICRandom(vector <Field>& etas)
{
  mt19937 rng(SEED); // seed the random number generator
  // create the generator itself - uniform integer distribution between 0 and GRID_X - 1
  // We use GRID_X - 1 because we only want GRID_X values
  uniform_real_distribution<double> dist(0.0, 1.0); 
  
  // Create the voronoi centers
  for (unsigned int i = 0; i < N_ETA; ++i)
  {
    int x = dist(rng) * MAX_X;
    int y = dist(rng) * MAX_Y;
    etas[i].center = Position (x,y);
  }
  
  // now assign each grid point to a specific eta.
  assignGridToEta(etas);
}

void generateICFluid(vector <Field>& etas)
{
  // random fluid
  mt19937 rng(SEED); // seed the random number generator
  uniform_real_distribution<double> dist(-1.0, 1.0);
  for (unsigned int i = 0; i < GRID_X; ++i)
  {
    for (unsigned int j = 0; j < GRID_Y; ++j)
    {
      for (unsigned int k = 0; k < N_ETA; ++k)
      {
        etas[k].values[i][j] = dist(rng) * 0.001;
      }
    }
  }
}

void generateICHex(vector <Field>& etas)
{
  unsigned int sqrt_num_centroids = sqrt(N_ETA);
  if (sqrt_num_centroids != sqrt(N_ETA))
  {
    cout << "Number of centroids must be a perfect square.\n";
    exit(INPUT_FORMAT_ERROR);
  }
  
  mt19937 rng(SEED); // seed the random number generator
  uniform_real_distribution<double> dist(-1.0, 1.0);
  int num_centroids = 0;
  for (unsigned int i = 0; i < sqrt_num_centroids; ++i)
  {
    for (unsigned int j = 0; j < sqrt_num_centroids; ++j)
    {
      for (unsigned int k = 0; k < 1; ++k)
      {
        double x = ((double)(i) / sqrt_num_centroids + (0.5 / sqrt_num_centroids * (j % 2)) + 0.5 / sqrt_num_centroids + dist(rng) * 0.1) * MAX_X;
        double y = ((double)(j) / sqrt_num_centroids + (0.5 / sqrt_num_centroids * (k % 2)) + dist(rng) * 0.1) * MAX_Y;
        etas[num_centroids].center = Position(x,y);
        ++num_centroids;
      }
    }
  }
  
  // now assign each grid point to a specific eta.
  assignGridToEta(etas);
}

void generateVoronoiIC(vector <Field>& etas, const vector <Position>& centroids)
{
  if (centroids.size() != etas.size())
  {
    cout << "Error: number of centroids must match the number of grains.\n";
    exit(INPUT_FORMAT_ERROR);
  }
  else
  {
    for (size_t i = 0; i < centroids.size(); ++i)
    {
      etas[i].center = centroids[i];
    }
  }
  
  // now assign each grid point to a specific eta.
  assignGridToEta(etas);
}

// Prints the eta field via the equation: phi(x,y) = sum(eta[i]^2)
// Should give a maximum value of 1 inside a grain, and less than 1 at the grain
// boundaries.  use_grain_num specifies whether the grain number will change the
// output values (useful for initial conditions that use a sharp interface)
void printField(const vector <Field> & etas, const string& filename, const bool& use_grain_num = false)
{
  ofstream fout(filename.c_str());
  checkFileStream(fout, filename);
  
  Field field; // the overall order parameter field
  
  
  for (unsigned int j = 0; j < GRID_X; ++j)
  {
    for (unsigned int k = 0; k < GRID_Y; ++k)
    {
      for (unsigned int i = 0; i < N_ETA; ++i)
      {
        if (!etas[i].is_active) {continue;}
        if (use_grain_num)
        {
          //cout << etas[i].values[j][k] << " is the value of eta " << i << " at " << j << "," << k << endl;
          field.values[j][k] += etas[i].values[j][k] * etas[i].values[j][k] * (i + 1);
        }
        else
        {
          // Eq. (13) from Fan and Chen.
          field.values[j][k] += etas[i].values[j][k] * etas[i].values[j][k];
        }
      }
      
      fout << field.values[j][k] << " ";
    }
    fout << "\n";
  }
  
  fout.close();
}

// Mostly for debugging purposes, this functions creates N_ETA files that
// show the values of the eta_{i} field over the whole domain
void printIndividualFields(const vector <Field>& etas, const int& step)
{
  string file;
  for (size_t i = 0; i < N_ETA; ++i)
  {
    stringstream ss;
    ss << "eta_" << i << "_step_" << step << ".txt";
    ss >> file;
    ofstream fout(file.c_str());
    checkFileStream(fout, file);
    for (unsigned int j = 0; j < GRID_X; ++j)
    {
      for (unsigned int k = 0; k < GRID_Y; ++k)
      {
        fout << etas[i].values[j][k] << " ";
      }
      fout << "\n";
    }
    fout.close();
  }
}

// find the updated center point of the order parameter field
Position findCenter(const Field& eta)
{
  // This is taken from Bai and Breen, Journal of Graphics Tools 13(4) (2008) 53-60
  struct Tube
  {
    vector <double> xs; 
    vector <double> ys; 
    vector <double> zs; 
    vector <double> thetas;
    vector <double> weights;
    double radius = (double)(GRID_X) / (2.0 * PI);
    vector <double> com;

    Tube()
    {
      xs.resize(GRID_X*GRID_Y, 0.0);
      ys.resize(GRID_X*GRID_Y, 0.0);
      zs.resize(GRID_X*GRID_Y, 0.0);
      thetas.resize(GRID_X*GRID_Y, 0.0);
      weights.resize(GRID_X*GRID_Y, 0.0);
      com.resize(3,0.0);
    }
  };
  
  Tube x;
  Tube y;
  
  for (int i = 0; i < GRID_X; ++i)
  {
    for (int j = 0; j < GRID_Y; ++j)
    {
      unsigned int index = GRID_X * i + j;
      x.ys[index] = j; // no change for y in tube x
      y.xs[index] = i; // similarly for x in tube y
      
      x.thetas[index] = (double)(i) * 2.0 * PI / GRID_X;
      y.thetas[index] = (double)(j) * 2.0 * PI / GRID_Y;
      x.xs[index] = x.radius * cos(x.thetas[index]);
      y.ys[index] = y.radius * cos(y.thetas[index]);
      
      x.zs[index] = x.radius * sin(x.thetas[index]);
      y.zs[index] = y.radius * sin(y.thetas[index]);
      
      x.weights[index] = eta.values[i][j];
      y.weights[index] = eta.values[i][j];
    }
  }
  
  for (size_t i = 0; i < x.xs.size(); ++i)
  {
    x.com[0] += x.xs[i] * x.weights[i];
    y.com[0] += y.xs[i] * y.weights[i];
    x.com[1] += x.ys[i] * x.weights[i];
    y.com[1] += y.ys[i] * y.weights[i];
    x.com[2] += x.zs[i] * x.weights[i];
    y.com[2] += y.zs[i] * y.weights[i];
  }
  
  for (int i = 0; i < 3; ++i)
  {
    x.com[i] /= accumulate(x.weights.begin(), x.weights.end(), 0.0);
    y.com[i] /= accumulate(y.weights.begin(), y.weights.end(), 0.0);
  }
  
  double th_avg = atan2(-x.com[2], -x.com[0]) + PI;
  double x_com = GRID_X * th_avg / (2.0 * PI);
  
  th_avg = atan2(-y.com[2], -y.com[1]) + PI;
  double y_com = GRID_Y * th_avg / (2.0 * PI);

  
  // WARNING: Viewing the field e.g. in python (using meshgrid for the x and y data) _switches_ the x and y axes
  // plotting these points on the same plot will _not_ match up unless you switch the x and y!
  return Position (x_com * DX, y_com * DX);
}

/*
 * multi_junctions is a vector of the MultiJunction struct, containing a list of 
 *   active etas for each junction, as well as the position of the junction (in length units)
 * etas is the order parameter field
 * Both are marked const - we do not want any changes occuring to these as we calculate the neighboring junctions.
 */
vector <JunctionNeighbors> findMultiJunctionNeighbors(const vector <MultiJunction>& multi_junctions, const vector <Field>& etas)
{
  vector <Position> eta_centers (N_ETA, Position()); // The list of centroids of each order parameter
  map <unsigned int, pair <Position, vector <MultiJunction> > > grain_junctions, unwrapped_grain_junctions; // map of the order parameter (key) to its center and associated multi junctions (value)
  vector <JunctionNeighbors> neighbors_list (multi_junctions.size(), JunctionNeighbors()); // A list of the neighbors for each multi junction
  
  for (size_t i = 0; i < etas.size(); ++i)
  {
    // Find the centers of the order parameters, in length units
    eta_centers[i] = findCenter(etas[i]);
    vector <MultiJunction> junctions, unwrapped_junctions; // the list of multi junction points
    // Find the junctions associated with this eta
    for (size_t j = 0; j < multi_junctions.size(); ++j)
    {
      for (vector <unsigned int>::const_iterator it = multi_junctions[j].active_etas.begin();
           it != multi_junctions[j].active_etas.end();
           ++it)
      {
        if (*it == i) 
        {
          // get the distance between this junction and the center of the order parameter
          Position diff = multi_junctions[j].position - eta_centers[i];
          double distance = diff * diff;
          
          Position pbc_diff = PBCs(diff); // apply pbcs, but keep original result
          double pbc_distance = pbc_diff * pbc_diff;
          if (abs(pbc_distance - distance) < 1.0e-8) // if the periodic distance is the same, do nothing
          {
            junctions.push_back(multi_junctions[j]);
            unwrapped_junctions.push_back(multi_junctions[j]);
          }
          else // otherwise, we need to figure out which direction(s) needed pbcs.
          {
            unwrapped_junctions.push_back(multi_junctions[j]);
            double tmp_x = multi_junctions[j].position.getX();
            double tmp_y = multi_junctions[j].position.getY();
            if (abs(pbc_diff.getX() - diff.getX()) < 1.0e-8)
            {  
              if (tmp_y - eta_centers[i].getY() > 0) {tmp_y -= MAX_Y;}
              else /*(tmp_y - eta_centers[i].getY() <= 0)*/ {tmp_y += MAX_Y;}
            }
            if (abs(pbc_diff.getY() - diff.getY()) < 1.0e-8)
            {
              if (tmp_x - eta_centers[i].getX() > 0) {tmp_x -= MAX_X;}
              else /*(tmp_x - eta_centers[i].getX() <= 0)*/ {tmp_x += MAX_X;}
            }
            MultiJunction tmp = multi_junctions[j];
            tmp.position = Position (tmp_x, tmp_y);
            junctions.push_back(tmp);
          }
          break;
        }
      }
    }
    // put the info in the map
    // This can be read as: the multi junctions for order parameter i with center eta_center are listed in vector junctions.
    grain_junctions[i] = make_pair(eta_centers[i], junctions);
    unwrapped_grain_junctions[i] = make_pair(eta_centers[i], unwrapped_junctions);
  }
  
  // sort the list of junctions clockwise, starting from 12 o'clock
  for (size_t i = 0; i < grain_junctions.size(); ++i)
  {
    vector <MultiJunction> junctions = grain_junctions[i].second; // grab the vector (grain_junctions[i])
    vector <MultiJunction> unwrapped_junctions = unwrapped_grain_junctions[i].second;
    
    vector <Position> junctions_positions (junctions.size());
    
    for (size_t j = 0; j < junctions.size(); ++j)
    {
      junctions_positions[j] = junctions[j].position;
    }
    
    for (size_t j = 0; j < junctions_positions.size(); ++j) {junctions_positions[j] = PBCs(junctions_positions[j] - grain_junctions[i].first);} // subtract the centers off the points 
    vector <size_t> junction_indices = sort_indices(junctions_positions); // indices will contain the original index for the map.
    for (size_t j = 0; j < junctions_positions.size(); ++j) {junctions_positions[j] = PBCs(junctions_positions[j] + grain_junctions[i].first);} // add the centers back to the positions of the junctions.
    
    vector <MultiJunction> unsorted_junctions = unwrapped_junctions;
    for (size_t j = 0; j < junction_indices.size(); ++j)
    {
      unwrapped_junctions[j] = unsorted_junctions[junction_indices[j]];
    }
    
    grain_junctions[i].second = unwrapped_junctions;
  }
  
  // Now determine the neighboring junctions
  for (size_t i = 0; i < multi_junctions.size(); ++i)
  {  
    neighbors_list[i].junction = multi_junctions[i];
    
    // loop through the active etas at this junction
    for (size_t j = 0; j < multi_junctions[i].active_etas.size(); ++j)
    {
      // for this eta, we need to find the next junction (clockwise) from the current one.
      // this is the sorted junction list
      vector <MultiJunction> tmp = grain_junctions[multi_junctions[i].active_etas[j]].second;
       
      auto it = find(tmp.begin(), tmp.end(), multi_junctions[i]);
      if (it != tmp.end())
      {
        if (++it == tmp.end()) {it = tmp.begin();}
        size_t index = distance(tmp.begin(), it);
        neighbors_list[i].neighbors.push_back(*it);
      }
      else
      {
        cout << "ERROR: Junction not found " << j << ".\n";
        exit(JUNCTION_NOT_FOUND);
      }
    }
  }
  
  return neighbors_list;
}

vector <pair <MultiJunction, vector <double> > > calculateJunctionAngles(const vector <JunctionNeighbors>& neighbors_list)
{
  vector <pair <MultiJunction, vector <double> > > junction_angles (neighbors_list.size());
  
  for (size_t i = 0; i < neighbors_list.size(); ++i)
  {
    Position junction = neighbors_list[i].junction.position;
    junction_angles[i].first = neighbors_list[i].junction;
    
    for (size_t j = 0; j < neighbors_list[i].neighbors.size() - 1; ++j) // because we can calculate the last angle more cheaply
    {
      Position line1 = neighbors_list[i].neighbors[j].position - junction;
      unsigned int j1 = (j + 1 >= neighbors_list[i].neighbors.size()) ? 0 : j + 1;
      Position line2 = neighbors_list[i].neighbors[j1].position - junction;
      
      junction_angles[i].second.push_back(acos(line1 * line2 / (sqrt(line1 * line1) * sqrt(line2 * line2))) * 180.0 / PI);
    }
    junction_angles[i].second.push_back(360.0 - accumulate(junction_angles[i].second.begin(), junction_angles[i].second.end(), 0.0));
  }
  
  return junction_angles;
}

void calculateInfo(const vector <Field>& etas, 
  ostream& fout1, /*for grain size*/
  ostream& fout2, /*for gb lengths*/
  ostream& fout3, /*for multi junctions*/
  const int& step)
{
  vector <double> sizes(N_ETA, 0.0);
  double total_gb_length = 0.0;
  double activity_threshold = 0.12; // if an eta value is greater than or equal to this cutoff, the order parameter is active at this point.
  // the value of 0.12 was calculated based on Moelans et al.'s paper - see page 172 of lab notebook Vol. 3
  //map <pair <unsigned int, unsigned int>, double> gb_lengths; // length of the grain boundary between etas i and j
  unsigned int num_junctions = 0;
  vector <MultiJunction> multi_junctions; // the list of all multi junctions
  Field checked_grid_points; // field for multi junction checking - we only need to check each grid point once // FIXME: This may break with triple junctions within the search radius of each other.

  // initialize the individual gb lengths
  /*for (unsigned int i = 0; i < N_ETA - 1; ++i)
  {
    for (unsigned int j = i + 1; j < N_ETA; ++j)
    {
      gb_lengths[make_pair(i,j)] = 0.0;
    }
  }*/
  
  #pragma omp parallel for collapse(2) num_threads(12)
  for (unsigned int i = 0; i < GRID_X; ++i)
  {
    for (unsigned int j = 0; j < GRID_Y; ++j)
    {
      unsigned int num_active = 0;
      vector <unsigned int> active_etas;
      bool bulk_value_added = false;

      for (unsigned int k = 0; k < N_ETA; ++k)
      {
        double val = etas[k].values[i][j];

        // check for whether the eta is active at this point
        if (val >= activity_threshold)
        {
          ++num_active;
          active_etas.push_back(k);
        }

        // Check if we are in the bulk of an order parameter
        if (val >= (1.0 - activity_threshold) && !bulk_value_added)
        {
          sizes[k] += DX * DX;
          bulk_value_added = true;
        }
      } // eta loop

      if (!bulk_value_added) {total_gb_length += DX;} // add to total GB length if no bulk value added.

      // calculate individual gb lengths
      /*if (num_active >= 2)
      {
        for (size_t k = 0; k < active_etas.size() - 1; ++k)
        {
          for (size_t l = k + 1; l < active_etas.size(); ++l)
          {
            unsigned int first, second;
            if (active_etas[k] < active_etas[l])
            {
              first = active_etas[k];
              second = active_etas[l];
            }
            else
            {  // TODO: This part of the code may never be called... I can check later, and possibly remove it.
              first = active_etas[l];
              second = active_etas[k];
            }
            gb_lengths[make_pair(first,second)] += DX;
          }
        }
      }*/
      
      // get positions of multi-junctions
      if (num_active >= 3 && checked_grid_points.values[i][j] < 1.0) // if the 'checked field' shows that this multijunction grid point has not been examined
      {
        multi_junctions.push_back(MultiJunction()); // add in a blank multijunction
        for (int ii = (int)(i) - (GRID_X / 20); ii < (int)(i) + (GRID_X / 20) + 1; ++ii) // check all grid points within GRID_X / 20
        {
          for (int jj = (int)(j) - (GRID_Y / 20); jj < (int)(j) + (GRID_Y / 20) + 1; ++jj) // similary for GRID_Y / 20
          {
            // apply PBCs
            int ii_before = ii;
            int jj_before = jj;
            int ii_after, jj_after;
            ii_after = (ii >= GRID_X) ? ii - GRID_X : ii;
            jj_after = (jj >= GRID_Y) ? jj - GRID_Y : jj;
            ii_after = (ii <= -1) ? ii + GRID_X : ii_after; // because we have already checked the right boundary
            jj_after = (jj <= -1) ? jj + GRID_Y : jj_after; // because we have already checked the top boundary
            
            if (checked_grid_points.values[ii_after][jj_after] > 0) {continue;}
            for (size_t a = 0; a < active_etas.size(); ++a)
            {
              checked_grid_points.values[ii_after][jj_after] = 1;
              if (!(etas[active_etas[a]].values[ii_after][jj_after] >= activity_threshold))
              {
                break; // this grid point is not a triple junction point
              }
              else
              {
                multi_junctions[num_junctions].active_etas = active_etas;
                multi_junctions[num_junctions].points.push_back(Position(ii_before, jj_before)); // use the before value so averaging works
              }
            }
          }
        }
        
        ++num_junctions;
      }
    } // GRID_Y loop
  } // GRID_X loop
  
  for (size_t a = 0; a < multi_junctions.size(); ++a)
  {
    for (size_t b = 0; b < multi_junctions[a].active_etas.size(); ++b)
    {
      multi_junctions[a].position = averagePositions(multi_junctions[a].points);
    }
  }
  
  vector <JunctionNeighbors> neighbors_list = findMultiJunctionNeighbors(multi_junctions, etas);
  
  // Calculate the angles at the multi junction
  vector <pair <MultiJunction, vector <double> > > junction_angles = calculateJunctionAngles(neighbors_list);

  // output the grain size data to fout1
  fout1 << step;
  for (unsigned int i = 0; i < N_ETA; ++i)
  {
    fout1 << " " << sizes[i];
  }
  fout1 << endl;

  // output the gb data to fout2
/*  double total2 = 0.0;
  pair <unsigned int, unsigned int> index;
  for (unsigned int i = 0; i < N_ETA - 1; ++i)
  {
    for (unsigned int j = i + 1; j < N_ETA; ++j)
    {
      index = make_pair(i,j);
      fout2 << gb_lengths[index] << " ";
      total2 += gb_lengths[index];
    }
  }
  if (abs(total2 - total_gb_length) < 1.0e-8)
  {
    fout2 << total2 << endl;
  }
  else
  {
    fout2 << total_gb_length << " " << total2 << endl;
  }*/

  // output the multi junction data to fout3
  // step junction_position active_etas junction_neighbors(positions) junction_angles
  // Assumes that multi_junctions, neighbors_list, and junction_angles all have the same number of elements
  fout3 << step << " ";
  for (size_t a = 0; a < multi_junctions.size(); ++a)
  {
    fout3 << multi_junctions[a].position;
    for (size_t b = 0; b < multi_junctions[a].active_etas.size(); ++b)
    {
      fout3 << " " << multi_junctions[a].active_etas[b];
    }
    
    for (size_t b = 0; b < neighbors_list[a].neighbors.size(); ++b)
    {
      fout3 << " " << neighbors_list[a].neighbors[b].position;
    }
    
    for (size_t b = 0; b < junction_angles[a].second.size(); ++b)
    {
      fout3 << " " << junction_angles[a].second[b];
    }
    fout3 << endl;
    if (a + 1 != multi_junctions.size()) {fout3 << "---- ";}
  }
}

// eq (11) from Fan and Chen
// Finite difference method
void calculateLaplacian(const vector <Field>& etas, vector <Field>& laplacian)
{
  for (unsigned int i = 0; i < N_ETA; ++i)
  {
    if (!etas[i].is_active) {continue;} // if the current eta is all zeros, move to the next one
    #pragma omp parallel for collapse(2)
    for (unsigned int j = 0; j < GRID_X; ++j)
    {
      for (unsigned int k = 0; k < GRID_Y; ++k)
      {
        // PBCs
        int j1 = (j + 1 == GRID_X) ? 0 : j + 1;
        int k1 = (k + 1 == GRID_Y) ? 0 : k + 1;
        int j_neg1 = ((int)(j) - 1 == -1) ? GRID_X - 1 : j - 1;
        int k_neg1 = ((int)(k) - 1 == -1) ? GRID_Y - 1 : k - 1;
        int j2 = (j + 2 >= GRID_X) ? j + 2 - GRID_X : j + 2;
        int k2 = (k + 2 >= GRID_Y) ? k + 2 - GRID_Y : k + 2;
        int j_neg2 = ((int)(j) - 2 <= -1) ? (int)(j) - 2 + GRID_X : j - 2;
        int k_neg2 = ((int)(k) - 2 <= -1) ? (int)(k) - 2 + GRID_Y : k - 2;

        // This value gets the contribution from first nearest neighbors as follows:
        //                                 j,k1
        //                                  ^
        //                                  |
        //                    j_neg1,k <-- j,k --> j1,k
        //                                  |
        //                                  V
        //                               j,k_neg1
        double first_nn_contribution = 0.5 * (etas[i].values[j1][k] + etas[i].values[j_neg1][k] + \
          etas[i].values[j][k1] + etas[i].values[j][k_neg1] - (4 * etas[i].values[j][k]));
          
        /* This the value for second nearest neighbors
        *                                 j,k2
        *                     j_neg1,k1    ^     j1,k1
        *                              \   |    /
        *                                \ |  /
        *              j_neg2,k <-------- j,k --------> j2,k
        *                               /  | \
        *                             /    |   \
        *                j_neg1,j_neg1     V    j1,k_neg1
        *                               j,k_neg2
        */
        double second_nn_contribution = 0.25 * (etas[i].values[j2][k] + etas[i].values[j1][k1] + \
          etas[i].values[j][k2] + etas[i].values[j_neg1][k1] + etas[i].values[j_neg2][k] + \
          etas[i].values[j_neg1][k_neg1] + etas[i].values[j][k_neg2] + etas[i].values[j1][k_neg1] - \
          (8 * etas[i].values[j][k]));
        
        // Grain growth only seems to occur if I have the factor of 10.0 here (as opposed to the usual factor of 1.0)
        // NOTE: is this where temperature plays a role?
        laplacian[i].values[j][k] = 10.0 / (DX * DX) * (first_nn_contribution + second_nn_contribution);
      }
    }
  }
}

double getMatrixValue(const vector <vector <double> >& mat, const vector <Field>& etas, const unsigned int& grid_x, const unsigned int& grid_y)
{
  double top = 0.0;
  double bott = 0.0;
  double etaij_sq = 0.0;
  
  for (unsigned int i = 0; i < N_ETA; ++i)
  {
    for (unsigned int j = i + 1; j < N_ETA; ++j)
    {
      etaij_sq = etas[i].values[grid_x][grid_y] * etas[i].values[grid_x][grid_y] * etas[j].values[grid_x][grid_y] * etas[j].values[grid_x][grid_y];
      top += mat[i][j] * etaij_sq;
      bott += etaij_sq;
    }
  }
  
  if (bott == 0)
  {
    return mat[0][0];
    // cout << "Matrix value is NaN, unable to continue.\n";
    // exit(BOUNDS_ERROR);
  }
  return top / bott;
}

void checkActiveEtas(vector <Field>& etas)
{
  bool zeros;
  for (unsigned int i = 0; i < N_ETA; ++i)
  {
    if (!etas[i].is_active) {continue;}
    zeros = all_of(etas[i].values.begin(), etas[i].values.end(),
      [](vector <double> v) {return all_of(v.begin(), v.end(), 
        [](double d) {return abs(d) < 1e-8;});});
    
    if (zeros) 
    {
      etas[i].is_active = false;
      //cout << "Setting eta " << i << " to inactive.\n";
    }
  }
}

void parseInput(const string& filename)
{
  int denominator = 1000;
  ifstream fin(filename.c_str());
  checkFileStream(fin, filename);
  
  string str; // junk variable
  // format of file is:
  // var_name = var_value
  getline(fin,str); // first line is a comment line
  fin >> str >> str >> GRID_X;
  fin >> str >> str >> GRID_Y;
  fin >> str >> str >> DX;
  fin >> str >> str >> SEED;
  fin >> str >> str >> N_ETA;
  fin >> str >> str >> NUMSTEPS;
  fin >> str >> str >> NSKIP;
  fin >> str >> str >> DT;
  fin >> str >> str >> mobility_file;
  
  fin.close();
  
  MAX_X = GRID_X * DX;
  MAX_Y = GRID_Y * DX;
  NCHECK = NUMSTEPS / denominator;
  while (NCHECK == 0)
  {
    denominator /= 2;
    NCHECK = NUMSTEPS / denominator;
  }
}

string createProgressString(const int& point_progress)
{
  stringstream ss;
  
  for (unsigned int i = 0; i < point_progress; ++i)
  {
    ss << ".";
  }
  
  for (unsigned int i = 0; i < 70 - point_progress; ++i)
  {
    ss << " ";
  }

  return ss.str();
}

int main(int argc, char** argv)
{
  if (!omp_get_cancellation()) // quietly handle setting up parallelism if not specified
  {
    char str [] = "OMP_CANCELLATION=true";
    putenv(str);
    execv(argv[0], argv);
  }
  string ic;
  if (argc < 2)
  {
    cout << "Please pass the input file on the command line.\n";
    return EXIT_SUCCESS;
  }
  else
  {
    parseInput(argv[1]);
    
    if (argc >= 3)
    {
      ic = argv[2];
      if (!(ic.compare("random")   == 0 ||
            ic.compare("fluid")    == 0 ||
            ic.compare("hex")      == 0 ||
            ic.compare("centroid") == 0 ||
            ic.compare("input")    == 0))
      {
        cout << "The initial condition should be specified by 'random,' 'fluid,' 'hex,' 'centroid,' or 'input'\n";
        return INPUT_FORMAT_ERROR;
      }
    }
    else
    {
      ic = "random";
    }
  }
  
  string next_filename;
  bool warn = false;
  int progress_indicator = 0;
  string progress = "";
  
  vector <Field> etas(N_ETA, Field()), etas_old(N_ETA, Field()); // initialize the etas
  vector <Field> laplacian(N_ETA, Field());
  vector <Field> dEta_dt(N_ETA, Field());
  vector <vector <double> > _L (N_ETA, vector <double> (N_ETA, 0.0)); // mobility coefficients
  vector <vector <double> > _kappa (N_ETA, vector <double> (N_ETA, 2.0)); // energy coefficients
  
  ifstream fin(mobility_file.c_str());
  checkFileStream(fin, mobility_file);
  
  for (unsigned int i = 0; i < N_ETA; ++i)
  {
    for (unsigned int j = 0; j < N_ETA; ++j)
    {
      fin >> _L[i][j];
    }
  }
  
  if (ic.compare("random") == 0)
  {
    // random centroids
    generateICRandom(etas_old); // voronoi tessellation (inefficient)
  }
  else if (ic.compare("fluid") == 0)
  {
    generateICFluid(etas_old);
  }
  else if (ic.compare("hex") == 0)
  {
    generateICHex(etas_old);
  }
  else if (ic.compare("centroid") == 0)
  {
    if (argc >= 4)
    {
      vector <Position> centroids = getCentroidsFromFile(argv[3]);
      generateVoronoiIC(etas_old, centroids);
    }
    else
    {
      cout << "Please specify the list of centroids in a file.\n";
      return INPUT_FORMAT_ERROR;
    }
  }
  else if (ic.compare("input") == 0)
  {
    if (argc >= 4)
    {
      etas_old = getFieldFromFile(argv[3]);
    }
    else
    {
      cout << "Please specify the input file containing the eta fields.\n";
      return INPUT_FORMAT_ERROR;
    }
  }
  
  ofstream fout_grain_size("grain_size_distribution.txt");
  checkFileStream(fout_grain_size, "grain_size_distribution.txt");
  
  fout_grain_size << "# Step ";
  for (unsigned int i = 0; i < N_ETA; ++i)
  {
    fout_grain_size << "gr" << i << " ";
  }
  fout_grain_size << "\n";
  
  ofstream fout_GB("grain_boundary_lengths.txt");
  checkFileStream(fout_GB, "grain_boundary_lengths.txt");
  
  fout_GB << "# ";
  for (unsigned int i = 0; i < N_ETA - 1; ++i)
  {
    for (unsigned int j = i + 1; j < N_ETA; ++j)
    {
      fout_GB << "gr" << i << "gr" << j << "_gb ";
    }
  }
  fout_GB << "total (total2 if different)\n";
  
  ofstream fout_multi_junctions("multi_junctions.txt");
  checkFileStream(fout_multi_junctions, "multi_junctions.txt");
  fout_multi_junctions << "# Step junction_position active_etas junction_neighbors junction_angles\n";
  
  printField(etas_old, "initial_structure.txt", true);
  
  cout << "DX = " << DX << endl 
       << "DT = " << DT << endl;
  
  for (unsigned int a = 1; a < NUMSTEPS + 1; ++a)
  {
    if (all_of(etas_old.begin(), etas_old.end(), [](Field f) {return f.is_active == false;}) || 
        accumulate(etas_old.begin(), etas_old.end(), 0, [](int a, Field f) {return a + f.is_active;}) == 1)
    {
      cout << "\nNumber of active order parameters <= 1.  Ending simulation loop.\n";
      break;
    }
    
    int count = 0;
    for (unsigned int i = 0; i < N_ETA; ++i)
    {
      count += etas_old[i].is_active;
    }
    if (count == 1)
    {
      cout << "\nOnly one active order parameter.  Ending simulation loop.\n";
      stringstream ss;
      ss << "out_" << a << ".txt";
      ss >> next_filename;
      printField(etas_old, next_filename);
      break;
    }
    calculateLaplacian(etas_old, laplacian);
    
    dEta_dt.assign(N_ETA, Field()); // reset the values to 0
    for (unsigned int i = 0; i < N_ETA; ++i)
    {
      if (!etas[i].is_active) {continue;}
      #pragma omp parallel for collapse(2) num_threads(12)
      for (unsigned int j = 0; j < GRID_X; ++j) 
        for (unsigned int k = 0; k < GRID_Y; ++k)
        {
          double eta3 = etas_old[i].values[j][k] * etas_old[i].values[j][k] * etas_old[i].values[j][k];
          double boundary_term = 0.0;
          for (unsigned int l = 0; l < N_ETA; ++l)
          {
            if (l == i) {continue;}
            boundary_term += etas_old[l].values[j][k] * etas_old[l].values[j][k];
          }
          
          // Note that the 3.0 comes from assuming gamma = 1.5 in eq. (10) from Fan and Chen, Acta Mat 45 (1997) 611-622
          // original equation says 2 * gamma * eta[i]
          double L = getMatrixValue(_L, etas_old, j, k);
          double kappa = getMatrixValue(_kappa, etas_old, j, k);
          dEta_dt[i].values[j][k] += -L * \
            (eta3 - etas_old[i].values[j][k] + 3.0 * etas_old[i].values[j][k] * boundary_term \
              - kappa * laplacian[i].values[j][k]);
          
          // if (isnan(dEta_dt[i].values[j][k]))
          // {
          //   #pragma omp cancel for
          // }
          
          // move forward one timestep - eq (12) from Fan and Chen
          etas[i].values[j][k] = etas_old[i].values[j][k] + dEta_dt[i].values[j][k] * DT;
        }
        // #pragma omp cancellation point for
        // {
        //   cout << "\nError: solution unstable.  Try reducing the timestep (or increasing DX).\n";
        //   return 10;
        // }
      
    }
    
    if (((int)(a) % NCHECK) == 0)
    {
      string progress = createProgressString((int)(a / (double)(NUMSTEPS) * 70.0));
      cout << "\r" << progress 
           << " (" << setprecision(1) << fixed 
           << (double)(a) / NUMSTEPS * 100.0 << "%)" 
           << flush;
      if (a == NUMSTEPS) 
      {
        cout << endl;
      }
    }
      
    if ((a % NSKIP) == 0)
    {
      stringstream ss;
      ss << "out_" << a << ".txt";
      ss >> next_filename;
      printField(etas, next_filename);
      printIndividualFields(etas, 0);
    }
    
    if ((a % 100) == 0 && a >= EQUILIBRIUM)
    {
      // calculateInfo(etas_old, fout_grain_size, fout_GB, fout_multi_junctions, a);
      // printIndividualFields(etas, a); // for debugging and individual field analysis
    }
    
    checkActiveEtas(etas);
    
    // update the old timestep to the current one
    etas_old = etas;
  }
  
  fout_grain_size.close();
  fout_GB.close();
  fout_multi_junctions.close();
  return EXIT_SUCCESS;
}
#ifndef CELL_LINKED_LIST_H
#define CELL_LINKED_LIST_H

#include <vector>
#include "atom.h"

class Cell_ll
{
private:
  int ncellx; // Number of cells in the x direction
  int ncelly; // Number of cells in the y direction
  int ncellz; // Number of cells in the z direction
  double lcellx; // Length of cells in the x direction
  double lcelly; // Length of cells in the y direction
  double lcellz; // Length of cells in the z direction
  double Lx; // Length of x boundary
  double Ly; // Length of y boundary
  double Lz; // Length of z boundary
  double cutoff; // cutoff distance.
  int N; // number of atoms

  int n_atoms_per_cell; // Number of atoms per cell

  std::vector <std::vector <int> > iatom; // The cell-linked list
  std::vector <std::vector <std::vector <int> > > icell; // cell index
  std::vector <std::vector <std::vector <std::vector <int> > > > pcell; // atom index in each cell

  void resizeVector();
  template <class T> int anInt(T x)
  {
    int temp;
    if (x > 0.0) x += 0.5;
    if (x < 0.0) x -= 0.5;
    temp = (int)(x);
    return (double)(temp);
  }

public:
  Cell_ll(); // Default constructor
  Cell_ll(double lx, double ly, double lz, double r_cut, unsigned int size);
  // Getters and setters
  const int getNcellx() {return ncellx;}
  const int getNcelly() {return ncelly;}
  const int getNcellz() {return ncellz;}
  const double getLcellx() {return lcellx;}
  const double getLcelly() {return lcelly;}
  const double getLcellz() {return lcellz;}
  const double getLx() {return Lx;}
  const double getLy() {return Ly;}
  const double getLz() {return Lz;}
  const double getCutoff() {return cutoff;}

  const int getNumNeighbors(int data_id) {return iatom[0][data_id];}
  const int getNeighborIndex(int l, int i) {return iatom[l][i];} // return the lth neighbor of the ith atom

  void setLx(double lx) {Lx = lx;}
  void setLy(double ly) {Ly = ly;}
  void setLz(double lz) {Lz = lz;}
  void setCutoff(double r_cut) {cutoff = r_cut;}

  void calculateNcells();
  void calculateLcells();
  void calculateNumberPerCell();
  void validateBounds();

  // Templated functions are best defined in the header file
  template <class T> void createCellLinkedList(const std::vector <T> &data)
  {
    for (unsigned int i = 0; i < data.size(); ++i)
    {
      int idx = (int)(data[i].getX() / lcellx);
      int idy = (int)(data[i].getY() / lcelly);
      int idz = (int)(data[i].getZ() / lcellz);

      if (idx >= ncellx) idx = ncellx - 1;
      if (idy >= ncelly) idy = ncelly - 1;
      if (idz >= ncellz) idz = ncellz - 1;

      ++icell[idx][idy][idz]; // increase the number in this cell
      // assign the number to this index.
      pcell[idx][idy][idz][icell[idx][idy][idz] - 1] = i;
    }

    for (int i = 0; i < ncellx; ++i) // For each x cell
    {
      for (int j = 0; j < ncelly; ++j) // For each y cell
      {
        for (int k = 0; k < ncellz; ++k) // For each z cell
        {
          for (int l = 0; l < icell[i][j][k]; ++l) // For each atom in this cell
          {
            int id = pcell[i][j][k][l]; // store this atom id
            // Now we check each sub cell around the current one
            for (int ii = -1; ii < 2; ++ii) // allowed values: -1, 0, and 1
            {
              for (int jj = -1; jj < 2; ++jj)
              {
                for (int kk = -1; kk < 2; ++kk)
                {
                  int ia = i + ii; // min value: -1.  Max value: number of cells in dimension
                  int ja = j + jj;
                  int ka = k + kk;
                  // Check to make sure we are still in bounds
                  // C++ indexes from 0, so we accomodate.
                  if (ia >= ncellx) ia = 0;
                  if (ja >= ncelly) ja = 0;
                  if (ka >= ncellz) ka = 0;
                  if (ia < 0) ia = ncellx - 1;
                  if (ja < 0) ja = ncelly - 1;
                  if (ka < 0) ka = ncellz - 1;

                  // Now check each atom in this cell
                  for (int m = 0; m < icell[ia][ja][ka]; ++m)
                  {
                    int jd = pcell[ia][ja][ka][m];
                    // If jd <= id, we've already dealt with this interaction
                    if (jd <= id)
                    {
                      continue;
                    }

                    // Now the actual calculations!
                    double rxij = data[id].getX() - data[jd].getX();
                    double ryij = data[id].getY() - data[jd].getY();
                    double rzij = data[id].getZ() - data[jd].getZ();

                    // Apply PBCs
                    rxij = rxij - anInt(rxij / Lx) * Lx;
                    ryij = ryij - anInt(ryij / Ly) * Ly;
                    rzij = rzij - anInt(rzij / Lz) * Lz;

                    // Now calculate the distance
                    double drij_sq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij);

                    if (drij_sq > cutoff * cutoff)
                    {
                      continue; // move to the next atom if we're too far away
                    }

                    if (drij_sq == 0.0) // This should never be hit, but just in case
                    {
                      continue; // This is the same atom!
                    }

                    // Create the neighbor list
                    iatom[0][id] += 1; //for atom id
                    iatom[(iatom[0][id])][id] = jd;
                    iatom[0][jd] += 1; // for atom jd
                    iatom[(iatom[0][jd])][jd] = id;
                  } // m
                } //kk
              } //jj
            } //ii
          } // l
        } // k
      } // j
    } // i
  }
};
#endif // CELL_LINKED_LIST_H

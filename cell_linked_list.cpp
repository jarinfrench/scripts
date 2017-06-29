#include "cell_linked_list.h"

// Default constructor
Cell_ll::Cell_ll() : Lx(0), Ly(0), Lz(0), cutoff(1.0) {};

// Constructor with important info
Cell_ll::Cell_ll(double lx, double ly, double lz, double r_cut, unsigned int size)
{
  Lx = lx;
  Ly = ly;
  Lz = lz;
  cutoff = r_cut;
  N = size;

  validateBounds();
  calculateNcells();
  calculateLcells();
  calculateNumberPerCell();
  resizeVector();
}

// This calculates the number of cells in each direction.
void Cell_ll::calculateNcells()
{
  ncellx = (int)(Lx / cutoff) + 1;
  ncelly = (int)(Ly / cutoff) + 1;
  ncellz = (int)(Lz / cutoff) + 1;
}

// This function calculates the length of the cells in each direction.
void Cell_ll::calculateLcells()
{
  lcellx = Lx / ncellx;
  lcelly = Ly / ncelly;
  lcellz = Lz / ncellz;
}

// Calculates the number of atoms per cell, with a minimum allowed of 200
void Cell_ll::calculateNumberPerCell()
{
  n_atoms_per_cell = std::max((int)(N / (double)(3 * ncellx * ncelly * ncellz)), 200);
}

void Cell_ll::resizeVector()
{
  icell.resize(ncellx, std::vector <std::vector <int> > // x dimension
              (ncelly, std::vector <int> // y dimension
              (ncellz, 0))); // z dimension

  pcell.resize(ncellx, std::vector <std::vector <std::vector <int> > > //x dimension
              (ncelly, std::vector <std::vector <int> > // y dimension
              (ncellz, std::vector <int> // z dimension
              (n_atoms_per_cell, 0)))); // atom number in cell

  iatom.resize(n_atoms_per_cell, std::vector <int> (N, 0)); // Cell-linked list
}

void Cell_ll::validateBounds()
{
  if (Lx < 0.0) Lx = 0.0;
  if (Ly < 0.0) Ly = 0.0;
  if (Lz < 0.0) Lz = 0.0;
  if (cutoff < 0) cutoff = 1.0;
  if (N < 1) N = 1;
}

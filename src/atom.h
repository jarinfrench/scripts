#ifndef ATOM_H
#define  ATOM_H

/******************************************************************************
* This class contains all of the relevant information pertaining to a LAMMPS
* atom.  id is the ID number of the atom (labeled from 1 to N); type is the type
* of atom (for the specific case of UO2, U is 1, O is 2); charge is the atomic
* charge (note that this in general is a fitted parameter, and not quantized
* amounts); x, y, and z are the Cartesian coordinates of the atom; mark
* specifies a marking number which can be used for various purposes. Only values
* of 0 are "unmarked".
******************************************************************************/

#include <iostream>
//#include "position.h"

class Atom
{
private:
  int id; // id number
  int type; // atomic type, given as a number
  double charge; // atomic charge
  double x; // x position
  double y; // y position
  double z; // z position
  // Position wrapped;
  double xu; // unwrapped x position
  double yu; // unwrapped y position
  double zu; // unwrapped z position
  //Position unwrapped;
  int mark; // the mark on the atom.

public:
  Atom(); // Default constructor
  // Constructor given the atom information
  Atom(int id, int type, double charge, double x, double y, double z);
  // All the getters
  int getId() const {return id;}
  int getType() const {return type;}
  double getCharge() const {return charge;}
  double getX() const {return x;}
  double getY() const {return y;}
  double getZ() const {return z;}
  //Position getWrapped const {return wrapped;}
  double getXu() const {return xu;}
  double getYu() const {return yu;}
  double getZu() const {return zu;}
  //Position getUnwrapped {return unwrapped;}
  int getMark() const {return mark;}

  // All the setters
  void setId(int id) {this->id = id;}
  void setType(int type) {this->type = type;}
  void setCharge(double charge) {this->charge = charge;}
  void setX(double x) {this->x = x;}
  void setY(double y) {this->y = y;}
  void setZ(double z) {this->z = z;}
  //void setWrapped(Position pos) {this->wrapped = pos;}
  void setXu(double x) {this->xu = x;}
  void setYu(double y) {this->yu = y;}
  void setZu(double z) {this->zu = z;}
  //void setUnwrapped(Position pos) {this->unwrapped = pos;}
  void setMark(int mark) {this->mark = mark;}

};
bool operator==(const Atom& lhs, const Atom& rhs);
std::ostream& operator << (std::ostream&, const Atom&);
#endif //  ATOM_H

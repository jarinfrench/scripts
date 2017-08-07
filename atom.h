#ifndef ATOM_H
#define  ATOM_H

/******************************************************************************
* This class contains all of the relevant information pertaining to a LAMMPS
* atom, specifically for U and O.  id is the ID number of the atom (labeled from
* 1 to N); type is the type of atom (for the specific case of UO2, U is 1, O
* is 2); charge is the atomic charge (note that this in general is a fitted
* parameter, and not quantized amounts); x, y, and z are the Cartesian coordinates
* of the atom; mark specifies whether or not the atom is going to be removed.
* Only values of 0 are "unmarked".
******************************************************************************/

class Atom
{
private:
  int id; // id number
  int type; // atomic type, given as a number
  double charge; // atomic charge
  double x; // x position
  double y; // y position
  double z; // z position
  int mark; // is the atom marked for removal?

public:
  Atom(); // Default constructor
  // Constructor given the atom information
  Atom(int id, int type, double charge, double x, double y, double z);
  // All the getters
  int getId() {return id;}
  int getType() {return type;}
  double getCharge() {return charge;}
  double getX() {return x;}
  double getY() {return y;}
  double getZ() {return z;}
  int getMark() {return mark;}

  // All the setters
  void setId(int id) {this->id = id;}
  void setType(int type) {this->type = type;}
  void setCharge(double charge) {this->charge = charge;}
  void setX(double x) {this->x = x;}
  void setY(double y) {this->y = y;}
  void setZ(double z) {this->z = z;}
  void setMark(int mark) {this->mark = mark;}

  //bool operator==(const Atom& rhs) const;
};
#endif //  ATOM_H

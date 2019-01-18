#ifndef POSITION_H
#define POSITION_H

/******************************************************************************
* This class defines a three dimensional position in Cartesian coordinates.
*******************************************************************************/

#include <iostream>

class Position
{
private:
  double x; // x position
  double y; // y position
  double z; // z position

public:
  Position(); // Default Constructor
  Position(double x, double y, double z); // Constructor given the coordinates

  double getX() const {return x;}
  double getY() const {return y;}
  double getZ() const {return z;}

  void setX(double x) {this->x = x;}
  void setY(double y) {this->y = y;}
  void setZ(double z) {this->z = z;}
};

Position operator -(const Position& rhs); // unary minus
Position operator +(const Position& lhs, const Position& rhs);
Position operator -(const Position& lhs, const Position& rhs);
std::ostream& operator << (std::ostream&, const Position&);
std::istream& operator >> (std::istream&, Position&);
bool operator==(const Position& lhs, const Position& rhs);
#endif // POSITION_H

#ifndef POSITION_H
#define POSITION_H

/******************************************************************************
* This class defines a three dimensional position in Cartesian coordinates.
*******************************************************************************/

class Position
{
private:
  double x; // x position
  double y; // y position
  double z; // z position

public:
  Position(); // Default Constructor
  Position(double x, double y, double z); // Constructor given the coordinates

  double getX() {return x;}
  double getY() {return y;}
  double getZ() {return z;}

  void setX(double x) {this->x = x;}
  void setY(double y) {this->y = y;}
  void setZ(double z) {this->z = z;}
};

Position operator +(const Position& lhs, const Position& rhs);
#endif // POSITION_H

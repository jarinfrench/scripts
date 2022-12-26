#ifndef POSITION_H
#define POSITION_H

/******************************************************************************
* This class defines a three dimensional position in Cartesian coordinates.
*******************************************************************************/

#include <iostream>
#include <vector>

class Position {
private:
  double x; // x position
  double y; // y position
  double z; // z position

public:
  Position(); // Default Constructor
  Position(double x, double y);
  Position(double x, double y, double z); // Constructor given the coordinates
  Position(std::vector <double> pos); // given the coordinates as a 1x3 vector

  double getX() const {return x;}
  double getY() const {return y;}
  double getZ() const {return z;}

  void setX(double x) {this->x = x;}
  void setY(double y) {this->y = y;}
  void setZ(double z) {this->z = z;}

  double distance(const Position& p); // distance between this point and another.

  void print2D(std::ostream& stream) {stream << "(" << this->x << "," << this->y << ")";}
  void print2DSpace(std::ostream& stream) {stream << this->x << " " << this->y;}
  void print3DSpace(std::ostream& stream) {stream << this->x << " " << this->y << " " << this->z;}

  // Member operator overloads - defined in the .cpp file
  // NOTE: Because the base type of each Point member is a double, all data types
  // are implicitly converted to doubles, thus the only operator overloads that
  // are required are the ones for doubles.
  Position& operator += (const Position& rhs);
  Position& operator -= (const Position& rhs);
  Position& operator *= (const double& rhs);
  Position& operator /= (const double& rhs);
  double& operator[](int index);
};

// Non member operator overloads
std::ostream& operator << (std::ostream&, const Position&);
std::istream& operator >> (std::istream&, Position&);
bool operator==(const Position& lhs, const Position& rhs);
bool operator!=(const Position& lhs, const Position& rhs);

// inline non member operator overloads - must be defined in the .h file
inline Position operator - (const Position& rhs) {
  return Position(-1.0 * rhs.getX(), -1.0 * rhs.getY(), -1.0 * rhs.getZ());
}

inline Position operator + (const Position& lhs, const Position& rhs) {
  return Position(lhs.getX() + rhs.getX(),
                  lhs.getY() + rhs.getY(),
                  lhs.getZ() + rhs.getZ());
}

inline Position operator - (const Position& lhs, const Position& rhs) {
  return lhs + (-rhs);
}

inline double operator * (const Position& lhs, const Position& rhs) {
  return (lhs.getX() * rhs.getX() + lhs.getY() * rhs.getY() + lhs.getZ() * rhs.getZ());
}

inline Position operator * (const Position& lhs, const double& rhs) {
  return Position(lhs.getX() * rhs, lhs.getY() * rhs, lhs.getZ() * rhs);
}

inline Position operator * (const double& lhs, const Position& rhs) {
  return rhs * lhs;
}

inline Position operator / (const Position& lhs, const double& rhs) {
  return (Position(lhs.getX() / rhs, lhs.getY() / rhs, lhs.getZ() / rhs));
}

#endif // POSITION_H

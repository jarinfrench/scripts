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
  Position(double x, double y);
  Position(double x, double y, double z); // Constructor given the coordinates

  double getX() const {return x;}
  double getY() const {return y;}
  double getZ() const {return z;}

  void setX(double x) {this->x = x;}
  void setY(double y) {this->y = y;}
  void setZ(double z) {this->z = z;}
  
  void print2D(std::ostream& stream) {stream << "(" << this->x << "," << this->y << ")";}
  
  Position& operator+=(const Position& rhs);
  Position& operator-=(const Position& rhs);
  double& operator[](int index);
};

// Non member operator overloads
std::ostream& operator << (std::ostream&, const Position&);
std::istream& operator >> (std::istream&, Position&);
bool operator==(const Position& lhs, const Position& rhs);
bool operator!=(const Position& lhs, const Position& rhs);

// inline non member operator overloads - must be defined in the .h file
inline Position operator - (const Position& rhs)
{
  return Position(-1.0 * rhs.getX(), -1.0 * rhs.getY(), -1.0 * rhs.getZ());
}

inline Position operator + (const Position& lhs, const Position& rhs)
{
  return Position(lhs.getX() + rhs.getX(),
                  lhs.getY() + rhs.getY(),
                  lhs.getZ() + rhs.getZ());
}

inline Position operator - (const Position& lhs, const Position& rhs)
{
  return lhs + (-rhs);
}

inline double operator * (const Position& lhs, const Position& rhs)
{
  return (lhs.getX() * rhs.getX() + lhs.getY() * rhs.getY() + lhs.getZ() * rhs.getZ());
}

template <typename T>
inline Position operator * (const Position& lhs, const T& rhs)
{
  static_assert(std::is_arithmetic<T>::value, "Cannot multiply non-numeric values");
  return Position(lhs.getX() * rhs, lhs.getY() * rhs, lhs.getZ() * rhs);
}

template <typename T>
inline Position operator * (const T& lhs, const Position& rhs)
{
  static_assert(std::is_arithmetic<T>::value, "Cannot multiply non-numeric values");
  return rhs * lhs;
}

inline Position operator * (const Position& lhs, const int& rhs)
{
  return Position(lhs.getX() * rhs, lhs.getY() * rhs, lhs.getZ() * rhs);
}

#endif // POSITION_H

#include "position.h"

// Default constructor sets position to the origin
Position::Position() : x(0), y(0), z(0) {}

// This constructor sets the position to the coordinates given
Position::Position(double x, double y, double z)
{
  this->x = x;
  this->y = y;
  this->z = z;
}

inline Position operator -(const Position& rhs)
{
  return Position(-1.0 * rhs.getX(), -1.0 * rhs.getY(), -1.0 * rhs.getZ());
}

inline Position operator +(const Position& lhs, const Position& rhs)
{
  return Position(lhs.getX() + rhs.getX(),
                  lhs.getY() + rhs.getY(),
                  lhs.getZ() + rhs.getZ());
}

inline Position operator -(const Position& lhs, const Position& rhs)
{
  return lhs + (-rhs);
}

std::ostream& operator << (std::ostream& os, const Position& rhs)
{
  os << "(" << rhs.getX() << ", " << rhs.getY() << ", "
     << rhs.getZ() << ")" << std::endl;
  return os;
}

std::istream& operator >> (std::istream& in, Position& rhs)
{
  float x, y, z;
  in >> x >> y >> z;
  rhs = Position(x,y,z);
  return in;
}

bool operator==(const Position& lhs, const Position& rhs)
{
  return (lhs.getX() == rhs.getX() &&
          lhs.getY() == rhs.getY() &&
          lhs.getZ() == rhs.getZ());
}

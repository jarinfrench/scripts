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

inline Position operator +(const Position& lhs, const Position& rhs)
{
  return Position(lhs.getX() + rhs.getX(),
                  lhs.getY() + rhs.getY(),
                  lhs.getZ() + rhs.getZ())
}

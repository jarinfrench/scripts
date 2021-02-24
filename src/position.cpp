#include "position.h"
#include <cmath>

// Default constructor sets position to the origin
Position::Position() : x(0), y(0), z(0) {}

Position::Position(double x, double y) {
  this->x = x;
  this->y = y;
  this->z = 0.0;
}
// This constructor sets the position to the coordinates given
Position::Position(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

Position::Position(std::vector <double> pos) {
  if (pos.size() != 3) {throw false;}
  this->x = pos[0];
  this->y = pos[1];
  this->z = pos[2];
}

double Position::distance(const Position& p) {
  Position pos;
  pos.setX(this->x - p.getX());
  pos.setY(this->y - p.getY());
  pos.setZ(this->z - p.getZ());
  
  return std::sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);
}

Position& Position::operator+=(const Position& rhs) {
  this->x += rhs.getX();
  this->y += rhs.getY();
  this->z += rhs.getZ();
  return *this;
}

Position& Position::operator-=(const Position& rhs) {
  this->x -= rhs.getX();
  this->y -= rhs.getY();
  this->z -= rhs.getZ();
  return *this;
}

Position& Position::operator*=(const double& rhs) {
  this->x *= rhs;
  this->y *= rhs;
  this->z *= rhs;
  return *this;
}

Position& Position::operator/=(const double& rhs) {
  this->x /= rhs;
  this->y /= rhs;
  this->z /= rhs;
  return *this;
}

double& Position::operator[](int index) {
  if (index < 0 || index >= 3) {throw false;}
  else if (index == 0) {return this->x;}
  else if (index == 1) {return this->y;}
  else if (index == 2) {return this->z;}
  else {throw false;}
}

std::ostream& operator << (std::ostream& os, const Position& rhs) {
  os << "(" << rhs.getX() << ", " << rhs.getY() << ", "
     << rhs.getZ() << ")";
  return os;
}

std::istream& operator >> (std::istream& in, Position& rhs) {
  float x, y, z;
  in >> x >> y >> z;
  rhs = Position(x,y,z);
  return in;
}

bool operator==(const Position& lhs, const Position& rhs) {
  return (lhs.getX() == rhs.getX() &&
          lhs.getY() == rhs.getY() &&
          lhs.getZ() == rhs.getZ());
}

bool operator!=(const Position& lhs, const Position& rhs) {
  return !(lhs == rhs);
}

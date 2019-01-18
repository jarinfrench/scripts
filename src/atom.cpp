#include "atom.h"

// Default constructor sets everything to 0
Atom::Atom() : id(0), type(0), charge(0.0), x(0.0), y(0.0), z(0.0), mark(0) {}

// This constructor sets everything as specified, and sets the mark to "unmarked"
Atom::Atom(int id, int type, double charge, double x, double y, double z)
{
  this->id = id;
  this->type = type;
  this->charge = charge;
  this->x = x;
  this->y = y;
  this->z = z;
  mark = 0;
}

std::ostream& operator << (std::ostream& os, const Atom& atom)
{
  os << "Atom " << atom.getId() << ":\n"
     << "  type: " << atom.getType()
     << "\n  charge: " << atom.getCharge()
     << "\n  x: " << atom.getX() << "(" << atom.getXu() << ")"
     << "\n  y: " << atom.getY() << "(" << atom.getYu() << ")"
     << "\n  z: " << atom.getZ() << "(" << atom.getZu() << ")"
     << "\n  mark: " << atom.getMark() << std::endl;

  return os;
}

bool operator==(const Atom& lhs, const Atom& rhs)
{
  return (lhs.getId() == rhs.getId() &&
          lhs.getType() == rhs.getType() &&
          lhs.getCharge() == rhs.getCharge() &&
          //lhs.getWrapped() == rhs.getWrapped()
          lhs.getX() == rhs.getX() &&
          lhs.getY() == rhs.getY() &&
          lhs.getZ() == rhs.getZ() &&
          lhs.getMark() == rhs.getMark());
}

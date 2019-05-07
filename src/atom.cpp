#include "atom.h"

// Default constructor sets everything to 0
Atom::Atom() : id(0), type(0), charge(0.0), wrapped(0.0, 0.0, 0.0), mark(0) {}

// This constructor sets everything as specified, and sets the mark to "unmarked"
Atom::Atom(int id, int type, double charge, Position p)
{
  this->id = id;
  this->type = type;
  this->charge = charge;
  this->wrapped = p;
  mark = 0;
}

void Atom::setExtraInfo(unsigned int index, double val)
{
  if (index > this->extra_info.size())
  {
    setExtraInfoSize(index);
  }
  this->extra_info[index] = val;
}

std::ostream& operator << (std::ostream& os, const Atom& atom)
{
  os << "Atom " << atom.getId() << ":\n"
     << "  type: " << atom.getType()
     << "\n  charge: " << atom.getCharge()
     << "\n  x: " << atom.getWrapped()[0] << "(" << atom.getUnwrapped()[0] << ")"
     << "\n  y: " << atom.getWrapped()[1] << "(" << atom.getUnwrapped()[1] << ")"
     << "\n  z: " << atom.getWrapped()[2] << "(" << atom.getUnwrapped()[2] << ")"
     << "\n  mark: " << atom.getMark() << std::endl;

  return os;
}

bool operator==(const Atom& lhs, const Atom& rhs)
{
  return (lhs.getId() == rhs.getId() &&
          lhs.getType() == rhs.getType() &&
          lhs.getCharge() == rhs.getCharge() &&
          lhs.getWrapped() == rhs.getWrapped() &&
          // lhs.getX() == rhs.getX() &&
          // lhs.getY() == rhs.getY() &&
          // lhs.getZ() == rhs.getZ() &&
          lhs.getMark() == rhs.getMark());
}

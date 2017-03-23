#include "atom.h"

Atom::Atom() : id(0), type(0), charge(0.0), x(0.0), y(0.0), z(0.0), mark(0) {}

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

int Atom::checkMark(int mark)
{
  if (this->mark >= 1)
    return 1;
  else //(this->mark <= 0)
    return 0;
}

/*bool Atom::operator==(const Atom& rhs) const
{
  return (id == rhs.id && type == rhs.type && charge == rhs.charge &&
          x == rhs.x && y == rhs.y && z == rhs.z && mark == rhs.mark);
}*/

#include "atom.h"

// Default constructor sets everything to 0
Atom::Atom() : id(0), type(0), charge(0.0), wrapped(0.0, 0.0, 0.0), mark(0) {}

// This constructor sets everything as specified, and sets the mark to "unmarked"
Atom::Atom(int id, int type, double charge, Position p) {
  this->id = id;
  this->type = type;
  this->charge = charge;
  this->wrapped = p;
  mark = 0;
}

void Atom::setExtraInfo(unsigned int index, double val) {
  if (index > this->extra_info.size()) {
    setExtraInfoSize(index);
  }
  this->extra_info[index] = val;
}

void Atom::setExtraInfo(std::vector <double> v) {
  this->extra_info = v;
}

void Atom::setExtraInfoSize(unsigned int index) {
  this->extra_info.resize(index, 0.0);
  this->extra_info_names.resize(index, "none");
}

void Atom::setExtraInfoNames(unsigned int index, std::string str) {
  if (index > this->extra_info.size()) {
    setExtraInfoSize(index);
  }
  this->extra_info_names[index] = str;
}

void Atom::setExtraInfoNames(std::vector<std::string> v) {
  this->extra_info_names = v;
}

std::ostream& operator << (std::ostream& os, const Atom& atom) {
  os << "Atom " << atom.getId() << ":\n"
     << "  type: " << atom.getType()
     << "\n  charge: " << atom.getCharge()
     << "\n  x: " << atom.getWrapped()[0] << "(" << atom.getUnwrapped()[0] << ")"
     << "\n  y: " << atom.getWrapped()[1] << "(" << atom.getUnwrapped()[1] << ")"
     << "\n  z: " << atom.getWrapped()[2] << "(" << atom.getUnwrapped()[2] << ")"
     << "\n  mark: " << atom.getMark() << "\n";
  for (size_t i = 0; i < atom.getExtraInfo().size(); ++i) {
    os << "  " << atom.getExtraInfoNames()[i] << ": " << atom.getExtraInfo()[i] << "\n";
  }

  return os;
}

bool operator==(const Atom& lhs, const Atom& rhs) {
  return (lhs.getId() == rhs.getId() &&
          lhs.getType() == rhs.getType() &&
          lhs.getCharge() == rhs.getCharge() &&
          lhs.getWrapped() == rhs.getWrapped() &&
          lhs.getMark() == rhs.getMark());
}

bool operator<(const Atom& lhs, const Atom& rhs) {
  if (lhs.getId() != rhs.getId()) {return lhs.getId() < rhs.getId();}
  else if (lhs.getType() != rhs.getType()) {return lhs.getType() < rhs.getType();}
  else if (lhs.getMark() != rhs.getMark()) {return lhs.getMark() < rhs.getMark();}
  else {return false;}
}

bool operator>(const Atom& lhs, const Atom& rhs) {
  return !(lhs < rhs);
}

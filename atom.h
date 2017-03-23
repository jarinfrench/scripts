#ifndef ATOM_H
#define  ATOM_H

class Atom
{
private:
  int id; // id number
  int type; // atomic type, given as a number
  double charge; // atomic charge
  double x; // x position
  double y; // y position
  double z; // z position
  int mark; // is the atom marked for removal?

public:
  Atom();
  Atom(int id, int type, double charge, double x, double y, double z);
  int getId() {return id;}
  int getType() {return type;}
  double getCharge() {return charge;}
  double getX() {return x;}
  double getY() {return y;}
  double getZ() {return z;}
  int getMark() {return mark;}

  void setId(int id) {this->id = id;}
  void setType(int type) {this->type = type;}
  void setCharge(double charge) {this->charge = charge;}
  void setX(double x) {this->x = x;}
  void setY(double y) {this->y = y;}
  void setZ(double z) {this->z = z;}
  void setMark(int mark) {this->mark = mark;}
  int checkMark(int mark);


  //bool operator==(const Atom& rhs) const;
};
#endif //  ATOM_H

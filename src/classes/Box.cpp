#include "Box.h"
#include "error_code_defines.h"

// default constructor
Box::Box() : xlo(0), xhi(0), ylo(0), yhi(0), zlo(0), zhi(0), xy(0), xz(0), yz(0), has_tilt(false) {}

// 2D, assumed 0 as lower bound
Box::Box(double x, double y) : xlo(0), ylo(0), zlo(0), zhi(0), xy(0), xz(0), yz(0), has_tilt(false) {
  setXHigh(x);
  setYHigh(y);
  validate();
}

// 3D, same lower bound assumption
Box::Box(double x, double y, double z) : xlo(0), ylo(0), zlo(0), xy(0), xz(0), yz(0), has_tilt(false) {
  setXHigh(x);
  setYHigh(y);
  setZHigh(z);
  validate();
}

// 2D, no lower bound assumption
Box::Box(double xlo, double xhi, double ylo, double yhi) : zlo(0), zhi(0), xy(0), xz(0), yz(0), has_tilt(false) {
  setXLow(xlo);
  setXHigh(xhi);
  setYLow(ylo);
  setYHigh(yhi);
  validate();
}

// 3D, no lower bound assumption
Box::Box(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi) : xy(0), xz(0), yz(0), has_tilt(false) {
  setXLow(xlo);
  setXHigh(xhi);
  setYLow(ylo);
  setYHigh(yhi);
  setZLow(zlo);
  setZHigh(zhi);
  validate();
}

bool Box::validate() const {
  if (xlo > xhi || ylo > yhi) {
    return false;
  }
  if (zhi != 0 && zlo > zhi) {
    return false;
  }
  return true;
}

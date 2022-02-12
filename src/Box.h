#ifndef BOX_H
#define BOX_H

/******************************************************************************
* This class holds all the information for the box bounds in a LAMMPS simulation
*******************************************************************************/

class Box {
private:
  double xlo, xhi, ylo, yhi, zlo, zhi; // lower and upper bounds of the box
  double xy, xz, yz; // tilt factors
  bool has_tilt;

public:
  Box();
  Box(double x, double y); // 2D, assume lower bound is 0
  Box(double x, double y, double z); // 3D, assume lower bound is 0
  Box(double xlo, double xhi, double ylo, double yhi); // 2D, lower bound not assumed
  Box(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi); // 3D, lower bound not assumed
  double getXBound() const {return this->xhi - this->xlo;}
  double getYBound() const {return this->yhi - this->ylo;}
  double getZBound() const {return this->zhi - this->zlo;}
  double getXLow() const {return this->xlo;}
  double getYLow() const {return this->ylo;}
  double getZLow() const {return this->zlo;}
  double getXHigh() const {return this->xhi;}
  double getYHigh() const {return this->yhi;}
  double getZHigh() const {return this->zhi;}
  double getXYTilt() const {return this->xy;}
  double getXZTilt() const {return this->xz;}
  double getYZTilt() const {return this->yz;}

  void setXLow(double xlo) {this->xlo = xlo;}
  void setYLow(double ylo) {this->ylo = ylo;}
  void setZLow(double zlo) {this->zlo = zlo;}
  void setXHigh(double xhi) {this->xhi = xhi;}
  void setYHigh(double yhi) {this->yhi = yhi;}
  void setZHigh(double zhi) {this->zhi = zhi;}
  void setXYTilt(double xy) {this->xy = xy; has_tilt = true;}
  void setXZTilt(double xz) {this->xz = xz; has_tilt = true;}
  void setYZTilt(double yz) {this->yz = yz; has_tilt = true;}

  bool validate() const;

  template <typename T>
  void writeBounds(T& stream) const {
    stream << xlo << " " << xhi << " xlo xhi\n"
           << ylo << " " << yhi << " ylo yhi\n"
           << zlo << " " << zhi << " zlo zhi";

    if (has_tilt) {
      stream << "\n" << xy << " " << xz << " " << yz << " xy xz yz";
    }
  }
};

#endif // BOX_H

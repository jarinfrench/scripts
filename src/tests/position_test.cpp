#include <iostream>
#include "../position.h"

using namespace std;

int main(int argc, char** argv)
{
  Position p1(1,1,1);
  Position p2(2,-1,-5);
  double d1 = 2.5;
  int i1 = 3;
  float f1 = 10.0;
  long double ld1 = 0.1;
  long int li1 = -10;

  cout << "Unary minus example: p1 = " << p1 << ", -p1 = " << -p1 << endl
       << "Binary plus example: p1 = " << p1 << ", p2 = " << p2 << ", p1 + p2 = " << p1 + p2 << endl
       << "Binary minus example: p1 - p2 = " << p1 - p2 << endl
       << "Dot product example: p1 * p2 = " << p1 * p2 << endl
       << "Scaling examples:\n\tDouble: d1 = " << d1 << ", d1 * p1 = " << d1 * p1 << ", p1 * d1 = " << p1 * d1 << endl
       << "\tInt: i1 = " << i1 << ", i1 * p1 = " << i1 * p1 << ", p1 * i1 = " << p1 * i1 << endl
       << "\tFloat: f1 = " << f1 << ", f1 * p1 = " << f1 * p1 << ", p1 * f1 = " << p1 * f1 << endl
       << "\tLong Double: ld1 = " << ld1 << ", ld1 * p1 = " << ld1 * p1 << ", p1 * ld1 = " << p1 * ld1 << endl
       << "\tLong int: li1 = " << li1 << ", li1 * p1 = " << li1 * p1 << ", p1 * li1 = " << p1 * li1 << endl;
  p1 += p2;
  cout << "+= operator: p1 += p2 --> p1 = " << p1 << endl;
  p1 -= p2;
  cout << "-= operator: p1 -= p2 --> p1 = " << p1 << endl;
  p1 *= d1;
  cout << "*= operator: p1 *= d1 --> p1 = " << p1 << endl;
  p1 /= d1;
  cout << "/= operator: p1 /= d1 --> p1 = " << p1 << endl;
  cout << "2D printing: p1 = ";
  p1.print2D(cout);
  cout << ", p2 = ";
  p2.print2D(cout);
  cout << "\n[] operator example: p1[0] = " << p1[0] << ", p1[1] = " << p1[1] << ", p1[2] = " << p1[2] << endl
       << "\n[] operator example: p2[0] = " << p2[0] << ", p2[1] = " << p2[1] << ", p2[2] = " << p2[2] << endl;

  try {
    cout << p1[-1];
  }
  catch (bool val) {
    cout << "Access to p1[-1]: " << val << endl;
  }

  try {
    cout << p1[3];
  }
  catch (bool val)   {
    cout << "Access to p1[3]: " << val << endl;
  }

}

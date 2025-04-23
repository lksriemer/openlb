/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Shota Ito, Timm Krüger
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

#ifndef PVECTEN_H
#define PVECTEN_H

#include <cmath>
#include <iostream>

namespace olb {

namespace membrane {

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////
////////////// PHYSICAL VECTOR (p_vec) /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

template <class T = double> struct p_vec {
  // Constructor.
  p_vec() {
    x = 0; y = 0; z = 0;
  }
  p_vec(T a, T b, T c) {
    x = a;  y = b;  z = c;
  }

  // Destructor.
  ~p_vec() {}

  // Components.
  T x, y, z;

  ////////////////////////////////////////////////////////////////////////////////////////////

  double length() const { /// Computes the length of the vector.
    return sqrt(x * x + y * y + z * z);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  void invert() { /// Inverts the direction of the vector.
    x = - x;
    y = - y;
    z = - z;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  void mult_by(double a) { /// Multiplies the vector by an arbitrary number.
    x *= a;
    y *= a;
    z *= a;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  void norm_to(double a) { /// Normalizes the vector to a given length.
    double l = length();
    if(l != 0) {
      x *= a / l;
      y *= a / l;
      z *= a / l;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  void set_to(T a, T b, T c) { /// Sets the values of the vector explicitly.
    x = a;
    y = b;
    z = c;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  void reset() { /// Resets the vector to zero.
    x = 0;
    y = 0;
    z = 0;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec<int> integer() { /// Casts vector to integer.
    p_vec temp(*this);
    p_vec<int> output;
    output.x = (int) temp.x;
    output.y = (int) temp.y;
    output.z = (int) temp.z;
    return output;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec operator-() const { /// Defines the operator - (unary negation).
    p_vec temp(*this);
    temp.x *= -1;
    temp.y *= -1;
    temp.z *= -1;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec<double> operator+(p_vec<int> op2) const { /// Defines the operator + (sum of two vectors).
    p_vec<double> temp(*this);
    temp.x += op2.x;
    temp.y += op2.y;
    temp.z += op2.z;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec<double> operator+(p_vec<double> op2) const { /// Defines the operator + (sum of two vectors).
    p_vec<double> temp(*this);
    temp.x += op2.x;
    temp.y += op2.y;
    temp.z += op2.z;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec& operator+=(const p_vec op2) { /// Defines the operator += (sum of two vectors).
    x += op2.x;
    y += op2.y;
    z += op2.z;
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec<double> operator-(p_vec<int> op2) const { /// Defines the operator - (difference of two vectors).
    p_vec<double> temp(*this);
    temp.x -= op2.x;
    temp.y -= op2.y;
    temp.z -= op2.z;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec<double> operator-(p_vec<double> op2) const { /// Defines the operator - (difference of two vectors).
    p_vec<double> temp(*this);
    temp.x -= op2.x;
    temp.y -= op2.y;
    temp.z -= op2.z;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec& operator-=(const p_vec op2) { /// Defines the operator -= (difference of two vectors).
    x -= op2.x;
    y -= op2.y;
    z -= op2.z;
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  template <class U> double operator*(p_vec<U> op2) const { /// Defines the operator * (scalar product of two vectors).
    p_vec temp(*this);
    return (temp.x * op2.x + temp.y * op2.y + temp.z * op2.z);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec operator%(p_vec op2) const { /// Defines the operatort % (cross product of two vectors).
    p_vec temp(*this);
    p_vec temp2(*this);
    temp.x = temp2.y * op2.z - temp2.z * op2.y;
    temp.y = temp2.z * op2.x - temp2.x * op2.z;
    temp.z = temp2.x * op2.y - temp2.y * op2.x;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec operator/(T op) const { /// Defines the operator / (divide vector by type).
    p_vec temp(*this);
    temp.x = temp.x / op;
    temp.y = temp.y / op;
    temp.z = temp.z / op;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec& operator/=(T op) { /// Defines the operator / (divide vector by type).
    x /= op;
    y /= op;
    z /= op;
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec<double> operator*(double op) const { /// Defines the operator * (multiplies vector by double).
    p_vec<double> temp(this->x, this->y, this->z);
    temp.x = temp.x * op;
    temp.y = temp.y * op;
    temp.z = temp.z * op;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec& operator*=(T op) { /// Defines the operator * (multiplies vector by type); must be p_vec * double and NOT double * p_vec.
    x *= op;
    y *= op;
    z *= op;
    return *this;
  }
};

template <class T> ostream& operator<<(ostream& os, const p_vec<T>& vec) { /// << operator
  os << vec.x << " " << vec.y << " " << vec.z;
  return os;
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////// PHYSICAL TENSOR RANK 2 (p_ten2) /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

struct p_ten2 {
  // Constructor.
  p_ten2() {
    xx = 0.;  xy = 0.;  xz = 0.;
    yx = 0.;  yy = 0.;  yz = 0.;
    zx = 0.;  zy = 0.;  zz = 0.;
  }
  p_ten2(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    xx = a; xy = b; xz = c;
    yx = d; yy = e; yz = f;
    zx = g; zy = h; zz = i;
  }

  // Destructor.
  ~p_ten2() {}

  // Components.
  double xx, xy, xz, yx, yy, yz, zx, zy, zz;

  ////////////////////////////////////////////////////////////////////////////////////////////

  double tr() { /// Computes the trace of the tensor.
    return xx + yy + zz;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  double det() { /// Computes the Determinant of the tensor.
    return xx * (yy * zz - yz * zy) + xy * (yz * zx - yx * zz) + xz * (yx * zy - yy * zx);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  double mag() { /// Computes the magnitude of the tensor.
    return sqrt(xx * xx + yy * yy + zz * zz + 2. * xy * xy + 2. * xz * xz + 2. * yz * yz);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2 invert() { /// Computes the inverse of the tensor.
    p_ten2 temp(*this);
    p_ten2 inverse(
      temp.yy * temp.zz - temp.yz * temp.zy,
      temp.xz * temp.zy - temp.xy * temp.zz,
      temp.xy * temp.yz - temp.xz * temp.yy,
      temp.yz * temp.zx - temp.yx * temp.zz,
      temp.xx * temp.zz - temp.xz * temp.zx,
      temp.xz * temp.yx - temp.xx * temp.yz,
      temp.yx * temp.zy - temp.yy * temp.zx,
      temp.xy * temp.zx - temp.xx * temp.zy,
      temp.xx * temp.yy - temp.xy * temp.yx
    );
    return inverse / temp.det();
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  void mult_by(double a) { /// Multiplies the tensor by an arbitrary number.
    xx *= a;  xy *= a;  xz *= a;
    yx *= a;  yy *= a;  yz *= a;
    yx *= a;  yy *= a;  yz *= a;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  void set_to(double a, double b, double c, double d, double e, double f, double g, double h, double i) { /// Sets the values of the tensor.
    xx = a; xy = b; xz = c;
    yx = d; yy = e; yz = f;
    zx = g; zy = h; zz = i;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  void reset() { /// Resets the values of the tensor.
    xx = 0.; xy = 0.; xz = 0.;
    yx = 0.; yy = 0.; yz = 0.;
    zx = 0.; zy = 0.; zz = 0.;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  bool is_symm() { /// Checks whether the tensor is symmetric.
    if(xy == yx && xz == zx && yz == zy) {
      return true;
    }
    else {
      return false;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  bool is_asymm() { /// Checks whether the tensor is asymmetric.
    if(xy == -yx && xz == -zx && yz == -zy && xx == 0 && yy == 0 && zz == 0) {
      return true;
    }
    else {
      return false;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2 operator-() const { /// Defines the operator - (unary negation).
    p_ten2 temp(*this);
    temp.xx *= -1;  temp.xy *= -1;  temp.xz *= -1;
    temp.yx *= -1;  temp.yy *= -1;  temp.yz *= -1;
    temp.zx *= -1;  temp.zy *= -1;  temp.zz *= -1;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2 operator+(p_ten2 op2) const { /// Defines the operator + (sum of two tensors).
    p_ten2 temp(*this);
    temp.xx += op2.xx;  temp.xy += op2.xy;  temp.xz += op2.xz;
    temp.yx += op2.yx;  temp.yy += op2.yy;  temp.yz += op2.yz;
    temp.zx += op2.zx;  temp.zy += op2.zy;  temp.zz += op2.zz;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2& operator+=(p_ten2 op2) { /// Defines the operator += (sum of two tensors).
    xx += op2.xx;  xy += op2.xy;  xz += op2.xz;
    yx += op2.yx;  yy += op2.yy;  yz += op2.yz;
    zx += op2.zx;  zy += op2.zy;  zz += op2.zz;
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2 operator-(p_ten2 op2) const { /// Defines the operator - (difference of two tensors).
    p_ten2 temp(*this);
    temp.xx -= op2.xx;  temp.xy -= op2.xy;  temp.xz -= op2.xz;
    temp.yx -= op2.yx;  temp.yy -= op2.yy;  temp.yz -= op2.yz;
    temp.zx -= op2.zx;  temp.zy -= op2.zy;  temp.zz -= op2.zz;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2& operator-=(p_ten2 op2) { // Defines the operator -= (difference of two tensors).
    xx -= op2.xx;  xy -= op2.xy;  xz -= op2.xz;
    yx -= op2.yx;  yy -= op2.yy;  yz -= op2.yz;
    zx -= op2.zx;  zy -= op2.zy;  zz -= op2.zz;
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2 operator*(p_ten2 op2) const { /// Defines the operator * (tensor product: tensor * tensor => tensor).
    p_ten2 temp(*this);
    p_ten2 temp2(*this);
    temp.xx = temp2.xx * op2.xx + temp2.xy * op2.yx + temp2.xz * op2.zx;
    temp.xy = temp2.xx * op2.xy + temp2.xy * op2.yy + temp2.xz * op2.zy;
    temp.xz = temp2.xx * op2.xz + temp2.xy * op2.yz + temp2.xz * op2.zz;
    temp.yx = temp2.yx * op2.xx + temp2.yy * op2.yx + temp2.yz * op2.zx;
    temp.yy = temp2.yx * op2.xy + temp2.yy * op2.yy + temp2.yz * op2.zy;
    temp.yz = temp2.yx * op2.xz + temp2.yy * op2.yz + temp2.yz * op2.zz;
    temp.zx = temp2.zx * op2.xx + temp2.zy * op2.yx + temp2.zz * op2.zx;
    temp.zy = temp2.zx * op2.xy + temp2.zy * op2.yy + temp2.zz * op2.zy;
    temp.zz = temp2.zx * op2.xz + temp2.zy * op2.yz + temp2.zz * op2.zz;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_vec<> operator*(p_vec<> op2) const { /// Defines the operator * (tensor product: tensor * vector => vector).
    p_ten2 temp2(*this);
    p_vec<> temp;
    temp.x = temp2.xx * op2.x + temp2.xy * op2.y + temp2.xz * op2.z;
    temp.y = temp2.yx * op2.x + temp2.yy * op2.y + temp2.yz * op2.z;
    temp.z = temp2.zx * op2.x + temp2.zy * op2.y + temp2.zz * op2.z;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  double operator%(p_ten2 op2) const { /// Defines the operator % (contraction of two tensors).
    p_ten2 temp(*this);
    return temp.xx * op2.xx + temp.xy * op2.xy + temp.xz * op2.xz +
           temp.yx * op2.yx + temp.yy * op2.yy + temp.yz * op2.yz +
           temp.zx * op2.zx + temp.zy * op2.zy + temp.zz * op2.zz;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2 operator*(double a) const { /// Defines the operator * (product of tensor with constant).
    p_ten2 temp(*this);
    temp.xx *= a;  temp.xy *= a;  temp.xz *= a;
    temp.yx *= a;  temp.yy *= a;  temp.yz *= a;
    temp.zx *= a;  temp.zy *= a;  temp.zz *= a;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2& operator*=(double a) { /// Defines the operator *= (product of tensor with constant).
    xx *= a;  xy *= a;  xz *= a;
    yx *= a;  yy *= a;  yz *= a;
    zx *= a;  zy *= a;  zz *= a;
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2 operator/(double a) const { /// Defines the operator / (quotient of tensor and constant).
    p_ten2 temp(*this);
    temp.xx /= a;  temp.xy /= a;  temp.xz /= a;
    temp.yx /= a;  temp.yy /= a;  temp.yz /= a;
    temp.zx /= a;  temp.zy /= a;  temp.zz /= a;
    return temp;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////

  p_ten2& operator/=(double a) { /// Defines the operator /= (quotient of tensor and constant).
    xx /= a;  xy /= a;  xz /= a;
    yx /= a;  yy /= a;  yz /= a;
    zx /= a;  zy /= a;  zz /= a;
    return *this;
  }
};

ostream &operator<< (ostream &ostr, const p_ten2 &tensor) { /// cout operator.
    ostr << tensor.xx << " " << tensor.xy << " " << tensor.xz << " "
         << tensor.yx << " " << tensor.yy << " " << tensor.yz << " "
         << tensor.zx << " " << tensor.zy << " " << tensor.zz;
    return ostr;
}

}

}

#endif

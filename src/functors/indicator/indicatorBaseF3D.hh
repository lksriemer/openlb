/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Albert Mink, Mathias J. Krause, Benjamin Förster
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

#ifndef INDICATOR_BASE_F_3D_HH
#define INDICATOR_BASE_F_3D_HH


#include<cmath>
#include "indicatorBaseF3D.h"
#include "math.h"

namespace olb {

template <typename S>
IndicatorF3D<S>::IndicatorF3D() : GenericF<bool,S>(1, 3)
{}

template <typename S>
Vector<S,3>& IndicatorF3D<S>::getMin()
{
  return _myMin;
}

template <typename S>
Vector<S,3>& IndicatorF3D<S>::getMax()
{
  return _myMax;
}

template <typename S>
bool IndicatorF3D<S>::distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC)
{
  bool originValue;
  bool currentValue;
  S precision = .0001;
  S pitch = 0.5;

  // start at origin and move into given direction
  Vector<S,3> currentPoint(origin);

  (*this)(&originValue, origin.data);
  (*this)(&currentValue, currentPoint.data);

  while (currentValue == originValue && isInsideBox(currentPoint)) {
    currentPoint += direction;
    // update currentValue until the first point on the other side (inside/outside) is found
    (*this)(&currentValue, currentPoint.data);
  }

  // return false if no point was found in given direction
  if (!isInsideBox(currentPoint) && !originValue) {
    return false;
  }


  while (pitch >= precision) {
    if (!isInsideBox(currentPoint) && originValue) {
      currentPoint -= pitch * direction;
      pitch /= 2.;
    } else {
      (*this)(&currentValue, currentPoint.data);
      if (currentValue == originValue) {
        currentPoint += pitch * direction;
        pitch /= 2.;
      } else {
        currentPoint -= pitch * direction;
        pitch /= 2.;
      }
    }
  }

  distance = (currentPoint - origin).norm();
  return true;
}

template <typename S>
bool IndicatorF3D<S>::isInsideBox(Vector<S,3> point)
{
  return point >= _myMin && point <= _myMax;
}



template <typename S>
SuperIndicatorF3D<S>::SuperIndicatorF3D (SuperF3D<S>& rhs) : _superF(rhs)
{ };



DiscreteIndicatorF3D::DiscreteIndicatorF3D() : GenericF<bool,int>(1, 3)
{ }


DiscreteIndicatorFalse3D::DiscreteIndicatorFalse3D()
{}

bool DiscreteIndicatorFalse3D::operator() (bool output[], const int input[])
{
  // always return false
  output[0] = false;
  return true;
}

DiscreteIndicatorTrue3D::DiscreteIndicatorTrue3D()
{}

bool DiscreteIndicatorTrue3D::operator() (bool output[], const int input[])
{
  // always return true
  output[0] = true;
  return true;
}


template <typename S>
DiscreteIndicatorMaterial3D<S>::DiscreteIndicatorMaterial3D (SuperGeometry3D<S>& rhs, std::vector<int> materialNumbers)
  : _superGeometry(rhs), _materialNumbers(materialNumbers)
{ }

template <typename S>
bool DiscreteIndicatorMaterial3D<S>::operator() (bool output[], const int input[])
{
  // iterate over material numbers and check if given point has that material number
  bool pointHasMaterial = false;
  for (int& m : _materialNumbers) {
    pointHasMaterial |= ( _superGeometry.get(input[0], input[1], input[2], input[3]) == m );
  }
  output[0] = pointHasMaterial;
  return true;
}


// \todo This method is correctly implemented but probably useless(?) Leaving it for the moment.
//template <typename T>
//void DiscreteIndicatorMaterial3D<T>::calculateMinMax() {
//  // Initialize `_myMin` with largest possible value
//  Vector<int,3> motherExtend = _superGeometry.getCuboidGeometry().getMotherCuboid().getExtend();
//  motherExtend -= 1;
//  _superGeometry.getCuboidGeometry().getMotherCuboid().getPhysR(this->_myMin.data, motherExtend.data);
//
//  // Iterate through all points and calculate min and max
//  int nCglob = _superGeometry.getCuboidGeometry().getNc();
//  for(int iCglob = 0; iCglob < nCglob; iCglob++) {
//    Vector<int,3> extend = _superGeometry.getCuboidGeometry().get(iCglob).getExtend();
//    for(int iX = 0; iX < extend[0]; iX++) {
//      for(int iY = 0; iY < extend[1]; iY++) {
//        for (int iZ = 0; iZ < extend[2]; iZ++) {
//
//          for(int& m: _materialNumbers) {
//            if(_superGeometry.get(iCglob, iX, iY, iZ) == m) {
//              Vector<T,3> physR;
//              _superGeometry.getCuboidGeometry().get(iCglob).getPhysR(physR.data, extend.data);
//              for(int iDim=0; iDim<3; iDim++) {
//                this->_myMin[iDim] = std::min(this->_myMin[iDim], physR[iDim]);
//                this->_myMax[iDim] = std::max(this->_myMax[iDim], physR[iDim]);
//              }
//            }
//          }
//
//        }
//      }
//    }
//  }
//}




// identity to "store results"
template <typename S>
IndicatorIdentity3D<S>::IndicatorIdentity3D(IndicatorF3D<S>& f) : _f(f)
{
  this->_myMin = _f.getMin();
  this->_myMax = _f.getMax();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename S>
bool IndicatorIdentity3D<S>::operator() (bool output[], const S input[])
{
  return _f(output, input);
}


template <typename T, typename S>
SmoothIndicatorF3D<T,S>::SmoothIndicatorF3D() : AnalyticalF3D<T,S>(1)
{}

template <typename T, typename S>
Vector<S,3> & SmoothIndicatorF3D<T,S>::getMin()
{
  return _myMin;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getMax()
{
  return _myMax;
}


// identity to "store results"
template <typename T, typename S>
SmoothIndicatorIdentity3D<T,S>::SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T,S>& f)
  : _f(f)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = _f.getMin()[i];
    this->_myMax[i] = _f.getMax()[i];
  }
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool SmoothIndicatorIdentity3D<T,S>::operator() (T output[], const S input[])
{
  _f(output, input);
  return true;
}

} // namespace olb

#endif

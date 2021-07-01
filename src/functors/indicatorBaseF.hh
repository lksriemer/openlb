/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Cyril Masquelier, Mathias J. Krause
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

#ifndef INDICATOR_BASE_F_HH
#define INDICATOR_BASE_F_HH

#include<vector>
#include<cmath>

#include "functors/indicatorBaseF.h"
#include "math.h"

namespace olb {

template <typename T, typename S>
bool IndicatorF2D<T,S>::distance(S& distance, const std::vector<S>& origin, const std::vector<S>& direction, int iC)
{
  T originValue = (*this)(origin)[0];
  std::vector<S> currentPoint(origin);

  S precision = .0001;

  while ((*this)(currentPoint)[0] == originValue && isInsideBox(currentPoint)) {
	for (int i = 0; i < 2; i++) {
      currentPoint[i] += direction[i];
    }
  }

  if (!isInsideBox(currentPoint) && !originValue) {
	  return false;
  }
  S pitch= 0.5;

  while (pitch >= precision) {
	  if (!isInsideBox(currentPoint) && originValue) {
	    for (int i = 0; i < 2; i++) {
        currentPoint[i] -= pitch*direction[i];
      }
      pitch = pitch/2.;
    }
	  else {
      if ((*this)(currentPoint)[0] == originValue) {
        for (int i = 0; i < 2; i++) {
          currentPoint[i] += pitch*direction[i];
        }
        pitch = pitch/2.;
      }
      else {
  	    for (int i = 0; i < 2; i++) {
            currentPoint[i] -= pitch*direction[i];
        }
        pitch = pitch/2.;
      }
	  }
  }

  S distance2 = 0.;
  for (int i = 0; i < 2; i++) {
	  distance2 += (currentPoint[i] - origin[i]) * (currentPoint[i] - origin[i]);
  }

  distance = sqrt(distance2); 
  return true;
}


template <typename T, typename S>
bool IndicatorF2D<T,S>::isInsideBox(std::vector<S> Point) {
  if (Point[0] >= _myMin[0] && Point[0] <= _myMax[0] &&
      Point[1] >= _myMin[1] && Point[1] <= _myMax[1] ) {
	  return true;
  }
  else {
	  return false;
  }
}



template <typename T, typename S>
bool IndicatorF3D<T,S>::distance(S& distance,const std::vector<S>& origin,
  const std::vector<S>& direction, int iC)
{
  T originValue = (*this)(origin)[0];
  std::vector<S> currentPoint(origin);

  S precision = .0001;

  while ((*this)(currentPoint)[0] == originValue && isInsideBox(currentPoint)) {
	  for (int i = 0; i < 3; i++) {
        currentPoint[i] += direction[i];
      }
  }

  if (!isInsideBox(currentPoint) && !originValue) {
	  return false;
  }
  S pitch = 0.5;

  while (pitch >= precision) {
	  if (!isInsideBox(currentPoint) && originValue) {
	    for (int i = 0; i < 3; i++) {
        currentPoint[i] -= pitch*direction[i];
      }
      pitch = pitch/2.;
    }
	  else {
      if ((*this)(currentPoint)[0] == originValue) {
        for (int i = 0; i < 3; i++) {
          currentPoint[i] += pitch*direction[i];
        }
        pitch = pitch/2.;
      }
      else {
  	    for (int i=0; i<3; i++) {
            currentPoint[i] -= pitch*direction[i];
        }
        pitch = pitch/2.;
      }
	  }
  }

  S distance2 = 0.;
  for (int i = 0; i < 3; i++) {
	  distance2 += (currentPoint[i] - origin[i]) * (currentPoint[i] - origin[i]);
  }

  distance = sqrt(distance2); 
  return true;
}

template <typename T, typename S>
bool IndicatorF3D<T,S>::isInsideBox(std::vector<S> Point) {
  if (Point[0]>= _myMin[0] && Point[0]<= _myMax[0] && Point[1]>= _myMin[1] &&
      Point[1]<= _myMax[1] && Point[2]>= _myMin[2] && Point[2]<= _myMax[2] ) {
  	return true;
  }
  else {
  	return false;
  }
}


// identity to "store results"
template <typename T, typename S>
IndicatorIdentity2D<T,S>::IndicatorIdentity2D(IndicatorF2D<T,S>& f)
  : IndicatorF2D<T,S>(f.getTargetDim()), _f(f) {
  for (int iDim = 0; iDim<2; iDim++) {
    this->_myMin.push_back(_f.getMin()[iDim]);
    this->_myMax.push_back(_f.getMax()[iDim]);
  }
  // add 'this' to father's child list to prevent father from being deleted
  _f.addChild(this);
}

template <typename T, typename S>
IndicatorIdentity2D<T,S>::~IndicatorIdentity2D() {
  // remove 'this' from father's child list
  _f.removeChild(this);
  // delete father from grandfather's child list
  _f.myErase(NULL);
}

template <typename T, typename S>
std::vector<T> IndicatorIdentity2D<T,S>::operator()(std::vector<S> input) 
{ return _f(input); }


// identity to "store results"
template <typename T, typename S>
IndicatorIdentity3D<T,S>::IndicatorIdentity3D(IndicatorF3D<T,S>& f)
  : IndicatorF3D<T,S>(f.getTargetDim()), _f(f) {
  for (int iDim = 0; iDim<3; iDim++) {
    this->_myMin.push_back(_f.getMin()[iDim]);
    this->_myMax.push_back(_f.getMax()[iDim]);
  }
  // add 'this' to father's child list to prevent father from being deleted
  _f.addChild(this);
}

template <typename T, typename S>
IndicatorIdentity3D<T,S>::~IndicatorIdentity3D() {
  // remove 'this' from father's child list
  _f.removeChild(this);
  // delete father from grandfather's child list
  _f.myErase(NULL);
}

//template <typename T, typename S>
//void IndicatorIdentity3D<T,S>::myErase(GenericF<T,S>* ptr) {
//  this->removeChild(ptr);
//  delete ptr;
//  if( this->_pointerVec.size() == 0 )  _f.myErase(this);
//};

template <typename T, typename S>
std::vector<T> IndicatorIdentity3D<T,S>::operator()(std::vector<S> input) 
{ return _f(input); }



// identity to "store results"
template <typename T, typename S>
SmoothIndicatorIdentity3D<T,S>::SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T,S>& f)
  : SmoothIndicatorF3D<T,S>(f.getTargetDim()), _f(f) {
  for (int iDim = 0; iDim<3; iDim++) {
    this->_myMin.push_back(_f.getMin()[iDim]);
    this->_myMax.push_back(_f.getMax()[iDim]);
  }
  // add 'this' to father's child list to prevent father from being deleted
  _f.addChild(this);
}

template <typename T, typename S>
SmoothIndicatorIdentity3D<T,S>::~SmoothIndicatorIdentity3D() {
  // remove 'this' from father's child list
  _f.removeChild(this);
  // delete father from grandfather's child list
  _f.myErase(NULL);
}

template <typename T, typename S>
std::vector<T> SmoothIndicatorIdentity3D<T,S>::operator()(std::vector<S> input)
{ return _f(input); }


} // namespace olb

#endif

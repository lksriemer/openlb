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

#ifndef INDICATOR_F_HH
#define INDICATOR_F_HH

#include<vector>
#include<cmath>

#include "functors/indicatorF.h"
#include "io/stlReader.h"
#include "io/xmlReader.h"
#include "utilities/vectorHelpers.h"
#include "math.h"

namespace olb {


// Warning : the cuboid is only defined parallel to the plans x=0 and y=0 !!!
// For other cuboids, please use the parallelepiped version
template <typename T, typename S>
IndicatorCuboid2D<T,S>::IndicatorCuboid2D(std::vector<S> extend, std::vector<S> origin)
  : IndicatorF2D<T,S>(1)
{
  for(int iDim = 0; iDim < 2; iDim++) {
    _center.push_back(origin[iDim] + .5*extend[iDim]);
    this->_myMin.push_back(origin[iDim]);
    this->_myMax.push_back(origin[iDim] + extend[iDim]);
  }
  _xLength = extend[0];
  _yLength = extend[1];
}

template <typename T, typename S>
IndicatorCuboid2D<T,S>::IndicatorCuboid2D(S xLength, S yLength, std::vector<S> center )
  : IndicatorF2D<T,S>(1), _center(center), _xLength(xLength), _yLength(yLength)
{
  S temp0[] = {_center[0] - _xLength/2.,
               _center[1] - _yLength/2.};
  this->_myMin.insert( this->_myMin.begin(), temp0, temp0 +2 );
  S temp1[] = {_center[0] + _xLength/2.,
               _center[1] + _yLength/2.};
  this->_myMax.insert( this->_myMax.begin(), temp1, temp1 +2 );
}


// returns true if x is inside the cuboid
template <typename T, typename S>
std::vector<T> IndicatorCuboid2D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  T result = T();
  if(    fabs(_center[0]-x[0]) <= _xLength/2.
      && fabs(_center[1]-x[1]) <= _yLength/2.) {
    result = T(1);
  }
  y.push_back(result);
  return y;
}


template <typename T, typename S>
IndicatorCircle2D<T,S>::IndicatorCircle2D(std::vector<S> center, S radius)
  :  IndicatorF2D<T,S>(1), _center(center), _radius2(radius*radius)
{
	S r = sqrt(_radius2);
	S temp0[] = {_center[0]-r, _center[1]-r};
	this->_myMin.insert(this->_myMin.begin(),temp0, temp0 +2);
	S temp1[] = {_center[0]+r, _center[1]+r};
	this->_myMax.insert(this->_myMax.begin(),temp1, temp1 +2);
}

// returns true if x is inside the circle
template <typename T, typename S>
std::vector<T> IndicatorCircle2D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  T result = T();
  if(  (_center[0]-x[0]) * (_center[0]-x[0])
      +(_center[1]-x[1]) * (_center[1]-x[1]) <= _radius2 ){
    result = T(1);
  }
  y.push_back(result);
  return y;
}


template <typename T, typename S>
bool IndicatorCircle2D<T,S>::distance(S& distance, const std::vector<S>& origin,
  const std::vector<S>& direction, int iC)
{
  S a = direction[0]*direction[0] + direction[1]*direction[1];   

  // returns 0 if point is at the boundary of the sphere 
  if ( a == _radius2 ) {
    distance = S();
    return true;
  }
  // norm of direction
  a = sqrt(a);

  S b = 2.*((origin[0] - _center[0])*direction[0] + 
            (origin[1] - _center[1])*direction[1])/a;
  S c = -_radius2 + (origin[0] - _center[0])*(origin[0] - _center[0])
                  + (origin[1] - _center[1])*(origin[1] - _center[1]);

  // discriminant
  S d = b*b - 4.*c;  
  if (d < 0) return false;

  S x1 = (- b + sqrt(d))/2.;
  S x2 = (- b - sqrt(d))/2.;

  // case if origin is inside the sphere
  if ((x1<0.) || (x2<0.)) {
    if (x1>0.) { distance = x1; return true; }
    if (x2>0.) { distance = x2; return true; }
  }
  // case if origin is ouside the sphere
  else {
    distance = std::min(x1,x2);
    return true;
  }

  return false;
}


template <typename T, typename S>
IndicatorCircle3D<T,S>::IndicatorCircle3D(std::vector<S> center,
  std::vector<S> normal, S radius)
  :  IndicatorF3D<T,S>(1), _center(center), _normal(util::normalize(normal) ),
    _radius2(radius*radius)
{
	S r = sqrt(_radius2);
	S temp0[] = {_center[0]-r, _center[1]-r, _center[2]-r};
	this->_myMin.insert(this->_myMin.begin(), temp0, temp0 +3);
	S temp1[] = {_center[0]+r, _center[1]+r, _center[2]+r};
	this->_myMax.insert(this->_myMax.begin(), temp1, temp1 +3);
}

template <typename T, typename S>
IndicatorCircle3D<T,S>::IndicatorCircle3D(S center0, S center1, S center2,
  S normal0, S normal1, S normal2, S radius)
  :  IndicatorF3D<T,S>(1), _radius2(radius*radius)
{
  _center.push_back(center0); _center.push_back(center1); _center.push_back(center2);
  _normal.push_back(normal0); _normal.push_back(normal1); _normal.push_back(normal2);
  _normal = util::normalize(_normal);
	S r = sqrt(_radius2);
	S temp0[] = {_center[0]-r, _center[1]-r, _center[2]-r};
	this->_myMin.insert(this->_myMin.begin(),temp0, temp0 +3);
	S temp1[] = {_center[0]+r, _center[1]+r, _center[2]+r};
	this->_myMax.insert(this->_myMax.begin(),temp1, temp1 +3);
}

// returns true if x is inside the circle
template <typename T, typename S>
std::vector<T> IndicatorCircle3D<T,S>::operator()(std::vector<S> x)
{
  S eps = std::numeric_limits<S>::epsilon();
  IndicatorCylinder3D<T,S> cylinder(_center, _normal, getRadius(), eps);
  return cylinder(x);
}


template <typename T, typename S>
IndicatorStl3D<T,S>::IndicatorStl3D(STLreader<S>& stlReader, S eps)
  : IndicatorF3D<T,S>(1), _stlReader(stlReader), _eps(eps) { }

template <typename T, typename S>
std::vector<T> IndicatorStl3D<T,S>::operator()(std::vector<S> x) {
  std::vector<T> result(1,false);
  std::vector<S> direction(3,S());
  std::vector<S> origin(3,S());
  S distance = 2.*_eps;
  
  direction[0] = 2.*_eps; direction[1] = 0; direction[2] = 0;
  origin[0] = x[0]-_eps; origin[1] = x[1]; origin[2] = x[2];
  _stlReader.distance( distance, origin, direction );
  if (distance<2.*_eps) { result[0] = true; return result;}
  direction[0] = 0; direction[1] = 2.*_eps; direction[2] = 0;
  origin[0] = x[0]; origin[1] = x[1]-_eps; origin[2] = x[2];
  _stlReader.distance( distance, origin, direction );
  if (distance<2.*_eps) { result[0] = true; return result;}
  direction[0] = 0; direction[1] = 0; direction[2] = 2.*_eps;
  origin[0] = x[0]; origin[1] = x[1]; origin[2] = x[2]-_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<2.*_eps) { result[0] = true; return result;}

  distance = 1.42*2.*_eps;

  direction[0] = 1.42*2.*_eps; direction[1] = 1.42*2.*_eps; direction[2] = 0.;
  origin[0] = x[0]-1.42*_eps; origin[1] = x[1]-1.42*_eps; origin[2] = x[2];
  _stlReader.distance( distance, origin, direction );
  if (distance<1.42*2.*_eps) { result[0] = true; return result;}
  direction[0] = 1.42*2.*_eps; direction[1] = -1.42*2.*_eps; direction[2] = 0.;
  origin[0] = x[0]-1.42*_eps; origin[1] = x[1]+1.42*_eps; origin[2] = x[2];
  _stlReader.distance( distance, origin, direction );
  if (distance<1.42*2.*_eps) { result[0] = true; return result;}
  direction[0] = 1.42*2.*_eps; direction[1] = 0.; direction[2] = 1.42*2.*_eps;
  origin[0] = x[0]-1.42*_eps; origin[1] = x[1]; origin[2] = x[2]-1.42*_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<1.42*2.*_eps) { result[0] = true; return result;}
  direction[0] = 1.42*2.*_eps; direction[1] = 0.; direction[2] = -1.42*2.*_eps;
  origin[0] = x[0]-1.42*_eps; origin[1] = x[1]; origin[2] = x[2]+1.42*_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<1.42*2.*_eps) { result[0] = true; return result;}
  direction[0] = 0.; direction[1] = 1.42*2.*_eps; direction[2] = 1.42*2.*_eps;
  origin[0] = x[0]; origin[1] = x[1]-1.42*_eps; origin[2] = x[2]-1.42*_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<1.42*2.*_eps) { result[0] = true; return result;}
  direction[0] = 0.; direction[1] = 1.42*2.*_eps; direction[2] = -1.42*2.*_eps;
  origin[0] = x[0]; origin[1] = x[1]-1.42*_eps; origin[2] = x[2]+1.42*_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<1.42*2.*_eps) { result[0] = true; return result;}

  distance = 1.7321*2.*_eps;

  direction[0] = 1.7321*2.*_eps; direction[1] = 1.7321*2.*_eps; direction[2] = 1.7321*2.*_eps;
  origin[0] = x[0]-1.7321*_eps; origin[1] = x[0]-1.7321*_eps; origin[2] = x[0]-1.7321*_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<1.7321*2.*_eps) { result[0] = true; return result;}
  direction[0] = 1.7321*2.*_eps; direction[1] = 1.7321*2.*_eps; direction[2] = -1.7321*2.*_eps;
  origin[0] = x[0]-1.7321*_eps; origin[1] = x[0]-1.7321*_eps; origin[2] = x[0]+1.7321*_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<1.7321*2.*_eps) { result[0] = true; return result;}
  direction[0] = 1.7321*2.*_eps; direction[1] = -1.7321*2.*_eps; direction[2] = 1.7321*2.*_eps;
  origin[0] = x[0]-1.7321*_eps; origin[1] = x[0]+1.7321*_eps; origin[2] = x[0]-1.7321*_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<1.7321*2.*_eps) { result[0] = true; return result;}
  direction[0] = 1.7321*2.*_eps; direction[1] = -1.7321*2.*_eps; direction[2] = -1.7321*2.*_eps;
  origin[0] = x[0]-1.7321*_eps; origin[1] = x[0]+1.7321*_eps; origin[2] = x[0]+1.7321*_eps;
  _stlReader.distance( distance, origin, direction );
  if (distance<1.7321*2.*_eps) { result[0] = true; return result;}

  return result;
}


template <typename T, typename S>
IndicatorSphere3D<T,S>::IndicatorSphere3D(std::vector<S> center, S radius)
  :  IndicatorF3D<T,S>(1), _center(center), _radius2(radius*radius)
{
	double r = sqrt(_radius2);
	S temp0[] = {_center[0]-r, _center[1]-r, _center[2]-r};
	this->_myMin.insert(this->_myMin.begin(),temp0, temp0 +3);
	S temp1[] = {_center[0]+r, _center[1]+r, _center[2]+r};
	this->_myMax.insert(this->_myMax.begin(),temp1, temp1 +3);
}

// returns true if x is inside the sphere
template <typename T, typename S>
std::vector<T> IndicatorSphere3D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  T result = T();
  if(  (_center[0]-x[0]) * (_center[0]-x[0])
      +(_center[1]-x[1]) * (_center[1]-x[1])
      +(_center[2]-x[2]) * (_center[2]-x[2]) <= _radius2 ){
    result = T(1);
  }
  y.push_back(result);
  return y;
}

template <typename T, typename S>
bool IndicatorSphere3D<T,S>::distance(S& distance, const std::vector<S>& origin,
  const std::vector<S>& direction, int iC)
{
  // computes pos. distance by solving quadratic equation by a-b-c-formula
  S a = direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2];   

  // returns 0 if point is at the boundary of the sphere 
  if (a ==_radius2) {
    distance = S();
    return true;
  }
  // norm of direction
  a = sqrt(a);

  S b = 2.*((origin[0] - _center[0])*direction[0] + 
            (origin[1] - _center[1])*direction[1] +
            (origin[2] - _center[2])*direction[2])/a;
  S c = -_radius2 + (origin[0] - _center[0])*(origin[0] - _center[0])
                  + (origin[1] - _center[1])*(origin[1] - _center[1])
                  + (origin[2] - _center[2])*(origin[2] - _center[2]);

  // discriminant
  S d = b*b - 4.*c;  
  if (d < 0) return false;

  S x1 = (- b + sqrt(d))/2.;
  S x2 = (- b - sqrt(d))/2.;

  // case if origin is inside the sphere
  if ((x1<0.) || (x2<0.)) {
    if (x1>0.) { distance = x1; return true; }
    if (x2>0.) { distance = x2; return true; }
  }
  // case if origin is ouside the sphere
  else {
    distance = std::min(x1,x2);
    return true;
  }

  return false;
}


template <typename T, typename S>
IndicatorLayer3D<T,S>::IndicatorLayer3D(IndicatorF3D<T,S>& indicatorF, S layerSize)
  :  IndicatorF3D<T,S>(1), _indicatorF(indicatorF), _layerSize(layerSize) {
  S temp0[] = {indicatorF.getMin()[0]-layerSize,
               indicatorF.getMin()[1]-layerSize,
               indicatorF.getMin()[2]-layerSize };
  this->_myMin.insert( this->_myMin.begin(),temp0, temp0 +3 );
  S temp1[] = {indicatorF.getMax()[0]+layerSize,
               indicatorF.getMax()[1]+layerSize,
               indicatorF.getMax()[2]+layerSize };
  this->_myMax.insert( this->_myMax.begin(),temp1, temp1 +3 );
}

// returns true if x is inside the layer
template <typename T, typename S>
std::vector<T> IndicatorLayer3D<T,S>::operator()(std::vector<S> x) {
  std::vector<T> y(1,false);
  for (int iX=-1;iX<2;iX++) {
    for (int iY=-1;iY<2;iY++) {
      for (int iZ=-1;iZ<2;iZ++) {
        std::vector<S> r(x);
        r[0]+=iX*_layerSize; r[1]+=iY*_layerSize; r[2]+=iZ*_layerSize;
        if (_indicatorF(r)[0]) { 
          y[0]=true;
          return y;
        }
      }
    }
  }
  return y;
}


/// Indicator function for a cylinder
template <typename T, typename S>
IndicatorCylinder3D<T,S>::IndicatorCylinder3D(std::vector<S> center1,
  std::vector<S> center2, S radius)
  :  IndicatorF3D<T,S>(1), _center1(center1), _center2(center2), _radius2(radius*radius)
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
	init();
}

/// Indicator function for a cylinder
template <typename T, typename S>
IndicatorCylinder3D<T,S>::IndicatorCylinder3D(std::vector<S> center1,
  std::vector<S> normal, S radius, S eps)
  :  IndicatorF3D<T,S>(1), _center1(center1), _radius2(radius*radius)
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
	_center1[0] -= .5*eps*normal[0];
	_center1[1] -= .5*eps*normal[1];
	_center1[2] -= .5*eps*normal[2];
	_center2.push_back(_center1[0] + eps*normal[0]);
	_center2.push_back(_center1[1] + eps*normal[1]);
	_center2.push_back(_center1[2] + eps*normal[2]);

	init();
}

/// Indicator function for a cylinder
template <typename T, typename S>
IndicatorCylinder3D<T,S>::IndicatorCylinder3D(IndicatorCircle3D<T,S> circleF, S eps)
  :  IndicatorF3D<T,S>(1), _center1(circleF.getCenter()),
    _radius2(circleF.getRadius()*circleF.getRadius())
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
	_center1[0] -= .5*eps*circleF.getNormal()[0];
	_center1[1] -= .5*eps*circleF.getNormal()[1];
	_center1[2] -= .5*eps*circleF.getNormal()[2];
	_center2.push_back(_center1[0] + eps*circleF.getNormal()[0]);
	_center2.push_back(_center1[1] + eps*circleF.getNormal()[1]);
	_center2.push_back(_center1[2] + eps*circleF.getNormal()[2]);

	init();
}

// returns true if x is inside the cylinder
template <typename T, typename S>
std::vector<T> IndicatorCylinder3D<T,S>::operator()(std::vector<S> x) {
  std::vector<T> y;
  T result = T();

  S X = _I[0]*(x[0]-_center1[0]) + _I[1]*(x[1]-_center1[1]) + _I[2]*(x[2]-_center1[2]);
  S Y = _J[0]*(x[0]-_center1[0]) + _J[1]*(x[1]-_center1[1]) + _J[2]*(x[2]-_center1[2]);
  S Z = _K[0]*(x[0]-_center1[0]) + _K[1]*(x[1]-_center1[1]) + _K[2]*(x[2]-_center1[2]);

  // X^2 + Y^2 <= _radius2
  if ( Z <= _length && Z >= 0 && X*X + Y*Y <= _radius2 ) {
    result = T(1);
  }
  y.push_back(result);
  return y;
}

template <typename T, typename S>
void IndicatorCylinder3D<T,S>::init() {
	  _length = sqrt( (_center2[0]-_center1[0]) * (_center2[0]-_center1[0])
	                 +(_center2[1]-_center1[1]) * (_center2[1]-_center1[1])
	                 +(_center2[2]-_center1[2]) * (_center2[2]-_center1[2]) );

	  // _K = centre2 - centre1 (normalized)
	  S tempk[] = {(_center2[0]-_center1[0]) / _length,
	               (_center2[1]-_center1[1]) / _length,
	               (_center2[2]-_center1[2]) / _length};
	  _K.insert(_K.begin(),tempk, tempk +3);

	  // _I and _J form an orthonormal base with _K
	  if ((_center2[1]-_center1[1]) == 0 && (_center2[0]-_center1[0]) == 0) {
	    if ((_center2[2]-_center1[2]) == 0) {
	      std::cout << "Warning: in the cylinder, the two centers have the same coordinates";
	    }
		  S tempi[] = {1,0,0};
		  S tempj[] = {0,1,0};
		  _I.insert(_I.begin(),tempi, tempi +3);
		  _J.insert(_J.begin(),tempj, tempj +3);
	  }
	  else
		{
		  S normi = sqrt (_K[1]*_K[1]+_K[0]*_K[0]);
		  S tempi[] = {-_K[1]/normi, _K[0]/normi,0};
		  _I.insert(_I.begin(),tempi, tempi +3);
		  S tempj[] = {_K[1]*_I[2] - _K[2]*_I[1],
	                 _K[2]*_I[0] - _K[0]*_I[2],
	                 _K[0]*_I[1] - _K[1]*_I[0] };
		  _J.insert(_J.begin(),tempj, tempj +3);
		}

	  double r = sqrt(_radius2);
	  S maxx, maxy, maxz;
	  S minx, miny, minz;

	  maxx= _center1[0] + sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + std::max(_K[0]*_length, 0.);
	  minx= _center1[0] - sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + std::min(_K[0]*_length, 0.);

	  maxy= _center1[1] + sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + std::max(_K[1]*_length, 0.);
	  miny= _center1[1] - sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + std::min(_K[1]*_length, 0.);

	  maxz= _center1[2] + sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + std::max(_K[2]*_length, 0.);
	  minz= _center1[2] - sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + std::min(_K[2]*_length, 0.);

	  S temp0[] = {minx, miny, minz};
	  this->_myMin.insert(this->_myMin.begin(),temp0, temp0 +3);
	  S temp1[] = {maxx, maxy, maxz};
	  this->_myMax.insert(this->_myMax.begin(),temp1, temp1 +3);
}


// cone defined by the centers of the two extremities and the radiuses of the two extremities
// the 2nd radius is optional: if it is not defined, the 2nd center is the vertex of the cone
template <typename T, typename S>
IndicatorCone3D<T,S>::IndicatorCone3D(std::vector<S> center1,
  std::vector<S> center2, S radius1, S radius2)
  :  IndicatorF3D<T,S>(1), _center1(center1), _center2(center2),
    _radius1(radius1), _radius2(radius2)
{
  // _I,_J,_K is the new base where _K is the axe of the cone

  // _K = centre2 - centre1 (normalized)
  _length = sqrt( (_center2[0]-_center1[0]) * (_center2[0]-_center1[0])
                 +(_center2[1]-_center1[1]) * (_center2[1]-_center1[1])
                 +(_center2[2]-_center1[2]) * (_center2[2]-_center1[2]) );
  S tempk[] = {(_center2[0]-_center1[0]) / _length,
               (_center2[1]-_center1[1]) / _length,
               (_center2[2]-_center1[2]) / _length};
  _K.insert(_K.begin(),tempk, tempk +3);  // _K = centre2 - centre1 (normalized)

  // _I and _J form an orthonormal base with _K
  if ( _center2[1] - _center1[1] == 0 && _center2[0] - _center1[0] == 0) {
    if ( _center2[2] - _center1[2] == 0) {
      std::cout << "Warning: in the cone, the two center have the same coordinates";
    }
	  S tempi[] = {1,0,0};
	  S tempj[] = {0,1,0};
	  _I.insert(_I.begin(),tempi, tempi +3);
	  _J.insert(_J.begin(),tempj, tempj +3);
  }
  else {
	  S normi = sqrt(_K[1]*_K[1] + _K[0]*_K[0]);
	  S tempi[] = {-_K[1]/normi, _K[0]/normi,0};
	  _I.insert(_I.begin(),tempi, tempi +3);
	  S tempj[] = {_K[1]*_I[2] - _K[2]*_I[1],
                 _K[2]*_I[0] - _K[0]*_I[2],
                 _K[0]*_I[1] - _K[1]*_I[0] };
	  _J.insert(_J.begin(),tempj, tempj +3);
  }

  S maxx, maxy, maxz;
  S minx, miny, minz;

  maxx= _center1[0] + std::max( sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                                sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);
  minx= _center1[0] + std::min(-sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                               -sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);

  maxy= _center1[1] + std::max( sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                                sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);
  miny= _center1[1] + std::min(-sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                               -sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);

  maxz= _center1[2] + std::max( sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                                sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);
  minz= _center1[2] + std::min(-sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                               -sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);

  S temp0[] = {minx, miny, minz};
  this->_myMin.insert(this->_myMin.begin(),temp0, temp0 +3);
  S temp1[] = {maxx, maxy, maxz};
  this->_myMax.insert(this->_myMax.begin(),temp1, temp1 +3);
}

// returns true if x is inside the cone(std::vector<S> center1, std::vector<S> center2, S radius1
template <typename T, typename S>
std::vector<T> IndicatorCone3D<T,S>::operator()(std::vector<S> x) {
  std::vector<T> y;
  T result = T();
  // radius: the radius of the cone at the point x

  S X = _I[0]*(x[0]-_center1[0]) + _I[1]*(x[1]-_center1[1]) + _I[2]*(x[2]-_center1[2]);
  S Y = _J[0]*(x[0]-_center1[0]) + _J[1]*(x[1]-_center1[1]) + _J[2]*(x[2]-_center1[2]);
  S Z = _K[0]*(x[0]-_center1[0]) + _K[1]*(x[1]-_center1[1]) + _K[2]*(x[2]-_center1[2]);
  S radius = _radius1 + (_radius2 - _radius1)*Z / _length;

  if ( Z <= _length && Z >= 0 && X*X + Y*Y <= radius*radius ) {
    result = T(1);
  }
  y.push_back(result);
  return y;
}


// Warning : the cuboid is only defined parallel to the plans x=0, y=0 and z=0 !!!
// For other cuboids, please use the parallelepiped version
template <typename T, typename S>
IndicatorCuboid3D<T,S>::IndicatorCuboid3D(std::vector<S> extend,
  std::vector<S> origin) : IndicatorF3D<T,S>(1)
{
  for(int iDim=0;iDim<3;iDim++) {
    _center.push_back(origin[iDim] + .5*extend[iDim]);
    this->_myMin.push_back(origin[iDim]);
    this->_myMax.push_back(origin[iDim] + extend[iDim]);
  }
  _xLength = extend[0];
  _yLength = extend[1];
  _zLength = extend[2];
}

template <typename T, typename S>
IndicatorCuboid3D<T,S>::IndicatorCuboid3D(S xLength, S yLength, S zLength, 
  std::vector<S> center)
  : IndicatorF3D<T,S>(1), _center(center), _xLength(xLength), _yLength(yLength),
    _zLength(zLength)
{
  S temp0[] = {_center[0] - _xLength/2.,
               _center[1] - _yLength/2.,
               _center[2] - _zLength/2.};
  this->_myMin.insert( this->_myMin.begin(),temp0, temp0 +3 );
  S temp1[] = {_center[0] + _xLength/2.,
               _center[1] + _yLength/2.,
               _center[2] + _zLength/2.};
  this->_myMax.insert( this->_myMax.begin(),temp1, temp1 +3 );
}

// returns true if x is inside the cuboid
template <typename T, typename S>
std::vector<T> IndicatorCuboid3D<T,S>::operator()(std::vector<S> x) {
  std::vector<T> y;
  T result = T();
  if(    fabs(_center[0]-x[0]) <= _xLength/2.
      && fabs(_center[1]-x[1]) <= _yLength/2.
      && fabs(_center[2]-x[2]) <= _zLength/2. ) {
    result = T(1);
  }
  y.push_back(result);
  return y;
}

template <typename T, typename S>
IndicatorCuboid3D<T,S>* createIndicatorCuboid3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCuboid3D");

  std::vector<S> origin(3,S());
  std::vector<S> extend(3,S(1));

  params.setWarningsOn(false);
  // params[parameter].read(value) sets the value or returns false if the parameter can not be found

  if (!params["originX"].read(origin[0], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: originX. Set default: originX=0" << std::endl;
  }
  if (!params["originY"].read(origin[1], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: originY. Set default: originY=0" << std::endl;
  }
  if (!params["originZ"].read(origin[2], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: originZ. Set default: originZ=0" << std::endl;
  }


  if (!params["extendX"].read(extend[0], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: extendX. Set default: extendX=1" << std::endl;
  }
  if (!params["extendY"].read(extend[1], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: extendY. Set default: extendY=1" << std::endl;
  }
  if (!params["extendZ"].read(extend[2], verbose)) {
    clout << "Warning: Cannot read parameter from Xml-file: extendZ. Set default: extendZ=1" << std::endl;
  }

  params.setWarningsOn(true);

  return new IndicatorCuboid3D<T,S>(extend, origin);
}


template <typename T, typename S>
IndicatorParallelepiped3D<T,S>::IndicatorParallelepiped3D(std::vector<S> origin,
  std::vector<S> corner1, std::vector<S> corner2, std::vector<S> corner3)
  :  IndicatorF3D<T,S>(1), _origin(origin), _corner1(corner1), _corner2(corner2),
    _corner3(corner3)
{ // _I,_J,_K is the new base , each vector is one of the side of the parallelepiped
  double tempI[] = {_corner1[0] - _origin[0],
                    _corner1[1] - _origin[1],
                    _corner1[2] - _origin[2]};
  _I.insert(_I.begin(),tempI, tempI +3);  // _I = corner1 - origin
  double tempJ[] = {_corner2[0] - _origin[0],
                    _corner2[1] - _origin[1],
                    _corner2[2] - _origin[2]};
  _J.insert(_J.begin(),tempJ, tempJ +3);  // _J = corner2 - origin
  double tempk[] = {_corner3[0] - _origin[0],
                    _corner3[1] - _origin[1],
                    _corner3[2] - _origin[2]};
  _K.insert(_K.begin(),tempk, tempk +3);  // _K = corner3 - origin

  _normI = sqrt( _I[0]*_I[0] + _I[1]*_I[1] + _I[2]*_I[2] );
  _normJ = sqrt( _J[0]*_J[0] + _J[1]*_J[1] + _J[2]*_J[2] );
  _normK = sqrt( _K[0]*_K[0] + _K[1]*_K[1] + _K[2]*_K[2] );

  for (int j = 0; j < 3; j++) {
    this->_myMax.push_back( _origin[j] +std::max(_corner1[j]-_origin[j],S(0))
                                       +std::max(_corner2[j]-_origin[j],S(0))
                                       +std::max(_corner3[j]-_origin[j],S(0)) );
    this->_myMin.push_back (_origin[j] +std::min(_corner1[j]-_origin[j],S(0))
                                       +std::min(_corner2[j]-_origin[j],S(0))
                                       +std::min(_corner3[j]-_origin[j],S(0)) );
  }
}

// returns true if x is inside the parallelepiped
template <typename T, typename S>
std::vector<T> IndicatorParallelepiped3D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  T result = T(0);
  S X = _I[0]*(x[0]-_origin[0]) + _I[1]*(x[1]-_origin[1]) + _I[2]*(x[2]-_origin[2]);
  S Y = _J[0]*(x[0]-_origin[0]) + _J[1]*(x[1]-_origin[1]) + _J[2]*(x[2]-_origin[2]);
  S Z = _K[0]*(x[0]-_origin[0]) + _K[1]*(x[1]-_origin[1]) + _K[2]*(x[2]-_origin[2]);

  if( X >= 0 && X <= _normI*_normI &&
	    Y >= 0 && Y <= _normJ*_normJ &&
	    Z >= 0 && Z <= _normK*_normK )
  {
    result = T(1);
  }
  y.push_back(result);
  return y;
}


template <typename T, typename S>
SmoothIndicatorCircle2D<T,S>::SmoothIndicatorCircle2D(std::vector<S> center,
  S radius, S epsilon)
  :  SmoothIndicatorF2D<T,S>(1), _center(center), _radius2(radius*radius), _epsilon(epsilon)
{
	double r = sqrt(_radius2);
	S temp0[] = {_center[0] - r,
               _center[1] - r};
	this->_myMin.insert(this->_myMin.begin(),temp0, temp0 +2);
	S temp1[] = {_center[0] + r,
               _center[1] + r};
	this->_myMax.insert(this->_myMax.begin(),temp1, temp1 +2);
}

// returns true if x is inside the sphere
template <typename T, typename S>
std::vector<T> SmoothIndicatorCircle2D<T,S>::operator()(std::vector<S> x) {
  std::vector<T> y;
  T result = T();
  double d;   // distance to the figure
  double distToCenter = sqrt( (_center[0]-x[0]) * (_center[0]-x[0])+
		                          (_center[1]-x[1]) * (_center[1]-x[1]) );
  double radiusOut = sqrt(_radius2) + _epsilon;

  if ( distToCenter <= sqrt(_radius2) ) {
	  y.push_back(T(1));
	  return y;
  } else if ( distToCenter >= radiusOut) {
	  y.push_back(T(0));
	  return y;
  } else if ( distToCenter >= sqrt(_radius2) && distToCenter <= radiusOut ) {
	  d = distToCenter - sqrt(_radius2);
    result = T( cos(M_PI*d/(2*_epsilon)) *cos(M_PI*d/(2*_epsilon)));
    y.push_back(result);
    return y;
  }
  return std::vector<T>(1,T());
}


template <typename T, typename S>
SmoothIndicatorSphere3D<T,S>::SmoothIndicatorSphere3D(std::vector<S> center,
  S radius, S epsilon)
  :  SmoothIndicatorF3D<T,S>(1), _center(center), _radius2(radius*radius), _epsilon(epsilon)
{
	double r = sqrt(_radius2);
	S temp0[] = {_center[0] - r,
               _center[1] - r,
               _center[2] - r};
	this->_myMin.insert(this->_myMin.begin(),temp0, temp0 +3);
	S temp1[] = {_center[0] + r,
               _center[1] + r,
               _center[2] + r};
	this->_myMax.insert(this->_myMax.begin(),temp1, temp1 +3);
}

// returns true if x is inside the sphere
template <typename T, typename S>
std::vector<T> SmoothIndicatorSphere3D<T,S>::operator()(std::vector<S> x) {
  std::vector<T> y;
  T result = T();
  double d;   // distance to the figure
  double distToCenter = sqrt( (_center[0]-x[0]) * (_center[0]-x[0])+
		                          (_center[1]-x[1]) * (_center[1]-x[1])+
		                          (_center[2]-x[2]) * (_center[2]-x[2]) );
  double radiusOut = sqrt(_radius2) + _epsilon;

  if ( distToCenter <=  sqrt(_radius2) ) {
	  y.push_back(T(1));
	  return y;
  } else if ( distToCenter >= radiusOut) {
	  y.push_back(T(0));
	  return y;
  } else if ( distToCenter >= sqrt(_radius2) && distToCenter <= radiusOut ) {
	  d = distToCenter - sqrt(_radius2);
    result = T( cos(M_PI*d/(2*_epsilon)) *cos(M_PI*d/(2*_epsilon)));
    y.push_back(result);
    return y;
  }
  return std::vector<T>(1,T());
}


template <typename T, typename S>
SmoothIndicatorCylinder3D<T,S>::SmoothIndicatorCylinder3D(std::vector<S> center1,
  std::vector<S> center2, S radius, S epsilon)
  : SmoothIndicatorF3D<T,S>(1), _center1(center1), _center2(center2),
    _radius2(radius*radius) , _epsilon(epsilon)
{ // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0

  _length = sqrt( (_center2[0]-_center1[0])*(_center2[0]-_center1[0])
                 +(_center2[1]-_center1[1])*(_center2[1]-_center1[1])
                 +(_center2[2]-_center1[2])*(_center2[2]-_center1[2]) );

  // _K = centre2 - centre1 (normalized)
  S tempk[] = {(_center2[0] - _center1[0])/_length,
               (_center2[1] - _center1[1])/_length,
               (_center2[2] - _center1[2])/_length};
  _K.insert(_K.begin(),tempk, tempk +3);

  // _I and _J form an orthonormal base with _K
  if( (_center2[1]-_center1[1]) ==0 && (_center2[0]-_center1[0]) == 0 ) {
    if ((_center2[2]-_center1[2])==0) {
      std::cout << "Warning: in the cylinder, the two centers have the same coordinates";
    }
	  S tempi[] = {1,0,0};
	  S tempj[] = {0,1,0};
	  _I.insert(_I.begin(),tempi, tempi +3);
	  _J.insert(_J.begin(),tempj, tempj +3);
  }
  else
	{
	  S normi = sqrt (_K[1]*_K[1] + _K[0]*_K[0]);
	  S tempi[] = {-_K[1]/normi, _K[0]/normi,0};
	  _I.insert(_I.begin(),tempi, tempi +3);
	  S tempj[] = {_K[1]*_I[2] - _K[2]*_I[1],
                 _K[2]*_I[0] - _K[0]*_I[2],
                 _K[0]*_I[1] - _K[1]*_I[0] };
	  _J.insert(_J.begin(),tempj, tempj +3);
	}

  double r = sqrt(_radius2);
  S maxx, maxy, maxz;
  S minx, miny, minz;

  maxx= _center1[0] + sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + std::max(_K[0]*_length, 0.);
  minx= _center1[0] - sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + std::min(_K[0]*_length, 0.);

  maxy= _center1[1] + sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + std::max(_K[1]*_length, 0.);
  miny= _center1[1] - sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + std::min(_K[1]*_length, 0.);

  maxz= _center1[2] + sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + std::max(_K[2]*_length, 0.);
  minz= _center1[2] - sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + std::min(_K[2]*_length, 0.);


  S temp0[] = {minx, miny, minz};
  this->_myMin.insert( this->_myMin.begin(),temp0, temp0 +3 );
  S temp1[] = {maxx, maxy, maxz};
  this->_myMax.insert( this->_myMax.begin(),temp1, temp1 +3 );
}

// returns true if x is inside the cylinder
template <typename T, typename S>
std::vector<T> SmoothIndicatorCylinder3D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  T result = T();

  double X = _I[0]*(x[0]-_center1[0]) + _I[1]*(x[1]-_center1[1]) + _I[2]*(x[2]-_center1[2]);
  double Y = _J[0]*(x[0]-_center1[0]) + _J[1]*(x[1]-_center1[1]) + _J[2]*(x[2]-_center1[2]);
  double Z = _K[0]*(x[0]-_center1[0]) + _K[1]*(x[1]-_center1[1]) + _K[2]*(x[2]-_center1[2]);
  double d = -1;  // distance from the point to the cylinder

  if (X*X + Y*Y  <= _radius2) {
    if (Z <= _length && Z >= 0) {
	    d = 0;
	  }
	  if (Z >= _length && Z <= _length + _epsilon){
	    d = Z-_length;
	  }
	  if (Z >= - _epsilon && Z <= 0){
	    d = - Z;
	  }
  }
  else if (X*X + Y*Y  <= std::pow(sqrt(_radius2) + _epsilon,2)) {
    if (Z <= _length && Z >= 0 ) {
      d = sqrt( X*X + Y*Y ) - sqrt(_radius2);
    }
    if (Z >= _length && Z <= _length + _epsilon){
      d = sqrt( std::pow(Z-_length,2) 
               +std::pow( (sqrt(X*X + Y*Y) - sqrt(_radius2) ) , 2) );
    }
    if (Z >= - _epsilon && Z <= 0){
      d = sqrt( Z*Z + std::pow((sqrt(X*X + Y*Y) - sqrt(_radius2)) , 2) );
    }
  }

  if (d>=0 && d<=_epsilon) {
    result = T( cos(M_PI*d/(2*_epsilon)) * cos(M_PI*d/(2*_epsilon)) );
  }
  y.push_back(result);

  return y;
}


// cone defined by the centers of the two extremities and the radiuses of the two extremities
// the 2nd radius is optional: if it is not defined, the 2nd center is the vertex of the cone
template <typename T, typename S>
SmoothIndicatorCone3D<T,S>::SmoothIndicatorCone3D(std::vector<S> center1,
  std::vector<S> center2, S radius1, S radius2, S epsilon)
  : SmoothIndicatorF3D<T,S>(1), _center1(center1), _center2(center2),
    _radius1(radius1), _radius2(radius2), _epsilon(epsilon)
{  // _I,_J,_K is the new base where _K is the axe of the cone

  // _K = centre2 - centre1 (normalized)
  _length = sqrt( (_center2[0]-_center1[0]) * (_center2[0]-_center1[0])
                 +(_center2[1]-_center1[1]) * (_center2[1]-_center1[1])
                 +(_center2[2]-_center1[2]) * (_center2[2]-_center1[2]) );
  S tempk[] = {(_center2[0] - _center1[0])/_length,
               (_center2[1] - _center1[1])/_length,
               (_center2[2] - _center1[2])/_length};
  _K.insert(_K.begin(),tempk, tempk +3);  // _K = centre2 - centre1 (normalized)

  // _I and _J form an orthonormal base with _K
  if ( _center2[1]-_center1[1] == 0 && _center2[0]-_center1[0] == 0 ) {
    if ( _center2[2]-_center1[2] == 0 ) {
      std::cout << "Warning: in the cone, the two center have the same coordinates";
    }
	  S tempi[] = {1,0,0};
	  S tempj[] = {0,1,0};
	  _I.insert(_I.begin(),tempi, tempi +3);
	  _J.insert(_J.begin(),tempj, tempj +3);
  }
  else {
	  S normi = sqrt (_K[1]*_K[1] + _K[0]*_K[0]);
	  S tempi[] = {-_K[1]/normi, _K[0]/normi,0};
	  _I.insert(_I.begin(),tempi, tempi +3);
	  S tempj[] = {_K[1]*_I[2] - _K[2]*_I[1],
                 _K[2]*_I[0] - _K[0]*_I[2],
                 _K[0]*_I[1] - _K[1]*_I[0] };
	  _J.insert(_J.begin(),tempj, tempj +3);
  }

  S maxx, maxy, maxz;
  S minx, miny, minz;

  maxx= _center1[0] + std::max( sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                                sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);
  minx= _center1[0] + std::min(-sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                               -sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);

  maxy= _center1[1] + std::max( sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                                sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);
  miny= _center1[1] + std::min(-sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                               -sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);

  maxz= _center1[2] + std::max( sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                                sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);
  minz= _center1[2] + std::min(-sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                               -sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);

  S temp0[] = {minx, miny, minz};
  this->_myMin.insert( this->_myMin.begin(),temp0, temp0 +3 );
  S temp1[] = {maxx, maxy, maxz};
  this->_myMax.insert( this->_myMax.begin(),temp1, temp1 +3 );
}

// returns true if x is inside the coneSmooth
template <typename T, typename S>
std::vector<T> SmoothIndicatorCone3D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  T result = T();
  // radius: the radius of the cone at the point x

  double X = _I[0]*(x[0]-_center1[0]) + _I[1]*(x[1]-_center1[1]) + _I[2]*(x[2]-_center1[2]);
  double Y = _J[0]*(x[0]-_center1[0]) + _J[1]*(x[1]-_center1[1]) + _J[2]*(x[2]-_center1[2]);
  double Z = _K[0]*(x[0]-_center1[0]) + _K[1]*(x[1]-_center1[1]) + _K[2]*(x[2]-_center1[2]);
  double radius = _radius1 + (_radius2 - _radius1)*Z/_length ;
  double d = -1;        // distance from the point to the cone
  double temp = X*X + Y*Y;

  // case 1
  if ( Z <= _length && Z >= 0 && temp <= radius*radius ) {
    d = 0;
  }
  // case 2
  if (Z > _length && Z <= _length + _epsilon && temp <= _radius2*_radius2) {
    d = Z - _length;
  }
  // case 3
  if (Z >= -_epsilon && Z < 0 && temp <= _radius1*_radius1) {
    d = -Z;
  }
  // case 4
  if( temp > radius*radius
      && (_radius1-_radius2)*sqrt(temp) <= Z*_length + _radius1*(_radius1-_radius2)
		  && (_radius1-_radius2)*sqrt(temp) >= (Z-_length)*_length + _radius2*(_radius1-_radius2) )
  {
	  d=(sqrt(temp) - radius) / sqrt( std::pow((_radius1-_radius2)/_length,2)+1 );
  }
  // case 5
  if( Z <= _length + _epsilon
		  && temp > _radius2*_radius2
		  && temp <= (_radius2+_epsilon)*(_radius2+_epsilon)
		  && (_radius1-_radius2)*sqrt(temp) < (Z-_length)*_length + _radius2*(_radius1-_radius2) )
  {
    d = sqrt(std::pow(Z-_length,2) + std::pow((sqrt(temp) - _radius2),2));
  }
  // case 6
  if( Z >= - _epsilon
		  && temp > _radius1*_radius1
		  && temp <= (_radius1+_epsilon)*(_radius1+_epsilon)
		  && (_radius1 -_radius2)*sqrt(temp) > Z*_length + _radius1*(_radius1-_radius2) )
  {
    d = sqrt( Z*Z + std::pow((sqrt(temp) - _radius1) , 2) );
  }

  if (d>=0 && d<=_epsilon) {
    result = S( cos(M_PI*d/(2*_epsilon)) * cos(M_PI*d/(2*_epsilon)) );
  }
  y.push_back(result);

  return y;
}


} // namespace olb

#endif

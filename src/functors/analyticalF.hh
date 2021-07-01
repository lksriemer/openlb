/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_F_HH
#define ANALYTICAL_F_HH

#include<vector>
#include<cmath>     // for lpnorm
#include<stdlib.h>  // for random
#include <iostream>

#include "functors/analyticalF.h"
#include "core/singleton.h"
#include "communication/mpiManager.h"

namespace olb {


//////////////////////////////////1D////////////////////////////////////////////
template <typename T, typename S>
AnalyticalConst1D<T,S>::AnalyticalConst1D(std::vector<T>& value)
  : AnalyticalF1D<T,S>(value.size()), _c(value)
{
  this->_name = "const";
}

template <typename T, typename S>
AnalyticalConst1D<T,S>::AnalyticalConst1D(T value) : AnalyticalF1D<T,S>(1)
{
  _c.push_back(value);
}

template <typename T, typename S>
std::vector<T> AnalyticalConst1D<T,S>::operator()(std::vector<S> x)
{
  return _c;
}


template <typename T, typename S>
AnalyticalLinear1D<T,S>::AnalyticalLinear1D(T a, T b) 
  : AnalyticalF1D<T,S>(1), _a(a), _b(b)
{
  this->_name = "linear";
}

template <typename T, typename S>
AnalyticalLinear1D<T,S>::AnalyticalLinear1D(S x0, T v0, S x1, T v1)
  : AnalyticalF1D<T,S>(1)
{
  if ( x1-x0 == 0 )
    std::cout << "Error: x1-x2=0" << std::endl;
  else {
    _a = ( v1-v0 ) / ( x1-x0 );
    _b = v0 - _a*x0;
  }
  this->_name = "linear";
}

template <typename T, typename S>
std::vector<T> AnalyticalLinear1D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  y.push_back( _a*x[0] + _b );
  return y;
}


template <typename T, typename S>
AnalyticalRandom1D<T,S>::AnalyticalRandom1D() : AnalyticalF1D<T,S>(1)
{
  this->_name = "random";
}

template <typename T, typename S>
std::vector<T> AnalyticalRandom1D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  #ifdef PARALLEL_MODE_MPI
    int rank = singleton::mpi().getRank();
    y.push_back( ((random()*(rank + 1))%RAND_MAX)/(T)RAND_MAX );
  #else
    y.push_back( (T)random()/(T)RAND_MAX );
  #endif
  return y;
}


template <typename T, typename S>
AnalyticalSquare1D<T,S>::AnalyticalSquare1D(S cp, S r, T maxi)
  : AnalyticalF1D<T,S>(1), _cp(cp), _r(r), _maxi(maxi)
{
  this->_name = "square";
}

template <typename T, typename S>
std::vector<T> AnalyticalSquare1D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  y.push_back( _maxi*(1-(x[0]-_cp) * (x[0]-_cp) / (T)_r / (T)_r) );
  return y;
}



/////////////////////////////////someOtherFunctors//////////////////////////////
template <typename T, typename S>
PolynomialStartScale<T,S>::PolynomialStartScale(S numTimeSteps, T maxValue)
  : AnalyticalF1D<T,S>(1), _numTimeSteps(numTimeSteps), _maxValue(maxValue)
{
  this->_name = "polyStartScale";
}

template <typename T, typename S>
std::vector<T> PolynomialStartScale<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  y.push_back( (T) x[0] / (T) _numTimeSteps);
  y[0] = _maxValue * y[0]*y[0]*y[0] * ( 10.0 + y[0] * (6.0*y[0] - 15.0) );
  return y;
}


template <typename T, typename S>
SinusStartScale<T,S>::SinusStartScale(int numTimeSteps, T maxValue)
  : AnalyticalF1D<T,S>(1), _numTimeSteps(numTimeSteps), _maxValue(maxValue),
    _pi(4.0*atan(1.0))
{
  this->_name = "sinusStartScale";
}

template <typename T, typename S>
std::vector<T> SinusStartScale<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  y.push_back( (_maxValue * (sin(-_pi / 2.0 + (T)x[0] / (T)_numTimeSteps * _pi) + 1.0)) / 2.0);
  return y;
}



///////////////////////////////////////2D///////////////////////////////////////
template <typename T, typename S>
AnalyticalComposed2D<T,S>::AnalyticalComposed2D( AnalyticalF2D<T,S>& f0,
  AnalyticalF2D<T,S>& f1)
  : AnalyticalF2D<T,S>(2), _f0(f0), _f1(f1)
{
  this->_name = "composed";
}

template <typename T, typename S>
std::vector<T> AnalyticalComposed2D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  y.push_back(_f0(x)[0]);
  y.push_back(_f1(x)[0]);
  return y;
}


template <typename T, typename S>
AnalyticalConst2D<T,S>::AnalyticalConst2D(std::vector<T>& value)
  : AnalyticalF2D<T,S>(value.size()), _c(value)
{
  this->_name = "const";
}

template <typename T, typename S>
AnalyticalConst2D<T,S>::AnalyticalConst2D(T value) : AnalyticalF2D<T,S>(1)
{
  _c.push_back(value);
}

template <typename T, typename S>
AnalyticalConst2D<T,S>::AnalyticalConst2D(T value0, T value1) : AnalyticalF2D<T,S>(2)
{ 
  _c.push_back(value0);
  _c.push_back(value1); 
}

template <typename T, typename S>
std::vector<T> AnalyticalConst2D<T,S>::operator()(std::vector<S> x)
{
  return _c;
}


template <typename T, typename S>
AnalyticalLinear2D<T,S>::AnalyticalLinear2D(T a, T b, T c)
  : AnalyticalF2D<T,S>(1), _a(a), _b(b), _c(c)
{
  this->_name = "linear";
}

template <typename T, typename S>
AnalyticalLinear2D<T,S>::AnalyticalLinear2D(S x0, S y0, T v0, S x1, S y1,
  T v1, S x2, S y2, T v2)
  : AnalyticalF2D<T,S>(1)
{
  this->_name = "linear";
  T n2= (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
  if ( n2 == 0 )
    std::cout << "Error function" << std::endl;
  else
  {
    T n0 = (y1-y0)*(v2-v0) - (v1-v0)*(y2-y0);
    T n1 = (v1-v0)*(x2-x0) - (x1-x0)*(v2-v0);
    _a = -n0 / n2;
    _b = -n1 / n2;
    _c = (x0*n0 + y0*n1 + v0*n2) / n2;
  }
}

template <typename T, typename S>
std::vector<T> AnalyticalLinear2D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  y.push_back( _a*x[0] + _b*x[1] + _c );
  return y;
}


template <typename T, typename S>
AnalyticalRandom2D<T,S>::AnalyticalRandom2D() : AnalyticalF2D<T,S>(1)
{
  this->_name = "random";
}

template <typename T, typename S>
std::vector<T> AnalyticalRandom2D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  #ifdef PARALLEL_MODE_MPI
    int rank = singleton::mpi().getRank();
    y.push_back( ((random()*(rank + 1))%RAND_MAX)/(T)RAND_MAX );
  #else
    y.push_back( (T)random()/(T)RAND_MAX );
  #endif
  return y;
}




///////////////////////////////////////3D///////////////////////////////////////
template <typename T, typename S>
AnalyticalComposed3D<T,S>::AnalyticalComposed3D(AnalyticalF3D<T,S>& f0,
  AnalyticalF3D<T,S>& f1, AnalyticalF3D<T,S>& f2)
  : AnalyticalF3D<T,S>(3), _f0(f0), _f1(f1), _f2(f2)
{
  this->_name = "composed";
}

template <typename T, typename S>
std::vector<T> AnalyticalComposed3D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  y.push_back(_f0(x)[0]);
  y.push_back(_f1(x)[0]);
  y.push_back(_f2(x)[0]);
  return y;
}


template <typename T, typename S>
AnalyticalConst3D<T,S>::AnalyticalConst3D(std::vector<T>& value)
  : AnalyticalF3D<T,S>(value.size()), _c(value) { }

template <typename T, typename S>
AnalyticalConst3D<T,S>::AnalyticalConst3D(T value) : AnalyticalF3D<T,S>(1)
{
  _c.push_back(value);
}

template <typename T, typename S>
AnalyticalConst3D<T,S>::AnalyticalConst3D(T value0, T value1) : AnalyticalF3D<T,S>(2)
{
  _c.push_back(value0);
  _c.push_back(value1); 
}

template <typename T, typename S>
AnalyticalConst3D<T,S>::AnalyticalConst3D(T value0, T value1, T value2)
  : AnalyticalF3D<T,S>(3)
{ 
  _c.push_back(value0);
  _c.push_back(value1); 
  _c.push_back(value2); 
}

template <typename T, typename S>
std::vector<T> AnalyticalConst3D<T,S>::operator()(std::vector<S> x)
{
  return _c;
}


template <typename T, typename S>
AnalyticalLinear3D<T,S>::AnalyticalLinear3D(T a, T b, T c, T d)
  : AnalyticalF3D<T,S>(1), _a(a), _b(b), _c(c), _d(d)
{
  this->_name = "linear";
}

template <typename T, typename S>
AnalyticalLinear3D<T,S>::AnalyticalLinear3D(S x0, S y0, S z0, T v0, S x1,
  S y1, S z1, T v1, S x2, S y2, S z2, T v2, S x3, S y3, S z3, T v3)
  : AnalyticalF3D<T,S>(1)
{
  this->_name = "linear";
  T n = ( (y3-y0)*(x1-x0)-(x3-x0)*(y1-y0) ) * ( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) )
       +( (z3-z0)*(x1-x0)-(x3-x0)*(z1-z0) ) * ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
  if ( n == 0 )
    std::cout << "Error function" << std::endl;
  else
  {
    T w = ( (y1-y0)*(x3-x0)-(x1-x0)*(y3-y0) ) * ( (v2-v0)-(x2-x0)*(v1-v0) / (x1-x0) )
         /( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) ) + (v3-v0) - (x3-x0)*(v1-v0) / (x1-x0);
    T zx = (y1-y0)*( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) )
          -(z1-z0)*( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T rx = (v1-v0)/(x1-x0) - (y1-y0)*(v2-v0) / ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) )
          +(y1-y0)*(x2-x0)*(v1-v0) / ( (y2-y0)*(x1-x0)*(x1-x0)-(x2-x0)*(y1-y0)*(x1-x0) );
    T zy = (x1-x0)*( (x2-x0)*(z1-z0)-(z2-z0)*(x1-x0) );
    T ry = ( (x1-x0)*(v2-v0)-(x2-x0)*(v1-v0) ) / ( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T zz = (x1-x0)*( (y2-y0)*(x1-x0)-(x2-x0)*(y1-y0) );
    T h = w/n;
    _a = rx + zx*h;
    _b = ry + zy*h;
    _c = zz*h;
    _d = v0 - x0*_a - y0*_b - z0*_c;
  }
}

template <typename T, typename S>
std::vector<T> AnalyticalLinear3D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  y.push_back( _a*x[0] + _b*x[1] + _c*x[2] + _d );
  return y;
}


template <typename T, typename S>
AnalyticalRandom3D<T,S>::AnalyticalRandom3D() : AnalyticalF3D<T,S>(1)
{
  this->_name = "random";
}

template <typename T, typename S>
std::vector<T> AnalyticalRandom3D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y;
  #ifdef PARALLEL_MODE_MPI
    int rank = singleton::mpi().getRank();
    y.push_back( ((random()*(rank + 1))%RAND_MAX)/(T)RAND_MAX );
  #else
    y.push_back( (T)random()/(T)RAND_MAX );
  #endif
  return y;
}


template <typename T, typename S>
AnalyticalScaled3D<T,S>::AnalyticalScaled3D(AnalyticalF3D<T,S>& f, T scale)
  : AnalyticalF3D<T,S>(f.getTargetDim()), _f(f), _scale(scale)
{
  this->_name = "scaled";
}

template <typename T, typename S>
std::vector<T> AnalyticalScaled3D<T,S>::operator()(std::vector<S> x)
{
  std::vector<T> y(this->getTargetDim(),T());
  for (int iDim = 0; iDim < this->getTargetDim(); iDim++)
    y[iDim] = _scale*_f(x)[iDim];
  return y;
}



} // end namespace olb

#endif

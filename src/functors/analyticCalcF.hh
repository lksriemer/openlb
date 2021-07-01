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

#ifndef ANALYTICAL_CALC_F_HH
#define ANALYTICAL_CALC_F_HH

#include<vector>

#include "functors/analyticCalcF.h"
#include "functors/analyticalBaseF.h"
#include "functors/genericF.h"


namespace olb {



//////////////////////////////// AnalyticCalc1D ////////////////////////////////
template <typename T, typename S>
AnalyticCalc1D<T,S>::AnalyticCalc1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g)
  : AnalyticalF1D<T,S>(f.getTargetDim()), _f(f), _g(g) { }

template <typename T, typename S>
void AnalyticCalc1D<T,S>::myErase(GenericF<T,S>* ptr) {
  // delete child...
  this->removeChild(ptr);
  delete ptr;
  // ... and delete myself from father's list
  if( this->_pointerVec.size() == 0 )  _f.myErase(this);
}


template <typename T, typename S>
AnalyticPlus1D<T,S>::AnalyticPlus1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g)
  : AnalyticCalc1D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticPlus1D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] + this->_g(input)[i] );
  }
  // start deleting PlusF and GenericF objects
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticMinus1D<T,S>::AnalyticMinus1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g)
  : AnalyticCalc1D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticMinus1D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] - this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0)  this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticMultiplication1D<T,S>::AnalyticMultiplication1D(AnalyticalF1D<T,S>& f,
  AnalyticalF1D<T,S>& g) : AnalyticCalc1D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticMultiplication1D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] * this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticDivision1D<T,S>::AnalyticDivision1D(AnalyticalF1D<T,S>& f,
  AnalyticalF1D<T,S>& g) : AnalyticCalc1D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticDivision1D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] / this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator+(AnalyticalF1D<T,S>& rhs) {
  AnalyticalF1D<T,S>* tmp = new AnalyticPlus1D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator-(AnalyticalF1D<T,S>& rhs) {
  AnalyticalF1D<T,S>* tmp = new AnalyticMinus1D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator*(AnalyticalF1D<T,S>& rhs) {
  AnalyticalF1D<T,S>* tmp = new AnalyticMultiplication1D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator/(AnalyticalF1D<T,S>& rhs) {
  AnalyticalF1D<T,S>* tmp = new AnalyticDivision1D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}




//////////////////////////////// AnalyticCalc2D ////////////////////////////////
template <typename T, typename S>
AnalyticCalc2D<T,S>::AnalyticCalc2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g)
  : AnalyticalF2D<T,S>(f.getTargetDim()), _f(f), _g(g) { }

template <typename T, typename S>
void AnalyticCalc2D<T,S>::myErase(GenericF<T,S>* ptr) {
  // delete child...
  this->removeChild(ptr);
  delete ptr;
  // ... and delete myself from father's list
  if( this->_pointerVec.size() == 0 )  _f.myErase(this);
}


template <typename T, typename S>
AnalyticPlus2D<T,S>::AnalyticPlus2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g)
  : AnalyticCalc2D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticPlus2D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] + this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticMinus2D<T,S>::AnalyticMinus2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g)
  : AnalyticCalc2D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticMinus2D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] - this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticMultiplication2D<T,S>::AnalyticMultiplication2D(AnalyticalF2D<T,S>& f,
  AnalyticalF2D<T,S>& g) : AnalyticCalc2D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticMultiplication2D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] * this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticDivision2D<T,S>::AnalyticDivision2D(AnalyticalF2D<T,S>& f,
  AnalyticalF2D<T,S>& g) : AnalyticCalc2D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticDivision2D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] / this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator+(AnalyticalF2D<T,S>& rhs) {
  AnalyticalF2D<T,S>* tmp = new AnalyticPlus2D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator-(AnalyticalF2D<T,S>& rhs) {
  AnalyticalF2D<T,S>* tmp = new AnalyticMinus2D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator*(AnalyticalF2D<T,S>& rhs) {
  AnalyticalF2D<T,S>* tmp = new AnalyticMultiplication2D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator/(AnalyticalF2D<T,S>& rhs) {
  AnalyticalF2D<T,S>* tmp = new AnalyticDivision2D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}




//////////////////////////////// AnalyticCalc3D ////////////////////////////////
template <typename T, typename S>
AnalyticCalc3D<T,S>::AnalyticCalc3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g)
  : AnalyticalF3D<T,S>(f.getTargetDim()), _f(f), _g(g) { }

template <typename T, typename S>
void AnalyticCalc3D<T,S>::myErase(GenericF<T,S>* ptr) {
  // delete child...
  this->removeChild(ptr);
  delete ptr;
  // ... and delete myself from father's list
  if( this->_pointerVec.size() == 0 )  _f.myErase(this);
}


template <typename T, typename S>
AnalyticPlus3D<T,S>::AnalyticPlus3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g)
  : AnalyticCalc3D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticPlus3D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] + this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticMinus3D<T,S>::AnalyticMinus3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g)
  : AnalyticCalc3D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticMinus3D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] - this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename S>
AnalyticMultiplication3D<T,S>::AnalyticMultiplication3D(AnalyticalF3D<T,S>& f,
  AnalyticalF3D<T,S>& g) : AnalyticCalc3D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticMultiplication3D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] * this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticDivision3D<T,S>::AnalyticDivision3D(AnalyticalF3D<T,S>& f,
  AnalyticalF3D<T,S>& g) : AnalyticCalc3D<T,S>(f,g) { }

template <typename T, typename S>
std::vector<T> AnalyticDivision3D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] / this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator/(AnalyticalF3D<T,S>& rhs) {
  AnalyticalF3D<T,S>* tmp = new AnalyticDivision3D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}
template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator*(AnalyticalF3D<T,S>& rhs) {
  AnalyticalF3D<T,S>* tmp = new AnalyticMultiplication3D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator-(AnalyticalF3D<T,S>& rhs) {
  AnalyticalF3D<T,S>* tmp = new AnalyticMinus3D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator+(AnalyticalF3D<T,S>& rhs) {
  AnalyticalF3D<T,S>* tmp = new AnalyticPlus3D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

} // end namespace olb

#endif

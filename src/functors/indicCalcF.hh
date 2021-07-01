/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause, Cyril Masquelier,
 *  Albert Mink
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

#ifndef INDIC_CALC_F_HH
#define INDIC_CALC_F_HH


#include<vector>

#include "functors/indicCalcF.h"
#include "functors/indicatorBaseF.h"
#include "functors/genericF.h"

namespace olb {



//////////////////////////////// IndicCalc1D ////////////////////////////////
template <typename T, typename S>
IndicCalc1D<T,S>::IndicCalc1D(IndicatorF1D<T,S>& f, IndicatorF1D<T,S>& g)
  : IndicatorF1D<T,S>(f.getTargetDim()), _f(f), _g(g) { }

template <typename T, typename S>
void IndicCalc1D<T,S>::myErase(GenericF<T,S>* ptr) {
  // delete child...
  this->removeChild(ptr);
  delete ptr;
  // ... and delete myself from father's list
  if( this->_pointerVec.size() == 0 ) _f.myErase(this);
}


template <typename T, typename S>
IndicPlus1D<T,S>::IndicPlus1D(IndicatorF1D<T,S>& f, IndicatorF1D<T,S>& g)
  : IndicCalc1D<T,S>(f,g) { }

// returns 1 if( f==1 || g==1 ) UNION
template <typename T, typename S>
std::vector<T> IndicPlus1D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
	  int tmp = (this->_f(input)[i]);
	  if (tmp == 0) {
      tmp = this->_g(input)[i];
    }
    output.push_back( tmp );
  }
  // start deleting PlusF and GenericF objects
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
IndicMinus1D<T,S>::IndicMinus1D(IndicatorF1D<T,S>& f, IndicatorF1D<T,S>& g)
  : IndicCalc1D<T,S>(f,g) { }

// returns 1 if( f==1 && g==0 ) WITHOUT
template <typename T, typename S>
std::vector<T> IndicMinus1D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
	  if (this->_f(input)[i] == 1) {
      output.push_back( this->_f(input)[i] - this->_g(input)[i] );
	  }
	  else {
      output.push_back(0);
    }
  }
  if (this->_pointerVec.size() == 0)  this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
IndicMultiplication1D<T,S>::IndicMultiplication1D(IndicatorF1D<T,S>& f, IndicatorF1D<T,S>& g)
  : IndicCalc1D<T,S>(f,g) { }

// returns 1 if( f==1 && g==1 ) INTERSECTION
template <typename T, typename S>
std::vector<T> IndicMultiplication1D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] * this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}



template <typename T, typename S>
IndicatorF1D<T,S>& IndicatorF1D<T,S>::operator+(IndicatorF1D<T,S>& rhs) {
  IndicatorF1D<T,S>* tmp = new IndicPlus1D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
IndicatorF1D<T,S>& IndicatorF1D<T,S>::operator-(IndicatorF1D<T,S>& rhs) {
  IndicatorF1D<T,S>* tmp = new IndicMinus1D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
IndicatorF1D<T,S>& IndicatorF1D<T,S>::operator*(IndicatorF1D<T,S>& rhs) {
  IndicatorF1D<T,S>* tmp = new IndicMultiplication1D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}



//////////////////////////////// IndicCalc2D ////////////////////////////////
template <typename T, typename S>
IndicCalc2D<T,S>::IndicCalc2D(IndicatorF2D<T,S>& f, IndicatorF2D<T,S>& g)
  : IndicatorF2D<T,S>(f.getTargetDim()), _f(f), _g(g) { }

template <typename T, typename S>
void IndicCalc2D<T,S>::myErase(GenericF<T,S>* ptr) {
  this->removeChild(ptr);
  delete ptr;
  if( this->_pointerVec.size() == 0 ) _f.myErase(this);
}


template <typename T, typename S>
IndicPlus2D<T,S>::IndicPlus2D(IndicatorF2D<T,S>& f, IndicatorF2D<T,S>& g)
  : IndicCalc2D<T,S>(f,g) { }

// returns 1 if( f==1 || g==1 ) UNION
template <typename T, typename S>
std::vector<T> IndicPlus2D<T,S>::operator()(std::vector<S> input) {
	std::vector<T> output;
	for(int i = 0; i < this->_f.getTargetDim(); i++) {
	  int tmp = (this->_f(input)[i]);
	  if (tmp == 0) {
      tmp = this->_g(input)[i];
    }
    output.push_back( tmp );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
IndicMinus2D<T,S>::IndicMinus2D(IndicatorF2D<T,S>& f, IndicatorF2D<T,S>& g)
  : IndicCalc2D<T,S>(f,g) { }

// returns 1 if( f==1 && g==0 ) WITHOUT
template <typename T, typename S>
std::vector<T> IndicMinus2D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
	  if (this->_f(input)[i] == 1) {
		  output.push_back( this->_f(input)[i] - this->_g(input)[i] );
	  }
	  else {
      output.push_back(0);
    }
  }

  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
IndicMultiplication2D<T,S>::IndicMultiplication2D(IndicatorF2D<T,S>& f, IndicatorF2D<T,S>& g)
  : IndicCalc2D<T,S>(f,g) { }

// returns 1 if( f==1 && g==1 ) INTERSECTION
template <typename T, typename S>
std::vector<T> IndicMultiplication2D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] * this->_g(input)[i] );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}



template <typename T, typename S>
IndicatorF2D<T,S>& IndicatorF2D<T,S>::operator+(IndicatorF2D<T,S>& rhs) {
  IndicatorF2D<T,S>* tmp = new IndicPlus2D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
IndicatorF2D<T,S>& IndicatorF2D<T,S>::operator-(IndicatorF2D<T,S>& rhs) {
  IndicatorF2D<T,S>* tmp = new IndicMinus2D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
IndicatorF2D<T,S>& IndicatorF2D<T,S>::operator*(IndicatorF2D<T,S>& rhs) {
  IndicatorF2D<T,S>* tmp = new IndicMultiplication2D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}



//////////////////////////////// IndicCalc3D ////////////////////////////////
template <typename T, typename S>
IndicCalc3D<T,S>::IndicCalc3D(IndicatorF3D<T,S>& f, IndicatorF3D<T,S>& g)
: IndicatorF3D<T,S>(f.getTargetDim()), _f(f), _g(g) { }


template <typename T, typename S>
void IndicCalc3D<T,S>::myErase(GenericF<T,S>* ptr) {
  this->removeChild(ptr);
  delete ptr;
  if( this->_pointerVec.size() == 0 ) _f.myErase(this);
}


template <typename T, typename S>
IndicPlus3D<T,S>::IndicPlus3D(IndicatorF3D<T,S>& f, IndicatorF3D<T,S>& g)
  : IndicCalc3D<T,S>(f,g) { 
  for (int iDim = 0; iDim<3; iDim++) {
    this->_myMin.push_back(std::min(f.getMin()[iDim],g.getMin()[iDim]) );
    this->_myMax.push_back(std::max(f.getMax()[iDim],g.getMax()[iDim]) );
  }
}

// returns 1 if( f==1 || g==1 ) UNION
template <typename T, typename S>
std::vector<T> IndicPlus3D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    int tmp = (this->_f(input)[i]);
	  if (tmp==0) {
      tmp= this->_g(input)[i];
    }
    output.push_back( tmp );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
IndicMinus3D<T,S>::IndicMinus3D(IndicatorF3D<T,S>& f, IndicatorF3D<T,S>& g)
  : IndicCalc3D<T,S>(f,g) {
	for (int iDim = 0; iDim<3; iDim++) {
	  this->_myMin.push_back(f.getMin()[iDim]);
	  this->_myMax.push_back(f.getMax()[iDim]);
	}
}

// returns 1 if( f==1 && g==0 ) WITHOUT
template <typename T, typename S>
std::vector<T> IndicMinus3D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    if (this->_f(input)[i] == 1) {
	    output.push_back( this->_f(input)[i] - this->_g(input)[i] );
    }
	  else {output.push_back(0);}
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


// returns 1 if( f==1 && g==1 ) INTERSECTION
template <typename T, typename S>
IndicMultiplication3D<T,S>::IndicMultiplication3D(IndicatorF3D<T,S>& f, IndicatorF3D<T,S>& g)
  : IndicCalc3D<T,S>(f,g) {
	for (int iDim = 0; iDim < 3; iDim++) {
	  this->_myMin.push_back(std::max(f.getMin()[iDim],g.getMin()[iDim]) );
	  this->_myMax.push_back(std::min(f.getMax()[iDim],g.getMax()[iDim]) );
	}
}

template <typename T, typename S>
std::vector<T> IndicMultiplication3D<T,S>::operator()(std::vector<S> input) {
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] * this->_g(input)[i] );
  }
  this->_myMin = std::max(this->_f.getMin(), this->_g.getMin());
  this->_myMax = std::min(this->_f.getMax(), this->_g.getMax());
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}



template <typename T, typename S>
IndicatorF3D<T,S>& IndicatorF3D<T,S>::operator*(IndicatorF3D<T,S>& rhs) {
  IndicatorF3D<T,S>* tmp = new IndicMultiplication3D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
IndicatorF3D<T,S>& IndicatorF3D<T,S>::operator-(IndicatorF3D<T,S>& rhs) {
  IndicatorF3D<T,S>* tmp = new IndicMinus3D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, typename S>
IndicatorF3D<T,S>& IndicatorF3D<T,S>::operator+(IndicatorF3D<T,S>& rhs) {
  IndicatorF3D<T,S>* tmp = new IndicPlus3D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}




//////////////////////////////// IndicSmoothCalc3D ////////////////////////////////
template <typename T, typename S>
SmoothIndicCalc3D<T,S>::SmoothIndicCalc3D(SmoothIndicatorF3D<T,S>& f,
  SmoothIndicatorF3D<T,S>& g)
  : 
SmoothIndicatorF3D<T,S>(f.getTargetDim()), _f(f), _g(g) { }

template <typename T, typename S>
void SmoothIndicCalc3D<T,S>::myErase(GenericF<T,S>* ptr) {
  this->removeChild(ptr);
  delete ptr;
  if( this->_pointerVec.size() == 0 )  _f.myErase(this);
}


template <typename T, typename S>
SmoothIndicPlus3D<T,S>::SmoothIndicPlus3D(SmoothIndicatorF3D<T,S>& f,
  SmoothIndicatorF3D<T,S>& g)
  : SmoothIndicCalc3D<T,S>(f,g)
{
  for (int iDim = 0; iDim < 3; iDim++) {
    this->_myMin.push_back(std::min(f.getMin()[iDim],g.getMin()[iDim]) );
    this->_myMax.push_back(std::max(f.getMax()[iDim],g.getMax()[iDim]) );
  }
}

// returns 1 if( f==1 || g==1 ) UNION
template <typename T, typename S>
std::vector<T> SmoothIndicPlus3D<T,S>::operator()(std::vector<S> input) {
  std::vector<S> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    int tmp = std::max(this->_f(input)[i], this->_g(input)[i]);
	  output.push_back( tmp );
  }
  if (this->_pointerVec.size() == 0) this->_f.myErase(this);
  return output;
}


template <typename T, typename S>
SmoothIndicatorF3D<T,S>& SmoothIndicatorF3D<T,S>::operator+(SmoothIndicatorF3D<T,S>& rhs) {
  SmoothIndicatorF3D<T,S>* tmp = new SmoothIndicPlus3D<T,S>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}



} // end namespace olb

#endif

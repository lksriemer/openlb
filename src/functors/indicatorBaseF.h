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

#ifndef INDICATOR_BASE_F_H
#define INDICATOR_BASE_F_H

#include<vector>
#include<cmath>

#include "functors/genericF.h"


namespace olb {


//////////////////////////////////////IndicatorF////////////////////////////////

template <typename T, typename S>
class IndicatorF1D : public GenericF<T,S> {
protected:
  IndicatorF1D(int n) : GenericF<T,S>(n,1) { }
public:
  IndicatorF1D<T,S>& operator-(IndicatorF1D<T,S>& rhs);
  IndicatorF1D<T,S>& operator+(IndicatorF1D<T,S>& rhs);
  IndicatorF1D<T,S>& operator*(IndicatorF1D<T,S>& rhs);
};


template <typename T, typename S>
class IndicatorF2D : public GenericF<T,S> {
protected:
  IndicatorF2D(int n) : GenericF<T,S>(n,2) { }
  std::vector<S> _myMin;
  std::vector<S> _myMax;
public:
  virtual std::vector<S>& getMin() {return _myMin; };
  virtual std::vector<S>& getMax() {return _myMax; };
  /// returns false or true and pos. distance if there was one found for an given origin and direction, mind that the default computation is done by an numerical approximation which searches .. [TODO]
  virtual bool distance(S& distance, const std::vector<S>& origin,
                        const std::vector<S>& direction, int iC = -1);
  bool isInsideBox(std::vector<S> origin);
  IndicatorF2D<T,S>& operator-(IndicatorF2D<T,S>& rhs);
  IndicatorF2D<T,S>& operator+(IndicatorF2D<T,S>& rhs);
  IndicatorF2D<T,S>& operator*(IndicatorF2D<T,S>& rhs);
};


template <typename T, typename S>
class IndicatorF3D : public GenericF<T,S> {
protected:
  IndicatorF3D(int n) : GenericF<T,S>(n,3) { };
  std::vector<S> _myMin;
  std::vector<S> _myMax;
public:
  virtual std::vector<T> operator()(std::vector<S> in) =0;
  virtual std::vector<S>& getMin() {return _myMin; };
  virtual std::vector<S>& getMax() {return _myMax; };
  /// returns false or true and pos. distance if there was one found for an given origin and direction, mind that the default computation is done by an numerical approximation which searches .. [TODO: CYRIL]
  virtual bool distance(S& distance,const std::vector<S>& origin,
                        const std::vector<S>& direction, int iC = -1);
  bool isInsideBox(std::vector<S> origin);
  IndicatorF3D<T,S>& operator-(IndicatorF3D<T,S>& rhs);
  IndicatorF3D<T,S>& operator+(IndicatorF3D<T,S>& rhs);
  IndicatorF3D<T,S>& operator*(IndicatorF3D<T,S>& rhs);
};



/////////////////////////////////////IdentityF//////////////////////////////////
template <typename T, typename S>
class IndicatorIdentity2D : public IndicatorF2D<T,S> {
protected:
  IndicatorF2D<T,S>& _f;
public:
  IndicatorIdentity2D<T,S>(IndicatorF2D<T,S>& f);
  ~IndicatorIdentity2D<T,S>();
  // access operator should not delete f, since f still has the identity as child
  std::vector<T> operator()(std::vector<S> input);
//  virtual std::string name() { return f.name(); }
};
/// identity functor
template <typename T, typename S>
class IndicatorIdentity3D : public IndicatorF3D<T,S> {
protected:
  IndicatorF3D<T,S>& _f;
public:
  IndicatorIdentity3D<T,S>(IndicatorF3D<T,S>& f);
  ~IndicatorIdentity3D<T,S>();
  // access operator should not delete f, since f still has the identity as child
  std::vector<T> operator()(std::vector<S> input);
  //virtual void myErase(GenericF<T,S>* ptr);
//  virtual std::string name() { return f.name(); }
};



//////////////////////////////////SmoothIndicatorF//////////////////////////////

/// SmoothIndicatorF2D has no jump between 0 and 1. Its values are continuous.
template <typename T, typename S>
class SmoothIndicatorF2D : public GenericF<T,S> {
protected:
  SmoothIndicatorF2D(int n) : GenericF<T,S>(n,2) { }
  std::vector<S> _myMin;
  std::vector<S> _myMax;
public:
  virtual std::vector<S>& getMin() {return _myMin; };
  virtual std::vector<S>& getMax() {return _myMax; };
  SmoothIndicatorF2D<T,S>& operator+(SmoothIndicatorF2D<T,S>& rhs);
};


/// SmoothIndicatorF3D has no jump between 0 and 1. Its values are continuous.
template <typename T, typename S>
class SmoothIndicatorF3D : public GenericF<T,S> {
protected:
  SmoothIndicatorF3D(int n) : GenericF<T,S>(n,3) { }
  std::vector<S> _myMin;
  std::vector<S> _myMax;
public:
  virtual std::vector<S>& getMin() {return _myMin; };
  virtual std::vector<S>& getMax() {return _myMax; };
  SmoothIndicatorF3D<T,S>& operator+(SmoothIndicatorF3D<T,S>& rhs);
};



/////////////////////////////////SmoothIdentityF////////////////////////////////
template <typename T, typename S>
class SmoothIndicatorIdentity3D : public SmoothIndicatorF3D<T,S> {
protected:
  SmoothIndicatorF3D<T,S>& _f;
public:
  SmoothIndicatorIdentity3D<T,S>(SmoothIndicatorF3D<T,S>& f);
  ~SmoothIndicatorIdentity3D<T,S>();
  // access operator should not delete f, since f still has the identity as child
  std::vector<T> operator()(std::vector<S> input);
//  virtual std::string name() { return f.name(); }
};

}

#endif


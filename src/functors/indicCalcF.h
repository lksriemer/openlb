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

#ifndef INDIC_CALC_F_H
#define INDIC_CALC_F_H


#include<vector>

#include "functors/indicatorBaseF.h"
#include "functors/genericF.h"


namespace olb {


/*
 *  arithmetic helper classes for IndicatorF1D, IndicatorF3D, IndicatorF3D
 *  smoothIndicator3D
 *  UNION         +
 *  WITHOUT       -
 *  INTERSECTION  *
*/



//////////////////////////////// IndicCalc1D ////////////////////////////////
/// arithmetic helper class for Indicator 1d functors
template <typename T, typename S>
class IndicCalc1D : public IndicatorF1D<T,S> {
protected:
  IndicatorF1D<T,S>& _f;
  IndicatorF1D<T,S>& _g;
public:
  // set dimensions as well
  IndicCalc1D(IndicatorF1D<T,S>& f, IndicatorF1D<T,S>& g);
  /// memory management
  virtual void myErase(GenericF<T,S>* ptr);
};

/// addition functor acts as union
template <typename T, typename S>
class IndicPlus1D : public IndicCalc1D<T,S> {
public:
  IndicPlus1D(IndicatorF1D<T,S>& f, IndicatorF1D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};

/// subtraction functor acts as without
template <typename T, typename S>
class IndicMinus1D : public IndicCalc1D<T,S> {
public:
  IndicMinus1D(IndicatorF1D<T,S>& f, IndicatorF1D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};

/// multiplication functor acts as intersection
template <typename T, typename S>
class IndicMultiplication1D : public IndicCalc1D<T,S> {
public:
  IndicMultiplication1D(IndicatorF1D<T,S>& f, IndicatorF1D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};



//////////////////////////////// IndicCalc2D ////////////////////////////////
/// arithmetic helper class for Indicator 2d functors
template <typename T, typename S>
class IndicCalc2D : public IndicatorF2D<T,S> {
protected:
  IndicatorF2D<T,S>& _f;
  IndicatorF2D<T,S>& _g;
public:
  IndicCalc2D(IndicatorF2D<T,S>& f, IndicatorF2D<T,S>& g);
  virtual void myErase(GenericF<T,S>* ptr);
};

/// addition functor acts as union
template <typename T, typename S>
class IndicPlus2D : public IndicCalc2D<T,S> {
public:
  IndicPlus2D(IndicatorF2D<T,S>& f, IndicatorF2D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};

/// subtraction functor acts as without
template <typename T, typename S>
class IndicMinus2D : public IndicCalc2D<T,S> {
public:
  IndicMinus2D(IndicatorF2D<T,S>& f, IndicatorF2D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};

/// multiplication functor acts as intersection
template <typename T, typename S>
class IndicMultiplication2D : public IndicCalc2D<T,S> {
public:
  IndicMultiplication2D(IndicatorF2D<T,S>& f, IndicatorF2D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};



//////////////////////////////// IndicCalc3D ////////////////////////////////
/// arithmetic helper class for Indicator 3d functors
template <typename T, typename S>
class IndicCalc3D : public IndicatorF3D<T,S> {
protected:
  IndicatorF3D<T,S>& _f;
  IndicatorF3D<T,S>& _g;
public:
  IndicCalc3D(IndicatorF3D<T,S>& f, IndicatorF3D<T,S>& g);
  virtual void myErase(GenericF<T,S>* ptr);
};

/// addition functor acts as union
template <typename T, typename S>
class IndicPlus3D : public IndicCalc3D<T,S> {
public:
  IndicPlus3D(IndicatorF3D<T,S>& f, IndicatorF3D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};

/// subtraction functor acts as without
template <typename T, typename S>
class IndicMinus3D : public IndicCalc3D<T,S> {
public:
  IndicMinus3D(IndicatorF3D<T,S>& f, IndicatorF3D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};

/// multiplication functor acts as intersection
template <typename T, typename S>
class IndicMultiplication3D : public IndicCalc3D<T,S> {
public:
  IndicMultiplication3D(IndicatorF3D<T,S>& f, IndicatorF3D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};




//////////////////////////////// IndicSmoothCalc3D ////////////////////////////////
/// arithmetic helper class for Indicator 3d functors
template <typename T, typename S>
class SmoothIndicCalc3D : public SmoothIndicatorF3D<T,S> {
protected:
  SmoothIndicatorF3D<T,S>& _f;
  SmoothIndicatorF3D<T,S>& _g;
public:
  SmoothIndicCalc3D(SmoothIndicatorF3D<T,S>& f, SmoothIndicatorF3D<T,S>& g);
  virtual void myErase(GenericF<T,S>* ptr);
};

/// addition functor acts as union
template <typename T, typename S>
class SmoothIndicPlus3D : public SmoothIndicCalc3D<T,S> {
public:
  SmoothIndicPlus3D(SmoothIndicatorF3D<T,S>& f, SmoothIndicatorF3D<T,S>& g);
  std::vector<T> operator()(std::vector<S> input);
};




} // end namespace olb

#endif


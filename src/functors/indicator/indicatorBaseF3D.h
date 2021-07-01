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

#ifndef INDICATOR_BASE_F_3D_H
#define INDICATOR_BASE_F_3D_H


#include <vector>

#include "core/vector.h"
#include "functors/analyticalBaseF.h"
#include "functors/genericF.h"
#include "functors/superBaseF3D.h"
#include "geometry/superGeometry3D.h"


namespace olb {

template<typename T> class SuperF3D;
template<typename T> class SuperGeometry3D;

/** IndicatorF3D is an application from \f$ \Omega \subset R^3 \to \{0,1\} \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  */
template <typename S>
class IndicatorF3D : public GenericF<bool,S> {
protected:
  IndicatorF3D();
  Vector<S,3> _myMin;
  Vector<S,3> _myMax;
public:
  virtual Vector<S,3>& getMin();
  virtual Vector<S,3>& getMax();
  /** \returns false or true and pos. distance if there was one found for a given origin and direction.
   * Mind that the default computation is done by a numerical approximation which searches .. [TODO: CYRIL]
   */
  virtual bool distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC=-1);
  /// Returns true if `point` is inside a cube with corners `_myMin` and `_myMax`
  bool isInsideBox(Vector<S,3> point);

  /// + Operator (Union)
  IndicatorF3D<S>& operator+(IndicatorF3D<S>& rhs);
  /// - Operator (Without)
  IndicatorF3D<S>& operator-(IndicatorF3D<S>& rhs);
  /// * Operator (Intersection)
  IndicatorF3D<S>& operator*(IndicatorF3D<S>& rhs);
};


/// Base indicator functor from SuperF3D
template <typename S>
class SuperIndicatorF3D : public IndicatorF3D<S> {
protected:
  SuperF3D<S>& _superF;
public:
  SuperIndicatorF3D (SuperF3D<S>& rhs);
};


/// Base indicator functor (discrete)
/**
 * Provides Union, Without and Intersection arithmetic.
 * _Note: `operator()` must be overloaded by child classes._
 */
class DiscreteIndicatorF3D : public GenericF<bool,int> {
public:
  DiscreteIndicatorF3D();

  /// + Operator (Union)
  DiscreteIndicatorF3D& operator+(DiscreteIndicatorF3D& rhs);
  /// - Operator (Without)
  DiscreteIndicatorF3D& operator-(DiscreteIndicatorF3D& rhs);
  /// * Operator (Intersection)
  DiscreteIndicatorF3D& operator*(DiscreteIndicatorF3D& rhs);
};



/// Indicator functor that returns false for all points
class DiscreteIndicatorFalse3D : public DiscreteIndicatorF3D {
public:
  DiscreteIndicatorFalse3D();
  virtual bool operator() (bool output[], const int input[]);
};


/// Indicator functor that returns true for all points
class DiscreteIndicatorTrue3D : public DiscreteIndicatorF3D {
public:
  DiscreteIndicatorTrue3D();
  virtual bool operator() (bool output[], const int input[]);
};


/// Indicator functor from material number
template <typename S>
class DiscreteIndicatorMaterial3D : public DiscreteIndicatorF3D {
protected:
  SuperGeometry3D<S> _superGeometry;
  std::vector<int> _materialNumbers;
public:
  DiscreteIndicatorMaterial3D (SuperGeometry3D<S>& rhs, std::vector<int> materialNumbers);
  virtual bool operator() (bool output[], const int input[]);

//protected:
  //void calculateMinMax();
};





template <typename S>
class IndicatorIdentity3D : public IndicatorF3D<S> {
protected:
  IndicatorF3D<S>& _f;
public:
  IndicatorIdentity3D(IndicatorF3D<S>& f);
  bool operator() (bool output[], const S input[]) override;
};



/** SmoothIndicatorF3D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  */
template <typename T, typename S>
class SmoothIndicatorF3D : public AnalyticalF3D<T,S> {
protected:
  SmoothIndicatorF3D();
  Vector<S,3> _myMin;
  Vector<S,3> _myMax;
public:
  virtual Vector<S,3>& getMin();
  virtual Vector<S,3>& getMax();
  SmoothIndicatorF3D<T,S>& operator+(SmoothIndicatorF3D<T,S>& rhs);
};


template <typename T, typename S>
class SmoothIndicatorIdentity3D : public SmoothIndicatorF3D<T,S> {
protected:
  SmoothIndicatorF3D<T,S>& _f;
public:
  SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T,S>& f);
  bool operator() (T output[], const S input[]) override;
};

}

#endif

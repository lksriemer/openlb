/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#ifndef SUPER_CALC_F_2D_H
#define SUPER_CALC_F_2D_H


#include "functors/superBaseF2D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// arithmetic helper class for SuperLatticeF2D functors
/** Warning: Allocation error possible in functors that have multiple functor evaluation like SuperSum2D */
template <typename T>
class SuperCalc2D : public SuperF2D<T> {
protected:
  SuperF2D<T>& _f;
  SuperF2D<T>& _g;
public:
  SuperCalc2D(SuperF2D<T>& f, SuperF2D<T>& g);
};

/// addition functor
template <typename T>
class SuperPlus2D : public SuperCalc2D<T> {
public:
  SuperPlus2D(SuperF2D<T>& f, SuperF2D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// subtraction functor
template <typename T>
class SuperMinus2D : public SuperCalc2D<T> {
public:
  SuperMinus2D(SuperF2D<T>& f, SuperF2D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// multiplication functor
template <typename T>
class SuperMultiplication2D : public SuperCalc2D<T> {
public:
  SuperMultiplication2D(SuperF2D<T>& f, SuperF2D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// division functor
template <typename T>
class SuperDivision2D : public SuperCalc2D<T> {
public:
  SuperDivision2D(SuperF2D<T>& f, SuperF2D<T>& g);
  bool operator() (T output[], const int input[]);
};


} // end namespace olb

#endif

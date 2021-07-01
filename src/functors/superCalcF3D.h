/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
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

#ifndef SUPER_CALC_F_3D_H
#define SUPER_CALC_F_3D_H


#include "superBaseF3D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// arithmetic helper class for SuperF3D functors
/** Warning: Allocation error possible in functors that have multiple functor evaluation like SuperSum3D */
template <typename T>
class SuperCalc3D : public SuperF3D< T > {
protected:
  SuperF3D<T>& _f;
  SuperF3D<T>& _g;
public:
  SuperCalc3D(SuperF3D<T>& f, SuperF3D<T>& g);
};

/// addition functor
template <typename T>
class SuperPlus3D : public SuperCalc3D<T> {
public:
  SuperPlus3D(SuperF3D<T>& f, SuperF3D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// subtraction functor
template <typename T>
class SuperMinus3D : public SuperCalc3D<T> {
public:
  SuperMinus3D(SuperF3D<T>& f, SuperF3D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// multiplication functor
template <typename T>
class SuperMultiplication3D : public SuperCalc3D<T> {
public:
  SuperMultiplication3D(SuperF3D<T>& f, SuperF3D<T>& g);
  bool operator() (T output[], const int input[]);
};

/// division functor
template <typename T>
class SuperDivision3D : public SuperCalc3D<T> {
public:
  SuperDivision3D(SuperF3D<T>& f, SuperF3D<T>& g);
  bool operator() (T output[], const int input[]);
};


} // end namespace olb

#endif

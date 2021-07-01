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

#ifndef SUPER_LATTICE_CALC_F_2D_H
#define SUPER_LATTICE_CALC_F_2D_H

#include<vector>
#include "functors/genericF.h"
#include "functors/superLatticeBaseF2D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// arithmetic helper class for SuperLatticeF2D functors
/** Warning: Allocation error possible in functors that have multiple functor evaluation like SuperSum2D */
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeCalc2D : public SuperLatticeF2D< T, DESCRIPTOR > {
protected:
  SuperLatticeF2D<T,DESCRIPTOR>& _f;
  SuperLatticeF2D<T,DESCRIPTOR>& _g;
public:
  SuperLatticeCalc2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
                     SuperLatticeF2D<T,DESCRIPTOR>& g);
  virtual void myErase(GenericF<T,int>* ptr);
};

/// addition functor
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePlus2D : public SuperLatticeCalc2D<T,DESCRIPTOR> {
public:
  SuperLatticePlus2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
                     SuperLatticeF2D<T,DESCRIPTOR>& g);
  std::vector<T> operator()(std::vector<int> input);
};

/// subtraction functor
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeMinus2D : public SuperLatticeCalc2D<T,DESCRIPTOR> {
public:
  SuperLatticeMinus2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
                      SuperLatticeF2D<T,DESCRIPTOR>& g);
  std::vector<T> operator()(std::vector<int> input);
};

/// multiplication functor
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeMultiplication2D : public SuperLatticeCalc2D<T,DESCRIPTOR> {
public:
  SuperLatticeMultiplication2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
                               SuperLatticeF2D<T,DESCRIPTOR>& g);
  std::vector<T> operator()(std::vector<int> input);
};

/// division functor
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeDivision2D : public SuperLatticeCalc2D<T,DESCRIPTOR> {
public:
  SuperLatticeDivision2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
                         SuperLatticeF2D<T,DESCRIPTOR>& g);
  std::vector<T> operator()(std::vector<int> input);
};


} // end namespace olb

#endif

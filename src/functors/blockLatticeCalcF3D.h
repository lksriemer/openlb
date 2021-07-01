/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_CALC_F_3D_H
#define BLOCK_LATTICE_CALC_F_3D_H

#include<vector>
#include "functors/genericF.h"
#include "functors/blockLatticeBaseF3D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// arithmetic helper class for BlockLatticeF3D functors
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeCalc3D : public BlockLatticeF3D< T, DESCRIPTOR > {
protected:
  BlockLatticeF3D<T,DESCRIPTOR>& _f;
  BlockLatticeF3D<T,DESCRIPTOR>& _g;
public:
  BlockLatticeCalc3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
                     BlockLatticeF3D<T,DESCRIPTOR>& g);
  virtual void myErase(GenericF<T,int>* ptr);
};

/// addition functor
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePlus3D : public BlockLatticeCalc3D<T,DESCRIPTOR> {
public:
  BlockLatticePlus3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
                     BlockLatticeF3D<T,DESCRIPTOR>& g);
  std::vector<T> operator()(std::vector<int> input);
};

/// subtraction functor
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeMinus3D : public BlockLatticeCalc3D<T,DESCRIPTOR> {
public:
  BlockLatticeMinus3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
                      BlockLatticeF3D<T,DESCRIPTOR>& g);
  std::vector<T> operator()(std::vector<int> input);
};

/// multiplication functor
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeMultiplication3D : public BlockLatticeCalc3D<T,DESCRIPTOR> {
public:
  BlockLatticeMultiplication3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
                               BlockLatticeF3D<T,DESCRIPTOR>& g);
  std::vector<T> operator()(std::vector<int> input);
};

/// division functor
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeDivision3D : public BlockLatticeCalc3D<T,DESCRIPTOR> {
public:
  BlockLatticeDivision3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
                         BlockLatticeF3D<T,DESCRIPTOR>& g);
  std::vector<T> operator()(std::vector<int> input);
};


} // end namespace olb

#endif

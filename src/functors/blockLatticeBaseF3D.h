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

#ifndef BLOCK_LATTICE_BASE_F_3D_H
#define BLOCK_LATTICE_BASE_F_3D_H

#include "functors/genericF.h"
#include "core/blockLattice3D.h"
#include "core/blockLatticeStructure3D.h"
#include "core/units.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// represents all functors that operate on a Lattice in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeF3D : public GenericF<T,int> {
protected:
  BlockLatticeStructure3D<T,DESCRIPTOR>& _blockLattice;
public:
  BlockLatticeF3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int targetDim);

  BlockLatticeF3D<T,DESCRIPTOR>& operator-(BlockLatticeF3D<T,DESCRIPTOR>& rhs);
  BlockLatticeF3D<T,DESCRIPTOR>& operator+(BlockLatticeF3D<T,DESCRIPTOR>& rhs);
  BlockLatticeF3D<T,DESCRIPTOR>& operator*(BlockLatticeF3D<T,DESCRIPTOR>& rhs);
  BlockLatticeF3D<T,DESCRIPTOR>& operator/(BlockLatticeF3D<T,DESCRIPTOR>& rhs);

  BlockLatticeStructure3D<T,DESCRIPTOR>& getBlockLattice3D() { return _blockLattice; }
};

/// represents all functors that operate on a Lattice with output in Phys, e.g. physVelocity(), physForce(), physPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysF3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const LBconverter<T>& _converter;
public:
  BlockLatticePhysF3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                      const LBconverter<T>& converter, int targetDim);
};


/// identity functor
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeIdentity3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  BlockLatticeF3D<T,DESCRIPTOR>& _f;
public:
  BlockLatticeIdentity3D<T,DESCRIPTOR>(BlockLatticeF3D<T,DESCRIPTOR>& f);
  ~BlockLatticeIdentity3D<T,DESCRIPTOR>();
  // access operator should not delete f, since f still has the identity as child
  std::vector<T> operator()(std::vector<int> input);
};


} // end namespace olb

#endif

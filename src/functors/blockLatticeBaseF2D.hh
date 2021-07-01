/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Albert Mink, Lukas Baron, Mathias J. Krause
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

#ifndef BLOCK_LATTICE_BASE_F_2D_HH
#define BLOCK_LATTICE_BASE_F_2D_HH

#include "functors/blockLatticeBaseF2D.h"


namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>::BlockLatticeF2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, int targetDim)
  : GenericF<T,int>(targetDim,2), _blockLattice(blockLattice) { }


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysF2D<T,DESCRIPTOR>::BlockLatticePhysF2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter,
  int targetDim)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice, targetDim), _converter(converter) { }


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeIdentity2D<T,DESCRIPTOR>::BlockLatticeIdentity2D
  (BlockLatticeF2D<T,DESCRIPTOR>& f)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice2D(),f.getTargetDim() ), _f(f)
{
  this->_name = _f.getName();
  // add 'this' to father's child list to prevent father from being deleted
  _f.addChild(this);
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeIdentity2D<T,DESCRIPTOR>::~BlockLatticeIdentity2D()
{
  // remove 'this' from father's child list
  _f.removeChild(this);
  // delete father from grandfather's child list
  _f.myErase(NULL);
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeIdentity2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{ return _f(input); }


} // end namespace olb

#endif

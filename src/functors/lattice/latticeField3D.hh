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

#ifndef LATTICE_FIELD_3D_HH
#define LATTICE_FIELD_3D_HH

#include "latticeField3D.h"
#include "core/fields.h"

namespace olb {

template<typename T, typename DESCRIPTOR, typename FIELD>
SuperLatticeField3D<T,DESCRIPTOR,FIELD>::SuperLatticeField3D(
  SuperLattice<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, DESCRIPTOR::template size<FIELD>())
{
  this->getName() = fields::name<FIELD>();
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeField3D<T,DESCRIPTOR,FIELD>(this->_sLattice.getBlock(iC)));
  }
}

template<typename T, typename DESCRIPTOR, typename FIELD>
BlockLatticeField3D<T,DESCRIPTOR,FIELD>::BlockLatticeField3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR::template size<FIELD>())
{
  this->getName() = fields::name<FIELD>();
}

template<typename T, typename DESCRIPTOR, typename FIELD>
bool BlockLatticeField3D<T,DESCRIPTOR,FIELD>::operator()(T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).template computeField<FIELD>(output);
  return true;
}

}
#endif

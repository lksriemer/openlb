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

#ifndef SUPER_LATTICE_BASE_F_2D_HH
#define SUPER_LATTICE_BASE_F_2D_HH

#include "superLatticeBaseF2D.h"


namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF2D<T,DESCRIPTOR>::SuperLatticeF2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, int targetDim)
  : GenericF<T,int>(targetDim,3), _sLattice(sLattice) { }


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysF2D<T,DESCRIPTOR>::SuperLatticePhysF2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter,
  int targetDim)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, targetDim), _converter(converter) { }


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeIdentity2D<T,DESCRIPTOR>::SuperLatticeIdentity2D
  (SuperLatticeF2D<T,DESCRIPTOR>& f)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),f.getTargetDim() ), _f(f)
{
  this->_name = _f.getName();
  // add 'this' to father's child list to prevent father from being deleted
  _f.addChild(this);
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeIdentity2D<T,DESCRIPTOR>::~SuperLatticeIdentity2D() {
  // remove 'this' from father's child list
  _f.removeChild(this);
  // delete father from grandfather's child list
  _f.myErase(NULL);
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeIdentity2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{ return _f(input); }


} // end namespace olb

#endif

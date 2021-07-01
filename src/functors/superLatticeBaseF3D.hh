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

#ifndef SUPER_LATTICE_BASE_F_3D_HH
#define SUPER_LATTICE_BASE_F_3D_HH


#include "superLatticeBaseF3D.h"

namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF3D<T,DESCRIPTOR>::SuperLatticeF3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
  int targetDim)
  : GenericF<T,int>(targetDim,4), _sLattice(sLattice) { }


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysF3D<T,DESCRIPTOR>::SuperLatticePhysF3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter,
  int targetDim)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, targetDim), _converter(converter) { }

template <typename T, template <typename U> class DESCRIPTOR>
LBconverter<T> const& SuperLatticePhysF3D<T,DESCRIPTOR>::getConverter() const {
  return this->_converter;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeIdentity3D<T,DESCRIPTOR>::SuperLatticeIdentity3D(SuperLatticeF3D<T,DESCRIPTOR>& f)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D() ,f.getTargetDim() ), _f(f) 
{
  this->_name = _f.getName();
  // add 'this' to father's child list to prevent father from being deleted
  _f.addChild(this);
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeIdentity3D<T,DESCRIPTOR>::~SuperLatticeIdentity3D() {
  // remove 'this' from father's child list
  _f.removeChild(this);
  // delete father from grandfather's child list
  _f.myErase(NULL);
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeIdentity3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{ return _f(input); }


template <typename T, template <typename U> class DESCRIPTOR>
ComposedSuperLatticeF3D<T,DESCRIPTOR>::ComposedSuperLatticeF3D
  (SuperLatticeF3D<T,DESCRIPTOR>& f0, SuperLatticeF3D<T,DESCRIPTOR>& f1,
  SuperLatticeF3D<T,DESCRIPTOR>& f2)
  : SuperLatticeF3D<T,DESCRIPTOR>(f0.getSuperLattice3D(), 3), _f0(f0), _f1(f1), _f2(f2)
{ this->_name = "composedSuperLatticeF3D"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> ComposedSuperLatticeF3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  std::vector<T> tmp(3,0);
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff0(_f0), ff1(_f1), ff2(_f2);
  tmp[0] = _f0(input)[0];
  tmp[1] = _f1(input)[0];
  tmp[2] = _f2(input)[0];
  return tmp;
}


} // end namespace olb

#endif

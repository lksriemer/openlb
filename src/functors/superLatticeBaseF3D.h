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

#ifndef SUPER_LATTICE_BASE_F_3D_H
#define SUPER_LATTICE_BASE_F_3D_H

#include<vector>

#include "functors/genericF.h"
#include "core/superLattice3D.h"
#include "core/units.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

template<typename T, template<typename U> class Lattice> class SuperLattice3D;


/// represents all functors that operate on a Lattice in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeF3D : public GenericF<T,int> {
protected:
  SuperLattice3D<T,DESCRIPTOR>& _sLattice;
public:
  SuperLatticeF3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, int targetDim);

  SuperLatticeF3D<T,DESCRIPTOR>& operator-(SuperLatticeF3D<T,DESCRIPTOR>& rhs);
  SuperLatticeF3D<T,DESCRIPTOR>& operator+(SuperLatticeF3D<T,DESCRIPTOR>& rhs);
  SuperLatticeF3D<T,DESCRIPTOR>& operator*(SuperLatticeF3D<T,DESCRIPTOR>& rhs);
  SuperLatticeF3D<T,DESCRIPTOR>& operator/(SuperLatticeF3D<T,DESCRIPTOR>& rhs);

  SuperLattice3D<T,DESCRIPTOR>& getSuperLattice3D() { return _sLattice; }
};

/// represents all functors that operate on a Lattice with output in Phys, e.g. physVelocity(), physForce(), physPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysF3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  const LBconverter<T>& _converter;
public:
  SuperLatticePhysF3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                      const LBconverter<T>& converter, int targetDim);
  LBconverter<T> const& getConverter() const;
};

/// identity functor for memory management
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeIdentity3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  SuperLatticeF3D<T,DESCRIPTOR>& _f;
public:
  SuperLatticeIdentity3D<T,DESCRIPTOR>(SuperLatticeF3D<T,DESCRIPTOR>& f);
  ~SuperLatticeIdentity3D<T,DESCRIPTOR>();
  // access operator should not delete f, since f still has the identity as child
  std::vector<T> operator()(std::vector<int> input);
};


template <typename T, template <typename U> class DESCRIPTOR>
class ComposedSuperLatticeF3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeF3D<T,DESCRIPTOR>& _f0;
  SuperLatticeF3D<T,DESCRIPTOR>& _f1;
  SuperLatticeF3D<T,DESCRIPTOR>& _f2;
public:
  ComposedSuperLatticeF3D(SuperLatticeF3D<T,DESCRIPTOR>& f0,
                          SuperLatticeF3D<T,DESCRIPTOR>& f1,
                          SuperLatticeF3D<T,DESCRIPTOR>& f2);
  std::vector<T> operator()(std::vector<int> x);
};

} // end namespace olb

#endif

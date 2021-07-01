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

#ifndef INTERPOLATION_F_3D_H
#define INTERPOLATION_F_3D_H


#include "functors/analyticalF.h"
#include "functors/blockBaseF3D.h"
#include "functors/superBaseF3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "geometry/blockGeometry3D.h"



namespace olb {

template< typename T, template <typename U> class DESCRIPTOR> class BlockLattice3D;
template< typename T, template <typename U> class DESCRIPTOR> class SuperLattice3D;

/// a class used to convert a block functior to analytical functor
template <typename T>
class AnalyticalFfromBlockF3D : public AnalyticalF3D<T,T> {
protected:
  BlockF3D<T>& _f;
  Cuboid3D<T>& _cuboid;
public:
  AnalyticalFfromBlockF3D(BlockF3D<T>& f, Cuboid3D<T>& cuboid);
  bool operator() (T output[], const T physC[]);
};

/// a class used to convert super lattice functions to analytical functions
template <typename T>
class AnalyticalFfromSuperF3D : public AnalyticalF3D<T,T> {
protected:
  SuperF3D<T>&                              _f;
  CuboidGeometry3D<T>&                      _cuboidGeometry;
  bool                                      _communicateToAll;
  int                                       _overlap;
  std::vector<AnalyticalFfromBlockF3D<T>* > _analyticalFfromBlockF;
public:
  AnalyticalFfromSuperF3D(SuperF3D<T>& f, bool communicateToAll=false, int overlap=-1);
  virtual ~AnalyticalFfromSuperF3D();
  bool operator() (T output[], const T physC[]);
};


/**
 *  class used to convert analytical functions to lattice functions
 *  input functions are interpreted as SI->SI units, the resulting lattice
 *  function will map lattice->lattice units
 */
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeFfromAnalyticalF3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  AnalyticalF3D<T,T>& _f;
public:
  SuperLatticeFfromAnalyticalF3D(AnalyticalF3D<T,T>& f,
                                 SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};


//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeFfromAnalyticalF3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  AnalyticalF3D<T,T>&  _f;
  BlockGeometry3D<T>&  _superGeometry;
  CuboidGeometry3D<T>& _cuboidGeometry;
public:
  BlockLatticeFfromAnalyticalF3D(AnalyticalF3D<T,T>& f,
                                 BlockLattice3D<T,DESCRIPTOR>& sLattice,
                                 BlockGeometry3D<T>& superGeometry,
                                 CuboidGeometry3D<T>& cuboidGeometry);
  bool operator() (T output[], const int input[]);
};



} // end namespace olb

#endif

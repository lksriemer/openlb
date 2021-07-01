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


#include<vector>

#include "functors/analyticalF.h"
#include "functors/blockLatticeBaseF3D.h"
#include "core/superLattice3D.h"
#include "geometry/superGeometry3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "core/blockLattice3D.h"


namespace olb {


/// a class used to convert lattice functions to analytical functions
template <typename T, template <typename U> class DESCRIPTOR>
class AnalyticalFfromSuperLatticeF3D : public AnalyticalF3D<T,T> {
protected:
  SuperLatticeF3D<T,DESCRIPTOR>&  _f;
  CuboidGeometry3D<T>&            _cuboidGeometry;
  bool                            _communicateToAll;
  int                             _overlap;
public:
  AnalyticalFfromSuperLatticeF3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
                                 bool communicateToAll=false, int overlap=2);
  std::vector<T> operator() (std::vector<T> physC);
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
  std::vector<T> operator() (std::vector<int> input);
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
    std::vector<T> operator() (std::vector<int> input);
};



} // end namespace olb

#endif

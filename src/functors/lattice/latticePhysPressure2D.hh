/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_PHYS_PRESSURE_2D_HH
#define LATTICE_PHYS_PRESSURE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePhysPressure2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticePhysPressure2D<T,DESCRIPTOR>::SuperLatticePhysPressure2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physPressure";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysPressure2D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysPressure2D<T,DESCRIPTOR>::BlockLatticePhysPressure2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physPressure";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysPressure2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  // lattice pressure = c_s^2 ( rho -1 )
  T latticePressure = ( this->_blockLattice.get( input[0], input[1] ).computeRho() - 1.0 ) / descriptors::invCs2<T,DESCRIPTOR>();
  output[0] = this->_converter.getPhysPressure(latticePressure);

  return true;
}


template<typename T,typename DESCRIPTOR>
SuperLatticePhysIncPressure2D<T,DESCRIPTOR>::SuperLatticePhysIncPressure2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physIncPressure";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysIncPressure2D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC), this->_converter));
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticeIncPressure2D<T,DESCRIPTOR>::SuperLatticeIncPressure2D(
  SuperLattice<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "incPressure";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeIncPressure2D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC)));
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysIncPressure2D<T,DESCRIPTOR>::BlockLatticePhysIncPressure2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physIncPressure";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysIncPressure2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T latticePressure = this->_blockLattice.get( input[0], input[1] ).computeRho();
  output[0] = this->_converter.getPhysPressure(latticePressure);

  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticeIncPressure2D<T,DESCRIPTOR>::BlockLatticeIncPressure2D
(BlockLattice<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{
  this->getName() = "incPressure";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticeIncPressure2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T latticePressure = this->_blockLattice.get( input[0], input[1] ).computeRho();
  output[0] = latticePressure;

  return true;
}

}
#endif

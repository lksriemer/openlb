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

#ifndef LATTICE_PHYS_STRAIN_RATE_3D_HH
#define LATTICE_PHYS_STRAIN_RATE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysStrainRate3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticePhysStrainRate3D<T, DESCRIPTOR>::SuperLatticePhysStrainRate3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 9)
{
  this->getName() = "physStrainRate";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysStrainRate3D<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysStrainRate3D<T, DESCRIPTOR>::BlockLatticePhysStrainRate3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 9)
{
  this->getName() = "strainRate";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysStrainRate3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(input).computeAllMomenta(rho, uTemp, pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();

  output[0] = -pi[0] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[1] = -pi[1] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[2] = -pi[2] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[3] = -pi[1] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[4] = -pi[3] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[5] = -pi[4] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[6] = -pi[2] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[7] = -pi[4] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[8] = -pi[5] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;

  return true;
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysStressTensor3D<T, DESCRIPTOR>::SuperLatticePhysStressTensor3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 9)
{
  this->getName() = "physStressTensor";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysStressTensor3D<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysStressTensor3D<T, DESCRIPTOR>::BlockLatticePhysStressTensor3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 9),
    _pressureF(blockLattice, converter)
{
  this->getName() = "stressTensor";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysStressTensor3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(input).computeAllMomenta(rho, uTemp, pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();

  output[0] = -pi[0] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[1] = -pi[1] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[2] = -pi[2] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[3] = -pi[1] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[4] = -pi[3] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[5] = -pi[4] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[6] = -pi[2] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[7] = -pi[4] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[8] = -pi[5] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;

  for (int i=0; i < 9; ++i) {
    output[i] *= 2*this->_converter.getPhysDensity()*this->_converter.getPhysViscosity();
  }

  T pressure{};
  _pressureF(&pressure, input);

  output[0] -= pressure;
  output[4] -= pressure;
  output[8] -= pressure;

  return true;
}

}
#endif

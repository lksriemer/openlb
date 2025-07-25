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

#ifndef LATTICE_PHYS_PRESSURE_2D_H
#define LATTICE_PHYS_PRESSURE_2D_H

#include <vector>

#include "superBaseF2D.h"
#include "core/superLattice2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "utilities/functorPtr.h"
#include "blockBaseF2D.h"
#include "geometry/blockGeometry.h"
#include "indicator/blockIndicatorF2D.h"
#include "dynamics/porousBGKdynamics.h"

namespace olb {

/// functor to get pointwise phys pressure from rho on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysPressure2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysPressure2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
};

template <typename T, typename DESCRIPTOR>
SuperLatticePhysPressure2D(SuperLattice<T,DESCRIPTOR>&,
                           const UnitConverter<T,DESCRIPTOR>&)
  -> SuperLatticePhysPressure2D<T,DESCRIPTOR>;

/// BlockLatticePhysPressure2D returns pointwise phys pressure from rho on local lattices.
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysPressure2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysPressure2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys pressure from incompressible model on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysIncPressure2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysIncPressure2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                                const UnitConverter<T,DESCRIPTOR>& converter);
};

template <typename T, typename DESCRIPTOR>
SuperLatticePhysIncPressure2D(SuperLattice<T,DESCRIPTOR>&,
                              const UnitConverter<T,DESCRIPTOR>&)
  -> SuperLatticePhysIncPressure2D<T,DESCRIPTOR>;

/// functor to get pointwise pressure from incompressible model on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticeIncPressure2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeIncPressure2D(SuperLattice<T,DESCRIPTOR>& sLattice);
};

template <typename T, typename DESCRIPTOR>
SuperLatticeIncPressure2D(SuperLattice<T,DESCRIPTOR>&)
  -> SuperLatticeIncPressure2D<T,DESCRIPTOR>;

  
/// BlockLatticePhysPressure2D returns pointwise phys pressure from incompressible model on local lattices.
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysIncPressure2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysIncPressure2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// BlockLatticePhysPressure2D returns pointwise pressure from incompressible model on local lattices.
template <typename T, typename DESCRIPTOR>
class BlockLatticeIncPressure2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
public:
  BlockLatticeIncPressure2D(BlockLattice<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

}
#endif

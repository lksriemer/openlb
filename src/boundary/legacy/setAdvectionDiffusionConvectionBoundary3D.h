/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

//This file contains the AdvectionDiffusionConvection Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_ADVECTION_DIFFUSION_CONVECTION_BOUNDARY_3D_H
#define SET_ADVECTION_DIFFUSION_CONVECTION_BOUNDARY_3D_H

#include <vector>
#include "utilities/functorPtr.h"
#include "geometry/superGeometry.h"
#include "core/superLattice.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/dynamics.h"
#include "descriptor/definition/rtlbm.h"
#include "io/ostreamManager.h"
#include "advectionDiffusionBoundaryPostProcessor3D.h"
#include "boundary/dynamics/advectionDiffusionBoundaries.h"
#include "setBoundary3D.h"


namespace olb {

///Initialising the AdvectionDiffusionConvectionBoundary function on the superLattice domain
///This is an Advection Diffusion Boundary therefore mostly--> MixinDynamics = AdvectionDiffusionRLBdynamics<T,DESCRIPTOR>
template<typename T, typename DESCRIPTOR>
void setAdvectionDiffusionConvectionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material);

///Initialising the AdvectionDiffusionConvectionBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setAdvectionDiffusionConvectionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,FunctorPtr<SuperIndicatorF3D<T>>&& indicator);


///Add AdvectionDiffusionConvection boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setAdvectionDiffusionConvectionBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, bool includeOuterCells=false);


}//namespace olb
#endif


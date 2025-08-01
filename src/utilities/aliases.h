/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
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

#ifndef ALIASES_H
#define ALIASES_H

#include <type_traits>

namespace olb {

// *INDENT-OFF*

template <typename T> class Communicator2D;
template <typename T> class Communicator3D;

template <typename T, unsigned D>
using Communicator = std::conditional_t<
  D == 2,
  Communicator2D<T>,
  Communicator3D<T>
>;

template <typename T, typename DESCRIPTOR> class PostProcessor2D;
template <typename T, typename DESCRIPTOR> class PostProcessor3D;

template <typename T, typename DESCRIPTOR>
using PostProcessor = std::conditional_t<
  DESCRIPTOR::d == 2,
  PostProcessor2D<T,DESCRIPTOR>,
  PostProcessor3D<T,DESCRIPTOR>
>;

template <typename T, typename DESCRIPTOR> class PostProcessorGenerator2D;
template <typename T, typename DESCRIPTOR> class PostProcessorGenerator3D;

template <typename T, typename DESCRIPTOR>
using PostProcessorGenerator = std::conditional_t<
  DESCRIPTOR::d == 2,
  PostProcessorGenerator2D<T,DESCRIPTOR>,
  PostProcessorGenerator3D<T,DESCRIPTOR>
>;


template <typename T> class SuperGeometryStatistics2D;
template <typename T> class SuperGeometryStatistics3D;

template <typename T, unsigned DIM>
using SuperGeometryStatistics = std::conditional_t<
                     DIM == 2,
                     SuperGeometryStatistics2D<T>,
                     SuperGeometryStatistics3D<T>
                     >;

template<typename T, typename OUT_T, typename W> class SuperVTMwriter2D;
template<typename T, typename OUT_T, typename W> class SuperVTMwriter3D;

template <typename T, unsigned DIM, typename OUT_T=float, typename W=T>
using SuperVTMwriter = std::conditional_t<
                     DIM == 2,
                     SuperVTMwriter2D<T,OUT_T,W>,
                     SuperVTMwriter3D<T,OUT_T,W>
                     >;

template <typename T> class SuperGeometryF2D;
template <typename T> class SuperGeometryF3D;

template <typename T, unsigned D>
using SuperGeometryF = std::conditional_t<
                     D == 2,
                     SuperGeometryF2D<T>,
                     SuperGeometryF3D<T>
                     >;

template <typename T, typename DESCRIPTOR> class SuperLatticeCuboid2D;
template <typename T, typename DESCRIPTOR> class SuperLatticeCuboid3D;

template <typename T, typename DESCRIPTOR>
using SuperLatticeCuboid = std::conditional_t<
                     DESCRIPTOR::d == 2,
                     SuperLatticeCuboid2D<T,DESCRIPTOR>,
                     SuperLatticeCuboid3D<T,DESCRIPTOR>
                     >;

template <typename T, typename DESCRIPTOR> class SuperLatticeRank2D;
template <typename T, typename DESCRIPTOR> class SuperLatticeRank3D;

template <typename T, typename DESCRIPTOR>
using SuperLatticeRank = std::conditional_t<
                     DESCRIPTOR::d == 2,
                     SuperLatticeRank2D<T,DESCRIPTOR>,
                     SuperLatticeRank3D<T,DESCRIPTOR>
                     >;

template <typename T, typename DESCRIPTOR> class SuperLatticeF2D;
template <typename T, typename DESCRIPTOR> class SuperLatticeF3D;

template <typename T, typename DESCRIPTOR>
using SuperLatticeF = std::conditional_t<
                     DESCRIPTOR::d == 2,
                     SuperLatticeF2D<T,DESCRIPTOR>,
                     SuperLatticeF3D<T,DESCRIPTOR>
                     >;

template <typename T, typename W> class SuperIntegral2D;
template <typename T, typename W> class SuperIntegral3D;
template <unsigned DIM, typename T, typename W = T>
using SuperIntegral = std::conditional_t<
                     DIM == 2,
                     SuperIntegral2D<T,W>,
                     SuperIntegral3D<T,W>
                     >;

template <typename T, typename DESCRIPTOR, typename FIELD> class SuperLatticeField2D;
template <typename T, typename DESCRIPTOR, typename FIELD> class SuperLatticeField3D;

template <typename T, typename DESCRIPTOR, typename FIELD>
using SuperLatticeField = std::conditional_t<
                          DESCRIPTOR::d == 2,
                          SuperLatticeField2D<T,DESCRIPTOR,FIELD>,
                          SuperLatticeField3D<T,DESCRIPTOR,FIELD>
                          >;

template <typename T, typename DESCRIPTOR> class BlockLatticeF2D;
template <typename T, typename DESCRIPTOR> class BlockLatticeF3D;

template <typename T, typename DESCRIPTOR>
using BlockLatticeF = std::conditional_t<
                     DESCRIPTOR::d == 2,
                     BlockLatticeF2D<T,DESCRIPTOR>,
                     BlockLatticeF3D<T,DESCRIPTOR>
                     >;

template <typename T, typename DESCRIPTOR> class LatticeCouplingGenerator2D;
template <typename T, typename DESCRIPTOR> class LatticeCouplingGenerator3D;

template <typename T, typename DESCRIPTOR>
using LatticeCouplingGenerator = std::conditional_t<
  DESCRIPTOR::d == 2,
  LatticeCouplingGenerator2D<T,DESCRIPTOR>,
  LatticeCouplingGenerator3D<T,DESCRIPTOR>
>;

template <typename T> class BlockGeometryStatistics2D;
template <typename T> class BlockGeometryStatistics3D;

template <typename T, unsigned D>
using BlockGeometryStatistics = std::conditional_t<
  D == 2,
  BlockGeometryStatistics2D<T>,
  BlockGeometryStatistics3D<T>
>;

template <typename T> class BlockF2D;
template <typename T> class BlockF3D;

template <typename T, unsigned D>
using BlockF = std::conditional_t<
  D == 2,
  BlockF2D<T>,
  BlockF3D<T>
>;

template <typename T, typename U> class SuperF2D;
template <typename T, typename U> class SuperF3D;

template <unsigned D, typename T, typename U=T>
using SuperF = std::conditional_t<
  D == 2,
  SuperF2D<T,U>,
  SuperF3D<T,U>
>;

template <typename T> class SuperIndicatorF2D;
template <typename T> class SuperIndicatorF3D;

template <typename T, unsigned D>
using SuperIndicatorF = std::conditional_t<
  D == 2,
  SuperIndicatorF2D<T>,
  SuperIndicatorF3D<T>
>;

template <typename T> class BlockIndicatorF2D;
template <typename T> class BlockIndicatorF3D;

template <typename T, unsigned D>
using BlockIndicatorF = std::conditional_t<
  D == 2,
  BlockIndicatorF2D<T>,
  BlockIndicatorF3D<T>
>;

template <typename T> class BlockIndicatorMaterial2D;
template <typename T> class BlockIndicatorMaterial3D;

template <typename T, unsigned D>
using BlockIndicatorMaterial = std::conditional_t<
  D == 2,
  BlockIndicatorMaterial2D<T>,
  BlockIndicatorMaterial3D<T>
>;

template <typename T> class BlockIndicatorBoundaryNeighbor2D;
template <typename T> class BlockIndicatorBoundaryNeighbor3D;

template <typename T, unsigned D>
using BlockIndicatorBoundaryNeighbor = std::conditional_t<
  D == 2,
  BlockIndicatorBoundaryNeighbor2D<T>,
  BlockIndicatorBoundaryNeighbor3D<T>
>;

template <typename T> class BlockIndicatorFfromIndicatorF2D;
template <typename T> class BlockIndicatorFfromIndicatorF3D;

template <typename T, unsigned D>
using BlockIndicatorFfromIndicatorF = std::conditional_t<
  D == 2,
  BlockIndicatorFfromIndicatorF2D<T>,
  BlockIndicatorFfromIndicatorF3D<T>
>;

template <typename T> class IndicatorF2D;
template <typename T> class IndicatorF3D;

template <typename T, unsigned D>
using IndicatorF = std::conditional_t<
  D == 2,
  IndicatorF2D<T>,
  IndicatorF3D<T>
>;

template <typename T> class IndicatorCuboid2D;
template <typename T> class IndicatorCuboid3D;

template <typename T, unsigned D>
using IndicatorCuboid = std::conditional_t<
  D == 2,
  IndicatorCuboid2D<T>,
  IndicatorCuboid3D<T>
>;

template <typename T, typename S, bool PARTICLE> class SmoothIndicatorF2D;
template <typename T, typename S, bool PARTICLE> class SmoothIndicatorF3D;

template <typename T, typename S, unsigned D, bool PARTICLE=false>
using SmoothIndicatorF = std::conditional_t<
  D == 2,
  SmoothIndicatorF2D<T,T,PARTICLE>,
  SmoothIndicatorF3D<T,T,PARTICLE>
>;

template <typename T> class SuperIndicatorFfromIndicatorF2D;
template <typename T> class SuperIndicatorFfromIndicatorF3D;

template <typename T, unsigned D>
using SuperIndicatorFfromIndicatorF = std::conditional_t<
  D == 2,
  SuperIndicatorFfromIndicatorF2D<T>,
  SuperIndicatorFfromIndicatorF3D<T>
>;

template <typename T> class SuperIndicatorMaterial2D;
template <typename T> class SuperIndicatorMaterial3D;

template <typename T, unsigned D>
using SuperIndicatorMaterial = std::conditional_t<
  D == 2,
  SuperIndicatorMaterial2D<T>,
  SuperIndicatorMaterial3D<T>
>;

template <typename T> class SuperIndicatorBoundaryNeighbor2D;
template <typename T> class SuperIndicatorBoundaryNeighbor3D;

template <typename T, unsigned D>
using SuperIndicatorBoundaryNeighbor = std::conditional_t<
  D == 2,
  SuperIndicatorBoundaryNeighbor2D<T>,
  SuperIndicatorBoundaryNeighbor3D<T>
>;

template <typename T, typename DESCRIPTOR, typename FIELD> class SuperField2D;
template <typename T, typename DESCRIPTOR, typename FIELD> class SuperField3D;

template <typename T, typename DESCRIPTOR, typename FIELD>
using SuperField = std::conditional_t<
  DESCRIPTOR::d == 2,
  SuperField2D<T,DESCRIPTOR,FIELD>,
  SuperField3D<T,DESCRIPTOR,FIELD>
>;


template <typename T, typename DESCRIPTOR> class SuperLatticePhysF2D;
template <typename T, typename DESCRIPTOR> class SuperLatticePhysF3D;

template <typename T, typename DESCRIPTOR>
using SuperLatticePhysF = std::conditional_t<
  DESCRIPTOR::d == 2,
  SuperLatticePhysF2D<T,DESCRIPTOR>,
  SuperLatticePhysF3D<T,DESCRIPTOR>
>;

template <typename T, typename DESCRIPTOR> class BlockLatticePhysF2D;
template <typename T, typename DESCRIPTOR> class BlockLatticePhysF3D;

template <typename T, typename DESCRIPTOR>
using BlockLatticePhysF = std::conditional_t<
  DESCRIPTOR::d == 2,
  BlockLatticePhysF2D<T,DESCRIPTOR>,
  BlockLatticePhysF3D<T,DESCRIPTOR>
>;

template <typename T, typename DESCRIPTOR> class SuperLatticeInterpPhysVelocity2D;
template <typename T, typename DESCRIPTOR> class SuperLatticeInterpPhysVelocity3D;

template <typename T, typename DESCRIPTOR>
using SuperLatticeInterpPhysVelocity = std::conditional_t<
  DESCRIPTOR::d == 2,
  SuperLatticeInterpPhysVelocity2D<T,DESCRIPTOR>,
  SuperLatticeInterpPhysVelocity3D<T,DESCRIPTOR>
>;


template <typename T, typename DESCRIPTOR> class BlockLatticeInterpPhysVelocity2D;
template <typename T, typename DESCRIPTOR> class BlockLatticeInterpPhysVelocity3D;

template <typename T, typename DESCRIPTOR>
using BlockLatticeInterpPhysVelocity = std::conditional_t<
  DESCRIPTOR::d == 2,
  BlockLatticeInterpPhysVelocity2D<T,DESCRIPTOR>,
  BlockLatticeInterpPhysVelocity3D<T,DESCRIPTOR>
>;

template <typename T, typename DESCRIPTOR> class SuperLatticeFpop2D;
template <typename T, typename DESCRIPTOR> class SuperLatticeFpop3D;

template <typename T, typename DESCRIPTOR>
using SuperLatticeFpop = std::conditional_t<
  DESCRIPTOR::d == 2,
  SuperLatticeFpop2D<T,DESCRIPTOR>,
  SuperLatticeFpop3D<T,DESCRIPTOR>
>;


template <typename T> class BlockIndicatorFfromCallableF2D;
template <typename T> class BlockIndicatorFfromCallableF3D;

template <typename T, unsigned D>
using BlockIndicatorFfromCallableF = std::conditional_t<
  D == 2,
  BlockIndicatorFfromCallableF2D<T>,
  BlockIndicatorFfromCallableF3D<T>
>;


template <typename T, typename BaseType> class BlockVTIreader2D;
template <typename T, typename BaseType> class BlockVTIreader3D;

template <typename T, typename BaseType, unsigned D >
using BlockVTIreader = std::conditional_t<
    D == 2,
    BlockVTIreader2D<T, BaseType>,
    BlockVTIreader3D<T, BaseType>
>;


template <typename T, typename BaseType> class BlockDataF2D;
template <typename T, typename BaseType> class BlockDataF3D;

template <typename T, typename BaseType, unsigned D >
using BlockDataF = std::conditional_t<
  D == 2,
  BlockDataF2D<T, BaseType>,
  BlockDataF3D<T, BaseType>
>;


template <typename T, typename W> class SpecialAnalyticalFfromBlockF2D;
template <typename T, typename W> class SpecialAnalyticalFfromBlockF3D;

template <typename T, typename W, unsigned D >
using SpecialAnalyticalFfromBlockF = std::conditional_t<
  D == 2,
  SpecialAnalyticalFfromBlockF2D<T, W>,
  SpecialAnalyticalFfromBlockF3D<T, W>
>;


template <typename T, typename DESCRIPTOR, typename FIELD> struct SuperLatticePhysField2D;
template <typename T, typename DESCRIPTOR, typename FIELD> struct SuperLatticePhysField3D;

template <typename T, typename DESCRIPTOR, typename FIELD>
using SuperLatticePhysField = std::conditional_t<
  DESCRIPTOR::d == 2,
  SuperLatticePhysField2D<T, DESCRIPTOR, FIELD>,
  SuperLatticePhysField3D< T, DESCRIPTOR, FIELD>
>;

// *INDENT-ON*

}

#endif

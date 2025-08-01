/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Eduardo Czelusniak
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

//This file contains all inlet and outlet boundary conditions of phase field models
//This boundary only contains free floating functions

#ifndef PHASE_FIELD_INLET_OUTLET_2D_H
#define PHASE_FIELD_INLET_OUTLET_2D_H

#include "setBoundary.h"
#include "phaseFieldInletOutlet.h"
#include "zouHeDynamics.h"
#include "dynamics/freeEnergyDynamics.h"
#include "phaseFieldBoundaryDynamics.h"

namespace olb {

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct IncompressibleZouHeVelocity<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      IncZouHeDynamics,MixinDynamics,momenta::IncDirichletVelocityBoundaryTuple
    >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    return meta::id<typename MixinDynamics::template exchange_momenta<momenta::FixedVelocityBoundaryTuple>>{};

  case DiscreteNormalType::InternalCorner:
    return boundaryhelper::PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple2D
    >::construct(n);

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::ExternalCorner:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,OuterVelocityCornerProcessor2D>(n);

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct IncompressibleZouHePressure<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      IncZouHeDynamics,MixinDynamics,momenta::IncDirichletPressureBoundaryTuple
    >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    return meta::id<typename MixinDynamics::template exchange_momenta<momenta::FixedVelocityBoundaryTuple>>{};

  case DiscreteNormalType::InternalCorner:
    return boundaryhelper::PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple2D
    >::construct(n);

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {

  case DiscreteNormalType::ExternalCorner:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,OuterVelocityCornerProcessor2D>(n);

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct IncompressibleConvective<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
          ConvectiveOutletDynamics,MixinDynamics,
          momenta::IncBulkTuple<momenta::ForcedMomentum<momenta::IncompressibleBulkMomentum>>
        >::construct(n);

  /*case DiscreteNormalType::ExternalCorner:
    return meta::id<typename MixinDynamics::template exchange_momenta<momenta::FixedVelocityBoundaryTuple>>{};

  case DiscreteNormalType::InternalCorner:
    return boundaryhelper::PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple2D
    >::construct(n);*/

  default:
    return std::nullopt;
  }
  //return DynamicsPromise(meta::id<MixinDynamics>{});
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  return meta::id<ConvectivePostProcessor2D<T,DESCRIPTOR>>();
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct PhaseFieldInlet<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
          PhaseFieldInletDynamics,MixinDynamics,momenta::EquilibriumBoundaryTuple
        >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    return meta::id<typename MixinDynamics::template exchange_momenta<momenta::FixedVelocityBoundaryTuple>>{};

  case DiscreteNormalType::InternalCorner:
    return boundaryhelper::PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple2D
    >::construct(n);

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::ExternalCorner:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,OuterVelocityCornerProcessor2D>(n);

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct PhaseFieldConvective<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
          PhaseFieldConvectiveOutletDynamics,MixinDynamics,momenta::ExternalVelocityTuple
        >::construct(n);

  /*case DiscreteNormalType::ExternalCorner:
    return meta::id<typename MixinDynamics::template exchange_momenta<momenta::FixedVelocityBoundaryTuple>>{};

  case DiscreteNormalType::InternalCorner:
    return boundaryhelper::PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple2D
    >::construct(n);*/

  default:
    return std::nullopt;
  }
  //return DynamicsPromise(meta::id<MixinDynamics>{});
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FlatConvectivePhaseFieldPostProcessor2D>(n);
  //return std::nullopt;
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct FreeEnergyVelocity<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedVelocityBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletMomentum2D>(n);

  case DiscreteNormalType::ExternalCorner:
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct FreeEnergyPressure<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletMomentum2D>(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 2)
struct FreeEnergyPressureConvective<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
      CombinedRLBdynamics,MixinDynamics,momenta::RegularizedPressureBoundaryTuple
      >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyOutletMomentum2D>(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyOutletMomentum2D>(n);

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyOutletMomentum2D>(n);

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
requires (DESCRIPTOR::d == 2)
struct FreeEnergyOrderParameter<T,DESCRIPTOR> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationDynamics<T,DESCRIPTOR,FreeEnergyInletOutletDynamics>
    ::construct(n);

  case DiscreteNormalType::ExternalCorner:
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyInletOrderParameter2D>(n);

  case DiscreteNormalType::ExternalCorner:
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

};

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
requires (DESCRIPTOR::d == 2)
struct FreeEnergyOrderParameterConvective<T,DESCRIPTOR> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::DirectionOrientationDynamics<T,DESCRIPTOR,FreeEnergyInletOutletDynamics>
    ::construct(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,FreeEnergyOutletOrderParameter2D>(n);

  case DiscreteNormalType::ExternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  case DiscreteNormalType::InternalCorner:
    throw std::runtime_error("No valid discrete normal found. This BC is not suited for curved walls.");
    return std::nullopt;

  default:
    return std::nullopt;
  }
}

};

}

}

#endif

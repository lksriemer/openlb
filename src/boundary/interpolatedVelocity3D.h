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

#ifndef INTERPOLATED_VELOCITY_3D_H
#define INTERPOLATED_VELOCITY_3D_H

#include "setBoundary.h"
#include "interpolatedVelocity.h"

namespace olb {

namespace boundary {

template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR, typename MixinDynamics>
requires (DESCRIPTOR::d == 3)
struct InterpolatedVelocity<T,DESCRIPTOR,MixinDynamics> {

using value_t = T;
using descriptor_t = DESCRIPTOR;

CellDistance getNeighborhoodRadius() {
  return 1;
}

std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                         DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::MixinDynamicsExchangeDirectionOrientationMomenta<T,DESCRIPTOR,
          MixinDynamics,momenta::BasicDirichletVelocityBoundaryTuple
          >::construct(n);

  case DiscreteNormalType::ExternalCorner:
    return meta::id<typename MixinDynamics::template exchange_momenta<momenta::FixedVelocityBoundaryTuple>>{};

  case DiscreteNormalType::InternalCorner:
    return boundaryhelper::PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
          CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple3D
          >::construct(n);

  case DiscreteNormalType::ExternalEdge:
    return meta::id<typename MixinDynamics::template exchange_momenta<momenta::FixedVelocityBoundaryTuple>>{};

  case DiscreteNormalType::InternalEdge:
    return boundaryhelper::PlainMixinDynamicsForNormalSpecialMomenta<T,DESCRIPTOR,
          CombinedRLBdynamics,MixinDynamics,momenta::InnerEdgeVelocityTuple3D
          >::construct(n);

  default:
    return std::nullopt;
  }
}

std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                   DiscreteNormal<DESCRIPTOR> n) {
  switch (type) {
  case DiscreteNormalType::Flat:
    return boundaryhelper::promisePostProcessorForDirectionOrientation<
          T,DESCRIPTOR,PlaneFdBoundaryProcessor3D>(n);

  case DiscreteNormalType::ExternalCorner:
    return boundaryhelper::promisePostProcessorForNormal<T,DESCRIPTOR,
        OuterVelocityCornerProcessor3D>(n);

  case DiscreteNormalType::ExternalEdge:
    return boundaryhelper::promisePostProcessorForNormalSpecial<T,DESCRIPTOR,
        OuterVelocityEdgeProcessor3D>(n);

  default:
    return std::nullopt;
  }
}

};

}

}

#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-24 Adrian Kummerlaender, Shota Ito
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

/*  ========================================================
 *  ==  WARNING: This is an automatically generated file, ==
 *  ==                  do not modify.                    ==
 *  ========================================================
 */

#pragma once


namespace olb {

namespace dynamics {

template <typename T, typename... FIELDS>
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedTemperatureMomentum<1, -1>, momenta::NoStress, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x5 = parameters.template get<descriptors::OMEGA>();
auto x6 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x7 = x6 + V{-1};
auto x8 = V{0.0555555555555556}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x9 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x10 = x9 + V{-1};
auto x11 = x6 + V{1};
auto x12 = x9 + V{1};
auto x13 = x5 + V{-1};
auto x14 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x15 = x13*(V{0.5}*cell[1] - V{0.5}*cell[3] + x11*x14 + x14*x7);
auto x16 = V{0.166666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0);
auto x17 = V{1.38777878078145e-17}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x12 + V{1}*cell[2];
auto x18 = x10*x16 + V{0.166666666666667};
cell[0] = V{0.0555555555555556}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x11 + V{0.0555555555555556}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x12 + V{0.111111111111111}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0) - x10*x8 - x7*x8 + V{-0.333333333333333};
cell[1] = -x15 - x16*x7 + V{-0.166666666666667};
cell[2] = x13*(-V{0.166666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x10 - x17 + V{-0.166666666666667}) - x18;
cell[3] = x11*x16 + x15 + V{-0.166666666666667};
cell[4] = V{0.166666666666667}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0)*x12 - x13*(-x17 - x18) + V{-0.166666666666667};
return { V{1}*cell.template getFieldComponent<momenta::FixedDensity::RHO>(0), cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1) };
}
};

}

}

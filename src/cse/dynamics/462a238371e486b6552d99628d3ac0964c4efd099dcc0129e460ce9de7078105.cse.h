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
struct CSE<dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::BGK, AdvectionDiffusionExternalVelocityCollision>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x5 = parameters.template get<descriptors::OMEGA>();
auto x6 = x5 + V{-1};
auto x7 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4];
auto x8 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(0);
auto x9 = x7 + V{1};
auto x10 = V{0.166666666666667}*x5;
auto x11 = V{3}*cell.template getFieldComponent<descriptors::VELOCITY>(1);
auto x0 = -cell[0]*x6 + V{0.333333333333333}*x5*x7;
auto x1 = -cell[1]*x6 - x10*(x9*(x8 + V{-1}) + V{1});
auto x2 = -cell[2]*x6 - x10*(x9*(x11 + V{-1}) + V{1});
auto x3 = -cell[3]*x6 + V{0.166666666666667}*x5*(x9*(x8 + V{1}) + V{-1});
auto x4 = -cell[4]*x6 + V{0.166666666666667}*x5*(x9*(x11 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
return { x7 + V{1}, cell.template getFieldComponent<descriptors::VELOCITY>(0)*cell.template getFieldComponent<descriptors::VELOCITY>(0) + cell.template getFieldComponent<descriptors::VELOCITY>(1)*cell.template getFieldComponent<descriptors::VELOCITY>(1) };
}
};

}

}

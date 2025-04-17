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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = (cell[0] + V{2}*cell[1] + cell[2] + cell[6] + V{2}*cell[7] + V{2}*cell[8] + V{1})/(cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x12 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x13 = V{1.5}*x12;
auto x14 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x15 = V{1.5}*x14;
auto x16 = x13 + x15 + V{-1};
auto x17 = V{0.0277777777777778}*x9;
auto x18 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x19 = -x18;
auto x20 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x21 = -x20;
auto x22 = V{1} - x15;
auto x23 = -x13;
auto x24 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x25 = x23 + x24;
auto x26 = x22 + x25;
auto x27 = V{0.111111111111111}*x9;
auto x28 = V{3}*x12 + x22;
auto x29 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x30 = x29*x29;
auto x31 = x23 - x24;
auto x32 = V{3}*x14 + V{1};
auto x0 = -cell[0]*x10 - V{0.444444444444444}*x9*(x11*x16 + V{1});
auto x1 = -(cell[1]*x10 + x17*(-x11*(x19 + x26 + V{4.5}*(x21*x21)) + V{1}));
auto x2 = -cell[2]*x10 - x27*(-x11*(x19 + x28) + V{1});
auto x3 = -cell[3]*x10 - x17*(-x11*(-x16 - x18 - x24 + V{4.5}*x30) + V{1});
auto x4 = -cell[4]*x10 - x27*(-x11*(x31 + x32) + V{1});
auto x5 = -(cell[5]*x10 + x17*(-x11*(x18 + x22 + x31 + V{4.5}*(x20*x20)) + V{1}));
auto x6 = -cell[6]*x10 - x27*(-x11*(x18 + x28) + V{1});
auto x7 = -cell[7]*x10 - x17*(-x11*(x18 + x26 + V{4.5}*x30) + V{1});
auto x8 = -cell[8]*x10 - x27*(-x11*(x25 + x32) + V{1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { V{1}*x11, x12 + x14 };
}
};

}

}

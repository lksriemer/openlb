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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = x11 + V{1};
auto x14 = V{1} / ((x13)*(x13));
auto x15 = V{1.5}*x14;
auto x16 = cell.template getFieldComponent<descriptors::POROSITY>(0)*cell.template getFieldComponent<descriptors::POROSITY>(0);
auto x17 = cell[1] - cell[5];
auto x18 = -cell[4] + cell[8];
auto x19 = -cell[3] + cell[7] + x17 + x18;
auto x20 = x19*x19;
auto x21 = x16*x20;
auto x22 = x15*x21;
auto x23 = cell[2] - cell[6];
auto x24 = cell[3] - cell[7] + x17 + x23;
auto x25 = -x24;
auto x26 = x15*x16*(x25*x25) + V{-1};
auto x27 = x22 + x26;
auto x28 = V{0.0277777777777778}*x9;
auto x29 = V{3}*cell.template getFieldComponent<descriptors::POROSITY>(0)/x13;
auto x30 = x19*x29;
auto x31 = V{4.5}*x14;
auto x32 = V{2}*cell[1] - V{2}*cell[5] + x18 + x23;
auto x33 = x25*x29;
auto x34 = x27 + x33;
auto x35 = V{0.111111111111111}*x9;
auto x36 = x24*x29;
auto x37 = V{3}*x14;
auto x38 = x24*x24;
auto x39 = x16*x38;
auto x40 = V{1} - x22;
auto x41 = x37*x39 + x40;
auto x42 = -V{2}*cell[3] - cell[4] + V{2}*cell[7] + cell[8] - x23;
auto x43 = x16*x31*(x42*x42);
auto x44 = x21*x37;
auto x45 = -x32;
auto x46 = -x36;
auto x47 = -x15*x39 + x30;
auto x0 = -cell[0]*x10 - V{0.444444444444444}*x9*(x12*x27 + V{1});
auto x1 = -(cell[1]*x10 + x28*(x12*(-x16*x31*x32*x32 - x30 + x34) + V{1}));
auto x2 = -cell[2]*x10 + x35*(x12*(x36 + x41) + V{-1});
auto x3 = -cell[3]*x10 - x28*(x12*(x30 + x34 - x43) + V{1});
auto x4 = -cell[4]*x10 - x35*(x12*(x26 + x30 - x44) + V{1});
auto x5 = -(cell[5]*x10 + x28*(x12*(-x16*x31*x45*x45 + x27 + x30 - x33) + V{1}));
auto x6 = -cell[6]*x10 + x35*(x12*(x41 + x46) + V{-1});
auto x7 = -cell[7]*x10 + x28*(x12*(x40 + x43 + x46 + x47) + V{-1});
auto x8 = -cell[8]*x10 + x35*(x12*(x44 + x47 + V{1}) + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x13, V{1}*x14*x16*(x20 + x38) };
}
};

}

}

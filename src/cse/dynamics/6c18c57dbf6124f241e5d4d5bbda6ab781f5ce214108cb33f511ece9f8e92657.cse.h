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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::fsi::HLBM>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x9 = parameters.template get<descriptors::OMEGA>();
auto x10 = x9 + V{-1};
auto x11 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x12 = x11 + V{1};
auto x13 = x11 + V{1};
auto x14 = V{1} / ((x13)*(x13));
auto x15 = V{1.5}*x14;
auto x16 = cell[1] - cell[5];
auto x17 = -cell[4] + cell[8];
auto x18 = -cell[3] + cell[7] + x16 + x17;
auto x19 = x18*x18;
auto x20 = x15*x19;
auto x21 = cell[2] - cell[6];
auto x22 = cell[3] - cell[7] + x16 + x21;
auto x23 = -x22;
auto x24 = x15*(x23*x23) + V{-1};
auto x25 = x20 + x24;
auto x26 = V{1} / (x13);
auto x27 = V{1}*cell[3];
auto x28 = V{1}*cell[1] - V{1}*cell[5];
auto x29 = x26*(-V{1}*cell[4] + V{1}*cell[7] + V{1}*cell[8] - x27 + x28);
auto x30 = cell.template getFieldComponent<descriptors::POROSITY>(0) + V{-1};
auto x31 = -x30;
auto x32 = -x29;
auto x33 = cell.template getFieldComponent<descriptors::VELOCITY>(1) + x32;
auto x34 = x31*x33;
auto x35 = x29 + x34;
auto x36 = V{1}*cell[2] - V{1}*cell[6] - V{1}*cell[7] + x27 + x28;
auto x37 = -x26*x36;
auto x38 = x31*(cell.template getFieldComponent<descriptors::VELOCITY>(0) - x37);
auto x39 = x37 + x38;
auto x40 = V{-1} + V{1.5}*(x35*x35) + V{1.5}*(x39*x39);
auto x41 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778};
auto x42 = V{3}*cell[3];
auto x43 = V{3}*cell[7];
auto x44 = V{3}*cell[1] - V{3}*cell[5];
auto x45 = V{3}*cell[2] - V{3}*cell[6] + x42 - x43 + x44;
auto x46 = x26*x45;
auto x47 = V{4.5}*x14;
auto x48 = V{2}*cell[1] - V{2}*cell[5] + x17 + x21;
auto x49 = -x48;
auto x50 = V{1} - x20;
auto x51 = x26*(-V{3}*cell[4] + V{3}*cell[8] - x42 + x43 + x44);
auto x52 = x22*x22;
auto x53 = x15*x52;
auto x54 = x51 - x53;
auto x55 = x50 + x54;
auto x56 = x46 + x47*(x49*x49) + x55;
auto x57 = x18*x26;
auto x58 = -x57;
auto x59 = x23*x26 + x38;
auto x60 = x31*x33 - x58 - x59;
auto x61 = -x51;
auto x62 = V{3}*x34;
auto x63 = -x26*x45;
auto x64 = V{3}*x38 + x40 + x63;
auto x65 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111};
auto x66 = V{3}*x14;
auto x67 = x50 + x52*x66;
auto x68 = x46 + x67;
auto x69 = V{0.0277777777777778}*x9;
auto x70 = -V{2}*cell[3] - cell[4] + V{2}*cell[7] + cell[8] - x21;
auto x71 = x47*(x70*x70);
auto x72 = x25 + x51 + x63 - x71;
auto x73 = x34 + x57;
auto x74 = x59 + x73;
auto x75 = x51 + x62;
auto x76 = V{0.111111111111111}*x9;
auto x77 = x19*x66;
auto x78 = x24 + x51 - x77;
auto x79 = x46 + V{-1};
auto x80 = x51 + x79;
auto x81 = -(-V{4.5}*x14*x48*x48 + x20 + x53 + x80);
auto x82 = x30*x33;
auto x83 = V{3}*x82;
auto x84 = -x83;
auto x85 = x26*x36;
auto x86 = x30*(cell.template getFieldComponent<descriptors::VELOCITY>(0) + x85);
auto x87 = x22*x26 + x86;
auto x88 = -x57 + x82 - x87;
auto x89 = x85 + x86;
auto x90 = V{1.5}*(x89*x89);
auto x91 = -x32 - x82;
auto x92 = V{1.5}*(x91*x91);
auto x93 = V{3}*x86 + x90 + x92;
auto x94 = -x46;
auto x95 = x67 + x94;
auto x96 = x79 + x93;
auto x97 = x55 + x71 + x94;
auto x98 = x58 + x82;
auto x99 = -x87 - x98;
auto x100 = x54 + x77 + V{1};
auto x101 = -x98;
auto x0 = -cell[0]*x10 + V{0.444444444444444}*x12*x25 - x40*(V{0.444444444444444}*cell[0] + V{0.444444444444444}*cell[1] + V{0.444444444444444}*cell[2] + V{0.444444444444444}*cell[3] + V{0.444444444444444}*cell[4] + V{0.444444444444444}*cell[5] + V{0.444444444444444}*cell[6] + V{0.444444444444444}*cell[7] + V{0.444444444444444}*cell[8] + V{0.444444444444444}) - V{0.444444444444444}*x9*(x12*x25 + V{1});
auto x1 = -(cell[1]*x10 + x41*x56 + x41*(x61 - x62 + x64 - V{4.5}*x60*x60) - V{0.0277777777777778}*x9*(x12*x56 + V{-1}));
auto x2 = -(cell[2]*x10 + x65*x68 + x65*(x64 - V{4.5}*x59*x59) - V{0.111111111111111}*x9*(x12*x68 + V{-1}));
auto x3 = -(cell[3]*x10 - V{0.0277777777777778}*x12*x72 + x41*(x64 + x75 - V{4.5}*x74*x74) + x69*(x12*x72 + V{1}));
auto x4 = -(cell[4]*x10 - V{0.111111111111111}*x12*x78 + x65*(x40 + x75 - V{4.5}*x73*x73) + x76*(x12*x78 + V{1}));
auto x5 = -cell[5]*x10 - x41*x81 - x41*(x80 + x84 + x93 - V{4.5}*x88*x88) + x69*(x12*x81 + V{-1});
auto x6 = -cell[6]*x10 - x65*x95 - x65*(x96 - V{4.5}*x87*x87) + x76*(x12*x95 + V{-1});
auto x7 = -cell[7]*x10 - x41*x97 - x41*(x61 + x83 + x96 - V{4.5}*x99*x99) + x69*(x12*x97 + V{-1});
auto x8 = -cell[8]*x10 - x100*x65 + x65*(x51 + x84 - x90 - x92 + V{1} + V{4.5}*(x101*x101)) + x76*(x100*x12 + V{-1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { x13, V{1}*x14*(x19 + x52) };
}
};

}

}

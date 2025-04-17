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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<2, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<2, 1, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{0.5}*cell[10];
auto x22 = V{0.5}*cell[11];
auto x23 = V{2}*cell[13];
auto x24 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{1});
auto x25 = V{0.25}*cell[0] + V{0.25}*cell[12] + V{0.5}*cell[13] + V{0.25}*cell[3] + V{0.25};
auto x26 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x27 = cell[0] + cell[12] + cell[3] + x23 + V{1};
auto x28 = V{2}*cell[10] + cell[11] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[2] + cell[8] + cell[9] + x27;
auto x29 = cell[10] + V{2}*cell[11] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + cell[1] + V{2}*cell[5] + cell[6] + cell[7] + x27;
auto x30 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x31 = V{1.5}*x30;
auto x32 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x33 = V{1.5}*x32;
auto x34 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x35 = V{1.5}*x34;
auto x36 = x33 + x35 + V{-1};
auto x37 = x31 + x36;
auto x38 = x37*(V{0.166666666666667}*x24*x28 + V{0.166666666666667}*x26*x29);
auto x39 = V{0.25}*x24*x28 + V{0.25}*x26*x29;
auto x40 = V{0.0138888888888889}*x24*x28 + V{0.0138888888888889}*x26*x29;
auto x41 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x42 = V{4.5}*(x41*x41);
auto x43 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x44 = -x31;
auto x45 = V{1} - x33;
auto x46 = x44 + x45;
auto x47 = x43 + x46;
auto x48 = -x35;
auto x49 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x50 = x48 + x49;
auto x51 = x42 + x47 + x50;
auto x52 = x37 + x43;
auto x53 = -x42 + x49 + x52;
auto x54 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x54;
auto x56 = -x55;
auto x57 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x58 = -x57;
auto x59 = x52 + x58;
auto x60 = x59 - V{4.5}*x56*x56;
auto x61 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x62 = V{4.5}*(x61*x61);
auto x63 = x52 + x57 - x62;
auto x64 = -x40*x63;
auto x65 = x48 + x57;
auto x66 = x47 + x62 + x65;
auto x67 = x40*x66;
auto x68 = -x43;
auto x69 = -V{4.5}*x55*x55;
auto x70 = x37 + x57;
auto x71 = x68 + x69 + x70;
auto x72 = -x40*x71;
auto x73 = V{0.0277777777777778}*x24*x28 + V{0.0277777777777778}*x26*x29;
auto x74 = V{3}*x34;
auto x75 = x47 + x74;
auto x76 = x73*x75;
auto x77 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x78 = V{4.5}*(x77*x77);
auto x79 = x46 + x49 + x65 + x78;
auto x80 = V{0.0555555555555556}*x24*x28 + V{0.0555555555555556}*x26*x29;
auto x81 = V{3}*x32;
auto x82 = x44 + x50 + x81 + V{1};
auto x83 = V{1.11022302462516e-16}*x24*x28 + V{1.11022302462516e-16}*x26*x29;
auto x84 = V{8.32667268468867e-17}*x24*x28 + V{8.32667268468867e-17}*x26*x29;
auto x85 = x31 + V{-1};
auto x86 = x33 + x43 - x74 + x85;
auto x87 = x73*x86;
auto x88 = -x49;
auto x89 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x90 = -V{4.5}*x89*x89;
auto x91 = x52 + x88 + x90;
auto x92 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x54;
auto x93 = -V{4.5}*x92*x92;
auto x94 = x70 + x88 + x93;
auto x95 = V{3}*x30;
auto x96 = x36 + x57 - x95;
auto x97 = -x73*x96;
auto x98 = x45 + x65 + x95;
auto x99 = x73*x98;
auto x100 = V{5.55111512312578e-17}*cell[0] + V{5.55111512312578e-17}*cell[12] + V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[3] + V{5.55111512312578e-17};
auto x101 = V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{1.11022302462516e-16}*cell[3] + V{3.33066907387547e-16}*cell[5] - x24*(V{1.11022302462516e-16}*cell[10] + V{5.55111512312578e-17}*cell[11] + V{1.11022302462516e-16}*cell[14] + V{1.11022302462516e-16}*cell[15] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{5.55111512312578e-17}*cell[2] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x100) - x26*(V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[11] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[17] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[1] + V{1.11022302462516e-16}*cell[5] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x100) - x30*(V{5.55111512312578e-17}*x24*x28 + V{5.55111512312578e-17}*x26*x29) - x38 + x51*(V{0.0277777777777778}*x24*x28 + V{0.0277777777777778}*x26*x29) - x53*(V{4.62592926927149e-18}*x24*x28 + V{4.62592926927149e-18}*x26*x29) + x97 + x99 + V{-2.22044604925031e-16};
auto x102 = V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x101 - x32*x83 - x34*x84 + x64 + x67 + x72 + x73*x79 - x73*x91 - x73*x94 + x76 + x80*x82 - x87;
auto x103 = -x92;
auto x104 = x37 + x49;
auto x105 = x104 + x58;
auto x106 = x105 - V{4.5}*x103*x103;
auto x107 = -x89;
auto x108 = x104 + x68;
auto x109 = x108 - V{4.5}*x107*x107;
auto x110 = x49 + x70 - x78;
auto x111 = -x110*x40;
auto x112 = x40*x79;
auto x113 = -x40*x94;
auto x114 = x73*x82;
auto x115 = x35 + x49 - x81 + x85;
auto x116 = x115*x73;
auto x117 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + V{2.22044604925031e-16}*cell[9] + x101 + x111 + x112 + x113 + x114 - x116 - x32*x84 - x34*x83 + x66*x73 - x71*x73 + x75*x80;
auto x118 = x24*(x102 - x40*x60) + x26*(-x106*x40 - x109*x73 + x117);
auto x119 = V{0.0833333333333333}*cell[11];
auto x120 = V{0.0833333333333333}*cell[2];
auto x121 = V{0.0833333333333334}*cell[15];
auto x122 = V{0.0833333333333334}*cell[16];
auto x123 = V{0.0833333333333334}*cell[6];
auto x124 = V{0.0833333333333334}*cell[7];
auto x125 = V{0.0833333333333333}*x24*x28 + V{0.0833333333333333}*x26*x29;
auto x126 = V{0.0416666666666667}*x24*x28 + V{0.0416666666666667}*x26*x29;
auto x127 = x126*x34;
auto x128 = V{3.46944695195361e-18}*cell[0] + V{3.46944695195361e-18}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[3] + V{3.46944695195361e-18};
auto x129 = x24*(V{6.93889390390723e-18}*cell[10] + V{3.46944695195361e-18}*cell[11] + V{6.93889390390723e-18}*cell[14] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{3.46944695195361e-18}*cell[2] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x128);
auto x130 = x26*(V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[11] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[1] + V{6.93889390390723e-18}*cell[5] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x128);
auto x131 = V{0.00115740740740741}*x24*x28 + V{0.00115740740740741}*x26*x29;
auto x132 = V{0.0833333333333333}*cell[12] - V{0.166666666666667}*cell[13] - V{0.0833333333333334}*cell[14] + V{0.0833333333333333}*cell[3] - V{0.0833333333333334}*cell[5] - x126*x30 + x129 + x130 + x131*x51 + x131*x53 + V{0.0555555555555555};
auto x133 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x119 + x120 - x121 - x122 - x123 - x124 + x125*x32 - x127 + x132;
auto x134 = V{0.0833333333333333}*cell[10];
auto x135 = V{0.0833333333333333}*cell[1];
auto x136 = V{0.0833333333333334}*cell[17];
auto x137 = V{0.0833333333333334}*cell[18];
auto x138 = V{0.0833333333333334}*cell[8];
auto x139 = V{0.0833333333333334}*cell[9];
auto x140 = x126*x32;
auto x141 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x125*x34 + x132 + x134 + x135 - x136 - x137 - x138 - x139 - x140;
auto x142 = V{0.00231481481481481}*x24*x28 + V{0.00231481481481481}*x26*x29;
auto x143 = V{0.166666666666667}*cell[12] - V{0.333333333333333}*cell[13] - V{0.166666666666667}*cell[14] + V{0.166666666666667}*cell[3] - V{0.166666666666667}*cell[5] - x119 - x120 + x121 + x122 + x123 + x124 - x125*x30 + x127 - x129 - x130 - x134 - x135 + x136 + x137 + x138 + x139 + x140 + x142*x51 + x142*x53 + V{-0.0555555555555555};
auto x144 = V{0.00578703703703704}*x24*x28 + V{0.00578703703703704}*x26*x29;
auto x145 = x24*x28 + x26*x29;
auto x146 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x145;
auto x147 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x146;
auto x148 = V{0.0208333333333333}*x24*x28 + V{0.0208333333333333}*x26*x29;
auto x149 = V{0.0208333333333333}*cell[0] + V{0.0208333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0208333333333333}*cell[3] + V{0.0208333333333333};
auto x150 = -x24*(V{0.0416666666666667}*cell[10] + V{0.0208333333333333}*cell[11] + V{0.0416666666666667}*cell[14] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0208333333333333}*cell[17] + V{0.0208333333333333}*cell[18] + V{0.0208333333333333}*cell[2] + V{0.0208333333333333}*cell[8] + V{0.0208333333333333}*cell[9] + x149);
auto x151 = -x26*(V{0.0208333333333333}*cell[10] + V{0.0416666666666667}*cell[11] + V{0.0208333333333333}*cell[15] + V{0.0208333333333333}*cell[16] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0208333333333333}*cell[1] + V{0.0416666666666667}*cell[5] + V{0.0208333333333333}*cell[6] + V{0.0208333333333333}*cell[7] + x149);
auto x152 = V{0.0416666666666667}*x24*x28 + V{0.0416666666666667}*x26*x29;
auto x153 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x150 + x151 - x152*x32 + V{0.0138888888888889};
auto x154 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] - x152*x34;
auto x155 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] + x148*x30 + x153 + x154;
auto x156 = V{0.833333333333333}*cell[13] - V{0.0833333333333333}*cell[14] - V{0.0833333333333333}*cell[5] - x147 + x155;
auto x157 = x108 + x90;
auto x158 = V{0.00115740740740741}*x24*x28 + V{0.00115740740740741}*x26*x29;
auto x159 = -V{0.166666666666667}*cell[13] + V{0.416666666666667}*cell[14] + V{0.416666666666667}*cell[5] + x147 + x155 + x158*x51 + x158*x53;
auto x160 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x146;
auto x161 = V{0.000578703703703704}*x24*x28 + V{0.000578703703703704}*x26*x29;
auto x162 = V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[5] - x152*x30 - x161*x51 - x161*x53;
auto x163 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] + x148*x34 + x153 + x162;
auto x164 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] - x160 + x163;
auto x165 = -x40*(x105 + x93);
auto x166 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] + x160 + x163;
auto x167 = V{0.125}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x145;
auto x168 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] + x148*x32 + x150 + x151 + x154 + x162 + V{0.0138888888888889};
auto x169 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] - x167 + x168;
auto x170 = -x40*(x59 + x69);
auto x171 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] + x167 + x168;
auto x172 = x24*(x102 + x170) + x26*(x117 - x157*x73 + x165);
auto x0 = -x19*(V{0.166666666666667}*x118*x37 + V{0.333333333333333}) + x20*(V{0.5}*cell[12] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21 + x22 + x23 - x24*(V{0.25}*cell[11] + V{0.5}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[16] + V{0.25}*cell[17] + V{0.25}*cell[18] + V{0.25}*cell[2] + V{0.25}*cell[8] + V{0.25}*cell[9] + x21 + x25) - x26*(V{0.25}*cell[10] + V{0.25}*cell[15] + V{0.25}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{0.25}*cell[1] + V{0.5}*cell[5] + V{0.25}*cell[6] + V{0.25}*cell[7] + x22 + x25) - x30*x39 - x32*x39 - x34*x39 + x38 - x40*x51 - x40*x53 + V{0.833333333333333});
auto x1 = -x19*(V{0.0277777777777778}*x115*x118 + V{0.0555555555555556}) + x20*(x116 + x133);
auto x2 = -x19*(V{0.0277777777777778}*x118*x86 + V{0.0555555555555556}) + x20*(x141 + x87);
auto x3 = -x19*(V{0.0277777777777778}*x118*x96 + V{0.0555555555555556}) - x20*(x143 + x97);
auto x4 = -x19*(V{0.0138888888888889}*x118*x53 + V{0.0277777777777778}) - x20*(-x144*x51 + x156 - x53*(V{0.0196759259259259}*x24*x28 + V{0.0196759259259259}*x26*x29));
auto x5 = -x19*(V{0.0138888888888889}*x109*x118 + V{0.0277777777777778}) - x20*(-x157*x40 + x159);
auto x6 = -x19*(V{0.0138888888888889}*x110*x118 + V{0.0277777777777778}) - x20*(x111 + x164);
auto x7 = -x19*(V{0.0138888888888889}*x106*x118 + V{0.0277777777777778}) - x20*(x165 + x166);
auto x8 = -x19*(V{0.0138888888888889}*x118*x63 + V{0.0277777777777778}) - x20*(x169 + x64);
auto x9 = -x19*(V{0.0138888888888889}*x118*x60 + V{0.0277777777777778}) - x20*(x170 + x171);
auto x10 = x19*(V{0.0277777777777778}*x172*x82 + V{-0.0555555555555556}) + x20*(-x114 + x133);
auto x11 = x19*(V{0.0277777777777778}*x172*x75 + V{-0.0555555555555556}) + x20*(x141 - x76);
auto x12 = x19*(V{0.0277777777777778}*x172*x98 + V{-0.0555555555555556}) - x20*(x143 + x99);
auto x13 = x19*(V{0.0138888888888889}*x172*x51 + V{-0.0277777777777778}) - x20*(-x144*x53 + x156 + x51*(V{0.00810185185185185}*x24*x28 + V{0.00810185185185185}*x26*x29));
auto x14 = -x19*(V{0.0138888888888889}*x118*x91 + V{0.0277777777777778}) - x20*(x159 - x40*x91);
auto x15 = x19*(V{0.0138888888888889}*x172*x79 + V{-0.0277777777777778}) - x20*(x112 + x164);
auto x16 = -x19*(V{0.0138888888888889}*x118*x94 + V{0.0277777777777778}) - x20*(x113 + x166);
auto x17 = x19*(V{0.0138888888888889}*x172*x66 + V{-0.0277777777777778}) - x20*(x169 + x67);
auto x18 = -x19*(V{0.0138888888888889}*x118*x71 + V{0.0277777777777778}) - x20*(x171 + x72);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { V{0.5}*x172, x30 + x32 + x34 };
}
};

}

}

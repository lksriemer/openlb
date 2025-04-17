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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<-1, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<-1, 1, -1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x20 = x19 + V{-1};
auto x21 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + V{-1});
auto x22 = V{0.166666666666667}*cell[0];
auto x23 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + V{0.333333333333333}*cell[6] + x22 + V{0.166666666666667};
auto x24 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + V{0.333333333333333}*cell[5];
auto x25 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2) + V{-1});
auto x26 = V{0.166666666666667}*cell[13];
auto x27 = V{0.166666666666667}*cell[14];
auto x28 = V{0.166666666666667}*cell[4];
auto x29 = V{0.166666666666667}*cell[5];
auto x30 = V{0.166666666666667}*cell[10] + V{0.333333333333333}*cell[18] + V{0.166666666666667}*cell[1];
auto x31 = V{1} / (cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + V{1});
auto x32 = V{0.166666666666667}*cell[15];
auto x33 = V{0.166666666666667}*cell[16];
auto x34 = V{0.166666666666667}*cell[6];
auto x35 = V{0.166666666666667}*cell[7];
auto x36 = cell[0] + cell[12] + cell[3] + V{2}*cell[5] + V{1};
auto x37 = cell[10] + V{2}*cell[18] + cell[1];
auto x38 = V{2}*cell[11] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + cell[6] + cell[7] + x36 + x37;
auto x39 = cell[11] + cell[2] + V{2}*cell[6];
auto x40 = cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[4] + V{2}*cell[7] + cell[8] + cell[9] + x36 + x39;
auto x41 = cell[0] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[8] + x37 + x39 + V{1};
auto x42 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x43 = V{1.5}*x42;
auto x44 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x45 = V{1.5}*x44;
auto x46 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x47 = V{1.5}*x46;
auto x48 = x45 + x47 + V{-1};
auto x49 = x43 + x48;
auto x50 = x49*(-V{0.111111111111111}*x21*x40 - V{0.111111111111111}*x25*x41 + V{0.111111111111111}*x31*x38);
auto x51 = -V{0.166666666666667}*x21*x40 - V{0.166666666666667}*x25*x41 + V{0.166666666666667}*x31*x38;
auto x52 = -V{0.00925925925925926}*x21*x40 - V{0.00925925925925926}*x25*x41 + V{0.00925925925925926}*x31*x38;
auto x53 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x54 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x55 = V{4.5}*(x54*x54);
auto x56 = -x43;
auto x57 = V{1} - x47;
auto x58 = x56 + x57;
auto x59 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x60 = -x45;
auto x61 = x59 + x60;
auto x62 = x53 + x55 + x58 + x61;
auto x63 = -V{5.55111512312578e-17}*x21*x40 - V{5.55111512312578e-17}*x25*x41 + V{5.55111512312578e-17}*x31*x38;
auto x64 = x49 + x59;
auto x65 = x53 - x55 + x64;
auto x66 = -cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x67 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1) + x66;
auto x68 = -x67;
auto x69 = -x59;
auto x70 = x49 + x53;
auto x71 = x69 + x70;
auto x72 = x71 - V{4.5}*x68*x68;
auto x73 = -x53;
auto x74 = -V{4.5}*x67*x67;
auto x75 = x64 + x73 + x74;
auto x76 = -V{0.0185185185185185}*x21*x40 - V{0.0185185185185185}*x25*x41 + V{0.0185185185185185}*x31*x38;
auto x77 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x78 = V{4.5}*(x77*x77);
auto x79 = V{3}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x80 = x49 + x79;
auto x81 = x53 - x78 + x80;
auto x82 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + x66;
auto x83 = -x82;
auto x84 = x69 + x80;
auto x85 = x84 - V{4.5}*x83*x83;
auto x86 = -V{0.037037037037037}*x21*x40 - V{0.037037037037037}*x25*x41 + V{0.037037037037037}*x31*x38;
auto x87 = V{3}*x44;
auto x88 = x43 + V{-1};
auto x89 = x47 + x79 - x87 + x88;
auto x90 = -V{7.40148683083438e-17}*x21*x40 - V{7.40148683083438e-17}*x25*x41 + V{7.40148683083438e-17}*x31*x38;
auto x91 = V{5.55111512312578e-17}*cell[0];
auto x92 = V{5.55111512312578e-17}*cell[11] + V{5.55111512312578e-17}*cell[2] + V{1.11022302462516e-16}*cell[6] + x91 + V{5.55111512312578e-17};
auto x93 = V{5.55111512312578e-17}*cell[12] + V{5.55111512312578e-17}*cell[3] + V{1.11022302462516e-16}*cell[5];
auto x94 = x21*(V{5.55111512312578e-17}*cell[17] + V{5.55111512312578e-17}*cell[18] + V{1.11022302462516e-16}*cell[1] + V{1.11022302462516e-16}*cell[4] + V{1.11022302462516e-16}*cell[7] + V{5.55111512312578e-17}*cell[8] + V{5.55111512312578e-17}*cell[9] + x92 + x93);
auto x95 = V{1.11022302462516e-16}*cell[3];
auto x96 = V{5.55111512312578e-17}*cell[10] + V{1.11022302462516e-16}*cell[18] + V{5.55111512312578e-17}*cell[1];
auto x97 = x25*(V{5.55111512312578e-17}*cell[13] + V{5.55111512312578e-17}*cell[14] + V{1.11022302462516e-16}*cell[16] + V{5.55111512312578e-17}*cell[4] + V{5.55111512312578e-17}*cell[5] + V{1.11022302462516e-16}*cell[8] + x92 + x95 + x96);
auto x98 = V{1.11022302462516e-16}*cell[11];
auto x99 = -x31*(V{1.11022302462516e-16}*cell[13] + V{5.55111512312578e-17}*cell[15] + V{5.55111512312578e-17}*cell[16] + V{1.11022302462516e-16}*cell[17] + V{5.55111512312578e-17}*cell[6] + V{5.55111512312578e-17}*cell[7] + x91 + x93 + x96 + x98 + V{5.55111512312578e-17});
auto x100 = V{3}*x46;
auto x101 = x53 + x60;
auto x102 = x100 + x101 + x56 + V{1};
auto x103 = -x100 + x45 + x53 + x88;
auto x104 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) + cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x105 = V{4.5}*(x104*x104);
auto x106 = -x105 + x59 + x80;
auto x107 = -x50;
auto x108 = x102*x76 - x103*x76 - x106*x76 + x107 + x94 + x97 + x99 + V{-2.22044604925031e-16};
auto x109 = V{3}*x42;
auto x110 = x109 + x57 + x61;
auto x111 = -x109 + x48 + x59;
auto x112 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0) - cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x113 = -x112;
auto x114 = x73 + x80;
auto x115 = x114 - V{4.5}*x113*x113;
auto x116 = -V{3.70074341541719e-17}*x21*x40 - V{3.70074341541719e-17}*x25*x41 + V{3.70074341541719e-17}*x31*x38;
auto x117 = V{1.11022302462516e-16}*cell[12] + V{3.33066907387547e-16}*cell[13] + V{3.33066907387547e-16}*cell[14] + V{3.33066907387547e-16}*cell[4] + V{3.33066907387547e-16}*cell[5] + x110*x76 - x111*x76 - x115*x76 - x116*x42 + x95;
auto x118 = x58 + x79;
auto x119 = x101 + x118 + x78;
auto x120 = -x79;
auto x121 = -V{4.5}*x112*x112;
auto x122 = x120 + x121 + x70;
auto x123 = -V{4.5}*x82*x82;
auto x124 = x120 + x123 + x64;
auto x125 = -x76*x89;
auto x126 = x118 + x87;
auto x127 = x126*x76;
auto x128 = V{1.66533453693773e-16}*cell[10] + V{2.22044604925031e-16}*cell[18] + V{1.66533453693773e-16}*cell[1] + V{2.22044604925031e-16}*cell[9] + x125 + x127 - x44*x63 - x75*x76;
auto x129 = x105 + x118 + x61;
auto x130 = x21*(V{2.22044604925031e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{4.44089209850063e-16}*cell[15] + V{2.22044604925031e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{2.22044604925031e-16}*cell[1] + V{1.66533453693773e-16}*cell[2] + V{4.44089209850063e-16}*cell[6] + V{2.22044604925031e-16}*cell[7] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + x108 + x117 - x44*x90 - x46*x63 + x52*x62 - x52*x65 - x52*x72 - x52*x75 - x76*x81 - x76*x85 - x86*x89) + x25*(V{2.22044604925031e-16}*cell[12] + V{4.71844785465692e-16}*cell[13] + V{5.27355936696949e-16}*cell[14] + V{3.33066907387547e-16}*cell[15] + V{3.33066907387547e-16}*cell[16] + V{1.11022302462516e-16}*cell[2] + V{2.22044604925031e-16}*cell[3] + V{4.71844785465692e-16}*cell[4] + V{5.27355936696949e-16}*cell[5] + V{3.33066907387547e-16}*cell[6] + V{3.33066907387547e-16}*cell[7] + x108 - x111*x86 - x115*x52 - x116*x46 + x119*x52 - x122*x52 - x124*x76 + x128 - x42*x90 - x52*x81 - x65*x76 + x98) - x31*(V{2.22044604925031e-16}*cell[11] + V{4.71844785465692e-16}*cell[15] + V{5.27355936696949e-16}*cell[16] + V{4.44089209850063e-16}*cell[17] + V{2.22044604925031e-16}*cell[2] + V{4.71844785465692e-16}*cell[6] + V{5.27355936696949e-16}*cell[7] + V{4.44089209850063e-16}*cell[8] + x102*x86 - x106*x52 + x107 + x117 + x119*x76 - x124*x52 + x128 + x129*x52 - x46*x90 - x52*x85 + x62*x76 + x94 + x97 + x99 + V{-2.22044604925031e-16});
auto x131 = -x130;
auto x132 = -V{0.0277777777777778}*x21*x40 - V{0.0277777777777778}*x25*x41 + V{0.0277777777777778}*x31*x38;
auto x133 = x21*x40 + x25*x41 - x31*x38;
auto x134 = V{3.46944695195361e-18}*cell[0];
auto x135 = V{3.46944695195361e-18}*cell[12] + V{3.46944695195361e-18}*cell[3] + V{6.93889390390723e-18}*cell[5] + x134 + V{3.46944695195361e-18};
auto x136 = V{3.46944695195361e-18}*cell[11] + V{3.46944695195361e-18}*cell[2] + V{6.93889390390723e-18}*cell[6];
auto x137 = x21*(V{3.46944695195361e-18}*cell[17] + V{3.46944695195361e-18}*cell[18] + V{6.93889390390723e-18}*cell[1] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[7] + V{3.46944695195361e-18}*cell[8] + V{3.46944695195361e-18}*cell[9] + x135 + x136);
auto x138 = V{3.46944695195361e-18}*cell[10] + V{6.93889390390723e-18}*cell[18] + V{3.46944695195361e-18}*cell[1];
auto x139 = x25*(V{3.46944695195361e-18}*cell[13] + V{3.46944695195361e-18}*cell[14] + V{6.93889390390723e-18}*cell[16] + V{6.93889390390723e-18}*cell[3] + V{3.46944695195361e-18}*cell[4] + V{3.46944695195361e-18}*cell[5] + V{6.93889390390723e-18}*cell[8] + x134 + x136 + x138 + V{3.46944695195361e-18});
auto x140 = -x31*(V{6.93889390390723e-18}*cell[11] + V{6.93889390390723e-18}*cell[13] + V{3.46944695195361e-18}*cell[15] + V{3.46944695195361e-18}*cell[16] + V{6.93889390390723e-18}*cell[17] + V{3.46944695195361e-18}*cell[6] + V{3.46944695195361e-18}*cell[7] + x135 + x138);
auto x141 = -V{0.0833333333333333}*cell[12] + V{0.0833333333333334}*cell[13] + V{0.0833333333333334}*cell[14] - V{0.0833333333333333}*cell[3] + V{0.0833333333333334}*cell[4] + V{0.0833333333333334}*cell[5] + x137 + x139 + x140 + V{-0.0555555555555555};
auto x142 = -V{0.0833333333333333}*cell[11] + V{0.0833333333333334}*cell[15] + V{0.0833333333333334}*cell[16] - V{0.0833333333333333}*cell[2] + V{0.0833333333333334}*cell[6] + V{0.0833333333333334}*cell[7];
auto x143 = V{0.166666666666667}*cell[10] - V{0.166666666666667}*cell[17] - V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[1] - V{0.166666666666667}*cell[8] - V{0.166666666666667}*cell[9] + x132*x42 + x132*x46 + V{0.0555555555555556}*x133*x44 + x141 + x142;
auto x144 = V{0.0185185185185185}*x21*x40 + V{0.0185185185185185}*x25*x41 - V{0.0185185185185185}*x31*x38;
auto x145 = V{0.0555555555555556}*x21*x40 + V{0.0555555555555556}*x25*x41 - V{0.0555555555555556}*x31*x38;
auto x146 = V{0.0277777777777778}*x21*x40 + V{0.0277777777777778}*x25*x41 - V{0.0277777777777778}*x31*x38;
auto x147 = -V{0.0833333333333333}*cell[10] + V{0.0833333333333334}*cell[17] + V{0.0833333333333334}*cell[18] - V{0.0833333333333333}*cell[1] + V{0.0833333333333334}*cell[8] + V{0.0833333333333334}*cell[9] - x146*x44;
auto x148 = V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[2] + x141 + x145*x46 - x146*x42 + x147 - x32 - x33 - x34 - x35;
auto x149 = V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[3] + x137 + x139 + x140 + x142 + x145*x42 - x146*x46 + x147 - x26 - x27 - x28 - x29 + V{-0.0555555555555555};
auto x150 = V{0.00925925925925926}*x21*x40 + V{0.00925925925925926}*x25*x41 - V{0.00925925925925926}*x31*x38;
auto x151 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(0)*x133;
auto x152 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*x151;
auto x153 = V{0.0138888888888889}*x21*x40 + V{0.0138888888888889}*x25*x41 - V{0.0138888888888889}*x31*x38;
auto x154 = V{0.0138888888888889}*cell[0];
auto x155 = V{0.0138888888888889}*cell[11] + V{0.0138888888888889}*cell[2] + V{0.0277777777777778}*cell[6] + x154 + V{0.0138888888888889};
auto x156 = V{0.0138888888888889}*cell[12] + V{0.0138888888888889}*cell[3] + V{0.0277777777777778}*cell[5];
auto x157 = x21*(V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[7] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + x155 + x156);
auto x158 = V{0.0138888888888889}*cell[10] + V{0.0277777777777778}*cell[18] + V{0.0138888888888889}*cell[1];
auto x159 = x25*(V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[3] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0277777777777778}*cell[8] + x155 + x158);
auto x160 = -x31*(V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[13] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + x154 + x156 + x158 + V{0.0138888888888889});
auto x161 = V{0.0277777777777778}*x21*x40 + V{0.0277777777777778}*x25*x41 - V{0.0277777777777778}*x31*x38;
auto x162 = V{0.0833333333333333}*cell[10] + V{0.0416666666666667}*cell[17] + V{0.0416666666666667}*cell[18] + V{0.0833333333333333}*cell[1] + V{0.0416666666666667}*cell[8] + V{0.0416666666666667}*cell[9] + x157 + x159 + x160 + x161*x44 + V{0.0138888888888889};
auto x163 = V{0.0833333333333333}*cell[11] + V{0.0416666666666667}*cell[15] + V{0.0416666666666667}*cell[16] + V{0.0833333333333333}*cell[2] + V{0.0416666666666667}*cell[6] + V{0.0416666666666667}*cell[7] + x161*x46;
auto x164 = -V{0.0416666666666667}*cell[12] - V{0.0416666666666667}*cell[3] - x153*x42 + x162 + x163;
auto x165 = V{0.416666666666667}*cell[13] - V{0.0833333333333333}*cell[14] + V{0.416666666666667}*cell[4] - V{0.0833333333333333}*cell[5] + x152 + x164;
auto x166 = -V{0.0833333333333333}*cell[13] + V{0.416666666666667}*cell[14] - V{0.0833333333333333}*cell[4] + V{0.416666666666667}*cell[5] - x152 + x164;
auto x167 = cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x151;
auto x168 = V{0.0833333333333333}*cell[12] + V{0.0416666666666667}*cell[13] + V{0.0416666666666667}*cell[14] + V{0.0833333333333333}*cell[3] + V{0.0416666666666667}*cell[4] + V{0.0416666666666667}*cell[5] + x161*x42;
auto x169 = -V{0.0416666666666667}*cell[11] - V{0.0416666666666667}*cell[2] - x153*x46 + x162 + x168;
auto x170 = V{0.416666666666667}*cell[15] - V{0.0833333333333333}*cell[16] + V{0.416666666666667}*cell[6] - V{0.0833333333333333}*cell[7] + x167 + x169;
auto x171 = -V{0.0833333333333333}*cell[15] + V{0.416666666666667}*cell[16] - V{0.0833333333333333}*cell[6] + V{0.416666666666667}*cell[7] - x167 + x169;
auto x172 = V{0.0833333333333333}*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(1)*cell.template getFieldComponent<momenta::FixedVelocityMomentumGeneric::VELOCITY>(2)*x133;
auto x173 = -V{0.0416666666666667}*cell[10] - V{0.0416666666666667}*cell[1] - x153*x44 + x157 + x159 + x160 + x163 + x168 + V{0.0138888888888889};
auto x174 = V{0.416666666666667}*cell[17] - V{0.0833333333333333}*cell[18] + V{0.416666666666667}*cell[8] - V{0.0833333333333333}*cell[9] + x172 + x173;
auto x175 = -V{0.0833333333333333}*cell[17] + V{0.416666666666667}*cell[18] - V{0.0833333333333333}*cell[8] + V{0.416666666666667}*cell[9] - x172 + x173;
auto x0 = -x19*(V{0.111111111111111}*x131*x49 + V{0.333333333333333}) + x20*(V{0.5}*cell[10] + V{0.5}*cell[11] + V{0.5}*cell[12] + V{1}*cell[13] + V{1}*cell[14] + V{1}*cell[15] + V{1}*cell[16] + V{1}*cell[17] + V{1}*cell[18] + V{0.5}*cell[1] + V{0.5}*cell[2] + V{0.5}*cell[3] + V{1}*cell[4] + V{1}*cell[5] + V{1}*cell[6] + V{1}*cell[7] + V{1}*cell[8] + V{1}*cell[9] + x21*(V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[7] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x23 + x24) + x25*(V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[8] + x23 + x26 + x27 + x28 + x29 + x30) - x31*(V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[17] + x22 + x24 + x30 + x32 + x33 + x34 + x35 + V{0.166666666666667}) - x42*x51 - x44*x51 - x46*x51 + x50 + V{0.833333333333333});
auto x1 = -x19*(V{0.0185185185185185}*x131*x89 + V{0.0555555555555556}) + x20*(-x125 - x143);
auto x2 = -x19*(V{0.0185185185185185}*x103*x131 + V{0.0555555555555556}) - x20*(x103*x144 + x148);
auto x3 = -x19*(V{0.0185185185185185}*x111*x131 + V{0.0555555555555556}) - x20*(x111*x144 + x149);
auto x4 = -x19*(V{0.00925925925925926}*x131*x81 + V{0.0277777777777778}) - x20*(x150*x81 + x165);
auto x5 = -x19*(V{0.00925925925925926}*x115*x131 + V{0.0277777777777778}) - x20*(x150*(x114 + x121) + x166);
auto x6 = -x19*(V{0.00925925925925926}*x106*x131 + V{0.0277777777777778}) - x20*(x106*x150 + x170);
auto x7 = -x19*(V{0.00925925925925926}*x131*x85 + V{0.0277777777777778}) - x20*(x150*(x123 + x84) + x171);
auto x8 = -x19*(V{0.00925925925925926}*x131*x65 + V{0.0277777777777778}) - x20*(x150*x65 + x174);
auto x9 = -x19*(V{0.00925925925925926}*x131*x72 + V{0.0277777777777778}) - x20*(x150*(x71 + x74) + x175);
auto x10 = -x19*(V{0.0185185185185185}*x126*x130 + V{0.0555555555555556}) + x20*(-x127 - x143);
auto x11 = -x19*(V{0.0185185185185185}*x102*x130 + V{0.0555555555555556}) - x20*(-x102*x144 + x148);
auto x12 = -x19*(V{0.0185185185185185}*x110*x130 + V{0.0555555555555556}) - x20*(-x110*x144 + x149);
auto x13 = -x19*(V{0.00925925925925926}*x119*x130 + V{0.0277777777777778}) - x20*(-x119*x150 + x165);
auto x14 = -x19*(V{0.00925925925925926}*x122*x131 + V{0.0277777777777778}) - x20*(x122*x150 + x166);
auto x15 = -x19*(V{0.00925925925925926}*x129*x130 + V{0.0277777777777778}) - x20*(-x129*x150 + x170);
auto x16 = -x19*(V{0.00925925925925926}*x124*x131 + V{0.0277777777777778}) - x20*(x124*x150 + x171);
auto x17 = -x19*(V{0.00925925925925926}*x130*x62 + V{0.0277777777777778}) - x20*(-x150*x62 + x174);
auto x18 = -x19*(V{0.00925925925925926}*x131*x75 + V{0.0277777777777778}) - x20*(x150*x75 + x175);
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
return { V{0.333333333333333}*x131, x42 + x44 + x46 };
}
};

}

}

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Fedor Bukreev
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

#ifndef STRAIN_RATE_TENSOR_FDM_3D_H
#define STRAIN_RATE_TENSOR_FDM_3D_H

namespace olb {

template <typename CELL, typename V = typename CELL::value_t>
Vector<V,6> strainRateTensorFDM3D(CELL& cell) any_platform
{
  using namespace olb::util::tensorIndices3D;
  V inv2dx = 1./(2.);
  V u_pXp[3], u_pXm[3], u_pYp[3], u_pYm[3], u_pZp[3], u_pZm[3] {0.};

  cell.neighbor({1,0,0}).computeU(u_pXp);
  cell.neighbor({-1,0,0}).computeU(u_pXm);
  cell.neighbor({0,1,0}).computeU(u_pYp);
  cell.neighbor({0,-1,0}).computeU(u_pYm);
  cell.neighbor({0,0,1}).computeU(u_pZp);
  cell.neighbor({0,0,-1}).computeU(u_pZm);

  Vector<V,6> strainRate;
  strainRate[xx] = (u_pXp[0] - u_pXm[0])*inv2dx;
  strainRate[xy] = 0.5*( (u_pYp[0] - u_pYm[0])*inv2dx + (u_pXp[1] - u_pXm[1])*inv2dx );
  strainRate[xz] = 0.5*( (u_pZp[0] - u_pZm[0])*inv2dx + (u_pXp[2] - u_pXm[2])*inv2dx );
  strainRate[yy] = (u_pYp[1] - u_pYm[1])*inv2dx;
  strainRate[yz] = 0.5*( (u_pZp[1] - u_pZm[1])*inv2dx + (u_pYp[2] - u_pYm[2])*inv2dx );
  strainRate[zz] = (u_pZp[2] - u_pZm[2])*inv2dx;
  return strainRate;
};
}

#endif

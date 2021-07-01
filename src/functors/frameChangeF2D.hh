/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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

#ifndef FRAME_CHANGE_F_2D_HH
#define FRAME_CHANGE_F_2D_HH

#include<vector>
#include<cmath>
#include<string>

#include "frameChangeF2D.h"
#include "functors/genericF.h"
#include "functors/analyticalF.h"
#include "functors/superLatticeBaseF2D.h"
#include "functors/superLatticeLocalF2D.h"
#include "core/superLattice2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "utilities/vectorHelpers.h"  // for normalize

namespace olb {

template <typename T>
Poiseuille2D<T>::Poiseuille2D(std::vector<T> axisPoint_, std::vector<T> axisDirection_,  T maxVelocity_, T radius_)
  : AnalyticalF2D<T,T>(2)
{
  axisPoint.resize(2);
  axisDirection.resize(2);
  for(int i = 0; i < 2; i++) {
    axisDirection[i] = axisDirection_[i];
    axisPoint[i] = axisPoint_[i];
  }
  maxVelocity=maxVelocity_;
  radius = radius_;
}


template <typename T>
Poiseuille2D<T>::Poiseuille2D(
  SuperGeometry2D<T>& superGeometry_, int material_, T maxVelocity_, T distance2Wall_)
  : AnalyticalF2D<T,T>(2)
{
  axisPoint = superGeometry_.getStatistics().getCenterPhysR(material_);

  axisPoint = superGeometry_.getStatistics().getCenterPhysR(material_);
  std::vector<int> discreteNormal = superGeometry_.getStatistics().computeDiscreteNormal(material_);
  axisDirection.push_back((T)(discreteNormal[0]));
  axisDirection.push_back((T)(discreteNormal[1]));

  radius = T(distance2Wall_);
  for (int iD = 0; iD < 2; iD++)
    radius += (superGeometry_.getStatistics().getPhysRadius(material_)[iD]);
  maxVelocity = maxVelocity_;
}


template <typename T>
std::vector<T> Poiseuille2D<T>::operator()(std::vector<T> x)
{
  std::vector<T> velocity;
  velocity.resize(2);
  velocity[0] = maxVelocity*axisDirection[0]*(1.-((x[1]-axisPoint[1])*(x[1]-axisPoint[1]))/radius/radius);
  velocity[1] = maxVelocity*axisDirection[1]*(1.-((x[0]-axisPoint[0])*(x[0]-axisPoint[0]))/radius/radius);

  return velocity;
}



template <typename T, typename S>
PoiseuilleStrainRate2D<T,S>::PoiseuilleStrainRate2D(LBconverter<S> const& converter, T ly)
  : AnalyticalF2D<T,S>(4)
{
  lengthY = ly;
  maxVelocity = converter.getCharU();
  this->_name = "PoiseuilleStrainRate2D";
}


template <typename T, typename S>
std::vector<T> PoiseuilleStrainRate2D<T,S>::operator()(std::vector<S> input)
{
  std::vector<T> output(4,T());
  T y = input[1];

  T DuDx = T();
  T DuDy = (T) maxVelocity*(-2.*(y-(lengthY/2.))/((lengthY/2.)*(lengthY/2.)));
  T DvDx = T();
  T DvDy = T();

  output[0] = (DuDx + DuDx)/2.;
  output[1] = (DuDy + DvDx)/2.;
  output[2] = (DvDx + DuDy)/2.;
  output[3] = (DvDy + DvDy)/2.;

  return output;
};



} // end namespace olb
#endif

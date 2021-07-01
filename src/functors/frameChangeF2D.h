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

#ifndef FRAME_CHANGE_F_2D_H
#define FRAME_CHANGE_F_2D_H

#include<vector>
#include<cmath>
#include<string>
#include"math.h"

#include "functors/genericF.h"
#include "functors/analyticalF.h"
#include "functors/superLatticeBaseF2D.h"
#include "functors/superLatticeLocalF2D.h"
#include "core/superLattice2D.h"

/** \file
  This file contains two different classes of functors, in the FIRST part
  - for simulations in a rotating frame
  - different functors for
      velocity (3d, RotatingLinear3D),
      pressure (1d, RotatingQuadratic1D) and
      force    (3d, RotatingForceField3D)
  The functors return the displacement of a point x in a fixed amount of time.

  The ones in the SECOND part are useful to set Poiseuille velocity profiles on
  - pipes with round cross-section and
  - pipes with square-shaped cross-section.
*/

/** To enable simulations in a rotating frame, the axis is set in the
  * constructor with axisPoint and axisDirection. The axisPoint can be the
  * coordinate of any point on the axis. The axisDirection has to be a normed to
  * 1. The pulse w is in rad/s. It determines the pulse speed by its norm while
  * the trigonometric or clockwise direction is determined by its sign: When the
  * axisDirection is pointing "towards you", a positive pulse makes it turn in
  * the trigonometric way. It has to be noticed that putting both axisDirection
  * into -axisDirection and w into -w yields an exactly identical situation.
  */


namespace olb {


template <typename T>
class Poiseuille2D : public AnalyticalF2D<T,T> {
protected:
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T maxVelocity;
  T radius;

public:
  Poiseuille2D(std::vector<T> axisPoint_, std::vector<T> axisDirection_,  T maxVelocity_, T radius_);
  /// construct from material number, note: untested
  Poiseuille2D(SuperGeometry2D<T>& superGeometry_, int material_, T maxVelocity_, T distance2Wall_);
  std::vector<T> operator()(std::vector<T> x);
};


template <typename T, typename S>
class PoiseuilleStrainRate2D : public AnalyticalF2D<T,S> {
protected:
	T lengthY;
  T maxVelocity;

public:
  PoiseuilleStrainRate2D(LBconverter<S> const& converter, T ly);
  std::vector<T> operator()(std::vector<S> input);
};


} // end namespace olb
#endif

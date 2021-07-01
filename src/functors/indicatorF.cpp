/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Lukas Baron, Mathias J. Krause
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

#include "io/stlReader.h"
#include "indicatorF.h"
#include "indicatorF.hh"

namespace olb {

// stl indicator functors
template class STLreader<double>;

template class IndicatorStl3D<bool,double>;


// indicator functors 2D
template class IndicatorCuboid2D<bool,double>;
template class IndicatorCircle2D<bool,double>;

// indicator functors 3D
template class IndicatorCircle3D<bool,double>;
template class IndicatorSphere3D<bool,double>;
template class IndicatorLayer3D<bool,double>;
template class IndicatorCylinder3D<bool,double>;
template class IndicatorCone3D<bool,double>;
template class IndicatorCuboid3D<bool,double>;
template IndicatorCuboid3D<bool,double>* createIndicatorCuboid3D(XMLreader const& params, bool verbose);
template class IndicatorParallelepiped3D<bool,double>;


// smoothIndicator functors
// double since they return values /in [0,1]
template class SmoothIndicatorCircle2D<double,double>;

template class SmoothIndicatorSphere3D<double,double>;
template class SmoothIndicatorCylinder3D<double,double>;
template class SmoothIndicatorCone3D<double,double>;

}

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause, Cyril Masquelier,
 *  Albert Mink
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


#include "functors/indicCalcF.h"
#include "functors/indicCalcF.hh"

namespace olb {


// arithmetic helper class for Indical 1d functors
template class IndicCalc1D<bool,double>;

template class IndicPlus1D<bool,double>;

template class IndicMinus1D<bool,double>;

template class IndicMultiplication1D<bool,double>;


// arithmetic helper class for Indical 2d functors
template class IndicCalc2D<bool,double>;

template class IndicPlus2D<bool,double>;

template class IndicMinus2D<bool,double>;

template class IndicMultiplication2D<bool,double>;



// arithmetic helper class for Indical 3d functors
template class IndicCalc3D<bool,double>;

template class IndicPlus3D<bool,double>;

template class IndicMinus3D<bool,double>;

template class IndicMultiplication3D<bool,double>;


// arithmetic helper class for Indical 3d functors Smooth
template class SmoothIndicCalc3D<double,double>;

template class SmoothIndicPlus3D<double,double>;

}


/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_BASE_F_HH
#define ANALYTICAL_BASE_F_HH

#include<vector>

#include "analyticalBaseF.h"

namespace olb {


// identity to "store results"
template <typename T, typename S>
AnalyticalIdentity2D<T,S>::AnalyticalIdentity2D(AnalyticalF2D<T,S>& f)
  : AnalyticalF2D<T,S>(f.getTargetDim()), _f(f)
{
  this->_name = _f.getName();
  // add 'this' to father's child list to prevent father from being deleted
  _f.addChild(this);
}

template <typename T, typename S>
AnalyticalIdentity2D<T,S>::~AnalyticalIdentity2D()
{
  // remove 'this' from father's child list
  _f.removeChild(this);
  // delete father from grandfather's child list
  _f.myErase(NULL);
}

template <typename T, typename S>
std::vector<T> AnalyticalIdentity2D<T,S>::operator()(std::vector<S> input) 
{
  return _f(input);
}

// identity to "store results"
template <typename T, typename S>
AnalyticalIdentity3D<T,S>::AnalyticalIdentity3D(AnalyticalF3D<T,S>& f)
  : AnalyticalF3D<T,S>(f.getTargetDim()), _f(f)
{
  this->_name = _f.getName();
  // add 'this' to father's child list to prevent father from being deleted
  _f.addChild(this);
}

template <typename T, typename S>
AnalyticalIdentity3D<T,S>::~AnalyticalIdentity3D()
{
  // remove 'this' from father's child list
  _f.removeChild(this);
  // delete father from grandfather's child list
  _f.myErase(NULL);
}

template <typename T, typename S>
std::vector<T> AnalyticalIdentity3D<T,S>::operator()(std::vector<S> input) 
{
  return _f(input);
}



} // end namespace olb

#endif

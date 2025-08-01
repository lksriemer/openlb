/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2018 Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_INDICATOR_F_2D_HH
#define SUPER_INDICATOR_F_2D_HH

#include <numeric>

#include "superIndicatorF2D.h"
#include "blockIndicatorF2D.h"
#include "core/util.h"

namespace olb {

template <typename T>
SuperIndicatorFfromIndicatorF2D<T>::SuperIndicatorFfromIndicatorF2D(
  FunctorPtr<IndicatorF2D<T>>&& indicatorF, SuperGeometry<T,2>& geometry)
  : SuperIndicatorF2D<T>(geometry),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "SuperIndicator_from_" + _indicatorF->getName();

  LoadBalancer<T>& load = this->getSuperStructure().getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFfromIndicatorF2D<T>(
        *_indicatorF, geometry.getBlockGeometry(iC))
    );
  }
}

template <typename T>
bool SuperIndicatorFfromIndicatorF2D<T>::operator() (bool output[], const int input[])
{
  LatticeR<3> latticeR(input);
  auto physR = this->_superStructure.getCuboidDecomposition().getPhysR(input);
  return _indicatorF(output, physR.data());
}


template <typename T, bool HLBM>
SuperIndicatorFfromSmoothIndicatorF2D<T, HLBM>::SuperIndicatorFfromSmoothIndicatorF2D(
  FunctorPtr<SmoothIndicatorF2D<T,T,HLBM>>&& indicatorF,
  SuperGeometry<T,2>&                     geometry)
  : SuperIndicatorF2D<T>(geometry),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "SuperIndicator_from_" + _indicatorF->getName();

  LoadBalancer<T>& load = this->getSuperStructure().getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFfromSmoothIndicatorF2D<T, HLBM>(
        *_indicatorF, geometry.getBlockGeometry(iC))
    );
  }
}

template <typename T, bool HLBM>
bool SuperIndicatorFfromSmoothIndicatorF2D<T, HLBM>::operator() (bool output[], const int input[])
{
  T physR[2];
  T inside[1];
  this->_superStructure.getCuboidDecomposition().getPhysR(physR, input);
  _indicatorF(inside, physR);
  return !util::nearZero(inside[0]);
}

template <typename T>
SuperIndicatorMaterial2D<T>::SuperIndicatorMaterial2D(
  SuperGeometry<T,2>& geometry, std::vector<int> materials)
  : SuperIndicatorF2D<T>(geometry)
{
  geometry.updateStatistics(false);
  const std::string matString = std::accumulate(
                                  materials.begin()+1,
                                  materials.end(),
                                  std::to_string(materials[0]),
  [](const std::string& a, int b) {
    return a + '_' + std::to_string(b);
  });
  this->getName() = "SuperIndicator_on_Material_" + matString;

  for (int iC = 0; iC < this->_superGeometry.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorMaterial2D<T>(this->_superGeometry.getBlockGeometry(iC),
                                      materials)
    );
  }
}

template <typename T>
SuperIndicatorMaterial2D<T>::SuperIndicatorMaterial2D(
  SuperGeometry<T,2>& geometry, std::list<int> materials)
  : SuperIndicatorMaterial2D(geometry,
                             std::vector<int>(materials.begin(), materials.end()))
{ }

template <typename T>
bool SuperIndicatorMaterial2D<T>::operator() (bool output[], const int input[])
{
  output[0] = false;

  LoadBalancer<T>& load = this->_superGeometry.getLoadBalancer();

  if (!this->_blockF.empty() && load.isLocal(input[0])) {
    // query material number of appropriate block indicator
    return this->getBlockF(load.loc(input[0]))(output, &input[1]);
  }
  else {
    return false;
  }
}

template <typename T>
SuperIndicatorIdentity2D<T>::SuperIndicatorIdentity2D(FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF)
  : SuperIndicatorF2D<T>(indicatorF->getSuperGeometry()),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = _indicatorF->getName();

  for (int iC = 0; iC < _indicatorF->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorIdentity2D<T>(_indicatorF->getBlockIndicatorF(iC)));
  }
}

template <typename T>
bool SuperIndicatorIdentity2D<T>::operator()(bool output[], const int input[])
{
  return _indicatorF(output, input);
}

template <typename T>
SuperIndicatorBoundaryNeighbor2D<T>::SuperIndicatorBoundaryNeighbor2D(FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF, int overlap)
  : SuperIndicatorF2D<T>(indicatorF->getSuperGeometry()),
    _indicatorF(std::move(indicatorF)),
    _overlap(overlap)
{
  this->getName() = "SuperIndicatorBoundaryNeighbor_for_BoundaryMaterial_" + _indicatorF->getName();

  for (int iC = 0; iC < _indicatorF->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorBoundaryNeighbor2D<T>(_indicatorF->getBlockIndicatorF(iC), _overlap));
  }
}

template <typename T>
bool SuperIndicatorBoundaryNeighbor2D<T>::operator() (bool output[], const int input[])
{
  output[0] = false;

  LoadBalancer<T>& load = this->_superGeometry.getLoadBalancer();

  if (load.isLocal(input[0])) {
    return this->getBlockF(load.loc(input[0]))(output, &input[1]);
  }
  else {
    return false;
  }
}

} // namespace olb

#endif

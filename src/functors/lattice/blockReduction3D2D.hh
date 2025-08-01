/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef BLOCK_REDUCTION_3D2D_HH
#define BLOCK_REDUCTION_3D2D_HH

#include "blockReduction3D2D.h"

#include <limits>
#include "utilities/omath.h"

#include "utilities/vectorHelpers.h"
#include "functors/analytical/interpolationF3D.h"
#include "communication/mpiManager.h"
#include "utilities/functorPtr.hh"

namespace olb {


template <typename T>
void BlockReduction3D2D<T>::updateBlockAnalytical(BlockData<2,T,T>& block)
{
  AnalyticalFfromSuperF3D<T> analyticalF(*_f);

  for ( std::tuple<int,int,int>& pos : _rankLocalSubplane ) {
    const int& iX = std::get<0>(pos);
    const int& iY = std::get<1>(pos);
    const Vector<T,3> physR = this->getPhysR(iX, iY);

    for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
      block.get({iX, iY}, iSize) = T();
    }

    T output[_f->getTargetDim()];
    const T input[3] { physR[0], physR[1], physR[2] };

    if (analyticalF(output, input)) {
      for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
        block.get({iX, iY}, iSize) += output[iSize];
      }
    }
  }
}

template <typename T>
void BlockReduction3D2D<T>::updateBlockDiscrete(BlockData<2,T,T>& block)
{
  auto& geometry = _f->getSuperStructure().getCuboidDecomposition();

  for ( std::tuple<int,int,int>& pos : _rankLocalSubplane ) {
    const int& iX = std::get<0>(pos);
    const int& iY = std::get<1>(pos);
    const int& iC = std::get<2>(pos);
    const Vector<T,3> physR = this->getPhysR(iX, iY);

    for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
      block.get({iX, iY}, iSize) = T();
    }

    T output[_f->getTargetDim()];
    auto input = geometry.get(iC).getLatticeR(physR).withPrefix(iC);

    if (_f(output, input.data())) {
      for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
        block.get({iX, iY}, iSize) += output[iSize];
      }
    }
  }
}

template <typename T>
BlockReduction3D2D<T>::BlockReduction3D2D(
  FunctorPtr<SuperF3D<T>>&& f,
  const HyperplaneLattice3D<T>& lattice,
  BlockDataSyncMode      syncMode,
  BlockDataReductionMode reductionMode)
  : HyperplaneLattice3D<T>(lattice),
    BlockDataF2D<T,T>(lattice.getNx(), lattice.getNy(), f->getTargetDim()),
    _f(std::move(f)),
    _syncMode(syncMode),
    _reductionMode(reductionMode)
{
  this->getName() = "planeReduction(" + _f->getName() + ")";

  if ( _reductionMode == BlockDataReductionMode::Discrete ) {
    const auto& geometry = _f->getSuperStructure().getCuboidDecomposition();
    const Hyperplane3D<T>& hyperplane   = this->getHyperplane();
    const bool spansAxisPlane = hyperplane.isXYPlane() ||
                                hyperplane.isXZPlane() ||
                                hyperplane.isYZPlane();
    // verify axes alignment and spacing of hyperplane parametrization
    if ( !spansAxisPlane ||
         lattice.getPhysSpacing() != geometry.getDeltaR() ) {
      // hyperplane lattice doesn't describe a trivially discretizable plane
      OstreamManager clerr(std::cerr, "BlockReduction3D2D");
      clerr << "Given hyperplane is not trivially discretizable. "
            << "Use BlockDataReductionMode::Analytical instead."
            << std::endl;
      exit(-1);
    }
  }

  // intialize list of relevant rank local points making up the reduced plane
  initialize();
  // first update of data
  update();
}

template <typename T>
BlockReduction3D2D<T>::BlockReduction3D2D(
  FunctorPtr<SuperF3D<T>>&& f,
  const Hyperplane3D<T>& hyperplane,
  BlockDataSyncMode      syncMode,
  BlockDataReductionMode reductionMode)
  : BlockReduction3D2D(
      std::forward<decltype(f)>(f),
      HyperplaneLattice3D<T>(f->getSuperStructure().getCuboidDecomposition(),
                             hyperplane),
      syncMode,
      reductionMode)
{ }

template <typename T>
BlockReduction3D2D<T>::BlockReduction3D2D(
  FunctorPtr<SuperF3D<T>>&& f,
  const Hyperplane3D<T>& hyperplane,
  int                    resolution,
  BlockDataSyncMode      syncMode)
  : BlockReduction3D2D(
      std::forward<decltype(f)>(f),
      HyperplaneLattice3D<T>(f->getSuperStructure().getCuboidDecomposition(),
                             hyperplane, resolution),
      syncMode)
{ }

template <typename T>
BlockReduction3D2D<T>::BlockReduction3D2D(
  FunctorPtr<SuperF3D<T>>&& f,
  const Vector<T,3>& origin, const Vector<T,3>& u, const Vector<T,3>& v,
  int resolution, BlockDataSyncMode syncMode)
  : BlockReduction3D2D(
      std::forward<decltype(f)>(f),
      Hyperplane3D<T>().originAt(origin).spannedBy(u, v),
      resolution, syncMode)
{ }

template <typename T>
BlockReduction3D2D<T>::BlockReduction3D2D(
  FunctorPtr<SuperF3D<T>>&& f,
  const Vector<T,3>& origin, const Vector<T,3>& normal,
  int resolution, BlockDataSyncMode syncMode)
  : BlockReduction3D2D(
      std::forward<decltype(f)>(f),
      Hyperplane3D<T>().originAt(origin).normalTo(normal),
      resolution, syncMode)
{ }

template <typename T>
BlockReduction3D2D<T>::BlockReduction3D2D(
  FunctorPtr<SuperF3D<T>>&& f,
  const Vector<T,3>& normal,
  int resolution, BlockDataSyncMode syncMode)
  : BlockReduction3D2D(
      std::forward<decltype(f)>(f),
      Hyperplane3D<T>()
      .centeredIn(f->getSuperStructure().getCuboidDecomposition().getMotherCuboid())
      .normalTo(normal),
      resolution, syncMode)
{ }

template <typename T>
void BlockReduction3D2D<T>::initialize()
{
  const auto& geometry = _f->getSuperStructure().getCuboidDecomposition();
  LoadBalancer<T>&           load     = _f->getSuperStructure().getLoadBalancer();

  _rankLocalSubplane.clear();

  for ( int iX = 0; iX < this->getNx(); ++iX ) {
    for ( int iY = 0; iY < this->getNy(); ++iY ) {
      const Vector<T,3> physR = this->getPhysR(iX, iY);

      // Schedule plane point for storage if its physical position intersects the
      // mother cuboid and the cuboid of the nearest lattice position is local to
      // the current rank:
      if (auto iC = geometry.getC(physR)) {
        if ( load.isLocal(*iC) ) {
          _rankLocalSubplane.emplace_back(iX, iY, *iC);
        }
      }
    }
  }
}

template <typename T>
void BlockReduction3D2D<T>::update()
{
  _f->getSuperStructure().communicate();

#ifdef PARALLEL_MODE_MPI
  std::unique_ptr<BlockData<2,T,T>> localBlockData(
    new BlockData<2,T,T>({{this->getNx(), this->getNy()}, 0}, _f->getTargetDim()));

  switch ( _reductionMode ) {
  case BlockDataReductionMode::Analytical:
    updateBlockAnalytical(*localBlockData);
    break;
  case BlockDataReductionMode::Discrete:
    updateBlockDiscrete(*localBlockData);
    break;
  }

  switch ( _syncMode ) {
  case BlockDataSyncMode::ReduceAndBcast:
    singleton::mpi().reduce(*localBlockData, this->getBlockData(), MPI_SUM);
    singleton::mpi().bCast(this->getBlockData());
    break;
  case BlockDataSyncMode::ReduceOnly:
    singleton::mpi().reduce(*localBlockData, this->getBlockData(), MPI_SUM);
    break;
  case BlockDataSyncMode::None:
    if (this->_owning) {
      delete this->_blockData;
    }
    this->_blockData = localBlockData.release();
    this->_owning = true;
    break;
  }
#else
  switch ( _reductionMode ) {
  case BlockDataReductionMode::Analytical:
    updateBlockAnalytical(this->getBlockData());
    break;
  case BlockDataReductionMode::Discrete:
    updateBlockDiscrete(this->getBlockData());
    break;
  }
#endif
}

template <typename T>
BlockStructureD<2>& BlockReduction3D2D<T>::getBlockStructure()
{
  return *this->_blockData;
}

template <typename T>
const std::vector<std::tuple<int,int,int>>& BlockReduction3D2D<T>::getRankLocalSubplane() const
{
  return this->_rankLocalSubplane;
}


} // end namespace olb

#endif

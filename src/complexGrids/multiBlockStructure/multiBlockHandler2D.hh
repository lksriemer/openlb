/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt and Bernd Stahl
 *  Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
 *  E-mail: jonas.latt@gmail.com
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

/** \file
 * Handler for 2D multiblock structure -- generic implementation.
 */

#ifndef MULTI_BLOCK_HANDLER_2D_HH
#define MULTI_BLOCK_HANDLER_2D_HH

#include "complexGrids/mpiManager/mpiManager.h"
#include "multiBlockHandler2D.h"
#include "parallelDynamics.h"
#include "core/util.h"
#include <algorithm>
#include <numeric>


namespace olb {

////////////////////// Class SerialMultiBlockHandler2D /////////////////////

template<typename T, template<typename U> class Lattice>
SerialMultiBlockHandler2D<T,Lattice>::SerialMultiBlockHandler2D (
        MultiDataDistribution2D const& dataDistribution_ )
    : dataDistribution(dataDistribution_)
{ }

template<typename T, template<typename U> class Lattice>
int SerialMultiBlockHandler2D<T,Lattice>::getNx() const {
    return dataDistribution.getNx();
}
template<typename T, template<typename U> class Lattice>
int SerialMultiBlockHandler2D<T,Lattice>::getNy() const {
    return dataDistribution.getNy();
}

template<typename T, template<typename U> class Lattice>
MultiDataDistribution2D const& SerialMultiBlockHandler2D<T,Lattice>::getMultiDataDistribution() const {
    return dataDistribution;
}

template<typename T, template<typename U> class Lattice>
bool SerialMultiBlockHandler2D<T,Lattice>::getLocalEnvelope(int iBlock, int& lx, int& ly) const {
    BlockParameters2D const& parameters = dataDistribution.getBlockParameters(iBlock);
    lx = parameters.getEnvelopeLx();
    ly = parameters.getEnvelopeLy();
    return true;
}

template<typename T, template<typename U> class Lattice>
T SerialMultiBlockHandler2D<T,Lattice>::reduceSum (T localSum) const {
    return localSum;
}

template<typename T, template<typename U> class Lattice>
T SerialMultiBlockHandler2D<T,Lattice>::reduceAverage (T localAverage, T localWeight) const {
    return localAverage;
}

template<typename T, template<typename U> class Lattice>
T SerialMultiBlockHandler2D<T,Lattice>::reduceMin (T localMin) const {
    return localMin;
}

template<typename T, template<typename U> class Lattice>
T SerialMultiBlockHandler2D<T,Lattice>::reduceMax (T localMax) const {
    return localMax;
}

template<typename T, template<typename U> class Lattice>
void SerialMultiBlockHandler2D<T,Lattice>::broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialMultiBlockHandler2D<T,Lattice>::broadCastScalar(T& scalar, int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialMultiBlockHandler2D<T,Lattice>::broadCastVector(T vect[Lattice<T>::d], int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialMultiBlockHandler2D<T,Lattice>::copyOverlap (
        Overlap2D const& overlap, SerialMultiBlockHandler2D<T,Lattice>::BlockVector2D& lattices ) const
{
    int originalId = overlap.getOriginalId();
    int overlapId  = overlap.getOverlapId();
    BlockParameters2D const& originalParameters = dataDistribution.getBlockParameters(originalId);
    BlockParameters2D const& overlapParameters  = dataDistribution.getBlockParameters(overlapId);

    BlockCoordinates2D originalCoords(originalParameters.toLocal(overlap.getOriginalCoordinates()));
    BlockCoordinates2D overlapCoords(overlapParameters.toLocal(overlap.getOverlapCoordinates()));

    OLB_PRECONDITION(originalCoords.x1-originalCoords.x0 == overlapCoords.x1-overlapCoords.x0);
    OLB_PRECONDITION(originalCoords.y1-originalCoords.y0 == overlapCoords.y1-overlapCoords.y0);

    int origX = originalCoords.x0;
    int overlapX = overlapCoords.x0;
    for (; origX<=originalCoords.x1; ++origX, ++overlapX) {
        int origY = originalCoords.y0;
        int overlapY = overlapCoords.y0;
        for (; origY<=originalCoords.y1; ++origY, ++overlapY) {
            lattices[overlapId] -> get(overlapX, overlapY).attributeValues (
                    lattices[originalId] -> get(origX, origY) );
        }
    }
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>& SerialMultiBlockHandler2D<T,Lattice>::
                     getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const
{
    OLB_PRECONDITION( baseCell );
    return *baseCell;
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice> const& SerialMultiBlockHandler2D<T,Lattice>::
                           getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const
{
    OLB_PRECONDITION( baseCell );
    return *baseCell;
}

#ifdef PARALLEL_MODE_MPI

////////////////////// Class ParallelMultiBlockHandler2D /////////////////////

template<typename T, template<typename U> class Lattice>
ParallelMultiBlockHandler2D<T,Lattice>::ParallelMultiBlockHandler2D (
        MultiDataDistribution2D const& dataDistribution_ )
    : dataDistribution(dataDistribution_),
      parallelDynamics( 0 )
{ }

template<typename T, template<typename U> class Lattice>
int ParallelMultiBlockHandler2D<T,Lattice>::getNx() const {
    return dataDistribution.getNx();
}
template<typename T, template<typename U> class Lattice>
int ParallelMultiBlockHandler2D<T,Lattice>::getNy() const {
    return dataDistribution.getNy();
}

template<typename T, template<typename U> class Lattice>
MultiDataDistribution2D const& ParallelMultiBlockHandler2D<T,Lattice>::getMultiDataDistribution() const {
    return dataDistribution;
}

template<typename T, template<typename U> class Lattice>
bool ParallelMultiBlockHandler2D<T,Lattice>::getLocalEnvelope(int iBlock, int& lx, int& ly) const {
    BlockParameters2D const& parameters = dataDistribution.getBlockParameters(iBlock);
    if ( parameters.getProcId() == singleton::mpi().getRank() ) {
        lx = parameters.getEnvelopeLx();
        ly = parameters.getEnvelopeLy();
        return true;
    }
    else {
        lx = ly = 0;
        return false;
    }
}

template<typename T, template<typename U> class Lattice>
T ParallelMultiBlockHandler2D<T,Lattice>::reduceSum(T localSum) const {
    T globalSum;
    singleton::mpi().reduce(localSum, globalSum, MPI_SUM);
    singleton::mpi().bCast(&globalSum, 1);
    return globalSum;
}

template<typename T, template<typename U> class Lattice>
T ParallelMultiBlockHandler2D<T,Lattice>::reduceAverage(T localAverage, T localWeight) const {
    T sumAverage, sumWeights;
    singleton::mpi().reduce(localAverage*localWeight, sumAverage, MPI_SUM);
    singleton::mpi().reduce(localWeight, sumWeights, MPI_SUM);
    if (singleton::mpi().isMainProcessor() && sumWeights>1.e-12) {
        sumAverage /= sumWeights;
    }
    singleton::mpi().bCast(&sumAverage, 1);
    return sumAverage;
}

template<typename T, template<typename U> class Lattice>
T ParallelMultiBlockHandler2D<T,Lattice>::reduceMin(T localMin) const {
    T globalMin;
    singleton::mpi().reduce(localMin, globalMin, MPI_MIN);
    singleton::mpi().bCast(&globalMin, 1);
    return globalMin;
}

template<typename T, template<typename U> class Lattice>
T ParallelMultiBlockHandler2D<T,Lattice>::reduceMax(T localMax) const {
    T globalMax;
    singleton::mpi().reduce(localMax, globalMax, MPI_MAX);
    singleton::mpi().bCast(&globalMax, 1);
    return globalMax;
}

template<typename T, template<typename U> class Lattice>
void ParallelMultiBlockHandler2D<T,Lattice>::broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const {
    const int sizeOfCell = Lattice<T>::q + Lattice<T>::ExternalField::numScalars;
    T* cellData = new T[sizeOfCell];
    int fromProc = dataDistribution.getBlockParameters(fromBlock).getProcId();
    if (singleton::mpi().getRank()==fromProc) {
        cell.serialize(cellData);
    }
    singleton::mpi().bCast(cellData, sizeOfCell, fromProc);
    cell.unSerialize(cellData);
    delete [] cellData;
}

template<typename T, template<typename U> class Lattice>
void ParallelMultiBlockHandler2D<T,Lattice>::broadCastScalar(T& scalar, int fromBlock) const {
    int fromProc = dataDistribution.getBlockParameters(fromBlock).getProcId();
    singleton::mpi().bCast(&scalar, 1, fromProc);
}

template<typename T, template<typename U> class Lattice>
void ParallelMultiBlockHandler2D<T,Lattice>::broadCastVector(T vect[Lattice<T>::d], int fromBlock) const {
    int fromProc = dataDistribution.getBlockParameters(fromBlock).getProcId();
    singleton::mpi().bCast(vect, Lattice<T>::d, fromProc);
}

template<typename T, template<typename U> class Lattice>
void ParallelMultiBlockHandler2D<T,Lattice>::copyOverlap (
        Overlap2D const& overlap, ParallelMultiBlockHandler2D<T,Lattice>::BlockVector2D& lattices ) const
{
    int originalId = overlap.getOriginalId();
    int overlapId  = overlap.getOverlapId();
    BlockParameters2D const& originalParameters = dataDistribution.getBlockParameters(originalId);
    BlockParameters2D const& overlapParameters  = dataDistribution.getBlockParameters(overlapId);
    int originalProc = originalParameters.getProcId();
    int overlapProc  = overlapParameters.getProcId();

    BlockCoordinates2D originalCoords(originalParameters.toLocal(overlap.getOriginalCoordinates()));
    BlockCoordinates2D overlapCoords(overlapParameters.toLocal(overlap.getOverlapCoordinates()));

    const int sizeOfCell = Lattice<T>::q + Lattice<T>::ExternalField::numScalars;
    int lx = originalCoords.x1-originalCoords.x0+1;
    int ly = originalCoords.y1-originalCoords.y0+1;
    OLB_PRECONDITION(lx == overlapCoords.x1-overlapCoords.x0+1);
    OLB_PRECONDITION(ly == overlapCoords.y1-overlapCoords.y0+1);

    // Case 1: both the original data and the region of overlap are on my processor
    if (originalProc==singleton::mpi().getRank() && overlapProc==singleton::mpi().getRank()) {
        int origX = originalCoords.x0;
        int overlapX = overlapCoords.x0;
        for (; origX<=originalCoords.x1; ++origX, ++overlapX) {
            int origY = originalCoords.y0;
            int overlapY = overlapCoords.y0;
            for (; origY<=originalCoords.y1; ++origY, ++overlapY) {
                lattices[overlapId] -> get(overlapX, overlapY).attributeValues (
                        lattices[originalId] -> get(origX, origY) );
            }
        }
    }
    // Case 2: only the original data is on my processor (I am sender)
    else if(originalProc==singleton::mpi().getRank()) {
        T* sendData = new T[lx*ly*sizeOfCell];
        int iData=0;
        for (int iX=originalCoords.x0; iX<=originalCoords.x1; ++iX) {
            for (int iY=originalCoords.y0; iY<=originalCoords.y1; ++iY) {
                lattices[originalId] -> get(iX,iY).serialize(sendData+iData);
                iData += sizeOfCell;
            }
        }
        singleton::mpi().send(sendData, lx*ly*sizeOfCell, overlapProc);
        delete [] sendData;
    }
    // Case 3: only overlap area is on my processor (I am receiver)
    else if (overlapProc==singleton::mpi().getRank()) {
        T* recvData = new T[lx*ly*sizeOfCell];
        singleton::mpi().receive(recvData, lx*ly*sizeOfCell, originalProc);
        int iData=0;
        for (int iX=overlapCoords.x0; iX<=overlapCoords.x1; ++iX) {
            for (int iY=overlapCoords.y0; iY<=overlapCoords.y1; ++iY) {
                lattices[overlapId] -> get(iX,iY).unSerialize(recvData+iData);
                iData += sizeOfCell;
            }
        }
        delete [] recvData;
    }
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>& ParallelMultiBlockHandler2D<T,Lattice>::
                     getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const
{
    int onProcessor = dataDistribution.getBlockParameters(onBlock).getProcId();
    delete parallelDynamics;
    parallelDynamics = new ParallelDynamics<T,Lattice>(baseCell, onProcessor);
    distributedCell.defineDynamics(parallelDynamics);
    return distributedCell;
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice> const& ParallelMultiBlockHandler2D<T,Lattice>::
                     getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const
{
    int onProcessor = dataDistribution.getBlockParameters(onBlock).getProcId();
    delete parallelDynamics;
    parallelDynamics = new ConstParallelDynamics<T,Lattice>(baseCell, onProcessor);
    distributedCell.defineDynamics(parallelDynamics);
    return distributedCell;
}

#endif  // PARALLEL_MODE_MPI

}  // namespace olb

#endif

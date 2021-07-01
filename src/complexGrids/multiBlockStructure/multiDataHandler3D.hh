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
 * Handler for 3D multiblock structure -- generic implementation.
 */

#ifndef MULTI_DATA_HANDLER_3D_HH
#define MULTI_DATA_HANDLER_3D_HH

#include "complexGrids/mpiManager/mpiManager.h"
#include "multiDataHandler3D.h"
#include "parallelDynamics.h"
#include "core/util.h"
#include <algorithm>
#include <numeric>


namespace olb {

////////////////////// Class SerialDataHandler3D /////////////////////

template<typename T, template<typename U> class Lattice>
SerialDataHandler3D<T,Lattice>::SerialDataHandler3D (
        MultiDataDistribution3D const& dataDistribution_ )
    : dataDistribution(dataDistribution_)
{ }

template<typename T, template<typename U> class Lattice>
int SerialDataHandler3D<T,Lattice>::getNx() const {
    return dataDistribution.getNx();
}
template<typename T, template<typename U> class Lattice>
int SerialDataHandler3D<T,Lattice>::getNy() const {
    return dataDistribution.getNy();
}
template<typename T, template<typename U> class Lattice>
int SerialDataHandler3D<T,Lattice>::getNz() const {
    return dataDistribution.getNz();
}

template<typename T, template<typename U> class Lattice>
MultiDataDistribution3D const& SerialDataHandler3D<T,Lattice>::getMultiDataDistribution() const {
    return dataDistribution;
}

template<typename T, template<typename U> class Lattice>
void SerialDataHandler3D<T,Lattice>::allocateBlockLattices (
        SerialDataHandler3D<T,Lattice>::BlockVector3D& lattices) const
{
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
        lattices.push_back(new BlockLattice3D<T,Lattice>(parameters.getEnvelopeLx(),
                                                         parameters.getEnvelopeLy(),
                                                         parameters.getEnvelopeLz()));
    }
}

template<typename T, template<typename U> class Lattice>
void SerialDataHandler3D<T,Lattice>::allocateVectorFields (
        SerialDataHandler3D<T,Lattice>::VectorFieldVector3D& fields) const
{
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
        fields.push_back(new TensorField3D<T,3>(parameters.getEnvelopeLx(),
                                                parameters.getEnvelopeLy(),
                                                parameters.getEnvelopeLz()));
    }
}

template<typename T, template<typename U> class Lattice>
void SerialDataHandler3D<T,Lattice>::allocateSymmetricMatrixFields (
        SerialDataHandler3D<T,Lattice>::SymmetricMatrixFieldVector3D& fields) const
{
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
        fields.push_back(new TensorField3D<T,6>(parameters.getEnvelopeLx(),
                                                parameters.getEnvelopeLy(),
                                                parameters.getEnvelopeLz()));
    }
}

template<typename T, template<typename U> class Lattice>
void SerialDataHandler3D<T,Lattice>::allocateScalarFields (
        SerialDataHandler3D<T,Lattice>::ScalarFieldVector3D& fields) const
{
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
        fields.push_back(new ScalarField3D<T>(parameters.getEnvelopeLx(),
                                              parameters.getEnvelopeLy(),
                                              parameters.getEnvelopeLz()));
    }
}

template<typename T, template<typename U> class Lattice>
T SerialDataHandler3D<T,Lattice>::reduceSum (T localSum) const {
    return localSum;
}

template<typename T, template<typename U> class Lattice>
T SerialDataHandler3D<T,Lattice>::reduceAverage (T localAverage, T localWeight) const {
    return localAverage;
}

template<typename T, template<typename U> class Lattice>
T SerialDataHandler3D<T,Lattice>::reduceMin (T localMin) const {
    return localMin;
}

template<typename T, template<typename U> class Lattice>
T SerialDataHandler3D<T,Lattice>::reduceMax (T localMax) const {
    return localMax;
}

template<typename T, template<typename U> class Lattice>
void SerialDataHandler3D<T,Lattice>::broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialDataHandler3D<T,Lattice>::broadCastScalar(T& scalar, int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialDataHandler3D<T,Lattice>::broadCastVector(T vect[Lattice<T>::d], int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialDataHandler3D<T,Lattice>::copyOverlap (
        Overlap3D const& overlap, SerialDataHandler3D<T,Lattice>::BlockVector3D& lattices ) const
{
    int originalId = overlap.getOriginalId();
    int overlapId  = overlap.getOverlapId();
    BlockParameters3D const& originalParameters = dataDistribution.getBlockParameters(originalId);
    BlockParameters3D const& overlapParameters  = dataDistribution.getBlockParameters(overlapId);

    BlockCoordinates3D originalCoords(originalParameters.toLocal(overlap.getOriginalCoordinates()));
    BlockCoordinates3D overlapCoords(overlapParameters.toLocal(overlap.getOverlapCoordinates()));

    assert(originalCoords.x1-originalCoords.x0 == overlapCoords.x1-overlapCoords.x0);
    assert(originalCoords.y1-originalCoords.y0 == overlapCoords.y1-overlapCoords.y0);
    assert(originalCoords.z1-originalCoords.z0 == overlapCoords.z1-overlapCoords.z0);

    int origX = originalCoords.x0;
    int overlapX = overlapCoords.x0;
    for (; origX<=originalCoords.x1; ++origX, ++overlapX) {
        int origY = originalCoords.y0;
        int overlapY = overlapCoords.y0;
        for (; origY<=originalCoords.y1; ++origY, ++overlapY) {
            int origZ = originalCoords.z0;
            int overlapZ = overlapCoords.z0;
            for (; origZ<=originalCoords.z1; ++origZ, ++overlapZ) {
                lattices[overlapId] -> get(overlapX, overlapY, overlapZ).attributeValues (
                        lattices[originalId] -> get(origX, origY, origZ) );
            }
        }
    }
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>& SerialDataHandler3D<T,Lattice>::
                     getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const
{
    OLB_PRECONDITION( baseCell );
    return *baseCell;
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice> const& SerialDataHandler3D<T,Lattice>::
                           getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const
{
    OLB_PRECONDITION( baseCell );
    return *baseCell;
}


#ifdef PARALLEL_MODE_MPI

////////////////////// Class ParallelDataHandler3D /////////////////////

template<typename T, template<typename U> class Lattice>
ParallelDataHandler3D<T,Lattice>::ParallelDataHandler3D (
        MultiDataDistribution3D const& dataDistribution_ )
    : dataDistribution(dataDistribution_),
      parallelDynamics( 0 )
{ }

template<typename T, template<typename U> class Lattice>
int ParallelDataHandler3D<T,Lattice>::getNx() const {
    return dataDistribution.getNx();
}

template<typename T, template<typename U> class Lattice>
int ParallelDataHandler3D<T,Lattice>::getNy() const {
    return dataDistribution.getNy();
}

template<typename T, template<typename U> class Lattice>
int ParallelDataHandler3D<T,Lattice>::getNz() const {
    return dataDistribution.getNz();
}

template<typename T, template<typename U> class Lattice>
MultiDataDistribution3D const& ParallelDataHandler3D<T,Lattice>::getMultiDataDistribution() const {
    return dataDistribution;
}

template<typename T, template<typename U> class Lattice>
void ParallelDataHandler3D<T,Lattice>::allocateBlockLattices (
        ParallelDataHandler3D<T,Lattice>::BlockVector3D& lattices) const
{
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
        if ( parameters.getProcId() == singleton::mpi().getRank() ) {
            lattices.push_back(new BlockLattice3D<T,Lattice>(parameters.getEnvelopeLx(),
                                                             parameters.getEnvelopeLy(),
                                                             parameters.getEnvelopeLz()));
        }
        else {
            lattices.push_back(0);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDataHandler3D<T,Lattice>::allocateVectorFields (
        ParallelDataHandler3D<T,Lattice>::VectorFieldVector3D& fields ) const
{
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
        if ( parameters.getProcId() == singleton::mpi().getRank() ) {
            fields.push_back(new TensorField3D<T,3>(parameters.getEnvelopeLx(),
                                                    parameters.getEnvelopeLy(),
                                                    parameters.getEnvelopeLz()));
        }
        else {
            fields.push_back( 0 );
        }
    }

}

template<typename T, template<typename U> class Lattice>
void ParallelDataHandler3D<T,Lattice>::allocateSymmetricMatrixFields (
        ParallelDataHandler3D<T,Lattice>::SymmetricMatrixFieldVector3D& fields) const
{
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
        if ( parameters.getProcId() == singleton::mpi().getRank() ) {
            fields.push_back(new TensorField3D<T,6>(parameters.getEnvelopeLx(),
                                                    parameters.getEnvelopeLy(),
                                                    parameters.getEnvelopeLz()));
        }
        else {
            fields.push_back( 0 );
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDataHandler3D<T,Lattice>::allocateScalarFields (
        ParallelDataHandler3D<T,Lattice>::ScalarFieldVector3D& fields) const
{
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
        if ( parameters.getProcId() == singleton::mpi().getRank() ) {
            fields.push_back(new ScalarField3D<T>(parameters.getEnvelopeLx(),
                                                  parameters.getEnvelopeLy(),
                                                  parameters.getEnvelopeLz()));
        }
        else {
            fields.push_back( 0 );
        }
    }
}

template<typename T, template<typename U> class Lattice>
T ParallelDataHandler3D<T,Lattice>::reduceSum(T localSum) const {
    T globalSum;
    singleton::mpi().reduce(localSum, globalSum, MPI_SUM);
    singleton::mpi().bCast(&globalSum, 1);
    return globalSum;
}

template<typename T, template<typename U> class Lattice>
T ParallelDataHandler3D<T,Lattice>::reduceAverage(T localAverage, T localWeight) const {
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
T ParallelDataHandler3D<T,Lattice>::reduceMin(T localMin) const {
    T globalMin;
    singleton::mpi().reduce(localMin, globalMin, MPI_MIN);
    singleton::mpi().bCast(&globalMin, 1);
    return globalMin;
}

template<typename T, template<typename U> class Lattice>
T ParallelDataHandler3D<T,Lattice>::reduceMax(T localMax) const {
    T globalMax;
    singleton::mpi().reduce(localMax, globalMax, MPI_MAX);
    singleton::mpi().bCast(&globalMax, 1);
    return globalMax;
}

template<typename T, template<typename U> class Lattice>
void ParallelDataHandler3D<T,Lattice>::broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const {
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
void ParallelDataHandler3D<T,Lattice>::broadCastScalar(T& scalar, int fromBlock) const {
    int fromProc = dataDistribution.getBlockParameters(fromBlock).getProcId();
    singleton::mpi().bCast(&scalar, 1, fromProc);
}

template<typename T, template<typename U> class Lattice>
void ParallelDataHandler3D<T,Lattice>::broadCastVector(T vect[Lattice<T>::d], int fromBlock) const {
    int fromProc = dataDistribution.getBlockParameters(fromBlock).getProcId();
    singleton::mpi().bCast(vect, Lattice<T>::d, fromProc);
}

template<typename T, template<typename U> class Lattice>
void ParallelDataHandler3D<T,Lattice>::copyOverlap (
        Overlap3D const& overlap, ParallelDataHandler3D<T,Lattice>::BlockVector3D& lattices ) const
{
    int originalId = overlap.getOriginalId();
    int overlapId  = overlap.getOverlapId();
    BlockParameters3D const& originalParameters = dataDistribution.getBlockParameters(originalId);
    BlockParameters3D const& overlapParameters  = dataDistribution.getBlockParameters(overlapId);
    int originalProc = originalParameters.getProcId();
    int overlapProc  = overlapParameters.getProcId();

    BlockCoordinates3D originalCoords(originalParameters.toLocal(overlap.getOriginalCoordinates()));
    BlockCoordinates3D overlapCoords(overlapParameters.toLocal(overlap.getOverlapCoordinates()));

    const int sizeOfCell = Lattice<T>::q + Lattice<T>::ExternalField::numScalars;
    int lx = originalCoords.x1-originalCoords.x0+1;
    int ly = originalCoords.y1-originalCoords.y0+1;
    int lz = originalCoords.z1-originalCoords.z0+1;
    OLB_PRECONDITION(lx == overlapCoords.x1-overlapCoords.x0+1);
    OLB_PRECONDITION(ly == overlapCoords.y1-overlapCoords.y0+1);
    OLB_PRECONDITION(lz == overlapCoords.z1-overlapCoords.z0+1);

    // Case 1: both the original data and the region of overlap are on my processor
    if (originalProc==singleton::mpi().getRank() && overlapProc==singleton::mpi().getRank()) {
        int origX = originalCoords.x0;
        int overlapX = overlapCoords.x0;
        for (; origX<=originalCoords.x1; ++origX, ++overlapX) {
            int origY = originalCoords.y0;
            int overlapY = overlapCoords.y0;
            for (; origY<=originalCoords.y1; ++origY, ++overlapY) {
                int origZ = originalCoords.z0;
                int overlapZ = overlapCoords.z0;
                for (; origZ<=originalCoords.z1; ++origZ, ++overlapZ) {
                    lattices[overlapId] -> get(overlapX, overlapY, overlapZ).attributeValues (
                            lattices[originalId] -> get(origX, origY, origZ) );
                }
            }
        }
    }
    // Case 2: only the original data is on my processor (I am sender)
    else if(originalProc==singleton::mpi().getRank()) {
        T* sendData = new T[lx*ly*lz*sizeOfCell];
        int iData=0;
        for (int iX=originalCoords.x0; iX<=originalCoords.x1; ++iX) {
            for (int iY=originalCoords.y0; iY<=originalCoords.y1; ++iY) {
                for (int iZ=originalCoords.z0; iZ<=originalCoords.z1; ++iZ) {
                    lattices[originalId] -> get(iX,iY,iZ).serialize(sendData+iData);
                    iData += sizeOfCell;
                }
            }
        }
        singleton::mpi().send(sendData, lx*ly*lz*sizeOfCell, overlapProc);
        delete [] sendData;
    }
    // Case 3: only overlap area is on my processor (I am receiver)
    else if (overlapProc==singleton::mpi().getRank()) {
        T* recvData = new T[lx*ly*lz*sizeOfCell];
        singleton::mpi().receive(recvData, lx*ly*lz*sizeOfCell, originalProc);
        int iData=0;
        for (int iX=overlapCoords.x0; iX<=overlapCoords.x1; ++iX) {
            for (int iY=overlapCoords.y0; iY<=overlapCoords.y1; ++iY) {
                for (int iZ=overlapCoords.z0; iZ<=overlapCoords.z1; ++iZ) {
                    lattices[overlapId] -> get(iX,iY,iZ).unSerialize(recvData+iData);
                    iData += sizeOfCell;
                }
            }
        }
        delete [] recvData;
    }
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>& ParallelDataHandler3D<T,Lattice>::
                     getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const
{
    int onProcessor = dataDistribution.getBlockParameters(onBlock).getProcId();
    delete parallelDynamics;
    parallelDynamics = new ParallelDynamics<T,Lattice>(baseCell, onProcessor);
    distributedCell.defineDynamics(parallelDynamics);
    return distributedCell;
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice> const& ParallelDataHandler3D<T,Lattice>::
                     getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const
{
    int onProcessor = dataDistribution.getBlockParameters(onBlock).getProcId();
    delete parallelDynamics;
    parallelDynamics = new ConstParallelDynamics<T,Lattice>(baseCell, onProcessor);
    distributedCell.defineDynamics(parallelDynamics);
    return distributedCell;
}

#endif

}  // namespace olb

#endif

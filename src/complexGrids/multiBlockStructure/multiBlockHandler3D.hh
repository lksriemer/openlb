/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2008 Jonas Latt and Bernd Stahl
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

#ifndef MULTI_BLOCK_HANDLER_3D_HH
#define MULTI_BLOCK_HANDLER_3D_HH

#include "complexGrids/mpiManager/mpiManager.h"
#include "multiBlockHandler3D.h"
#include "parallelDynamics.h"
#include "core/util.h"
#include <algorithm>
#include <numeric>


namespace olb {

////////////////////// Class SerialMultiBlockHandler3D /////////////////////

template<typename T, template<typename U> class Lattice>
SerialMultiBlockHandler3D<T,Lattice>::SerialMultiBlockHandler3D (
        MultiDataDistribution3D const& dataDistribution_ )
    : dataDistribution(dataDistribution_)
{ }

template<typename T, template<typename U> class Lattice>
int SerialMultiBlockHandler3D<T,Lattice>::getNx() const {
    return dataDistribution.getNx();
}

template<typename T, template<typename U> class Lattice>
int SerialMultiBlockHandler3D<T,Lattice>::getNy() const {
    return dataDistribution.getNy();
}

template<typename T, template<typename U> class Lattice>
int SerialMultiBlockHandler3D<T,Lattice>::getNz() const {
    return dataDistribution.getNz();
}

template<typename T, template<typename U> class Lattice>
MultiDataDistribution3D const& SerialMultiBlockHandler3D<T,Lattice>::getMultiDataDistribution() const {
    return dataDistribution;
}

template<typename T, template<typename U> class Lattice>
bool SerialMultiBlockHandler3D<T,Lattice>::getLocalEnvelope(int iBlock, int& lx, int& ly, int& lz) const {
    BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
    lx = parameters.getEnvelopeLx();
    ly = parameters.getEnvelopeLy();
    lz = parameters.getEnvelopeLz();
    return true;
}

template<typename T, template<typename U> class Lattice>
T SerialMultiBlockHandler3D<T,Lattice>::reduceSum (T localSum) const {
    return localSum;
}

template<typename T, template<typename U> class Lattice>
T SerialMultiBlockHandler3D<T,Lattice>::reduceAverage (T localAverage, T localWeight) const {
    return localAverage;
}

template<typename T, template<typename U> class Lattice>
T SerialMultiBlockHandler3D<T,Lattice>::reduceMin (T localMin) const {
    return localMin;
}

template<typename T, template<typename U> class Lattice>
T SerialMultiBlockHandler3D<T,Lattice>::reduceMax (T localMax) const {
    return localMax;
}

template<typename T, template<typename U> class Lattice>
void SerialMultiBlockHandler3D<T,Lattice>::broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialMultiBlockHandler3D<T,Lattice>::broadCastScalar(T& scalar, int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialMultiBlockHandler3D<T,Lattice>::broadCastVector(T vect[Lattice<T>::d], int fromBlock) const {
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Lattice>
void SerialMultiBlockHandler3D<T,Lattice>::copyOverlap (
        Overlap3D const& overlap, SerialMultiBlockHandler3D<T,Lattice>::BlockVector3D& lattices ) const
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
void SerialMultiBlockHandler3D<T,Lattice>::connectBoundaries (
        SerialMultiBlockHandler3D<T,Lattice>::BlockVector3D& lattices, bool periodicCommunication ) const
{
    for (int iOverlap=0; iOverlap<dataDistribution.getNumNormalOverlaps(); ++iOverlap) {
        copyOverlap(dataDistribution.getNormalOverlap(iOverlap), lattices);
    }
    if (periodicCommunication) {
        for (int iOverlap=0; iOverlap<dataDistribution.getNumPeriodicOverlaps(); ++iOverlap) {
            copyOverlap(dataDistribution.getPeriodicOverlap(iOverlap), lattices);
        }
    }
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>& SerialMultiBlockHandler3D<T,Lattice>::
    getDistributedCell(std::vector<Cell<T,Lattice>*>& baseCell, bool hasBulkCell) const
{
    OLB_PRECONDITION( baseCell.size()>0 && baseCell[0] );
    return *baseCell[0];
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice> const& SerialMultiBlockHandler3D<T,Lattice>::
    getDistributedCell(std::vector<Cell<T,Lattice> const*>& baseCell, bool hasBulkCell) const
{
    OLB_PRECONDITION( baseCell.size()>0 && baseCell[0] );
    return *baseCell[0];
}

template<typename T, template<typename U> class Lattice>
int SerialMultiBlockHandler3D<T,Lattice>::locateLocally(int iX, int iY, int iZ,
                                                        std::vector<int>& foundId,
                                                        bool& hasBulkCell, int guess) const
{
    hasBulkCell = true;
    return dataDistribution.locateInEnvelopes(iX,iY,iZ, foundId, guess);
}


#ifdef PARALLEL_MODE_MPI

////////////////////// Class DataTransmittor3D /////////////////////

template<typename T, template<typename U> class Lattice>
DataTransmittor3D<T,Lattice>::DataTransmittor3D (
        Overlap3D const& overlap, MultiDataDistribution3D const& dataDistribution )
{
    originalId = overlap.getOriginalId();
    overlapId  = overlap.getOverlapId();

    BlockParameters3D const& originalParameters = dataDistribution.getBlockParameters(originalId);
    BlockParameters3D const& overlapParameters  = dataDistribution.getBlockParameters(overlapId);

    originalProc = originalParameters.getProcId();
    overlapProc  = overlapParameters.getProcId();
    int myProc   = singleton::mpi().getRank();

    originalCoords = originalParameters.toLocal(overlap.getOriginalCoordinates());
    overlapCoords  = overlapParameters.toLocal(overlap.getOverlapCoordinates());

    int lx = originalCoords.x1-originalCoords.x0+1;
    int ly = originalCoords.y1-originalCoords.y0+1;
    int lz = originalCoords.z1-originalCoords.z0+1;
    OLB_PRECONDITION(lx == overlapCoords.x1-overlapCoords.x0+1);
    OLB_PRECONDITION(ly == overlapCoords.y1-overlapCoords.y0+1);
    OLB_PRECONDITION(lz == overlapCoords.z1-overlapCoords.z0+1);

    sizeOfCell = Lattice<T>::q + Lattice<T>::ExternalField::numScalars;
    bufferSize = lx*ly*lz*sizeOfCell;

    buffer.resize(bufferSize);

    if (originalProc==myProc && overlapProc==myProc) {
        myRole = senderAndReceiver;
    }
    else if (originalProc==myProc) {
        myRole = sender;
    }
    else if (overlapProc==myProc) {
        myRole = receiver;
    }
    else {
        myRole = nothing;
    }
}

template<typename T, template<typename U> class Lattice>
void DataTransmittor3D<T,Lattice>::prepareTransmission(BlockVector3D& lattices) {
    if (myRole==sender) {
        T* bufferP = &buffer[0];
        BlockLattice3D<T,Lattice>* lattice = lattices[originalId];
        int iData=0;
        for (int iX=originalCoords.x0; iX<=originalCoords.x1; ++iX) {
            for (int iY=originalCoords.y0; iY<=originalCoords.y1; ++iY) {
                for (int iZ=originalCoords.z0; iZ<=originalCoords.z1; ++iZ) {
                    lattice -> get(iX,iY,iZ).serialize(bufferP+iData);
                    iData += sizeOfCell;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Lattice>
void DataTransmittor3D<T,Lattice>::executeTransmission(BlockVector3D& lattices) {
    switch(myRole) {
        case senderAndReceiver:
        {
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
            break;
        }
        case sender:
            singleton::mpi().send(&buffer[0], bufferSize, overlapProc);
            break;
        case receiver:
            singleton::mpi().receive(&buffer[0], bufferSize, originalProc);
            break;
        case nothing:
            break;
    }
}

template<typename T, template<typename U> class Lattice>
void DataTransmittor3D<T,Lattice>::finalizeTransmission(BlockVector3D& lattices) {
    if (myRole==receiver) {
        T* bufferP = &buffer[0];
        BlockLattice3D<T,Lattice>* lattice = lattices[overlapId];
        int iData=0;
        for (int iX=overlapCoords.x0; iX<=overlapCoords.x1; ++iX) {
            for (int iY=overlapCoords.y0; iY<=overlapCoords.y1; ++iY) {
                for (int iZ=overlapCoords.z0; iZ<=overlapCoords.z1; ++iZ) {
                    lattice -> get(iX,iY,iZ).unSerialize(bufferP+iData);
                    iData += sizeOfCell;
                }
            }
        }
    }
}


////////////////////// Class ParallelMultiBlockHandler3D /////////////////////

template<typename T, template<typename U> class Lattice>
ParallelMultiBlockHandler3D<T,Lattice>::ParallelMultiBlockHandler3D (
        MultiDataDistribution3D const& dataDistribution_ )
    : dataDistribution(dataDistribution_),
      parallelDynamics( 0 )
{
    computeLocallyRelevantBlocks();
    normalTransmittors.resize(relevantNormalOverlaps.size());
    for (unsigned iRelevant=0; iRelevant<relevantNormalOverlaps.size(); ++iRelevant) {
        int iOverlap = relevantNormalOverlaps[iRelevant];
        normalTransmittors[iRelevant] =
                new DataTransmittor3D<T,Lattice>(dataDistribution.getNormalOverlap(iOverlap), dataDistribution);
    }
    periodicTransmittors.resize(relevantPeriodicOverlaps.size());
    for (unsigned iRelevant=0; iRelevant<relevantPeriodicOverlaps.size(); ++iRelevant) {
        int iOverlap = relevantPeriodicOverlaps[iRelevant];
        periodicTransmittors[iRelevant] =
                new DataTransmittor3D<T,Lattice>(dataDistribution.getPeriodicOverlap(iOverlap), dataDistribution);
    }
}

template<typename T, template<typename U> class Lattice>
ParallelMultiBlockHandler3D<T,Lattice>::~ParallelMultiBlockHandler3D() {
    for (unsigned iRelevant=0; iRelevant<relevantNormalOverlaps.size(); ++iRelevant) {
        delete normalTransmittors[iRelevant];
    }
    for (unsigned iRelevant=0; iRelevant<relevantPeriodicOverlaps.size(); ++iRelevant) {
        delete periodicTransmittors[iRelevant];
    }
}

template<typename T, template<typename U> class Lattice>
int ParallelMultiBlockHandler3D<T,Lattice>::getNx() const {
    return dataDistribution.getNx();
}

template<typename T, template<typename U> class Lattice>
int ParallelMultiBlockHandler3D<T,Lattice>::getNy() const {
    return dataDistribution.getNy();
}

template<typename T, template<typename U> class Lattice>
int ParallelMultiBlockHandler3D<T,Lattice>::getNz() const {
    return dataDistribution.getNz();
}

template<typename T, template<typename U> class Lattice>
MultiDataDistribution3D const& ParallelMultiBlockHandler3D<T,Lattice>::getMultiDataDistribution() const {
    return dataDistribution;
}

template<typename T, template<typename U> class Lattice>
bool ParallelMultiBlockHandler3D<T,Lattice>::getLocalEnvelope(int iBlock, int& lx, int& ly, int& lz) const {
    BlockParameters3D const& parameters = dataDistribution.getBlockParameters(iBlock);
    lx = parameters.getEnvelopeLx();
    ly = parameters.getEnvelopeLy();
    lz = parameters.getEnvelopeLz();
    if ( parameters.getProcId() == singleton::mpi().getRank() ) {
        return true;
    }
    else {
        return false;
    }
}

template<typename T, template<typename U> class Lattice>
T ParallelMultiBlockHandler3D<T,Lattice>::reduceSum(T localSum) const {
    T globalSum;
    singleton::mpi().reduce(localSum, globalSum, MPI_SUM);
    singleton::mpi().bCast(&globalSum, 1);
    return globalSum;
}

template<typename T, template<typename U> class Lattice>
T ParallelMultiBlockHandler3D<T,Lattice>::reduceAverage(T localAverage, T localWeight) const {
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
T ParallelMultiBlockHandler3D<T,Lattice>::reduceMin(T localMin) const {
    T globalMin;
    singleton::mpi().reduce(localMin, globalMin, MPI_MIN);
    singleton::mpi().bCast(&globalMin, 1);
    return globalMin;
}

template<typename T, template<typename U> class Lattice>
T ParallelMultiBlockHandler3D<T,Lattice>::reduceMax(T localMax) const {
    T globalMax;
    singleton::mpi().reduce(localMax, globalMax, MPI_MAX);
    singleton::mpi().bCast(&globalMax, 1);
    return globalMax;
}

template<typename T, template<typename U> class Lattice>
void ParallelMultiBlockHandler3D<T,Lattice>::broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const {
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
void ParallelMultiBlockHandler3D<T,Lattice>::broadCastScalar(T& scalar, int fromBlock) const {
    int fromProc = dataDistribution.getBlockParameters(fromBlock).getProcId();
    singleton::mpi().bCast(&scalar, 1, fromProc);
}

template<typename T, template<typename U> class Lattice>
void ParallelMultiBlockHandler3D<T,Lattice>::broadCastVector(T vect[Lattice<T>::d], int fromBlock) const {
    int fromProc = dataDistribution.getBlockParameters(fromBlock).getProcId();
    singleton::mpi().bCast(vect, Lattice<T>::d, fromProc);
}

template<typename T, template<typename U> class Lattice>
void ParallelMultiBlockHandler3D<T,Lattice>::connectBoundaries (
        ParallelMultiBlockHandler3D<T,Lattice>::BlockVector3D& lattices, bool periodicCommunication ) const
{
    for (unsigned iRelevant=0; iRelevant<relevantNormalOverlaps.size(); ++iRelevant) {
        normalTransmittors[iRelevant]->prepareTransmission(lattices);
    }
    if (periodicCommunication) {
        for (unsigned iRelevant=0; iRelevant<relevantPeriodicOverlaps.size(); ++iRelevant) {
            periodicTransmittors[iRelevant]->prepareTransmission(lattices);
        }
    }
    for (unsigned iRelevant=0; iRelevant<relevantNormalOverlaps.size(); ++iRelevant) {
        normalTransmittors[iRelevant]->executeTransmission(lattices);
    }
    if (periodicCommunication) {
        for (unsigned iRelevant=0; iRelevant<relevantPeriodicOverlaps.size(); ++iRelevant) {
            periodicTransmittors[iRelevant]->executeTransmission(lattices);
        }
    }
    for (unsigned iRelevant=0; iRelevant<relevantNormalOverlaps.size(); ++iRelevant) {
        normalTransmittors[iRelevant]->finalizeTransmission(lattices);
    }
    if (periodicCommunication) {
        for (unsigned iRelevant=0; iRelevant<relevantPeriodicOverlaps.size(); ++iRelevant) {
            periodicTransmittors[iRelevant]->finalizeTransmission(lattices);
        }
    }
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>& ParallelMultiBlockHandler3D<T,Lattice>::
    getDistributedCell(std::vector<Cell<T,Lattice>*>& baseCell, bool hasBulkCell) const
{
    delete parallelDynamics;
    parallelDynamics = new ParallelDynamics<T,Lattice>(baseCell, hasBulkCell);
    distributedCell.defineDynamics(parallelDynamics);
    return distributedCell;
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice> const& ParallelMultiBlockHandler3D<T,Lattice>::
    getDistributedCell(std::vector<Cell<T,Lattice> const*>& baseCell, bool hasBulkCell) const
{
    delete parallelDynamics;
    parallelDynamics = new ConstParallelDynamics<T,Lattice>(baseCell, hasBulkCell);
    distributedCell.defineDynamics(parallelDynamics);
    return distributedCell;
}

template<typename T, template<typename U> class Lattice>
int ParallelMultiBlockHandler3D<T,Lattice>::locateLocally(int iX, int iY, int iZ, std::vector<int>& foundId,
                                                          bool& hasBulkCell, int guess) const
{
    hasBulkCell = false;
    if (!util::contained(iX,iY,iZ,
                boundingBox.x0, boundingBox.x1, boundingBox.y0, boundingBox.y1,
                boundingBox.z0, boundingBox.z1) )
    {
        return -1;
    }
    for (int iBlock=0; iBlock < (int)myBlocks.size(); ++iBlock) {
        BlockCoordinates3D const& coord = dataDistribution.getBlockParameters(myBlocks[iBlock]).getEnvelope();
        if (util::contained(iX, iY, iZ, coord.x0, coord.x1, coord.y0, coord.y1, coord.z0, coord.z1)) {
            BlockCoordinates3D const& bulk =
                dataDistribution.getBlockParameters(myBlocks[iBlock]).getBulk();
            if (util::contained(iX, iY, iZ, bulk.x0, bulk.x1, bulk.y0, bulk.y1, bulk.z0, bulk.z1)) {
                hasBulkCell = true;
                foundId.insert(foundId.begin(),myBlocks[iBlock]);
            }
            foundId.push_back(myBlocks[iBlock]);
        }
    }
    return -1;
}


template<typename T, template<typename U> class Lattice>
void ParallelMultiBlockHandler3D<T,Lattice>::computeLocallyRelevantBlocks() {
    for (int iBlock=0; iBlock<dataDistribution.getNumBlocks(); ++iBlock) {
        if (dataDistribution.getBlockParameters(iBlock).getProcId() == singleton::mpi().getRank()) {
            BlockCoordinates3D const& newBlock = dataDistribution.getBlockParameters(iBlock).getEnvelope();
            if (myBlocks.empty()) {
                boundingBox = newBlock;
            }
            else {
                if (newBlock.x0 < boundingBox.x0) boundingBox.x0 = newBlock.x0;
                if (newBlock.x1 > boundingBox.x1) boundingBox.x1 = newBlock.x1;
                if (newBlock.y0 < boundingBox.y0) boundingBox.y0 = newBlock.y0;
                if (newBlock.y1 > boundingBox.y1) boundingBox.y1 = newBlock.y1;
                if (newBlock.z0 < boundingBox.z0) boundingBox.z0 = newBlock.z0;
                if (newBlock.z1 > boundingBox.z1) boundingBox.z1 = newBlock.z1;
            }
            myBlocks.push_back(iBlock);
            nearbyBlocks.push_back(iBlock);
        }
    }
    int myRank = singleton::mpi().getRank();
    for (int iOverlap=0; iOverlap<dataDistribution.getNumNormalOverlaps(); ++iOverlap) {
        Overlap3D const& overlap = dataDistribution.getNormalOverlap(iOverlap);
        int originalProc = dataDistribution.getBlockParameters(overlap.getOriginalId()).getProcId();
        int overlapProc = dataDistribution.getBlockParameters(overlap.getOverlapId()).getProcId();
        if (originalProc == myRank) {
            nearbyBlocks.push_back( overlap.getOverlapId() );
        }
        if (originalProc == myRank || overlapProc == myRank) {
            relevantNormalOverlaps.push_back(iOverlap);
        }
    }
    for (int iOverlap=0; iOverlap<dataDistribution.getNumPeriodicOverlaps(); ++iOverlap) {
        Overlap3D const& overlap = dataDistribution.getPeriodicOverlap(iOverlap);
        int originalProc = dataDistribution.getBlockParameters(overlap.getOriginalId()).getProcId();
        int overlapProc = dataDistribution.getBlockParameters(overlap.getOverlapId()).getProcId();
        if (originalProc == myRank) {
            nearbyBlocks.push_back( overlap.getOverlapId() );
        }
        if (originalProc == myRank || overlapProc == myRank) {
            relevantPeriodicOverlaps.push_back(iOverlap);
        }
    }
    // Erase duplicates in nearbyBlocks
    std::sort(nearbyBlocks.begin(), nearbyBlocks.end());
    std::vector<int>::iterator newEnd = unique(nearbyBlocks.begin(), nearbyBlocks.end());
    nearbyBlocks.erase(newEnd, nearbyBlocks.end());
}

#endif  // PARALLEL_MODE_MPI

}  // namespace olb

#endif

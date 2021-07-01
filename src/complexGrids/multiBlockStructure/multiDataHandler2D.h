/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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
 * Handler for 2D multiblock structure -- header file.
 */
#ifndef MULTI_DATA_HANDLER_2D_H
#define MULTI_DATA_HANDLER_2D_H

#include <vector>
#include "core/blockLattice2D.h"
#include "core/dataFields2D.h"
#include "core/cell.h"
#include "multiDataGeometry2D.h"


/*
 * STATUS OF THE CODE
 *
 * One should disable the feature that lattice.get() lets you change the dynamics, because this will
 * not work in parallel. Everybody should use defineDynamics() instead.
 *
 * Local boundary conditions work fine, because
 * - in the cache-optimized version, overlapping regions point twice to the same dynamics, which is consistent
 * - in the parallel version, overlapping dynamics is duplicated, and it is updated consistently through the
 *   data-parallel model
 *
 * Non-local boundary conditions work fine, because
 * - they instantiate local dynamics, which works as explained above
 * - the postprocessing step acts on the bulk (without envelope), which works just fine
 *
 *
 * DISCUSSION
 *
 * BlockLatticeView should be able to have periodic boundaries
 *
 */


namespace olb {

template<typename T, template<typename U> class Lattice>
struct MultiBlockHandler2D {
    typedef std::vector<BlockLattice2D<T,Lattice>*> BlockVector2D;

    virtual ~MultiBlockHandler2D() { }
    virtual int getNx() const =0;
    virtual int getNy() const =0;
    virtual MultiDataDistribution2D const& getMultiDataDistribution() const =0;
    virtual bool getLocalEnvelope(int iBlock, int& lx, int& ly) const =0;
    virtual T reduceSum(T localSum) const =0;
    virtual T reduceAverage(T localAverage, T localWeight) const =0;
    virtual T reduceMin(T localMin) const =0;
    virtual T reduceMax(T localMax) const =0;
    virtual void broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const =0;
    virtual void broadCastScalar(T& scalar, int fromBlock) const =0;
    virtual void broadCastVector(T vect[Lattice<T>::d], int fromBlock) const =0;
    virtual void copyOverlap(Overlap2D const& overlap, BlockVector2D& lattices) const =0;
    virtual Cell<T,Lattice>& getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const =0;
    virtual Cell<T,Lattice> const& getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const =0;
};

template<typename T>
struct MultiDataFieldHandler2D {
    virtual ~MultiDataFieldHandler2D() { }
    virtual int getNx() const =0;
    virtual int getNy() const =0;
    virtual MultiDataDistribution2D const& getMultiDataDistribution() const =0;
    virtual bool getLocalEnvelope(int iBlock, int& lx, int& ly) const =0;
    virtual T reduceSum(T localSum) const =0;
    virtual T reduceAverage(T localAverage, T localWeight) const =0;
    virtual T reduceMin(T localMin) const =0;
    virtual T reduceMax(T localMax) const =0;
    virtual void broadCastScalar(T& scalar, int fromBlock) const =0;
    virtual void broadCastVector(T* vect, int size, int fromBlock) const =0;
};

template<typename T, template<typename U> class Lattice>
class SerialMultiBlockHandler2D : public MultiBlockHandler2D<T,Lattice> {
public:
    typedef std::vector<BlockLattice2D<T,Lattice>*> BlockVector2D;
public:
    SerialMultiBlockHandler2D(MultiDataDistribution2D const& dataDistribution_);
    virtual int getNx() const;
    virtual int getNy() const;
    virtual MultiDataDistribution2D const& getMultiDataDistribution() const;
    virtual bool getLocalEnvelope(int iBlock, int& lx, int& ly) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const;
    virtual void broadCastScalar(T& scalar, int fromBlock) const;
    virtual void broadCastVector(T vect[Lattice<T>::d], int fromBlock) const;
    virtual void copyOverlap(Overlap2D const& overlap, BlockVector2D& lattices) const;
    virtual Cell<T,Lattice>& getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const;
    virtual Cell<T,Lattice> const& getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const;
private:
    MultiDataDistribution2D dataDistribution;
};

template<typename T>
class SerialMultiDataFieldHandler2D : public MultiDataFieldHandler2D<T> {
public:
    SerialMultiDataFieldHandler2D(MultiDataDistribution2D const& dataDistribution_);
    virtual int getNx() const;
    virtual int getNy() const;
    virtual MultiDataDistribution2D const& getMultiDataDistribution() const;
    virtual bool getLocalEnvelope(int iBlock, int& lx, int& ly) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastScalar(T& scalar, int fromBlock) const;
    virtual void broadCastVector(T* vect, int size, int fromBlock) const;
private:
    MultiDataDistribution2D dataDistribution;
};


#ifdef PARALLEL_MODE_MPI

template<typename T, template<typename U> class Lattice>
class ParallelMultiBlockHandler2D : public MultiBlockHandler2D<T,Lattice> {
public:
    typedef std::vector<BlockLattice2D<T,Lattice>*> BlockVector2D;
public:
    ParallelMultiBlockHandler2D(MultiDataDistribution2D const& dataDistribution_);
    virtual int getNx() const;
    virtual int getNy() const;
    virtual MultiDataDistribution2D const& getMultiDataDistribution() const;
    virtual bool getLocalEnvelope(int iBlock, int& lx, int& ly) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const;
    virtual void broadCastScalar(T& scalar, int fromBlock) const;
    virtual void broadCastVector(T vect[Lattice<T>::d], int fromBlock) const;
    virtual void copyOverlap(Overlap2D const& overlap, BlockVector2D& lattices) const;
    virtual Cell<T,Lattice>& getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const;
    virtual Cell<T,Lattice> const& getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const;
private:
    MultiDataDistribution2D dataDistribution;
    mutable Cell<T,Lattice> distributedCell;
    mutable Dynamics<T,Lattice>* parallelDynamics;
};

template<typename T>
class ParallelMultiDataFieldHandler2D : public MultiDataFieldHandler2D<T> {
public:
    ParallelMultiDataFieldHandler2D(MultiDataDistribution2D const& dataDistribution_);
    virtual int getNx() const;
    virtual int getNy() const;
    virtual MultiDataDistribution2D const& getMultiDataDistribution() const;
    virtual bool getLocalEnvelope(int iBlock, int& lx, int& ly) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastScalar(T& scalar, int fromBlock) const;
    virtual void broadCastVector(T* vect, int size, int fromBlock) const;
private:
    MultiDataDistribution2D dataDistribution;
};

#endif

}  // namespace olb

#endif

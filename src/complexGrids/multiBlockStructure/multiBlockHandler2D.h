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
#ifndef MULTI_BLOCK_HANDLER_2D_H
#define MULTI_BLOCK_HANDLER_2D_H

#include <vector>
#include "core/blockLattice2D.h"
#include "core/cell.h"
#include "multiDataGeometry2D.h"


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
#endif

}  // namespace olb

#endif

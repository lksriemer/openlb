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
 * Handler for 3D multiblock structure -- header file.
 */
#ifndef MULTI_DATA_HANDLER_3D_H
#define MULTI_DATA_HANDLER_3D_H

#include <vector>
#include "core/blockLattice3D.h"
#include "core/dataFields3D.h"
#include "core/cell.h"
#include "multiDataGeometry3D.h"


namespace olb {

template<typename T, template<typename U> class Lattice>
struct MultiDataHandler3D {
    typedef std::vector<BlockLattice3D<T,Lattice>*> BlockVector3D;
    typedef std::vector<ScalarField3D<T>*>          ScalarFieldVector3D;
    typedef std::vector<TensorField3D<T,3>*>        VectorFieldVector3D;
    typedef std::vector<TensorField3D<T,6>*>        SymmetricMatrixFieldVector3D;

    virtual ~MultiDataHandler3D() { }
    virtual int getNx() const =0;
    virtual int getNy() const =0;
    virtual int getNz() const =0;
    virtual MultiDataDistribution3D const& getMultiDataDistribution() const =0;
    virtual void allocateBlockLattices(BlockVector3D& lattices) const =0;
    virtual void allocateScalarFields(ScalarFieldVector3D& fields) const =0;
    virtual void allocateVectorFields(VectorFieldVector3D& fields) const =0;
    virtual void allocateSymmetricMatrixFields(SymmetricMatrixFieldVector3D& fields) const =0;
    virtual T reduceSum(T localSum) const =0;
    virtual T reduceAverage(T localAverage, T localWeight) const =0;
    virtual T reduceMin(T localMin) const =0;
    virtual T reduceMax(T localMax) const =0;
    virtual void broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const =0;
    virtual void broadCastScalar(T& scalar, int fromBlock) const =0;
    virtual void broadCastVector(T vect[Lattice<T>::d], int fromBlock) const =0;
    virtual void copyOverlap(Overlap3D const& overlap, BlockVector3D& lattices) const =0;
    virtual Cell<T,Lattice>& getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const =0;
    virtual Cell<T,Lattice> const& getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const =0;
};

template<typename T, template<typename U> class Lattice>
class SerialDataHandler3D : public MultiDataHandler3D<T,Lattice> {
public:
    typedef std::vector<BlockLattice3D<T,Lattice>*> BlockVector3D;
    typedef std::vector<ScalarField3D<T>*>          ScalarFieldVector3D;
    typedef std::vector<TensorField3D<T,3>*>        VectorFieldVector3D;
    typedef std::vector<TensorField3D<T,6>*>        SymmetricMatrixFieldVector3D;
public:
    SerialDataHandler3D(MultiDataDistribution3D const& dataDistribution_);
    virtual int getNx() const;
    virtual int getNy() const;
    virtual int getNz() const;
    virtual MultiDataDistribution3D const& getMultiDataDistribution() const;
    virtual void allocateBlockLattices(BlockVector3D& lattices) const;
    virtual void allocateScalarFields(ScalarFieldVector3D& fields) const;
    virtual void allocateVectorFields(VectorFieldVector3D& fields) const;
    virtual void allocateSymmetricMatrixFields(SymmetricMatrixFieldVector3D& fields) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const;
    virtual void broadCastScalar(T& scalar, int fromBlock) const;
    virtual void broadCastVector(T vect[Lattice<T>::d], int fromBlock) const;
    virtual void copyOverlap(Overlap3D const& overlap, BlockVector3D& lattices) const;
    virtual Cell<T,Lattice>& getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const;
    virtual Cell<T,Lattice> const& getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const;
private:
    MultiDataDistribution3D dataDistribution;
};


#ifdef PARALLEL_MODE_MPI
template<typename T, template<typename U> class Lattice>
class ParallelDataHandler3D : public MultiDataHandler3D<T,Lattice> {
public:
    typedef std::vector<BlockLattice3D<T,Lattice>*> BlockVector3D;
    typedef std::vector<ScalarField3D<T>*>          ScalarFieldVector3D;
    typedef std::vector<TensorField3D<T,3>*>        VectorFieldVector3D;
    typedef std::vector<TensorField3D<T,6>*>        SymmetricMatrixFieldVector3D;
public:
    ParallelDataHandler3D(MultiDataDistribution3D const& dataDistribution_);
    virtual int getNx() const;
    virtual int getNy() const;
    virtual int getNz() const;
    virtual MultiDataDistribution3D const& getMultiDataDistribution() const;
    virtual void allocateBlockLattices(BlockVector3D& lattices) const;
    virtual void allocateScalarFields(ScalarFieldVector3D& fields) const;
    virtual void allocateVectorFields(VectorFieldVector3D& fields) const;
    virtual void allocateSymmetricMatrixFields(SymmetricMatrixFieldVector3D& fields) const;
    virtual T reduceSum(T localSum) const;
    virtual T reduceAverage(T localAverage, T localWeight) const;
    virtual T reduceMin(T localMin) const;
    virtual T reduceMax(T localMax) const;
    virtual void broadCastCell(Cell<T,Lattice>& cell, int fromBlock) const;
    virtual void broadCastScalar(T& scalar, int fromBlock) const;
    virtual void broadCastVector(T vect[Lattice<T>::d], int fromBlock) const;
    virtual void copyOverlap(Overlap3D const& overlap, BlockVector3D& lattices) const;
    virtual Cell<T,Lattice>& getDistributedCell(Cell<T,Lattice>* baseCell, int onBlock) const;
    virtual Cell<T,Lattice> const& getDistributedCell(Cell<T,Lattice> const* baseCell, int onBlock) const;
private:
    MultiDataDistribution3D dataDistribution;
    mutable Cell<T,Lattice> distributedCell;
    mutable Dynamics<T,Lattice>* parallelDynamics;
};
#endif

}  // namespace olb

#endif

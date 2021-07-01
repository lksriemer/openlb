/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Orestis Malaspinas, Jonas Latt
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
 * Interface for post-processing steps -- header file.
 */
#ifndef N_BLOCK_LATTICE_POST_PROCESSING_H
#define N_BLOCK_LATTICE_POST_PROCESSING_H

#include <vector>
#include "core/ompManager.h"
#include "core/postProcessing.h"

namespace olb {

/////////////////// Forward Declarations /////////////////////////////

template<typename T, template<typename U> class Lattice>
class BlockLattice2D;

// template<typename T, template<typename U> class Lattice>
// class BlockLattice3D;


/////////////////// 2D Postprocessing ///////////////////////////////

/// Interface of 2D post-processing steps.
template<typename T, template<typename U> class Lattice>
struct NblockLatticePostProcessor2D {
    virtual ~NblockLatticePostProcessor2D() { }
    /// Execute post-processing step
    virtual void process(BlockLattice2D<T,Lattice> &blockLattice) =0;
    /// Execute post-processing step on a sublattice
    virtual void processSubDomain(std::vector<BlockLattice2D<T,Lattice>* > vecBlockLattice,
                                  int x0_, int x1_, int y0_, int y1_
                                  /*std::vector<BlockCoordinates2D> vecCoordinates*/) =0;
    /// Extent of application area (0 for purely local operations)
    virtual int extent() const =0;
    /// Extent of application area along a direction (0 or 1)
    virtual int extent(int direction) const =0;
    virtual bool hasReductions() const =0;
    virtual void subscribeReductions(std::vector<BlockLattice2D<T,Lattice>* > vecBlockLattice,
                                     Reductor<T>* reductor) =0;
};

template<typename T, template<typename U> class Lattice>
class NblockLatticePostProcessorGenerator2D : public LatticeCouplingGenerator2D
{
public:
    NblockLatticePostProcessorGenerator2D(/*std::vector<BlockCoordinates2D> vecCoordinates_*/);
    virtual ~NblockLatticePostProcessorGenerator2D() { }
    void shift(int deltaX0, int deltaX1, int deltaY0, int deltaY1/*, std::vector<int> deltas*/);
    bool extract(int x0_, int x1_, int y0_, int y1_/*,std::vector<BlockCoordinates2D> vecCoordinates_*/);
    virtual NblockLatticePostProcessor2D<T,Lattice>* generate() const =0;
    virtual NblockLatticePostProcessorGenerator2D<T,Lattice>* clone() const =0;
protected:
    int x0, x1, y0, y1;
//     std::vector<BlockCoordinates2D> vecCoordinates;
};


template<typename T, template<typename U> class Lattice>
struct NblockLatticeLocalPostProcessor2D : public NblockLatticePostProcessor2D<T,Lattice> {
    virtual bool hasReductions() const { return false; }
    virtual void subscribeReductions(std::vector<BlockLattice2D<T,Lattice>* > vecBlockLattice,
                                     Reductor<T>* reductor)
    { }
};

template<typename T, template<typename U> class Lattice>
struct NblockLatticeGlobalPostProcessor2D : public NblockLatticePostProcessor2D<T,Lattice> {
    virtual bool hasReductions() const { return true; }
    virtual void process(std::vector<BlockLattice2D<T,Lattice>* > vecBlockLattice) =0;
    virtual void processSubDomain(std::vector<BlockLattice2D<T,Lattice>* > vecBlockLattice,
                                  int x0_, int x1_, int y0_, int y1_/*,
                                  std::vector<BlockCoordinates2D> vecCoordinates*/)
    {
        this -> process(vecBlockLattice);
    }
    virtual int extent() const {
        return 0;
    }
    virtual int extent(int direction) const {
        return 0;
    }
};


// /////////////////// 3D Postprocessing ///////////////////////////////

}  // namespace olb

#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
#ifndef TWO_BLOCK_LATTICE_POST_PROCESSING_H
#define TWO_BLOCK_LATTICE_POST_PROCESSING_H

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
struct TwoBlockLatticePostProcessor2D {
    virtual ~TwoBlockLatticePostProcessor2D() { }
    /// Execute post-processing step
    virtual void process(BlockLattice2D<T,Lattice>& blockLatticeOne,
						 BlockLattice2D<T,Lattice>& blockLatticeTwo) =0;
    /// Execute post-processing step on a sublattice
    virtual void processSubDomain(BlockLattice2D<T,Lattice>& blockLatticeOne,
								  int x0One_, int x1One_, int y0One_, int y1One_,
								  BlockLattice2D<T,Lattice>& blockLatticeOne,
                                  int x0Two_, int x1Two_, int y0Two_, int y1Two_) =0;
    /// Extent of application area (0 for purely local operations)
    virtual int extent() const =0;
    /// Extent of application area along a direction (0 or 1)
    virtual int extent(int direction) const =0;
    virtual bool hasReductions() const =0;
    virtual void subscribeReductions(BlockLattice2D<T,Lattice>& blockLatticeOne,
									 BlockLattice2D<T,Lattice>& blockLatticeTwo,
                                     Reductor<T>* reductor) =0;
};

template<typename T, template<typename U> class Lattice>
class TwoBlockLatticePostProcessorGenerator2D {
public:
    TwoBlockLatticePostProcessorGenerator2D(int x0One_, int x1One_, int y0One_, int y1One_,
											int x0Two_, int x1Two_, int y0Two_, int y1Two_);
    virtual ~TwoBlockLatticePostProcessorGenerator2D() { }
    void shift(int deltaXOne, int deltaYOne,
			   int deltaXTwo, int deltaYTwo);
    bool extract(int x0One_, int x1One_, int y0One_, int y1One_,
				 int x0Two_, int x1Two_, int y0Two_, int y1Two_);
    virtual TwoBlockLatticePostProcessor2D<T,Lattice>* generate() const =0;
    virtual TwoBlockLatticePostProcessorGenerator2D<T,Lattice>* clone() const =0;
protected:
    int x0One, x1One, y0One, y1One;
	int x0Two, x1Two, y0Two, y1Two;
};


template<typename T, template<typename U> class Lattice>
struct TwoBlockLatticeLocalPostProcessor2D : public TwoBlockLatticePostProcessor2D<T,Lattice> {
    virtual bool hasReductions() const { return false; }
    virtual void subscribeReductions(BlockLattice2D<T,Lattice>& blockLatticeOne,
                                     BlockLattice2D<T,Lattice>& blockLatticeTwo,
                                     Reductor<T>* reductor)
    { }
};

template<typename T, template<typename U> class Lattice>
struct TwoBlockLatticeGlobalPostProcessor2D : public TwoBlockLatticePostProcessor2D<T,Lattice> {
    virtual bool hasReductions() const { return true; }
    virtual void process(BlockLattice2D<T,Lattice>& blockLatticeOne,
                         BlockLattice2D<T,Lattice>& blockLatticeTwo) =0;
    virtual void processSubDomain(BlockLattice2D<T,Lattice>& blockLatticeOne,
                                  int x0One_, int x1One_, int y0One_, int y1One_,
                                  BlockLattice2D<T,Lattice>& blockLatticeTwo,
                                  int x0Two_, int x1Two_, int y0Two_, int y1Two_ )
    {
        this -> process(blockLatticeOne,blockLatticeTwo);
    }
    virtual int extent() const {
        return 0;
    }
    virtual int extent(int direction) const {
        return 0;
    }
};


// /////////////////// 3D Postprocessing ///////////////////////////////
// 
// template<typename T, template<typename U> class Lattice>
// struct PostProcessor3D {
//     virtual ~PostProcessor3D() { }
//     /// Execute post-processing step
//     virtual void process(BlockLattice3D<T,Lattice>& blockLattice) =0;
//         /// Execute post-processing step on a sublattice
//     virtual void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
//                                   int x0_, int x1_, int y0_, int y1_,
//                          int z0_, int z1_ ) =0;
//     /// Extent of application area (0 for purely local operations)
//     virtual int extent() const =0;
//     /// Extent of application area along a direction (0 or 1)
//     virtual int extent(int direction) const =0;
//     virtual bool hasReductions() const =0;
//     virtual void subscribeReductions(BlockLattice3D<T,Lattice>& blockLattice,
//                                      Reductor<T>* reductor) =0;
// };
// 
// template<typename T, template<typename U> class Lattice>
// class PostProcessorGenerator3D {
// public:
//     PostProcessorGenerator3D( int x0_, int x1_, int y0_, int y1_,
//                               int z0_, int z1_ );
//     virtual ~PostProcessorGenerator3D() { }
//     void shift(int deltaX, int deltaY, int deltaZ);
//     bool extract(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
//     virtual PostProcessor3D<T,Lattice>* generate() const =0;
//     virtual PostProcessorGenerator3D<T,Lattice>* clone() const =0;
// protected:
//     int x0, x1, y0, y1, z0, z1;
// };
// 
// 
// template<typename T, template<typename U> class Lattice>
// struct LocalPostProcessor3D : public PostProcessor3D<T,Lattice> {
//     virtual bool hasReductions() const { return false; }
//     virtual void subscribeReductions(BlockLattice3D<T,Lattice>& blockLattice,
//                                      Reductor<T>* reductor)
//     { }
// };
// 
// template<typename T, template<typename U> class Lattice>
// struct GlobalPostProcessor3D : public PostProcessor3D<T,Lattice> {
//     virtual bool hasReductions() const { return true; }
//     virtual void process(BlockLattice3D<T,Lattice>& blockLattice) =0;
//     virtual void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
//                                   int x0_, int x1_, int y0_, int y1_,
//                                   int z0_, int z1_ )
//     {
//         this -> process(blockLattice);
//     }
//     virtual int extent() const {
//         return 0;
//     }
//     virtual int extent(int direction) const {
//         return 0;
//     }
// };

}  // namespace olb

#endif

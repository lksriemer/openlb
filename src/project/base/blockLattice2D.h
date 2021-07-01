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
 * The dynamics of a 2D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_2D_H
#define BLOCK_LATTICE_2D_H

#include <vector>
#include "olbDebug.h"
#include "postProcessing.h"
#include "dataFields2D.h"
#include "blockStructure2D.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T, template<typename U> class Lattice> struct Dynamics;
template<typename T, template<typename U> class Lattice> class Cell;

/// A regular lattice for highly efficient 2D LB dynamics.
/** A block lattice contains a regular array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Lattice>
class BlockLattice2D : public BlockStructure2D<T,Lattice> {
public:
    typedef std::vector<PostProcessor2D<T,Lattice>*> PostProcVector;
public:
    /// Construction of an nx_ by ny_ lattice
    BlockLattice2D(int nx_, int ny_);
    /// Destruction of the lattice
    ~BlockLattice2D();
    /// Copy construction
    BlockLattice2D(BlockLattice2D<T,Lattice> const& rhs);
    /// Copy assignment
    BlockLattice2D& operator=(BlockLattice2D<T,Lattice> const& rhs);
    /// Swap the content of two BlockLattices
    void swap(BlockLattice2D& rhs);
public:
    /// Read access to lattice width
    virtual int getNx() const { return nx; }
    /// Read access to lattice height
    virtual int getNy() const { return ny; }
    /// Read/write access to lattice cells
    virtual Cell<T,Lattice>& get(int iX, int iY) {
        OLB_PRECONDITION(iX<nx);
        OLB_PRECONDITION(iY<ny);
        return grid[iX][iY];
    }
    /// Read only access to lattice cells
    virtual Cell<T,Lattice> const& get(int iX, int iY) const {
        OLB_PRECONDITION(iX<nx);
        OLB_PRECONDITION(iY<ny);
        return grid[iX][iY];
    }
    /// Initialize the lattice cells to become ready for simulation
    virtual void initialize();
    /// Define the dynamics on a rectangular domain
    virtual void defineDynamics (int x0, int x1, int y0, int y1,
                                 Dynamics<T,Lattice>* dynamics );
    /// Specify wheter statistics measurements are done on given rect. domain
    virtual void specifyStatisticsStatus (int x0, int x1, int y0, int y1,
                                          bool status );
    /// Apply collision step to a rectangular domain
    virtual void collide(int x0, int x1, int y0, int y1);
    /// Apply collision step to the whole domain
    virtual void collide();
    /// Apply collision step to a rectangular domain, with fixed velocity
    virtual void staticCollide (int x0, int x1, int y0, int y1,
                                TensorField2D<T,2> const& u);
    /// Apply collision step to the whole domain, with fixed velocity
    virtual void staticCollide(TensorField2D<T,2> const& u);
    /// Apply streaming step to a rectangular domain
    virtual void stream(int x0, int x1, int y0, int y1, bool periodic=false);
    /// Apply streaming step to the whole domain
    virtual void stream(bool periodic=false);
    /// Apply first collision, then streaming step to a rectangular domain
    virtual void collideAndStream (
                int x0, int x1, int y0, int y1, bool periodic=false );
    /// Apply first collision, then streaming step to the whole domain
    virtual void collideAndStream(bool periodic=false);
    /// Compute the average density within a rectangular domain
    virtual T computeAverageDensity(int x0, int x1, int y0, int y1) const;
    /// Compute the average density within the whole domain
    virtual T computeAverageDensity() const;
    /// Subtract a constant offset from the density within the whole domain
    virtual void stripeOffDensityOffset (
                int x0, int x1, int y0, int y1, T offset );
    /// Subtract a constant offset from the density within a rect. domain
    virtual void stripeOffDensityOffset(T offset);
    /// Add a non-local post-processing step
    virtual void addPostProcessor (
                PostProcessorGenerator2D<T,Lattice> const& ppGen );
    /// Clean up all non-local post-processing steps
    virtual void resetPostProcessors();
    /// Execute post-processing on a sub-lattice
    virtual void postProcess(int x0_, int x1_, int y0_, int y1_);
    /// Execute post-processing steps
    virtual void postProcess();
    /// Subscribe postProcessors for reduction operations
    virtual void subscribeReductions(Reductor<T>& reductor);
    /// Return a handle to the LatticeStatistics object
    virtual LatticeStatistics<T>& getStatistics();
    /// Return a constant handle to the LatticeStatistics object
    virtual LatticeStatistics<T> const& getStatistics() const;
private:
    /// Helper method for memory allocation
    void allocateMemory();
    /// Helper method for memory de-allocation
    void releaseMemory();
    /// Release memory for post processors
    void clearPostProcessors();
    /// Apply streaming step to boundary cells
    void boundaryStream (
            int lim_x0, int lim_x1, int lim_y0, int lim_y1,
            int x0, int x1, int y0, int y1 );
    /// Apply streaming step to bulk (non-boundary) cells
    void bulkStream(int x0, int x1, int y0, int y1);
    /// Apply collision and streaming step to bulk (non-boundary) cells
    void bulkCollideAndStream(int x0, int x1, int y0, int y1);
    template<int normalX, int normalY> void periodicEdge(int from, int to);
    void makePeriodic();
private:
    int                  nx, ny;
    Cell<T,Lattice>      *rawData;
    Cell<T,Lattice>      **grid;
    PostProcVector       postProcessors;
    LatticeStatistics<T> statistics;
};

}  // namespace olb

#endif

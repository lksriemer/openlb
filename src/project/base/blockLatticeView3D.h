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
 * Dynamics for a generic 3D block lattice view -- header file.
 */

#ifndef BLOCK_LATTICE_VIEW_3D_H
#define BLOCK_LATTICE_VIEW_3D_H

#include <vector>
#include "blockLattice3D.h"

namespace olb {

template<typename T, template<typename U> class Lattice>
class BlockLatticeView3D : public BlockStructure3D<T,Lattice> {
public:
    BlockLatticeView3D(BlockLattice3D<T,Lattice>& originalLattice_);
    BlockLatticeView3D(BlockLattice3D<T,Lattice>& originalLattice_,
                       int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
    ~BlockLatticeView3D();
    BlockLatticeView3D(BlockLatticeView3D const& rhs);
    BlockLatticeView3D<T,Lattice>& operator= 
      (BlockLatticeView3D<T,Lattice> const& rhs);
    void swap(BlockLatticeView3D<T,Lattice>& rhs);

    virtual int getNx() const { return nx; }
    virtual int getNy() const { return ny; }
    virtual int getNz() const { return nz; }
    virtual Cell<T,Lattice>& get(int iX, int iY, int iZ);
    virtual Cell<T,Lattice> const& get(int iX, int iY, int iZ) const;
    virtual void initialize();
    virtual void defineDynamics (
        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
        Dynamics<T,Lattice>* dynamics );
    virtual void specifyStatisticsStatus (
                int x0_, int x1_, int y0_, int y1_,
                int z0_, int z1_, bool status );
    virtual void collide (
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ );
    virtual void collide();
    virtual void staticCollide (int x0, int x1, int y0, int y1,
                                int z0_, int z1_,
                                TensorField3D<T,3> const& u);
    virtual void staticCollide (TensorField3D<T,3> const& u);
    virtual void stream (
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
            bool periodic=false );
    virtual void stream(bool periodic=false);
    virtual void collideAndStream (
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
            bool periodic=false );
    virtual void collideAndStream(bool periodic=false);
    virtual T computeAverageDensity (
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) const;
    virtual T computeAverageDensity() const;
    virtual void stripeOffDensityOffset (
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T offset );
    virtual void stripeOffDensityOffset(T offset);
    virtual void addPostProcessor (
                     PostProcessorGenerator3D<T,Lattice> const& ppGen);
    virtual void resetPostProcessors();
    virtual void postProcess(int x0_, int x1_, int y0_, int y1_,
                             int z0_, int z1_);
    virtual void postProcess();
    virtual void subscribeReductions(Reductor<T>& reductor);
    virtual LatticeStatistics<T>& getStatistics();
    virtual LatticeStatistics<T> const& getStatistics() const;
private:
    BlockLattice3D<T,Lattice>  *originalLattice;
    int                        x0, y0, z0;
    int                        nx,ny,nz;
};

}  // namespace olb

#endif

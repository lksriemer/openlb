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

/** \file A helper for initialising 2D boundaries -- header file.  */
#ifndef DYNAMICS_HELPERS_2D_H
#define DYNAMICS_HELPERS_2D_H

#include "blockStructure2D.h"
#include "boundaries2D.h"
#include "fdBoundaries2D.h"
#include <vector>
#include <list>

namespace olb {

template<typename T, template<typename U> class Lattice>
class BoundaryCondition2D {
public:
    BoundaryCondition2D(BlockStructure2D<T,Lattice>& block_);

    virtual ~BoundaryCondition2D() { }

    virtual void addVelocityBoundary0N(int x0, int x1, int y0, int y1,
                                       T omega) =0;
    virtual void addVelocityBoundary0P(int x0, int x1, int y0, int y1,
                                       T omega) =0;
    virtual void addVelocityBoundary1N(int x0, int x1, int y0, int y1,
                                       T omega) =0;
    virtual void addVelocityBoundary1P(int x0, int x1, int y0, int y1,
                                       T omega) =0;

    virtual void addPressureBoundary0N(int x0, int x1, int y0, int y1,
                                       T omega) =0;
    virtual void addPressureBoundary0P(int x0, int x1, int y0, int y1,
                                       T omega) =0;
    virtual void addPressureBoundary1N(int x0, int x1, int y0, int y1,
                                       T omega) =0;
    virtual void addPressureBoundary1P(int x0, int x1, int y0, int y1,
                                       T omega) =0;

    virtual void addExternalVelocityCornerNN(int x, int y, T omega) =0;
    virtual void addExternalVelocityCornerNP(int x, int y, T omega) =0;
    virtual void addExternalVelocityCornerPN(int x, int y, T omega) =0;
    virtual void addExternalVelocityCornerPP(int x, int y, T omega) =0;

    virtual void addInternalVelocityCornerNN(int x, int y, T omega) =0;
    virtual void addInternalVelocityCornerNP(int x, int y, T omega) =0;
    virtual void addInternalVelocityCornerPN(int x, int y, T omega) =0;
    virtual void addInternalVelocityCornerPP(int x, int y, T omega) =0;

    BlockStructure2D<T,Lattice>& getBlock();
    BlockStructure2D<T,Lattice> const& getBlock() const;
private:
    BlockStructure2D<T,Lattice>& block;
};

template<typename T, template<typename U> class Lattice, typename Dynamics>
BoundaryCondition2D<T,Lattice>*
    createLocalBoundaryCondition2D(BlockStructure2D<T,Lattice>& lattice);

template<typename T, template<typename U> class Lattice, typename Dynamics>
BoundaryCondition2D<T,Lattice>*
createInterpBoundaryCondition2D(BlockStructure2D<T,Lattice>& lattice);

template<typename T, template<typename U> class Lattice>
BoundaryCondition2D<T,Lattice>*
createInterpBoundaryCondition2D(BlockStructure2D<T,Lattice>& lattice) {
    return createInterpBoundaryCondition2D<T,Lattice,RLBdynamics<T,Lattice> >(lattice);
}

template<typename T, template<typename U> class Lattice>
BoundaryCondition2D<T,Lattice>*
createLocalBoundaryCondition2D(BlockStructure2D<T,Lattice>& lattice) {
    return createLocalBoundaryCondition2D<T,Lattice,RLBdynamics<T,Lattice> >(lattice);
}

}  // namespace olb

#endif

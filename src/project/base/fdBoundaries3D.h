/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
 *  E-mail: jonas.latt@gmail.com
 *
 *  Generic collision, which modifies the particle distribution
 *  functions, implemented by Orestis Malaspinas, 2007
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

#ifndef FD_BOUNDARIES_3D_H
#define FD_BOUNDARIES_3D_H

#include "postProcessing.h"
#include "boundaries.h"
#include "blockLattice3D.h"

namespace olb {

template<typename T, template<typename U> class Lattice,
         int direction, int orientation, int deriveDirection,
         bool orthogonal>
struct DirectedGradients3D {
    static void interpolate(T rhoU[Lattice<T>::d],
                            BlockLattice3D<T,Lattice> const& blockLattice,
                            int iX, int iY, int iZ, bool multiplyRho);
};


/**
* This class computes the skordos BC
* on a plane wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
class PlaneFdBoundary3D : public LocalPostProcessor3D<T,Lattice>
{
public:
    PlaneFdBoundary3D (int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                       T omega_, BlockLattice3D<T,Lattice>& blockLattice);
    ~PlaneFdBoundary3D();
    virtual int extent() const { return 1; }
    virtual int extent(int whichDirection) const { return 1; }
    virtual void process(BlockLattice3D<T,Lattice>& blockLattice);
    virtual void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_,
                                  int z0_, int z1_ );
private:
    void initialize(BlockLattice3D<T,Lattice>& blockLattice);
    template<int deriveDirection>
    void interpolateGradients (
            BlockLattice3D<T,Lattice> const& blockLattice,
            T rhoU[Lattice<T>::d], int iX, int iY, int iZ ) const;
private:
    int x0, x1, y0, y1, z0, z1;
    T   omega;
    StressInterpolationDirichletBM<T,Lattice,HydroBM,direction,orientation>
                            *boundaryMomenta;
    Dynamics  **boundaryDynamics;
};

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
class PlaneFdBoundaryGenerator3D
    : public PostProcessorGenerator3D<T,Lattice>
{
public:
    PlaneFdBoundaryGenerator3D(int x0_, int x1_, int y0_, int y1_,
                               int z0_, int z1_, T omega_);
    virtual PostProcessor3D<T,Lattice>* generate (
               BlockLattice3D<T,Lattice>& blockLattice ) const; 
    virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
    T omega;
};


/**
* This class computes the skordos BC
* on a convex edge wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
class ConvexVelocityEdge3D : public LocalPostProcessor3D<T,Lattice> {
public:
    enum { direction1 = (plane+1)%3, direction2 = (plane+2)%3 };
public:
    ConvexVelocityEdge3D (
        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
        T omega_, BlockLattice3D<T,Lattice>& blockLattice);
    ~ConvexVelocityEdge3D();
    virtual int extent() const { return 2; }
    virtual int extent(int whichDirection) const { return 2; }
    virtual void process(BlockLattice3D<T,Lattice>& blockLattice);
    virtual void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_,
                                  int z0_, int z1_ );
private:
    void initialize(BlockLattice3D<T,Lattice>& blockLattice);
    T getNeighborRho(int x, int y, int z, int step1, int step2,
                     BlockLattice3D<T,Lattice> const& blockLattice);
    template<int deriveDirection, int orientation>
    void interpolateGradients (
            BlockLattice3D<T,Lattice> const& blockLattice,
            T rhoU[Lattice<T>::d], int iX, int iY, int iZ ) const;
private:
    int x0, x1, y0, y1, z0, z1;
    T   omega;
    GenericBoundaryMomenta<T,Lattice> *boundaryMomenta;
    Dynamics           **boundaryDynamics;
};

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
class ConvexVelocityEdgeGenerator3D
    : public PostProcessorGenerator3D<T,Lattice>
{
public:
    ConvexVelocityEdgeGenerator3D(int x0_, int x1_, int y0_, int y1_,
                                  int z0_, int z1_, T omega_);
    virtual PostProcessor3D<T,Lattice>* generate (
               BlockLattice3D<T,Lattice>& blockLattice ) const; 
    virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
    T omega;
};


template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
class ConvexVelocityCorner3D : public LocalPostProcessor3D<T,Lattice> {
public:
    ConvexVelocityCorner3D(int x_, int y_, int z_, T omega_,
                           BlockLattice3D<T,Lattice>& blockLattice);
    ~ConvexVelocityCorner3D();
    virtual int extent() const { return 2; }
    virtual int extent(int whichDirection) const { return 2; }
    virtual void process(BlockLattice3D<T,Lattice>& blockLattice);
    virtual void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_,
                                  int z0_, int z1_ );
private:
    void initialize(BlockLattice3D<T,Lattice>& blockLattice);
    int x,y,z;
    T   omega;
    GenericBoundaryMomenta<T,Lattice> *boundaryMomenta;
    Dynamics           *boundaryDynamics;
};

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
class ConvexVelocityCornerGenerator3D
    : public PostProcessorGenerator3D<T,Lattice>
{
public:
    ConvexVelocityCornerGenerator3D(int x_, int y_, int z_, T omega_);
    virtual PostProcessor3D<T,Lattice>* generate (
               BlockLattice3D<T,Lattice>& blockLattice ) const; 
    virtual PostProcessorGenerator3D<T,Lattice>*  clone() const;
private:
    T omega;
};


}

#endif

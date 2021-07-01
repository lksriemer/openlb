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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef DYNAMICS_HELPERS_2D_HH
#define DYNAMICS_HELPERS_2D_HH

#include "dynamicsHelpers2D.h"
#include "cell.h"
#include "blockLattice2D.h"

namespace olb {

///////// class InterpBoundaryCondition2D (declaration)  //////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
class InterpBoundaryCondition2D:
    public BoundaryCondition2D<T,Lattice>
{
public:
    typedef StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,VelocityBM, 0,-1>
                LeftVelBdGen;
    typedef StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,VelocityBM, 0, 1>
                RightVelBdGen;
    typedef StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,VelocityBM, 1,-1>
                LowerVelBdGen;
    typedef StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,VelocityBM, 1, 1>
                UpperVelBdGen;

    typedef StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,PressureBM, 0,-1>
                LeftPresBdGen;
    typedef StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,PressureBM, 0, 1>
                RightPresBdGen;
    typedef StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,PressureBM, 1,-1>
                LowerPresBdGen;
    typedef StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,PressureBM, 1, 1>
                UpperPresBdGen;

    typedef ConvexVelocityCornerGenerator2D<T,Lattice,Dynamics, -1,-1>
                NN_ExternalCornerGen;
    typedef ConvexVelocityCornerGenerator2D<T,Lattice,Dynamics, -1, 1>
                NP_ExternalCornerGen;
    typedef ConvexVelocityCornerGenerator2D<T,Lattice,Dynamics,  1,-1>
                PN_ExternalCornerGen;
    typedef ConvexVelocityCornerGenerator2D<T,Lattice,Dynamics,  1, 1>
                PP_ExternalCornerGen;

    typedef ConcaveCornerVelBM2D<T,Lattice, -1,-1> NN_InternalCorner;
    typedef ConcaveCornerVelBM2D<T,Lattice, -1, 1> NP_InternalCorner;
    typedef ConcaveCornerVelBM2D<T,Lattice,  1,-1> PN_InternalCorner;
    typedef ConcaveCornerVelBM2D<T,Lattice,  1, 1> PP_InternalCorner;

    typedef Dynamics BdDynamics;

public:
    InterpBoundaryCondition2D(BlockStructure2D<T,Lattice>& block);
    ~InterpBoundaryCondition2D();

    virtual void addVelocityBoundary0N(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addVelocityBoundary0P(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addVelocityBoundary1N(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addVelocityBoundary1P(int x0, int x1, int y0, int y1,
                                       T omega);

    virtual void addPressureBoundary0N(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addPressureBoundary0P(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addPressureBoundary1N(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addPressureBoundary1P(int x0, int x1, int y0, int y1,
                                       T omega);

    virtual void addExternalVelocityCornerNN(int x, int y, T omega);
    virtual void addExternalVelocityCornerNP(int x, int y, T omega);
    virtual void addExternalVelocityCornerPN(int x, int y, T omega);
    virtual void addExternalVelocityCornerPP(int x, int y, T omega);

    virtual void addInternalVelocityCornerNN(int x, int y, T omega);
    virtual void addInternalVelocityCornerNP(int x, int y, T omega);
    virtual void addInternalVelocityCornerPN(int x, int y, T omega);
    virtual void addInternalVelocityCornerPP(int x, int y, T omega);

private:
    template<typename ListT>
      static void cleanList(ListT& list);
    template<typename BoundaryGenerator>
      void addStraightBoundary(int x0, int x1, int y0, int y1, T omega);
    template<typename BoundaryGenerator>
      void addNonLocalCorner(int x, int y, T omega);
    template<typename BoundaryType>
      void addLocalCorner(int x, int y, T omega,
                          std::list<BoundaryType*>& bdList);
private:
    std::list<NN_InternalCorner*> nnInternalBdList;
    std::list<NP_InternalCorner*> npInternalBdList;
    std::list<PN_InternalCorner*> pnInternalBdList;
    std::list<PP_InternalCorner*> ppInternalBdList;

    std::list<BdDynamics*> bdDynamicsArrayList;
};


///////// class LocalBoundaryCondition2D (declaration)  //////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
class LocalBoundaryCondition2D : public BoundaryCondition2D<T,Lattice> {
public:
    typedef RegularizedVelocityBM<T,Lattice, 0,-1> LeftVelBdMomenta;
    typedef RegularizedVelocityBM<T,Lattice, 0,1>  RightVelBdMomenta;
    typedef RegularizedVelocityBM<T,Lattice, 1,-1> LowerVelBdMomenta;
    typedef RegularizedVelocityBM<T,Lattice, 1,1>  UpperVelBdMomenta;

    typedef RegularizedPressureBM<T,Lattice, 0,-1> LeftPresBdMomenta;
    typedef RegularizedPressureBM<T,Lattice, 0,1>  RightPresBdMomenta;
    typedef RegularizedPressureBM<T,Lattice, 1,-1> LowerPresBdMomenta;
    typedef RegularizedPressureBM<T,Lattice, 1,1>  UpperPresBdMomenta;

    typedef std::vector<LeftVelBdMomenta>  LeftVelBdArray;
    typedef std::vector<RightVelBdMomenta> RightVelBdArray;
    typedef std::vector<LowerVelBdMomenta> LowerVelBdArray;
    typedef std::vector<UpperVelBdMomenta> UpperVelBdArray;

    typedef std::vector<LeftPresBdMomenta>  LeftPresBdArray;
    typedef std::vector<RightPresBdMomenta> RightPresBdArray;
    typedef std::vector<LowerPresBdMomenta> LowerPresBdArray;
    typedef std::vector<UpperPresBdMomenta> UpperPresBdArray;

    typedef CombinedRLBdynamics<T, Lattice, Dynamics>  BdDynamics;
    typedef std::vector<BdDynamics*> BdDynamicsArray;
public:
    LocalBoundaryCondition2D(BlockStructure2D<T,Lattice>& block);
    ~LocalBoundaryCondition2D();

    virtual void addVelocityBoundary0N(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addVelocityBoundary0P(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addVelocityBoundary1N(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addVelocityBoundary1P(int x0, int x1, int y0, int y1,
                                       T omega);

    virtual void addPressureBoundary0N(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addPressureBoundary0P(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addPressureBoundary1N(int x0, int x1, int y0, int y1,
                                       T omega);
    virtual void addPressureBoundary1P(int x0, int x1, int y0, int y1,
                                       T omega);

    virtual void addExternalVelocityCornerNN(int x, int y, T omega);
    virtual void addExternalVelocityCornerNP(int x, int y, T omega);
    virtual void addExternalVelocityCornerPN(int x, int y, T omega);
    virtual void addExternalVelocityCornerPP(int x, int y, T omega);

    virtual void addInternalVelocityCornerNN(int x, int y, T omega);
    virtual void addInternalVelocityCornerNP(int x, int y, T omega);
    virtual void addInternalVelocityCornerPN(int x, int y, T omega);
    virtual void addInternalVelocityCornerPP(int x, int y, T omega);

private:
    template<typename ListT>
      static void cleanList(ListT& list);
    template<typename ListT>
      static void cleanArrayList(ListT& list);
    template<typename ArrayType>
      void addGenericBoundary(int x0, int x1, int y0, int y1, T omega,
                              std::list<ArrayType*>& bdMomentaList);
private:
    std::list<LeftVelBdArray*>  leftVelBdList;
    std::list<RightVelBdArray*> rightVelBdList;
    std::list<LowerVelBdArray*> lowerVelBdList;
    std::list<UpperVelBdArray*> upperVelBdList;

    std::list<LeftPresBdArray*>  leftPresBdList;
    std::list<RightPresBdArray*> rightPresBdList;
    std::list<LowerPresBdArray*> lowerPresBdList;
    std::list<UpperPresBdArray*> upperPresBdList;

    std::list<BdDynamicsArray*>  bdDynamicsArrayList;

    InterpBoundaryCondition2D<T,Lattice,Dynamics > interpBoundary;
};


///////// class BoundaryCondition2D ///////////////////////////////////

template<typename T, template<typename U> class Lattice>
BoundaryCondition2D<T,Lattice>::BoundaryCondition2D (
        BlockStructure2D<T,Lattice>& block_ )
    : block(block_)
{ }

template<typename T, template<typename U> class Lattice>
BlockStructure2D<T,Lattice>& BoundaryCondition2D<T,Lattice>::getBlock() {
    return block;
}

template<typename T, template<typename U> class Lattice>
BlockStructure2D<T,Lattice> const& BoundaryCondition2D<T,Lattice>::
    getBlock() const
{
    return block;
}


///////// class InterpBoundaryCondition2D ////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
InterpBoundaryCondition2D<T,Lattice,Dynamics>::InterpBoundaryCondition2D (
        BlockStructure2D<T,Lattice>& block )
    : BoundaryCondition2D<T,Lattice>(block)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics>
InterpBoundaryCondition2D<T,Lattice,Dynamics>::~InterpBoundaryCondition2D() {
    cleanList(bdDynamicsArrayList);

    cleanList(nnInternalBdList);
    cleanList(npInternalBdList);
    cleanList(pnInternalBdList);
    cleanList(ppInternalBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename ListT>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::cleanList(ListT& list) {
    typename ListT::iterator iL = list.begin();
    for (; iL != list.end(); ++iL) {
        delete *iL;
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename BoundaryGenerator>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::
    addStraightBoundary(int x0, int x1, int y0, int y1, T omega)
{
    OLB_PRECONDITION(x0==x1 || y0==y1);

    BoundaryGenerator bGenerator(x0,x1,y0,y1, omega);
    (this->getBlock()).addPostProcessor(bGenerator);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename BoundaryGenerator>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::
    addNonLocalCorner(int x, int y, T omega)
{
    BoundaryGenerator bGenerator(x,y, omega);
    (this->getBlock()).addPostProcessor(bGenerator);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename BoundaryType>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::
    addLocalCorner(int x, int y, T omega,
                   std::list<BoundaryType*>& bdMomentaList)
{
    BoundaryType* boundaryMomenta = new BoundaryType;

    BdDynamics* dynamics = new BdDynamics (
            omega,
            *boundaryMomenta,
            this->getBlock().getStatistics()
    );
    this->getBlock().defineDynamics(x,x,y,y, dynamics);

    bdMomentaList.push_back(boundaryMomenta);
    bdDynamicsArrayList.push_back(dynamics);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addVelocityBoundary0N (
    int x0, int x1, int y0, int y1, T omega )
{
    addStraightBoundary<LeftVelBdGen>(x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addVelocityBoundary0P (
    int x0, int x1, int y0, int y1, T omega )
{
    addStraightBoundary<RightVelBdGen>(x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addVelocityBoundary1N (
    int x0, int x1, int y0, int y1, T omega )
{
    addStraightBoundary<LowerVelBdGen>(x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addVelocityBoundary1P (
    int x0, int x1, int y0, int y1, T omega )
{
    addStraightBoundary<UpperVelBdGen>(x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addPressureBoundary0N (
    int x0, int x1, int y0, int y1, T omega )
{
    addStraightBoundary<LeftPresBdGen>(x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addPressureBoundary0P (
    int x0, int x1, int y0, int y1, T omega )
{
    addStraightBoundary<RightPresBdGen>(x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addPressureBoundary1N (
    int x0, int x1, int y0, int y1, T omega )
{
    addStraightBoundary<LowerPresBdGen>(x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addPressureBoundary1P (
    int x0, int x1, int y0, int y1, T omega )
{
    addStraightBoundary<UpperPresBdGen>(x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNN(int x, int y, T omega)
{
    addNonLocalCorner<NN_ExternalCornerGen>(x, y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::
   addExternalVelocityCornerNP(int x, int y, T omega)
{
    addNonLocalCorner<NP_ExternalCornerGen>(x, y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPN(int x, int y, T omega)
{
    addNonLocalCorner<PN_ExternalCornerGen>(x, y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPP(int x, int y, T omega)
{
    addNonLocalCorner<PP_ExternalCornerGen>(x, y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNN(int x, int y, T omega)
{
    addLocalCorner(x, y, omega, nnInternalBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addInternalVelocityCornerNP(
    int x, int y, T omega)
{
    addLocalCorner(x, y, omega, npInternalBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addInternalVelocityCornerPN(
    int x, int y, T omega)
{
    addLocalCorner(x, y, omega, pnInternalBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition2D<T,Lattice,Dynamics>::addInternalVelocityCornerPP(
    int x, int y, T omega)
{
    addLocalCorner(x, y, omega, ppInternalBdList);
}


///////// class LocalBoundaryCondition2D ////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
LocalBoundaryCondition2D<T,Lattice,Dynamics>::LocalBoundaryCondition2D (
        BlockStructure2D<T,Lattice>& lattice)
    : BoundaryCondition2D<T,Lattice>(lattice),
      interpBoundary(lattice)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics>
LocalBoundaryCondition2D<T,Lattice,Dynamics>::~LocalBoundaryCondition2D()
{
    cleanArrayList(bdDynamicsArrayList);

    cleanList(leftVelBdList);
    cleanList(rightVelBdList);
    cleanList(lowerVelBdList);
    cleanList(upperVelBdList);
    cleanList(leftPresBdList);
    cleanList(rightPresBdList);
    cleanList(lowerPresBdList);
    cleanList(upperPresBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename ListT>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::cleanList(ListT& list) {
    typename ListT::iterator iL = list.begin();
    for (; iL != list.end(); ++iL) {
        delete *iL;
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename ListT>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::cleanArrayList(ListT& list) {
    typename ListT::iterator iOuter = list.begin();
    for (; iOuter != list.end(); ++iOuter) {
        for (unsigned iInner=0; iInner < (**iOuter).size(); ++iInner) {
           delete (**iOuter)[iInner];
        }
        delete *iOuter;
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename ArrayType>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addGenericBoundary(int x0, int x1, int y0, int y1, T omega,
                       std::list<ArrayType*>& bdMomentaList)
{
    OLB_PRECONDITION(x0==x1 || y0==y1);
    int l = std::max(x1-x0, y1-y0) + 1;

    ArrayType* bdMomenta = new ArrayType;
    bdMomenta -> resize(l);

    BdDynamicsArray* dynamicsArray = new BdDynamicsArray;
    dynamicsArray -> resize(l);

    int iDyn=0;
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            BdDynamics* newDyn = new BdDynamics (
                omega,
                (*bdMomenta)[iDyn],
                this->getBlock().getStatistics()
            );
            this->getBlock().defineDynamics(iX,iX,iY,iY, newDyn);
            (*dynamicsArray)[iDyn] = newDyn;
            ++iDyn;
        }
    }

    bdMomentaList.push_back(bdMomenta);
    bdDynamicsArrayList.push_back(dynamicsArray);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addVelocityBoundary0N(int x0, int x1, int y0, int y1, T omega)
{
    addGenericBoundary(x0, x1, y0, y1, omega, leftVelBdList);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addVelocityBoundary0P(int x0, int x1, int y0, int y1, T omega)
{
    addGenericBoundary(x0, x1, y0, y1, omega, rightVelBdList);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addVelocityBoundary1N(int x0, int x1, int y0, int y1, T omega)
{
    addGenericBoundary(x0, x1, y0, y1, omega, lowerVelBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addVelocityBoundary1P(int x0, int x1, int y0, int y1, T omega)
{
    addGenericBoundary(x0, x1, y0, y1, omega, upperVelBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addPressureBoundary0N(int x0, int x1, int y0, int y1, T omega)
{
    addGenericBoundary(x0, x1, y0, y1, omega, leftPresBdList);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addPressureBoundary0P(int x0, int x1, int y0, int y1, T omega)
{
    addGenericBoundary(x0, x1, y0, y1, omega, rightPresBdList);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addPressureBoundary1N(int x0, int x1, int y0, int y1, T omega)
{
    addGenericBoundary(x0, x1, y0, y1, omega, lowerPresBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addPressureBoundary1P(int x0, int x1, int y0, int y1, T omega)
{
    addGenericBoundary(x0, x1, y0, y1, omega, upperPresBdList);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNN(int x, int y, T omega)
{
    interpBoundary.addExternalVelocityCornerNN(x,y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNP(int x, int y, T omega)
{
    interpBoundary.addExternalVelocityCornerNP(x,y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPN(int x, int y, T omega)
{
    interpBoundary.addExternalVelocityCornerPN(x,y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPP(int x, int y, T omega)
{
    interpBoundary.addExternalVelocityCornerPP(x,y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNN(int x, int y, T omega)
{
    interpBoundary.addInternalVelocityCornerNN(x,y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNP(int x, int y, T omega)
{
    interpBoundary.addInternalVelocityCornerNP(x,y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPN(int x, int y, T omega)
{
    interpBoundary.addInternalVelocityCornerPN(x,y, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition2D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPP(int x, int y, T omega)
{
    interpBoundary.addInternalVelocityCornerPP(x,y, omega);
}


////////////////// Factory functions //////////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
BoundaryCondition2D<T,Lattice>*
    createLocalBoundaryCondition2D(BlockStructure2D<T,Lattice>& block)
{
    return new LocalBoundaryCondition2D<T,Lattice,Dynamics>(block);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
BoundaryCondition2D<T,Lattice>*
createInterpBoundaryCondition2D(BlockStructure2D<T,Lattice>& block)
{
    return new InterpBoundaryCondition2D<T,Lattice,Dynamics>(block);
}


}  // namespace olb

#endif

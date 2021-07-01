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
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef DYNAMICS_HELPERS_3D_HH
#define DYNAMICS_HELPERS_3D_HH

#include "dynamicsHelpers3D.h"
#include "cell.h"
#include "blockLattice3D.h"

namespace olb {

///////// class InterpBoundaryCondition3D (declaration)  //////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
class InterpBoundaryCondition3D:
    public BoundaryCondition3D<T,Lattice>
{
public:
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,VelocityBM, 0,-1>
        VelBd_0N;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,VelocityBM, 0, 1>
        VelBd_0P;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,VelocityBM, 1,-1>
        VelBd_1N;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,VelocityBM, 1, 1>
        VelBd_1P;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,VelocityBM, 2,-1>
        VelBd_2N;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,VelocityBM, 2, 1>
        VelBd_2P;

    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,PressureBM, 0,-1>
        PressureBd_0N;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,PressureBM, 0, 1>
        PressureBd_0P;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,PressureBM, 1,-1>
        PressureBd_1N;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,PressureBM, 1, 1>
        PressureBd_1P;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,PressureBM, 2,-1>
        PressureBd_2N;
    typedef PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,PressureBM, 2, 1>
        PressureBd_2P;

    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 0, -1,-1>
        ExternalEdge_0NN;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 0, -1, 1>
        ExternalEdge_0NP;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 0,  1,-1>
        ExternalEdge_0PN;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 0,  1, 1>
        ExternalEdge_0PP;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 1, -1,-1>
        ExternalEdge_1NN;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 1, -1, 1>
        ExternalEdge_1NP;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 1,  1,-1>
        ExternalEdge_1PN;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 1,  1, 1>
        ExternalEdge_1PP;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 2, -1,-1>
        ExternalEdge_2NN;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 2, -1, 1>
        ExternalEdge_2NP;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 2,  1,-1>
        ExternalEdge_2PN;
    typedef ConvexVelocityEdgeGenerator3D<T, Lattice,Dynamics, 2,  1, 1>
        ExternalEdge_2PP;

    typedef ConcaveEdgeVelBM3D<T, Lattice, 0, -1,-1> InternalEdge_0NN;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 0, -1, 1> InternalEdge_0NP;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 0,  1,-1> InternalEdge_0PN;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 0,  1, 1> InternalEdge_0PP;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 1, -1,-1> InternalEdge_1NN;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 1, -1, 1> InternalEdge_1NP;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 1,  1,-1> InternalEdge_1PN;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 1,  1, 1> InternalEdge_1PP;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 2, -1,-1> InternalEdge_2NN;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 2, -1, 1> InternalEdge_2NP;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 2,  1,-1> InternalEdge_2PN;
    typedef ConcaveEdgeVelBM3D<T, Lattice, 2,  1, 1> InternalEdge_2PP;

    typedef std::vector<InternalEdge_0NN> InternalEdgeArray_0NN;
    typedef std::vector<InternalEdge_0NP> InternalEdgeArray_0NP;
    typedef std::vector<InternalEdge_0PN> InternalEdgeArray_0PN;
    typedef std::vector<InternalEdge_0PP> InternalEdgeArray_0PP;
    typedef std::vector<InternalEdge_1NN> InternalEdgeArray_1NN;
    typedef std::vector<InternalEdge_1NP> InternalEdgeArray_1NP;
    typedef std::vector<InternalEdge_1PN> InternalEdgeArray_1PN;
    typedef std::vector<InternalEdge_1PP> InternalEdgeArray_1PP;
    typedef std::vector<InternalEdge_2NN> InternalEdgeArray_2NN;
    typedef std::vector<InternalEdge_2NP> InternalEdgeArray_2NP;
    typedef std::vector<InternalEdge_2PN> InternalEdgeArray_2PN;
    typedef std::vector<InternalEdge_2PP> InternalEdgeArray_2PP;

    typedef ConvexVelocityCornerGenerator3D<T, Lattice,Dynamics, -1,-1,-1>
        NNN_ExternalCorner;
    typedef ConvexVelocityCornerGenerator3D<T, Lattice,Dynamics, -1,-1, 1>
        NNP_ExternalCorner;
    typedef ConvexVelocityCornerGenerator3D<T, Lattice,Dynamics, -1, 1,-1>
        NPN_ExternalCorner;
    typedef ConvexVelocityCornerGenerator3D<T, Lattice,Dynamics, -1, 1, 1>
        NPP_ExternalCorner;
    typedef ConvexVelocityCornerGenerator3D<T, Lattice,Dynamics,  1,-1,-1>
        PNN_ExternalCorner;
    typedef ConvexVelocityCornerGenerator3D<T, Lattice,Dynamics,  1,-1, 1>
        PNP_ExternalCorner;
    typedef ConvexVelocityCornerGenerator3D<T, Lattice,Dynamics,  1, 1,-1>
        PPN_ExternalCorner;
    typedef ConvexVelocityCornerGenerator3D<T, Lattice,Dynamics,  1, 1, 1>
        PPP_ExternalCorner;

    typedef ConcaveCornerVelBM3D<T, Lattice, -1,-1,-1> NNN_InternalCorner;
    typedef ConcaveCornerVelBM3D<T, Lattice, -1,-1, 1> NNP_InternalCorner;
    typedef ConcaveCornerVelBM3D<T, Lattice, -1, 1,-1> NPN_InternalCorner;
    typedef ConcaveCornerVelBM3D<T, Lattice, -1, 1, 1> NPP_InternalCorner;
    typedef ConcaveCornerVelBM3D<T, Lattice,  1,-1,-1> PNN_InternalCorner;
    typedef ConcaveCornerVelBM3D<T, Lattice,  1,-1, 1> PNP_InternalCorner;
    typedef ConcaveCornerVelBM3D<T, Lattice,  1, 1,-1> PPN_InternalCorner;
    typedef ConcaveCornerVelBM3D<T, Lattice,  1, 1, 1> PPP_InternalCorner;

    typedef Dynamics                 BdDynamics;
    typedef std::vector<BdDynamics*> BdDynamicsArray;
public:
    InterpBoundaryCondition3D(BlockStructure3D<T,Lattice>& block);
    ~InterpBoundaryCondition3D();

    virtual void addVelocityBoundary0N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary0P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary1N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary1P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary2N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary2P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);

    virtual void addPressureBoundary0N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary0P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary1N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary1P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary2N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary2P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);

    virtual void addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);

    virtual void addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);

    virtual void addExternalVelocityCornerNNN(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerNNP(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerNPN(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerNPP(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerPNN(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerPNP(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerPPN(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerPPP(int x, int y, int z, T omega);

    virtual void addInternalVelocityCornerNNN(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerNNP(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerNPN(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerNPP(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerPNN(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerPNP(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerPPN(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerPPP(int x, int y, int z, T omega);

private:
    template<typename ListT>
      static void cleanList(ListT& list);

    template<typename ListT>
      static void cleanArrayList(ListT& list);

    template<typename BoundaryGenerator>
      void addNonLocalBoundary(int x0, int x1, int y0, int y1,
                               int z0, int z1, T omega);
    template<typename BoundaryGenerator>
      void addNonLocalCorner(int x, int y, int z, T omega);
    template<typename BoundaryArrayType>
      void addLocalEdge(int x0, int x1, int y0, int y1, int z0, int z1,
                        T omega, std::list<BoundaryArrayType*>& bdArrayList);
    template<typename BoundaryType>
      void addLocalCorner(int x, int y, int z, T omega,
                          std::list<BoundaryType*>& bdList);
private:
    std::list<InternalEdgeArray_0NN*> internalEdgeBdList0NN;
    std::list<InternalEdgeArray_0NP*> internalEdgeBdList0NP;
    std::list<InternalEdgeArray_0PN*> internalEdgeBdList0PN;
    std::list<InternalEdgeArray_0PP*> internalEdgeBdList0PP;
    std::list<InternalEdgeArray_1NN*> internalEdgeBdList1NN;
    std::list<InternalEdgeArray_1NP*> internalEdgeBdList1NP;
    std::list<InternalEdgeArray_1PN*> internalEdgeBdList1PN;
    std::list<InternalEdgeArray_1PP*> internalEdgeBdList1PP;
    std::list<InternalEdgeArray_2NN*> internalEdgeBdList2NN;
    std::list<InternalEdgeArray_2NP*> internalEdgeBdList2NP;
    std::list<InternalEdgeArray_2PN*> internalEdgeBdList2PN;
    std::list<InternalEdgeArray_2PP*> internalEdgeBdList2PP;

    std::list<BdDynamicsArray*>       bdDynamicsArrayList;

    std::list<NNN_InternalCorner*> internalCornerBdListNNN;
    std::list<NNP_InternalCorner*> internalCornerBdListNNP;
    std::list<NPN_InternalCorner*> internalCornerBdListNPN;
    std::list<NPP_InternalCorner*> internalCornerBdListNPP;
    std::list<PNN_InternalCorner*> internalCornerBdListPNN;
    std::list<PNP_InternalCorner*> internalCornerBdListPNP;
    std::list<PPN_InternalCorner*> internalCornerBdListPPN;
    std::list<PPP_InternalCorner*> internalCornerBdListPPP;

    std::list<BdDynamics*>         bdDynamicsList;
};


///////// class LocalBoundaryCondition3D (declaration)  //////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
class LocalBoundaryCondition3D : public BoundaryCondition3D<T,Lattice> {
public:
    typedef RegularizedVelocityBM<T,Lattice, 0,-1> VelBdMomenta_0N;
    typedef RegularizedVelocityBM<T,Lattice, 0, 1> VelBdMomenta_0P;
    typedef RegularizedVelocityBM<T,Lattice, 1,-1> VelBdMomenta_1N;
    typedef RegularizedVelocityBM<T,Lattice, 1, 1> VelBdMomenta_1P;
    typedef RegularizedVelocityBM<T,Lattice, 2,-1> VelBdMomenta_2N;
    typedef RegularizedVelocityBM<T,Lattice, 2, 1> VelBdMomenta_2P;

    typedef RegularizedPressureBM<T,Lattice, 0,-1> PressureBdMomenta_0N;
    typedef RegularizedPressureBM<T,Lattice, 0, 1> PressureBdMomenta_0P;
    typedef RegularizedPressureBM<T,Lattice, 1,-1> PressureBdMomenta_1N;
    typedef RegularizedPressureBM<T,Lattice, 1, 1> PressureBdMomenta_1P;
    typedef RegularizedPressureBM<T,Lattice, 2,-1> PressureBdMomenta_2N;
    typedef RegularizedPressureBM<T,Lattice, 2, 1> PressureBdMomenta_2P;

    typedef std::vector<VelBdMomenta_0N> VelBdMomentaArray_0N;
    typedef std::vector<VelBdMomenta_0P> VelBdMomentaArray_0P;
    typedef std::vector<VelBdMomenta_1N> VelBdMomentaArray_1N;
    typedef std::vector<VelBdMomenta_1P> VelBdMomentaArray_1P;
    typedef std::vector<VelBdMomenta_2N> VelBdMomentaArray_2N;
    typedef std::vector<VelBdMomenta_2P> VelBdMomentaArray_2P;

    typedef std::vector<PressureBdMomenta_0N> PressureBdMomentaArray_0N;
    typedef std::vector<PressureBdMomenta_0P> PressureBdMomentaArray_0P;
    typedef std::vector<PressureBdMomenta_1N> PressureBdMomentaArray_1N;
    typedef std::vector<PressureBdMomenta_1P> PressureBdMomentaArray_1P;
    typedef std::vector<PressureBdMomenta_2N> PressureBdMomentaArray_2N;
    typedef std::vector<PressureBdMomenta_2P> PressureBdMomentaArray_2P;

    typedef CombinedRLBdynamics<T, Lattice, Dynamics>  BdDynamics;
    typedef std::vector<BdDynamics*>    BdDynamicsArray;
public:
    LocalBoundaryCondition3D(BlockStructure3D<T,Lattice>& block);
    ~LocalBoundaryCondition3D();

    virtual void addVelocityBoundary0N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary0P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary1N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary1P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary2N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addVelocityBoundary2P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);

    virtual void addPressureBoundary0N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary0P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary1N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary1P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary2N(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);
    virtual void addPressureBoundary2P(int x0, int x1, int y0, int y1,
                                       int z0, int z1, T omega);

    virtual void addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);

    virtual void addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);
    virtual void addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1,
                                            int z0, int z1, T omega);

    virtual void addExternalVelocityCornerNNN(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerNNP(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerNPN(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerNPP(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerPNN(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerPNP(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerPPN(int x, int y, int z, T omega);
    virtual void addExternalVelocityCornerPPP(int x, int y, int z, T omega);

    virtual void addInternalVelocityCornerNNN(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerNNP(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerNPN(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerNPP(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerPNN(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerPNP(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerPPN(int x, int y, int z, T omega);
    virtual void addInternalVelocityCornerPPP(int x, int y, int z, T omega);

private:
    template<typename ListT>
      static void cleanList(ListT& list);
    template<typename ArrayType>
      void addGenericBoundary(int x0, int x1, int y0, int y1, int z0, int z1,
                              T omega, std::list<ArrayType*>& bdMomentaList);
private:
    std::list<VelBdMomentaArray_0N*> velBdMomentaList_0N;
    std::list<VelBdMomentaArray_0P*> velBdMomentaList_0P;
    std::list<VelBdMomentaArray_1N*> velBdMomentaList_1N;
    std::list<VelBdMomentaArray_1P*> velBdMomentaList_1P;
    std::list<VelBdMomentaArray_2N*> velBdMomentaList_2N;
    std::list<VelBdMomentaArray_2P*> velBdMomentaList_2P;

    std::list<PressureBdMomentaArray_0N*> pressureBdMomentaList_0N;
    std::list<PressureBdMomentaArray_0P*> pressureBdMomentaList_0P;
    std::list<PressureBdMomentaArray_1N*> pressureBdMomentaList_1N;
    std::list<PressureBdMomentaArray_1P*> pressureBdMomentaList_1P;
    std::list<PressureBdMomentaArray_2N*> pressureBdMomentaList_2N;
    std::list<PressureBdMomentaArray_2P*> pressureBdMomentaList_2P;

    std::list<BdDynamicsArray*>  bdDynamicsList;

    InterpBoundaryCondition3D<T,Lattice, Dynamics > interpBoundary;
};


///////// class BoundaryCondition3D ///////////////////////////////////

template<typename T, template<typename U> class Lattice>
BoundaryCondition3D<T,Lattice>::BoundaryCondition3D (
        BlockStructure3D<T,Lattice>& block_ )
    : block(block_)
{ }

template<typename T, template<typename U> class Lattice>
BlockStructure3D<T,Lattice>& BoundaryCondition3D<T,Lattice>::getBlock() {
    return block;
}

template<typename T, template<typename U> class Lattice>
BlockStructure3D<T,Lattice> const& BoundaryCondition3D<T,Lattice>::
    getBlock() const
{
    return block;
}


///////// class InterpBoundaryCondition3D ////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
InterpBoundaryCondition3D<T,Lattice,Dynamics>::InterpBoundaryCondition3D (
        BlockStructure3D<T,Lattice>& block )
    : BoundaryCondition3D<T,Lattice>(block)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics>
InterpBoundaryCondition3D<T,Lattice,Dynamics>::~InterpBoundaryCondition3D()
{
    // Clean the list of dynamics objects for internal (local) edges
    cleanArrayList(bdDynamicsArrayList);

    // Clean the list of dynamics objects for internal (local) corners
    cleanList(bdDynamicsList);

    // Clean the momenta objects for internal (local) edges.
    // Given that the momenta objects are static (not dynamically allocated),
    // it is sufficient to delete the containing vectors, and thus to invoke
    // cleanList instead of cleanArrayList
    cleanList(internalEdgeBdList0NN);
    cleanList(internalEdgeBdList0NP);
    cleanList(internalEdgeBdList0PN);
    cleanList(internalEdgeBdList0PP);
    cleanList(internalEdgeBdList1NN);
    cleanList(internalEdgeBdList1NP);
    cleanList(internalEdgeBdList1PN);
    cleanList(internalEdgeBdList1PP);
    cleanList(internalEdgeBdList2NN);
    cleanList(internalEdgeBdList2NP);
    cleanList(internalEdgeBdList2PN);
    cleanList(internalEdgeBdList2PP);

    // Clean the momenta objects for internal (local) edges.
    // Given that the momenta objects are static (not dynamically allocated),
    // it is sufficient to delete the containing vectors, and thus to invoke
    // cleanList instead of cleanArrayList
    cleanList(internalCornerBdListNNN);
    cleanList(internalCornerBdListNNP);
    cleanList(internalCornerBdListNPN);
    cleanList(internalCornerBdListNPP);
    cleanList(internalCornerBdListPNN);
    cleanList(internalCornerBdListPNP);
    cleanList(internalCornerBdListPPN);
    cleanList(internalCornerBdListPPP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename ListT>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::cleanList(ListT& list) {
    typename ListT::iterator iL = list.begin();
    for (; iL != list.end(); ++iL) {
        delete *iL;
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename ListT>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::cleanArrayList(ListT& list) {
    typename ListT::iterator iOuter = list.begin();
    for (; iOuter != list.end(); ++iOuter) {
        for (unsigned iInner=0; iInner < (**iOuter).size(); ++iInner) {
           delete (**iOuter)[iInner];
        }
        delete *iOuter;
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename BoundaryGenerator>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addNonLocalBoundary(int x0, int x1, int y0, int y1,
                        int z0, int z1, T omega)
{
    BoundaryGenerator bGenerator(x0,x1,y0,y1,z0,z1, omega);
    (this->getBlock()).addPostProcessor(bGenerator);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename BoundaryGenerator>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addNonLocalCorner(int x, int y, int z, T omega)
{
    BoundaryGenerator bGenerator(x,y,z, omega);
    (this->getBlock()).addPostProcessor(bGenerator);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename BoundaryArrayType>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addLocalEdge(int x0, int x1, int y0, int y1, int z0, int z1,
                 T omega, std::list<BoundaryArrayType*>& bdMomentaList)
{
    OLB_PRECONDITION(
            ( x0==x1 && y0==y1 ) ||
            ( x0==x1 && z0==z1 ) ||
            ( y0==y1 && z0==z1 ) );
    int l = std::max( std::max(x1-x0, y1-y0), z1-z0 ) + 1;
    BoundaryArrayType* bdMomenta = new BoundaryArrayType;
    bdMomenta -> resize(l);

    BdDynamicsArray* dynamicsArray = new BdDynamicsArray;
    dynamicsArray -> resize(l);

    int iDyn=0;
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
                BdDynamics* newDyn = new BdDynamics (
                    omega,
                    (*bdMomenta)[iDyn],
                    this->getBlock().getStatistics()
                );
                this->getBlock().defineDynamics(iX,iX,iY,iY,iZ,iZ, newDyn);
                (*dynamicsArray)[iDyn] = newDyn;
                ++iDyn;
            }
        }
    }

    bdMomentaList.push_back(bdMomenta);
    bdDynamicsArrayList.push_back(dynamicsArray);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename BoundaryType>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addLocalCorner(int x, int y, int z, T omega,
                   std::list<BoundaryType*>& bdList)
{
    BoundaryType* boundary = new BoundaryType;
    BdDynamics* dynamics = new BdDynamics (
            omega,
            *boundary,
            this->getBlock().getStatistics()
    );
    this->getBlock().defineDynamics(x,x,y,y,z,z, dynamics);

    bdList.push_back(boundary);
    bdDynamicsList.push_back(dynamics);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<VelBd_0N>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<VelBd_0P>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<VelBd_1N>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<VelBd_1P>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<VelBd_2N>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<VelBd_2P>(x0,x1,y0,y1,z0,z1, omega);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<PressureBd_0N>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<PressureBd_0P>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<PressureBd_1N>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<PressureBd_1P>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<PressureBd_2N>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addNonLocalBoundary<PressureBd_2P>(x0,x1,y0,y1,z0,z1, omega);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_0NN>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_0NP>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_0PN>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_0PP>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_1NN>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_1NP>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_1PN>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_1PP>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_2NN>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_2NP>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_2PN>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addNonLocalBoundary<ExternalEdge_2PP>(x0,x1,y0,y1,z0,z1, omega);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList0NN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList0NP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList0PN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList0PP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList1NN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList1NP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList1PN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList1PP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList2NN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList2NP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList2PN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    addLocalEdge(x0,x1,y0,y1,z0,z1, omega, internalEdgeBdList2PP);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNNN(int x, int y, int z, T omega)
{
    addNonLocalCorner<NNN_ExternalCorner>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNNP(int x, int y, int z, T omega)
{
    addNonLocalCorner<NNP_ExternalCorner>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNPN(int x, int y, int z, T omega)
{
    addNonLocalCorner<NPN_ExternalCorner>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNPP(int x, int y, int z, T omega)
{
    addNonLocalCorner<NPP_ExternalCorner>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPNN(int x, int y, int z, T omega)
{
    addNonLocalCorner<PNN_ExternalCorner>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPNP(int x, int y, int z, T omega)
{
    addNonLocalCorner<PNP_ExternalCorner>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPPN(int x, int y, int z, T omega)
{
    addNonLocalCorner<PPN_ExternalCorner>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPPP(int x, int y, int z, T omega)
{
    addNonLocalCorner<PPP_ExternalCorner>(x,y,z, omega);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNNN(int x, int y, int z, T omega)
{
    addLocalCorner(x,y,z, omega, internalCornerBdListNNN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNNP(int x, int y, int z, T omega)
{
    addLocalCorner(x,y,z, omega, internalCornerBdListNNP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNPN(int x, int y, int z, T omega)
{
    addLocalCorner(x,y,z, omega, internalCornerBdListNPN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNPP(int x, int y, int z, T omega)
{
    addLocalCorner(x,y,z, omega, internalCornerBdListNPP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPNN(int x, int y, int z, T omega)
{
    addLocalCorner(x,y,z, omega, internalCornerBdListPNN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPNP(int x, int y, int z, T omega)
{
    addLocalCorner(x,y,z, omega, internalCornerBdListPNP);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPPN(int x, int y, int z, T omega)
{
    addLocalCorner(x,y,z, omega, internalCornerBdListPPN);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void InterpBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPPP(int x, int y, int z, T omega)
{
    addLocalCorner(x,y,z, omega, internalCornerBdListPPP);
}


///////// class LocalBoundaryCondition3D ////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
LocalBoundaryCondition3D<T,Lattice,Dynamics>::LocalBoundaryCondition3D (
        BlockStructure3D<T,Lattice>& block )
    : BoundaryCondition3D<T,Lattice>(block),
      interpBoundary(block)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics>
LocalBoundaryCondition3D<T,Lattice,Dynamics>::~LocalBoundaryCondition3D()
{
    typename std::list<BdDynamicsArray*>::iterator iL
        = bdDynamicsList.begin();
    for (; iL != bdDynamicsList.end(); ++iL) {
        BdDynamicsArray& dArray = **iL;
        for (unsigned iA=0; iA<dArray.size(); ++iA) {
            delete dArray[iA];
        }
        delete *iL;
    }

    cleanList(velBdMomentaList_0N);
    cleanList(velBdMomentaList_0P);
    cleanList(velBdMomentaList_1N);
    cleanList(velBdMomentaList_1P);
    cleanList(velBdMomentaList_2N);
    cleanList(velBdMomentaList_2P);

    cleanList(pressureBdMomentaList_0N);
    cleanList(pressureBdMomentaList_0P);
    cleanList(pressureBdMomentaList_1N);
    cleanList(pressureBdMomentaList_1P);
    cleanList(pressureBdMomentaList_2N);
    cleanList(pressureBdMomentaList_2P);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename ListT>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::cleanList(ListT& list) {
    typename ListT::iterator iL = list.begin();
    for (; iL != list.end(); ++iL) {
        delete *iL;
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
template<typename BoundaryArrayType>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addGenericBoundary(int x0, int x1, int y0, int y1, int z0, int z1,
                       T omega, std::list<BoundaryArrayType*>& bdMomentaList)
{
    OLB_PRECONDITION( x0==x1 || y0==y1 || z0==z1 );

    int size = (x1-x0+1) * (y1-y0+1) * (z1-z0+1);
    BoundaryArrayType* bdMomenta = new BoundaryArrayType;
    bdMomenta -> resize(size);

    BdDynamicsArray* dynamicsArray = new BdDynamicsArray;
    dynamicsArray -> resize(size);

    int iDyn=0;
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
                BdDynamics* newDyn = new BdDynamics (
                    omega,
                    (*bdMomenta)[iDyn],
                    this->getBlock().getStatistics()
                );
                this->getBlock().defineDynamics(iX,iX,iY,iY,iZ,iZ, newDyn);
                (*dynamicsArray)[iDyn] = newDyn;
                ++iDyn;
            }
        }
    }

    bdMomentaList.push_back(bdMomenta);
    bdDynamicsList.push_back(dynamicsArray);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, velBdMomentaList_0N);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, velBdMomentaList_0P);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, velBdMomentaList_1N);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, velBdMomentaList_1P);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, velBdMomentaList_2N);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addVelocityBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, velBdMomentaList_2P);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, pressureBdMomentaList_0N);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, pressureBdMomentaList_0P);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, pressureBdMomentaList_1N);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, pressureBdMomentaList_1P);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, pressureBdMomentaList_2N);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addPressureBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1,
                          T omega)
{
    addGenericBoundary(x0,x1,y0,y1,z0,z1, omega, pressureBdMomentaList_2P);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge0NN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge0NP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge0PN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge0PP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge1NN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge1NP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge1PN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge1PP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge2NN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge2NP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge2PN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addExternalVelocityEdge2PP(x0,x1,y0,y1,z0,z1, omega);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge0NN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge0NP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge0PN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge0PP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge1NN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge1NP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge1PN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge1PP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge2NN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge2NP(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge2PN(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1,
                               T omega)
{
    interpBoundary.addInternalVelocityEdge2PP(x0,x1,y0,y1,z0,z1, omega);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNNN(int x, int y, int z, T omega)
{
    interpBoundary.addExternalVelocityCornerNNN(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNNP(int x, int y, int z, T omega)
{
    interpBoundary.addExternalVelocityCornerNNP(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNPN(int x, int y, int z, T omega)
{
    interpBoundary.addExternalVelocityCornerNPN(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerNPP(int x, int y, int z, T omega)
{
    interpBoundary.addExternalVelocityCornerNPP(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPNN(int x, int y, int z, T omega)
{
    interpBoundary.addExternalVelocityCornerPNN(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPNP(int x, int y, int z, T omega)
{
    interpBoundary.addExternalVelocityCornerPNP(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPPN(int x, int y, int z, T omega)
{
    interpBoundary.addExternalVelocityCornerPPN(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addExternalVelocityCornerPPP(int x, int y, int z, T omega)
{
    interpBoundary.addExternalVelocityCornerPPP(x,y,z, omega);
}


template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNNN(int x, int y, int z, T omega)
{
    interpBoundary.addInternalVelocityCornerNNN(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNNP(int x, int y, int z, T omega)
{
    interpBoundary.addInternalVelocityCornerNNP(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNPN(int x, int y, int z, T omega)
{
    interpBoundary.addInternalVelocityCornerNPN(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerNPP(int x, int y, int z, T omega)
{
    interpBoundary.addInternalVelocityCornerNPP(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPNN(int x, int y, int z, T omega)
{
    interpBoundary.addInternalVelocityCornerPNN(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPNP(int x, int y, int z, T omega)
{
    interpBoundary.addInternalVelocityCornerPNP(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPPN(int x, int y, int z, T omega)
{
    interpBoundary.addInternalVelocityCornerPPN(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void LocalBoundaryCondition3D<T,Lattice,Dynamics>::
    addInternalVelocityCornerPPP(int x, int y, int z, T omega)
{
    interpBoundary.addInternalVelocityCornerPPP(x,y,z, omega);
}


///////// Factory functions ///////////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
BoundaryCondition3D<T,Lattice>* createLocalBoundaryCondition3D (
        BlockStructure3D<T,Lattice>& block)
{
    return new LocalBoundaryCondition3D<T,Lattice,Dynamics>(block);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
BoundaryCondition3D<T,Lattice>* createInterpBoundaryCondition3D (
        BlockStructure3D<T,Lattice>& block)
{
    return new InterpBoundaryCondition3D<T,Lattice,Dynamics>(block);
}


}  // namespace olb

#endif

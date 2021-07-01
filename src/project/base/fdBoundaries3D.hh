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

#ifndef FD_BOUNDARIES_3D_HH
#define FD_BOUNDARIES_3D_HH

#include "fdBoundaries3D.h"
#include "finiteDifference.h"
#include "blockLattice3D.h"
#include "util.h"

namespace olb {

////////  PlaneFdBoundary3D ///////////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
PlaneFdBoundary3D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    PlaneFdBoundary3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                      T omega_, BlockLattice3D<T,Lattice>& blockLattice)
    : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
      omega(omega_),
      boundaryMomenta(0),
      boundaryDynamics(0)
{
    OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
    initialize(blockLattice);
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
PlaneFdBoundary3D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    ~PlaneFdBoundary3D()
{
    int index=0; 
    if (boundaryDynamics) {
        for (int iX=x0; iX<=x1; ++iX) {
            for (int iY=y0; iY<=y1; ++iY) {
                for (int iZ=z0; iZ<=z1; ++iZ) {
                    delete boundaryDynamics[index];
                    ++index;
                }
            }
        }
    }
    delete [] boundaryDynamics;
    delete [] boundaryMomenta;
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
void PlaneFdBoundary3D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    initialize(BlockLattice3D<T,Lattice>& blockLattice)
{
    int length=0;
    if (x0==x1) {
        int dim1 = y1-y0+1;
        int dim2 = z1-z0+1;
        length = dim1*dim2;
    }
    else if (y0==y1) {
        int dim1 = x1-x0+1;
        int dim2 = z1-z0+1;
        length = dim1*dim2;
    }
    else {
        int dim1 = x1-x0+1;
        int dim2 = y1-y0+1;
        length = dim1*dim2;
    }

    boundaryMomenta = new StressInterpolationDirichletBM <
        T,Lattice, HydroBM, direction,orientation > [length];
    boundaryDynamics = new Dynamics*[length];
    int index=0; 
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
               boundaryDynamics[index] =
                   new Dynamics (
                       omega,
                       boundaryMomenta[index],
                       blockLattice.getStatistics()
                   );
               blockLattice.defineDynamics (
                    iX,iX,iY,iY,iZ,iZ, boundaryDynamics[index] );
                ++index;
            }
        }
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
void PlaneFdBoundary3D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    process(BlockLattice3D<T,Lattice>& blockLattice,
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
    using namespace olb::util::tensorIndices3D;

    int newX0, newX1, newY0, newY1, newZ0, newZ1;
    if ( util::intersect (
                x0, x1, y0, y1, z0, z1,
                x0_, x1_, y0_, y1_, z0_, z1_,
                newX0, newX1, newY0, newY1, newZ0, newZ1 ) )
    {

        T dx_rhoU[Lattice<T>::d], dy_rhoU[Lattice<T>::d],
          dz_rhoU[Lattice<T>::d];
        // in the following sum, two of the three terms are zero.
        int index=newX0-x0 + newY0-y0 + newZ0-z0;
        for (int iX=newX0; iX<=newX1; ++iX) {
            for (int iY=newY0; iY<=newY1; ++iY) {
                for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
                    interpolateGradients<0> (
                            blockLattice, dx_rhoU, iX, iY, iZ  );
                    interpolateGradients<1> (
                            blockLattice, dy_rhoU, iX, iY, iZ );
                    interpolateGradients<2> (
                            blockLattice, dz_rhoU, iX, iY, iZ );
                    T dx_rhoUx = dx_rhoU[0];
                    T dy_rhoUx = dy_rhoU[0];
                    T dz_rhoUx = dz_rhoU[0];
                    T dx_rhoUy = dx_rhoU[1];
                    T dy_rhoUy = dy_rhoU[1];
                    T dz_rhoUy = dz_rhoU[1];
                    T dx_rhoUz = dx_rhoU[2];
                    T dy_rhoUz = dy_rhoU[2];
                    T dz_rhoUz = dz_rhoU[2];
                    T omega = boundaryDynamics[index]->getOmega();
                    T sToPi = - (T)1 / Lattice<T>::invCs2 / omega;
                    boundaryMomenta[index][xx] = (T)2 * dx_rhoUx * sToPi;
                    boundaryMomenta[index][yy] = (T)2 * dy_rhoUy * sToPi;
                    boundaryMomenta[index][zz] = (T)2 * dz_rhoUz * sToPi;
                    boundaryMomenta[index][xy] = (dx_rhoUy + dy_rhoUx)
                                                   * sToPi;
                    boundaryMomenta[index][xz] = (dx_rhoUz + dz_rhoUx)
                                                   * sToPi;
                    boundaryMomenta[index][yz] = (dy_rhoUz + dz_rhoUy)
                                                   * sToPi;

                    // Computation of the particle distribution functions
                    // according to the regularized formula; implementation
                    // by Orestis Malaspinas
                    T rho, u[Lattice<T>:: d], pi[util::TensorVal<Lattice<T> >::n];
                    blockLattice.get(iX,iY,iZ).computeRhoU(rho,u);
                    for (int iPi = 0; iPi < util::TensorVal<Lattice<T> >::n; ++iPi)
                        pi[iPi] = boundaryMomenta[index][iPi];


                    T uSqr = util::normSqr<T,Lattice<T>::d>(u);

                    for (int iPop = 0; iPop < Lattice<T>::q; ++iPop)
                        blockLattice.get(iX,iY,iZ)[iPop] =
                            lbHelpers<T,Lattice>::equilibrium(iPop,rho,u,uSqr) +
                            firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);

                    ++index;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
void PlaneFdBoundary3D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    process(BlockLattice3D<T,Lattice>& blockLattice)
{
    process(blockLattice, x0, x1, y0, y1, z0, z1);
}


template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
template<int deriveDirection>
void PlaneFdBoundary3D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    interpolateGradients(BlockLattice3D<T,Lattice> const& blockLattice,
                         T rhoUDeriv[Lattice<T>::d],
                         int iX, int iY, int iZ) const
{
    DirectedGradients3D<T, Lattice, direction, orientation, deriveDirection,
                      direction==deriveDirection>::
        interpolate(rhoUDeriv, blockLattice, iX, iY, iZ, true);
}


////////  PlaneFdBoundaryGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,
                           HydroBM,direction,orientation>::
    PlaneFdBoundaryGenerator3D(int x0_, int x1_, int y0_, int y1_,
                               int z0_, int z1_, T omega_)
    : PostProcessorGenerator3D<T,Lattice>(x0_, x1_, y0_, y1_, z0_, z1_),
      omega(omega_)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
PostProcessor3D<T,Lattice>*
    PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,
                               HydroBM,direction,orientation>::
        generate(BlockLattice3D<T,Lattice>& blockLattice) const
{
    return new PlaneFdBoundary3D<T,Lattice,Dynamics,
                                 HydroBM,direction,orientation>
                   ( this->x0, this->x1, this->y0, this->y1,
                     this->z0, this->z1, omega,
                     blockLattice );
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
    PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,
                               HydroBM,direction,orientation>::
        clone() const
{
    return new PlaneFdBoundaryGenerator3D<T,Lattice,Dynamics,
                                             HydroBM,direction,orientation>
                   (this->x0, this->x1, this->y0, this->y1,
                    this->z0, this->z1, omega);
}
 


////////  ConvexVelocityEdge3D ///////////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
ConvexVelocityEdge3D<T,Lattice,Dynamics, plane, normal1,normal2>::
    ConvexVelocityEdge3D (
        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T omega_,
        BlockLattice3D<T,Lattice>& blockLattice)
    : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), omega(omega_) 
{
    OLB_PRECONDITION (
            (plane==2 && x0==x1 && y0==y1) ||
            (plane==1 && x0==x1 && z0==z1) ||
            (plane==0 && y0==y1 && z0==z1)     );
    initialize(blockLattice);

}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
ConvexVelocityEdge3D<T,Lattice,Dynamics, plane, normal1,normal2>::
    ~ConvexVelocityEdge3D()
{
    int index=0; 
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
                delete boundaryDynamics[index];
                ++index;
            }
        }
    }
    delete [] boundaryDynamics;
    delete [] boundaryMomenta;
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
void ConvexVelocityEdge3D<T,Lattice,Dynamics, plane, normal1,normal2>::
    initialize(BlockLattice3D<T,Lattice>& blockLattice)
{
    int length=0;
    switch(plane) {
        case 0: length = x1-x0+1; break;
        case 1: length = y1-y0+1; break;
        case 2: length = z1-z0+1; break;
        default: OLB_ASSERT(false, "Plane must be within [0,2]");
    }

    boundaryMomenta  = new GenericBoundaryMomenta<T,Lattice> [length];
    boundaryDynamics = new Dynamics* [length];
    int index=0; 
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
               boundaryDynamics[index] =
                   new Dynamics (
                       omega,
                       boundaryMomenta[index],
                       blockLattice.getStatistics()
                   );
               blockLattice.defineDynamics (
                    iX,iX, iY,iY, iZ,iZ, boundaryDynamics[index] );
                ++index;
            }
        }
    }
}


template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
void ConvexVelocityEdge3D<T,Lattice,Dynamics, plane, normal1,normal2>::
    process(BlockLattice3D<T,Lattice>& blockLattice,
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
    using namespace olb::util::tensorIndices3D;

    int newX0, newX1, newY0, newY1, newZ0, newZ1;
    if ( util::intersect (
                x0, x1, y0, y1, z0, z1,
                x0_, x1_, y0_, y1_, z0_, z1_,
                newX0, newX1, newY0, newY1, newZ0, newZ1 ) )
    {

        // in the following sum, two of the three terms are zero.
        int index=newX0-x0 + newY0-y0 + newZ0-z0;
        for (int iX=newX0; iX<=newX1; ++iX) {
            for (int iY=newY0; iY<=newY1; ++iY) {
                for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
                    T dA_rhouB_[3][3];
                    interpolateGradients<plane,0> (
                            blockLattice, dA_rhouB_[0], iX, iY, iZ );
                    interpolateGradients<direction1,normal1> (
                            blockLattice, dA_rhouB_[1], iX, iY, iZ );
                    interpolateGradients<direction2,normal2> (
                            blockLattice, dA_rhouB_[2], iX, iY, iZ );
                    T dA_rhouB[3][3];
                    for (int iBeta=0; iBeta<3; ++iBeta) {
                        dA_rhouB[plane][iBeta]
                            = dA_rhouB_[0][iBeta];
                        dA_rhouB[direction1][iBeta]
                            = dA_rhouB_[1][iBeta];
                        dA_rhouB[direction2][iBeta]
                            = dA_rhouB_[2][iBeta];
                    }
                    T omega = boundaryDynamics[index]->getOmega();
                    T sToPi = - (T)1 / Lattice<T>::invCs2 / omega;
                    T pi[util::TensorVal<Lattice<T> >::n];
                    pi[xx] = (T)2 * dA_rhouB[0][0] * sToPi;
                    pi[yy] = (T)2 * dA_rhouB[1][1] * sToPi;
                    pi[zz] = (T)2 * dA_rhouB[2][2] * sToPi;
                    pi[xy] = (dA_rhouB[0][1]+dA_rhouB[1][0]) * sToPi;
                    pi[xz] = (dA_rhouB[0][2]+dA_rhouB[2][0]) * sToPi;
                    pi[yz] = (dA_rhouB[1][2]+dA_rhouB[2][1]) * sToPi;
                    boundaryMomenta[index].defineStress(pi);

                    T rho10 = getNeighborRho(iX,iY,iZ,1,0, blockLattice);
                    T rho01 = getNeighborRho(iX,iY,iZ,0,1, blockLattice);
                    T rho20 = getNeighborRho(iX,iY,iZ,2,0, blockLattice);
                    T rho02 = getNeighborRho(iX,iY,iZ,0,2, blockLattice);
                    T rho = (T)2/(T)3*(rho01+rho10)-(T)1/(T)6*(rho02+rho20);
                    boundaryMomenta[index].defineRho(rho);


                    // Computation of the particle distribution functions
                    // according to the regularized formula; implementation
                    // by Orestis Malaspinas
                    T u[Lattice<T>::d];
                    boundaryMomenta[index].computeU(u);
                    T uSqr = util::normSqr<T,Lattice<T>::d>(u);

                    for (int iPop = 0; iPop < Lattice<T>::q; ++iPop)
                        blockLattice.get(iX,iY,iZ)[iPop] =
                            lbHelpers<T,Lattice>::equilibrium(iPop,rho,u,uSqr) +
                            firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);

                    ++index;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
void ConvexVelocityEdge3D<T,Lattice,Dynamics, plane, normal1,normal2>::
    process(BlockLattice3D<T,Lattice>& blockLattice)
{
    process(blockLattice, x0, x1, y0, y1, z0, z1);
}


template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
T ConvexVelocityEdge3D<T,Lattice,Dynamics, plane, normal1,normal2>::
    getNeighborRho(int x, int y, int z, int step1, int step2,
                   BlockLattice3D<T,Lattice> const& blockLattice)
{
    int coords[3] = {x, y, z};
    coords[direction1] += -normal1*step1;
    coords[direction2] += -normal2*step2;
    return blockLattice.get(coords[0], coords[1], coords[2]).computeRho();
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
template<int deriveDirection, int orientation>
void ConvexVelocityEdge3D<T,Lattice,Dynamics, plane, normal1,normal2>::
    interpolateGradients(BlockLattice3D<T,Lattice> const& blockLattice,
                         T rhoUDeriv[Lattice<T>::d],
                         int iX, int iY, int iZ) const
{
    DirectedGradients3D<T, Lattice, deriveDirection, orientation,
                        deriveDirection, deriveDirection!=plane>::
        interpolate(rhoUDeriv, blockLattice, iX, iY, iZ, true);
}

////////  ConvexVelocityEdgeGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
ConvexVelocityEdgeGenerator3D<T,Lattice,Dynamics, plane,normal1,normal2>::
    ConvexVelocityEdgeGenerator3D(int x0_, int x1_, int y0_, int y1_,
                                  int z0_, int z1_, T omega_)
    : PostProcessorGenerator3D<T,Lattice>(x0_, x1_, y0_, y1_, z0_, z1_),
      omega(omega_)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
PostProcessor3D<T,Lattice>*
ConvexVelocityEdgeGenerator3D<T,Lattice,Dynamics, plane,normal1,normal2>::
        generate(BlockLattice3D<T,Lattice>& blockLattice) const
{
    return new ConvexVelocityEdge3D <
                   T,Lattice,Dynamics, plane,normal1,normal2 >
                       ( this->x0, this->x1, this->y0, this->y1,
                         this->z0, this->z1, omega,
                         blockLattice );
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>*
ConvexVelocityEdgeGenerator3D<T,Lattice,Dynamics, plane,normal1,normal2>::
        clone() const
{
    return new ConvexVelocityEdgeGenerator3D<T,Lattice,Dynamics,
                                             plane, normal1, normal2 >
                   (this->x0, this->x1, this->y0, this->y1,
                    this->z0, this->z1, omega);
}


////////// DirectedGradients3D /////////////////////////////////////////////

// Implementation for orthogonal==true; i.e. the derivative is along
// the boundary normal.
template<typename T, template<typename U> class Lattice,
         int direction, int orientation, int deriveDirection>
struct DirectedGradients3D<T, Lattice, direction, orientation,
                           deriveDirection, true>
{
    static void interpolate(T velDeriv[Lattice<T>::d],
                            BlockLattice3D<T,Lattice> const& blockLattice,
                            int iX, int iY, int iZ, bool multiplyRho)
    {
        using namespace fd;

        T rho0, rho1, rho2;
        T u0[Lattice<T>::d], u1[Lattice<T>::d], u2[Lattice<T>::d];
        
        // note that the derivative runs along direction.
        blockLattice.get(iX,iY,iZ).computeRhoU(rho0, u0);
        blockLattice.get (
            iX+(direction==0 ? (-orientation):0),
            iY+(direction==1 ? (-orientation):0),
            iZ+(direction==2 ? (-orientation):0)  ).computeRhoU(rho1, u1);
        blockLattice.get (
            iX+(direction==0 ? (-2*orientation):0),
            iY+(direction==1 ? (-2*orientation):0),
            iZ+(direction==2 ? (-2*orientation):0) ).computeRhoU(rho2, u2);

        for (int iD=0; iD<Lattice<T>::d; ++iD) {
            if (multiplyRho) {
                velDeriv[iD] = -orientation *
                    boundaryGradient(rho0*u0[iD], rho1*u1[iD], rho2*u2[iD]);
            }
            else {
                velDeriv[iD] = -orientation *
                    boundaryGradient(u0[iD], u1[iD], u2[iD]);
            }
        }
    }
};

// Implementation for orthogonal==false; i.e. the derivative is aligned
// with the boundary.
template<typename T, template<typename U> class Lattice,
         int direction, int orientation, int deriveDirection>
struct DirectedGradients3D<T, Lattice, direction, orientation,
                           deriveDirection, false>
{
    static void  interpolate(T velDeriv[Lattice<T>::d],
                             BlockLattice3D<T,Lattice> const& blockLattice,
                             int iX, int iY, int iZ, bool multiplyRho)
    {
        using namespace fd;

        T rho_p1, rho_m1;
        T u_p1[Lattice<T>::d], u_m1[Lattice<T>::d];
        
        blockLattice.get (
            iX+(deriveDirection==0 ? 1:0),
            iY+(deriveDirection==1 ? 1:0),
            iZ+(deriveDirection==2 ? 1:0) ).computeRhoU(rho_p1, u_p1);

        blockLattice.get (
            iX+(deriveDirection==0 ? (-1):0),
            iY+(deriveDirection==1 ? (-1):0),
            iZ+(deriveDirection==2 ? (-1):0) ).computeRhoU(rho_m1, u_m1);

        for (int iD=0; iD<Lattice<T>::d; ++iD) {
            if (multiplyRho) {
                velDeriv[iD] = fd::centralGradient(rho_p1*u_p1[iD],
                                                   rho_m1*u_m1[iD]);
            }
            else {
                velDeriv[iD] = fd::centralGradient(u_p1[iD],u_m1[iD]);
            }
        }
    }
};

/////////// ConvexVelocityCorner3D /////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
ConvexVelocityCorner3D<T, Lattice,Dynamics, xNormal, yNormal, zNormal>::
    ConvexVelocityCorner3D ( int x_, int y_, int z_, T omega_,
                             BlockLattice3D<T,Lattice>& blockLattice )
    : x(x_), y(y_), z(z_),
      omega(omega_),
      boundaryMomenta(0),
      boundaryDynamics(0)
{
    initialize(blockLattice);
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
ConvexVelocityCorner3D<T, Lattice,Dynamics, xNormal, yNormal, zNormal>::
    ~ConvexVelocityCorner3D()
{
    delete boundaryDynamics;
    delete boundaryMomenta;
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
void ConvexVelocityCorner3D<T, Lattice,Dynamics, xNormal, yNormal, zNormal>::
    initialize(BlockLattice3D<T,Lattice>& blockLattice_)
{
     boundaryMomenta = new GenericBoundaryMomenta<T,Lattice>;
     boundaryDynamics = new Dynamics (
          omega,
          *boundaryMomenta,
          blockLattice_.getStatistics() );
    blockLattice_.defineDynamics (
        x, x, y, y, z,z, boundaryDynamics
    );
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
void ConvexVelocityCorner3D<T, Lattice,Dynamics, xNormal, yNormal, zNormal>::
    process(BlockLattice3D<T,Lattice>& blockLattice)
{
    using namespace olb::util::tensorIndices3D;
    T dx_rhoU[Lattice<T>::d], dy_rhoU[Lattice<T>::d], dz_rhoU[Lattice<T>::d];
    DirectedGradients3D<T, Lattice, 0, xNormal, 0, true>::
        interpolate(dx_rhoU, blockLattice, x,y,z, true);
    DirectedGradients3D<T, Lattice, 1, yNormal, 0, true>::
        interpolate(dy_rhoU, blockLattice, x,y,z, true);
    DirectedGradients3D<T, Lattice, 2, zNormal, 0, true>::
        interpolate(dz_rhoU, blockLattice, x,y,z, true);

    T dx_rhoUx = dx_rhoU[0];
    T dy_rhoUx = dy_rhoU[0];
    T dz_rhoUx = dz_rhoU[0];
    T dx_rhoUy = dx_rhoU[1];
    T dy_rhoUy = dy_rhoU[1];
    T dz_rhoUy = dz_rhoU[1];
    T dx_rhoUz = dx_rhoU[2];
    T dy_rhoUz = dy_rhoU[2];
    T dz_rhoUz = dz_rhoU[2];
    T omega = boundaryDynamics->getOmega();
    T sToPi = - (T)1 / Lattice<T>::invCs2 / omega;
    T pi[util::TensorVal<Lattice<T> >::n];
    pi[xx] = (T)2 * dx_rhoUx * sToPi;
    pi[yy] = (T)2 * dy_rhoUy * sToPi;
    pi[zz] = (T)2 * dz_rhoUz * sToPi;
    pi[xy] = (dx_rhoUy + dy_rhoUx) * sToPi;
    pi[xz] = (dx_rhoUz + dz_rhoUx) * sToPi;
    pi[yz] = (dy_rhoUz + dz_rhoUy) * sToPi;
    boundaryMomenta->defineStress(pi);

    T rho100 = blockLattice.get(x - 1*xNormal,
                                y - 0*yNormal,
                                z - 0*zNormal).computeRho();
    T rho010 = blockLattice.get(x - 0*xNormal,
                                y - 1*yNormal,
                                z - 0*zNormal).computeRho();
    T rho001 = blockLattice.get(x - 0*xNormal,
                                y - 0*yNormal,
                                z - 1*zNormal).computeRho();
    T rho200 = blockLattice.get(x - 2*xNormal,
                                y - 0*yNormal,
                                z - 0*zNormal).computeRho();
    T rho020 = blockLattice.get(x - 0*xNormal,
                                y - 2*yNormal,
                                z - 0*zNormal).computeRho();
    T rho002 = blockLattice.get(x - 0*xNormal,
                                y - 0*yNormal,
                                z - 2*zNormal).computeRho();

    T rho = (T)4/(T)9 * (rho001 + rho010 + rho100) -
            (T)1/(T)9 * (rho002 + rho020 + rho200);

    boundaryMomenta->defineRho(rho);

    // Computation of the particle distribution functions
    // according to the regularized formula; implementation
    // by Orestis Malaspinas
    T u[Lattice<T>::d];
    boundaryMomenta->computeU(u);
    T uSqr = util::normSqr<T,Lattice<T>::d>(u);

    for (int iPop = 0; iPop < Lattice<T>::q; ++iPop)
        blockLattice.get(x,y,z)[iPop] =
            lbHelpers<T,Lattice>::equilibrium(iPop,rho,u,uSqr) +
            firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);

}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
void ConvexVelocityCorner3D<T, Lattice,Dynamics, xNormal, yNormal, zNormal>::
    process(BlockLattice3D<T,Lattice>& blockLattice,
            int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
    if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
        process(blockLattice);
    }
}

////////  ConvexVelocityCornerGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
ConvexVelocityCornerGenerator3D<T,Lattice,Dynamics,
                                xNormal,yNormal,zNormal>::
    ConvexVelocityCornerGenerator3D(int x_, int y_, int z_, T omega_)
    : PostProcessorGenerator3D<T,Lattice>(x_,x_, y_,y_, z_,z_),
      omega(omega_)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
PostProcessor3D<T,Lattice>*
    ConvexVelocityCornerGenerator3D<T,Lattice,Dynamics,
                                    xNormal,yNormal,zNormal>::
        generate(BlockLattice3D<T,Lattice>& blockLattice) const
{
    return new ConvexVelocityCorner3D<T,Lattice,Dynamics,
                                      xNormal,yNormal,zNormal>
                   ( this->x0, this->y0, this->z0, omega, blockLattice );
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>*
    ConvexVelocityCornerGenerator3D<T,Lattice,Dynamics,
                                    xNormal,yNormal,zNormal>::
        clone() const
{
    return new ConvexVelocityCornerGenerator3D<T,Lattice,Dynamics,
                                               xNormal, yNormal, zNormal>
                   (this->x0, this->y0, this->z0, omega);
}


}  // namespace olb

#endif

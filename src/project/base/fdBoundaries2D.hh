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

#ifndef FD_BOUNDARIES_2D_HH
#define FD_BOUNDARIES_2D_HH

#include "fdBoundaries2D.h"
#include "finiteDifference.h"
#include "blockLattice2D.h"
#include "util.h"
#include "lbHelpers.h"
#include "firstOrderLbHelpers.h"

namespace olb {

///////////  StraightFdBoundary2D ///////////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
StraightFdBoundary2D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    StraightFdBoundary2D(int x0_, int x1_, int y0_, int y1_, T omega_,
                         BlockLattice2D<T,Lattice>& blockLattice)
    : x0(x0_), x1(x1_), y0(y0_), y1(y1_), 
      omega(omega_),
      boundaryMomenta(0),
      boundaryDynamics(0)
 
{
    OLB_PRECONDITION(x0==x1 || y0==y1);
    initialize(blockLattice);
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
void StraightFdBoundary2D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    initialize(BlockLattice2D<T,Lattice>& blockLattice)
{
    int length = x0==x1 ? y1-y0+1 : x1-x0+1;
    boundaryMomenta =
        new StressInterpolationDirichletBM< T,Lattice, HydroBM,
                                            direction,orientation > [length];
    boundaryDynamics = new Dynamics*[length];
    int index=0; 
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
             boundaryDynamics[index] =
                 new Dynamics (
                     omega,
                     boundaryMomenta[index],
                     blockLattice.getStatistics()
                 );
             blockLattice.defineDynamics (
                  iX,iX,iY,iY, boundaryDynamics[index] );
              ++index;
        }
    }
}


template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
StraightFdBoundary2D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    ~StraightFdBoundary2D()
{
    int index=0; 
    if (boundaryDynamics) {
        for (int iX=x0; iX<=x1; ++iX) {
            for (int iY=y0; iY<=y1; ++iY) {
                delete boundaryDynamics[index];
                ++index;
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
void StraightFdBoundary2D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                     int x0_, int x1_, int y0_, int y1_)
{
    using namespace olb::util::tensorIndices2D;

    int newX0, newX1, newY0, newY1;
    if ( util::intersect (
                x0, x1, y0, y1,
                x0_, x1_, y0_, y1_,
                newX0, newX1, newY0, newY1 ) )
    {
        T dx_rhoU[Lattice<T>::d], dy_rhoU[Lattice<T>::d];
        // in the following sum, one of the two terms is zero.
        int index=newX0-x0 + newY0-y0;
        for (int iX=newX0; iX<=newX1; ++iX) {
            for (int iY=newY0; iY<=newY1; ++iY) {
                interpolateGradients<0>(blockLattice, dx_rhoU, iX, iY);
                interpolateGradients<1>(blockLattice, dy_rhoU, iX, iY);
                T dx_rhoUx = dx_rhoU[0];
                T dy_rhoUx = dy_rhoU[0];
                T dx_rhoUy = dx_rhoU[1];
                T dy_rhoUy = dy_rhoU[1];
                T omega = boundaryDynamics[index]->getOmega();
                T sToPi = - (T)1 / Lattice<T>::invCs2 / omega;
                boundaryMomenta[index][xx] = (T)2 * dx_rhoUx * sToPi;
                boundaryMomenta[index][yy] = (T)2 * dy_rhoUy * sToPi;
                boundaryMomenta[index][xy] = (dx_rhoUy + dy_rhoUx) * sToPi;

                // Computation of the particle distribution functions
                // according to the regularized formula; implementation
                // by Orestis Malaspinas
                T rho, u[Lattice<T>::d],pi[util::TensorVal<Lattice<T> >::n];
                blockLattice.get(iX,iY).computeRhoU(rho,u);
                for (int iPi = 0; iPi < util::TensorVal<Lattice<T> >::n; ++iPi)
                    pi[iPi] = boundaryMomenta[index][iPi];


                T uSqr = util::normSqr<T,2>(u);

                for (int iPop = 0; iPop < Lattice<T>::q; ++iPop)
                {
                    blockLattice.get(iX,iY)[iPop] =
                        lbHelpers<T,Lattice>::equilibrium(iPop,rho,u,uSqr) +
                        firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
                }
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
void StraightFdBoundary2D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    process(BlockLattice2D<T,Lattice>& blockLattice)
{
    processSubDomain(blockLattice, x0, x1, y0, y1);
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
template<int deriveDirection>
void StraightFdBoundary2D<T,Lattice,Dynamics,HydroBM,direction,orientation>::
    interpolateGradients(BlockLattice2D<T,Lattice> const& blockLattice,
                         T rhoUDeriv[Lattice<T>::d], int iX, int iY) const
{
    DirectedGradients2D<T, Lattice, direction, orientation,
                      direction==deriveDirection>::
        interpolate(rhoUDeriv, blockLattice, iX, iY, true);
}

////////  StraightFdBoundaryGenerator2D ////////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,
                              HydroBM,direction,orientation>::
    StraightFdBoundaryGenerator2D(int x0_, int x1_, int y0_, int y1_,
                                  T omega_)
    : PostProcessorGenerator2D<T,Lattice>(x0_, x1_, y0_, y1_),
      omega(omega_)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
PostProcessor2D<T,Lattice>*
    StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,
                                  HydroBM,direction,orientation>::
        generate(BlockLattice2D<T,Lattice>& blockLattice) const
{
    return new StraightFdBoundary2D<T,Lattice,Dynamics,
                                    HydroBM,direction,orientation>
                   ( this->x0, this->x1, this->y0, this->y1, omega,
                     blockLattice );
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         template <
             typename T_, template<typename U_> class Lattice_,
             int direction_, int orientation_ >
         class HydroBM,
         int direction, int orientation>
PostProcessorGenerator2D<T,Lattice>*
    StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,
                                  HydroBM,direction,orientation>::
        clone() const
{
    return new StraightFdBoundaryGenerator2D<T,Lattice,Dynamics,
                                             HydroBM,direction,orientation>
                   (this->x0, this->x1, this->y0, this->y1, omega);
}
 

////////// DirectedGradients2D /////////////////////////////////////////////

// Implementation for orthogonal==true; i.e. the derivative is along
// the boundary normal.
template<typename T, template<typename U> class Lattice,
         int direction, int orientation>
struct DirectedGradients2D<T, Lattice, direction, orientation, true> {
    static void interpolate(T velDeriv[Lattice<T>::d],
                            BlockLattice2D<T,Lattice> const& blockLattice,
                            int iX, int iY, bool multiplyRho)
    {
        using namespace fd;

        T rho0, rho1, rho2;
        T u0[Lattice<T>::d], u1[Lattice<T>::d], u2[Lattice<T>::d];
        
        // note that the derivative runs along direction.
        blockLattice.get(iX,iY).computeRhoU(rho0, u0);
        blockLattice.get (
            iX+(direction==0 ? (-orientation):0),
            iY+(direction==1 ? (-orientation):0) ).computeRhoU(rho1, u1);
        blockLattice.get (
            iX+(direction==0 ? (-2*orientation):0),
            iY+(direction==1 ? (-2*orientation):0) ).computeRhoU(rho2, u2);

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
         int direction, int orientation>
struct DirectedGradients2D<T, Lattice, direction, orientation, false> {
    static void  interpolate(T velDeriv[Lattice<T>::d],
                             BlockLattice2D<T,Lattice> const& blockLattice,
                             int iX, int iY, bool multiplyRho)
    {
        using namespace fd;

        T rho_p1, rho_m1;
        T u_p1[Lattice<T>::d], u_m1[Lattice<T>::d];
        
        int deriveDirection = 1-direction;
        blockLattice.get (
            iX+(deriveDirection==0 ? 1:0),
            iY+(deriveDirection==1 ? 1:0) ).computeRhoU(rho_p1, u_p1);
        blockLattice.get (
            iX+(deriveDirection==0 ? (-1):0),
            iY+(deriveDirection==1 ? (-1):0) ).computeRhoU(rho_m1, u_m1);

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


/////////// ConvexVelocityCorner2D /////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal>
ConvexVelocityCorner2D<T, Lattice, Dynamics, xNormal, yNormal>::
    ConvexVelocityCorner2D(int x_, int y_, T omega_,
                           BlockLattice2D<T,Lattice>& blockLattice)
    : x(x_), y(y_), omega(omega_),
      boundaryMomenta(0),
      boundaryDynamics(0)
{
    initialize(blockLattice);
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal>
ConvexVelocityCorner2D<T, Lattice, Dynamics, xNormal, yNormal>::
    ~ConvexVelocityCorner2D()
{
    delete boundaryMomenta;
    delete boundaryDynamics;
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal>
void ConvexVelocityCorner2D<T, Lattice, Dynamics, xNormal, yNormal>::
    initialize(BlockLattice2D<T,Lattice>& blockLattice_)
{
    boundaryMomenta = new GenericBoundaryMomenta<T,Lattice>();
    boundaryDynamics = new Dynamics (
          omega,
          *boundaryMomenta,
          blockLattice_.getStatistics() );
    blockLattice_.defineDynamics (
        x,x, y,y, boundaryDynamics
    );
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal>
void ConvexVelocityCorner2D<T, Lattice, Dynamics, xNormal, yNormal>::
    process(BlockLattice2D<T,Lattice>& blockLattice)
{
    using namespace olb::util::tensorIndices2D;
    T dx_rhoU[Lattice<T>::d], dy_rhoU[Lattice<T>::d];
    DirectedGradients2D<T, Lattice, 0, xNormal, true>::
        interpolate(dx_rhoU, blockLattice, x,y, true);
    DirectedGradients2D<T, Lattice, 1, yNormal, true>::
        interpolate(dy_rhoU, blockLattice, x,y, true);
    T dx_rhoUx = dx_rhoU[0];
    T dy_rhoUx = dy_rhoU[0];
    T dx_rhoUy = dx_rhoU[1];
    T dy_rhoUy = dy_rhoU[1];

    T omega = boundaryDynamics->getOmega();
    T sToPi = - (T)1 / Lattice<T>::invCs2 / omega;
    T pi[util::TensorVal<Lattice<T> >::n];
    pi[xx] = (T)2 * dx_rhoUx * sToPi;
    pi[yy] = (T)2 * dy_rhoUy * sToPi;
    pi[xy] = (dx_rhoUy + dy_rhoUx) * sToPi;
    boundaryMomenta->defineStress(pi);

    T rho10 = blockLattice.get(x-1*xNormal, y-0*yNormal).computeRho();
    T rho01 = blockLattice.get(x-0*xNormal, y-1*yNormal).computeRho();

    T rho20 = blockLattice.get(x-2*xNormal, y-0*yNormal).computeRho();
    T rho02 = blockLattice.get(x-0*xNormal, y-2*yNormal).computeRho();

    T rho = (T)2/(T)3*(rho01+rho10) - (T)1/(T)6*(rho02+rho20);
    boundaryMomenta->defineRho(rho);
    
    // Computation of the particle distribution functions
    // according to the regularized formula; implementation
    // by Orestis Malaspinas
    T u[Lattice<T>:: d];
    blockLattice.get(x,y).computeU(u);

    T uSqr = util::normSqr<T,2>(u);

    for (int iPop = 0; iPop < Lattice<T>::q; ++iPop)
        blockLattice.get(x,y)[iPop] =
            lbHelpers<T,Lattice>::equilibrium(iPop,rho,u,uSqr) +
            firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);

}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal>
void ConvexVelocityCorner2D<T, Lattice, Dynamics, xNormal, yNormal>::
    processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                     int x0_, int x1_, int y0_, int y1_ )
{
    if (util::contained(x, y, x0_, x1_, y0_, y1_)) {
        process(blockLattice);
    }
}


////////  ConvexVelocityCornerGenerator2D ////////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal>
ConvexVelocityCornerGenerator2D<T, Lattice, Dynamics, xNormal, yNormal>::
    ConvexVelocityCornerGenerator2D(int x_, int y_, T omega_)
        : PostProcessorGenerator2D<T,Lattice>(x_, x_, y_, y_),
          omega(omega_)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal>
PostProcessor2D<T,Lattice>*
ConvexVelocityCornerGenerator2D<T, Lattice, Dynamics, xNormal, yNormal>::
    generate( BlockLattice2D<T,Lattice>& blockLattice) const
{
    return new ConvexVelocityCorner2D<T, Lattice, Dynamics,
                                      xNormal, yNormal>
               ( this->x0, this->y0, omega, blockLattice );
}

template<typename T, template<typename U> class Lattice, typename Dynamics,
         int xNormal, int yNormal>
PostProcessorGenerator2D<T,Lattice>*
ConvexVelocityCornerGenerator2D<T, Lattice, Dynamics, xNormal, yNormal>::
    clone() const
{
    return new ConvexVelocityCornerGenerator2D<T, Lattice, Dynamics,
                                      xNormal, yNormal>
               ( this->x0, this->y0, omega );
}

}  // namespace olb

#endif

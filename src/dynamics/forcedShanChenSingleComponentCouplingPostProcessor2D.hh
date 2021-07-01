/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
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

#ifndef FORCED_SHAN_CHEN_SINGLE_COMPONENT_COUPLING_POST_PROCESSOR_2D_HH
#define FORCED_SHAN_CHEN_SINGLE_COMPONENT_COUPLING_POST_PROCESSOR_2D_HH

#include "forcedShanChenSingleComponentCouplingPostProcessor2D.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "core/finiteDifference2D.h"

namespace olb {

////////  ForcedShanChenSingleComponentCouplingPostProcessor2D ///////////////////////////////////

// Interaction Potentials
template<typename T>
T BoundedPsi(T rho) {
  return 1-exp(-rho);
};

template<typename T>
T ShanChen(T rho, T psiZero=4, T rhoZero=200) {
  return psiZero*exp(-rhoZero/rho);
};


template<typename T, template<typename U> class Lattice>
ForcedShanChenSingleComponentCouplingPostProcessor2D <T,Lattice>::
ForcedShanChenSingleComponentCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_, T G_,
                                      std::vector<T> rho0_,
                                      std::vector<SpatiallyExtendedObject2D*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), G(G_), rho0(rho0_), partners(partners_)
{ }

template<typename T, template<typename U> class Lattice>
void ForcedShanChenSingleComponentCouplingPostProcessor2D<T,Lattice>::
processSubDomain( BlockLattice2D<T,Lattice>& blockLattice,
                  int x0_, int x1_, int y0_, int y1_ )
{
  typedef Lattice<T> L;
  enum {
    uOffset     = L::ExternalField::velocityBeginsAt,
    forceOffset = L::ExternalField::forceBeginsAt
  };

  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( x0, x1, y0, y1,
                         x0_, x1_, y0_, y1_,
                         newX0, newX1, newY0, newY1 ) )
  {
    int nx = newX1-newX0+3; // include a one-cell boundary
    int ny = newY1-newY0+3; // include a one-cell boundary
    int offsetX = newX0-1;
    int offsetY = newY0-1;
    ScalarField2D<T> rhoField1(nx, ny);
    rhoField1.construct();
    // Compute density and velocity on every site of first lattice, and store result
    //   in external scalars; envelope cells are included, because they are needed
    //   to compute the interaction potential in what follows.
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        Cell<T,Lattice>& cell = blockLattice.get(iX,iY);
        rhoField1.get(iX-offsetX, iY-offsetY) = cell.computeRho()*rho0[0];
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,Lattice>& blockCell   = blockLattice.get(iX,iY);

        T* j = blockCell.getExternal(uOffset);
        lbHelpers<T,Lattice>::computeJ(blockCell,j);

        T blockOmega   = blockCell.getDynamics()->getOmega();
        // Computation of the common velocity, shared among the two populations
        T rhoTot = rhoField1.get(iX-offsetX, iY-offsetY)*blockOmega;

        T uTot[Lattice<T>::d];
        T *blockU = blockCell.getExternal(uOffset);      // contains precomputed value rho*u
        for (int iD = 0; iD < Lattice<T>::d; ++iD) {
          uTot[iD] = (blockU[iD]*rho0[0]*blockOmega) / rhoTot;
        }

        // Computation of the interaction potential
        T rhoBlockContribution[L::d]   = {T(), T()};
        for (int iPop = 0; iPop < L::q; ++iPop) {
          int nextX = iX + L::c[iPop][0];
          int nextY = iY + L::c[iPop][1];
          T blockRho   = ShanChen(rhoField1.get(nextX-offsetX, nextY-offsetY));//rho0[0]; // 1-exp(-rhoField1.get(nextX-offsetX, nextY-offsetY))
          for (int iD = 0; iD < L::d; ++iD) {
            rhoBlockContribution[iD]   += ShanChen(rhoField1.get(iX-offsetX, iY-offsetY)) * blockRho * L::c[iPop][iD]* L::t[iPop]; // 1-exp(-rhoField2.get(iX-offsetX, iY-offsetY))
          }
        }

        // Computation and storage of the final velocity, consisting
        //   of u and the momentum difference due to interaction
        //   potential plus external force
        T *blockForce   = blockCell.getExternal(forceOffset);
        for (int iD = 0; iD < L::d; ++iD) {
          blockU[iD] = uTot[iD] + 1./blockOmega *
                       (blockForce[iD] - G*rhoBlockContribution[iD]/rhoField1.get(iX-offsetX, iY-offsetY));
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void ForcedShanChenSingleComponentCouplingPostProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


/// LatticeCouplingGenerator for NS coupling

template<typename T, template<typename U> class Lattice>
ForcedShanChenSingleComponentCouplingGenerator2D<T,Lattice>::ForcedShanChenSingleComponentCouplingGenerator2D (
  int x0_, int x1_, int y0_, int y1_, T G_, std::vector<T> rho0_ )
  : LatticeCouplingGenerator2D<T,Lattice>(x0_, x1_, y0_, y1_), G(G_), rho0(rho0_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>* ForcedShanChenSingleComponentCouplingGenerator2D<T,Lattice>::generate (
  std::vector<SpatiallyExtendedObject2D*> partners) const
{
  return new ForcedShanChenSingleComponentCouplingPostProcessor2D<T,Lattice>(
           this->x0,this->x1,this->y0,this->y1,G, rho0, partners);
}

template<typename T, template<typename U> class Lattice>
LatticeCouplingGenerator2D<T,Lattice>* ForcedShanChenSingleComponentCouplingGenerator2D<T,Lattice>::clone() const {
  return new ForcedShanChenSingleComponentCouplingGenerator2D<T,Lattice>(*this);
}



}  // namespace olb

#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006,2007 Orestis Malaspinas and Jonas Latt
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

#ifndef KEEP_INCOMING_DYNAMICS_HH
#define KEEP_INCOMING_DYNAMICS_HH

#include "keepIncomingDynamics.h"
#include "core/latticeDescriptors.h"
#include "core/util.h"
#include "core/lbHelpers.h"
#include <cmath>


namespace olb {

using namespace descriptors;

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
KeepIncomingDynamics<T,Lattice,Dynamics,direction,orientation>::KeepIncomingDynamics (
        T omega_, Momenta<T,Lattice>& momenta_)
    : BasicDynamics<T,Lattice>(momenta_),
      boundaryDynamics(omega_, momenta_)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
KeepIncomingDynamics<T,Lattice,Dynamics,direction,orientation>*
    KeepIncomingDynamics<T,Lattice, Dynamics, direction, orientation>::clone() const
{
    return new KeepIncomingDynamics<T,Lattice,Dynamics,direction,orientation>(*this);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void KeepIncomingDynamics<T,Lattice,Dynamics,direction,orientation>::collide (
        Cell<T,Lattice>& cell,
        LatticeStatistics<T>& statistics )
{
    typedef lbHelpers<T,Lattice> lbH;
    typedef Lattice<T> L;

    int dir1 = direction;
    int dir2 = (direction+1)%2;
    int leftOrientation = orientation;
    int rightOrientation = orientation;
    if (direction==0) {
        rightOrientation *= -1;
    }
    else {
        leftOrientation *= -1;
    }
    std::vector<int> indexes(L::q);
    
    for (int iPop=0; iPop<L::q; ++iPop) {
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==leftOrientation) {
           indexes[1] = iPop;
       }
       if (L::c[iPop][dir1]==0 && L::c[iPop][dir2]==leftOrientation) {
           indexes[2] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==leftOrientation) {
           indexes[3] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==0) {
           indexes[4] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==rightOrientation) {
           indexes[5] = iPop;
       }
       if (L::c[iPop][dir1]==0 && L::c[iPop][dir2]==rightOrientation) {
           indexes[6] = iPop;
       }
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==rightOrientation) {
           indexes[7] = iPop;
       }
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==0) {
           indexes[8] = iPop;
       }
    }

    T rho, u[L::d], fNeq[L::q];
    cell.computeRhoU(rho,u);
    T uSqr = util::normSqr<T,L::d>(u);
    lbH::computeFneq(cell, fNeq, rho, u);
    
    T piValue = -(T)2/L::invCs2*rho/boundaryDynamics.getOmega() * StrainValue<T>::value;
    
    cell[indexes[4]] = lbH::equilibrium(indexes[4], rho, u, uSqr) + 
                       (T)2*(fNeq[indexes[1]]+fNeq[indexes[7]]) +
                       fNeq[indexes[8]] + fNeq[indexes[2]] + fNeq[indexes[6]] - piValue;
    cell[indexes[3]] = lbH::equilibrium(indexes[3], rho, u, uSqr) +
                       piValue/(T)2 - fNeq[indexes[1]] - fNeq[indexes[2]];
    cell[indexes[5]] = lbH::equilibrium(indexes[5], rho, u, uSqr) +
                       piValue/(T)2 - fNeq[indexes[6]] - fNeq[indexes[7]];

    boundaryDynamics.collide(cell, statistics);

    statistics.gatherStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void KeepIncomingDynamics<T,Lattice,Dynamics,direction,orientation>::staticCollide (
        Cell<T,Lattice>& cell,
        const T u[Lattice<T>::d],
        LatticeStatistics<T>& statistics )
{
    typedef lbHelpers<T,Lattice> lbH;
    typedef Lattice<T> L;

    int dir1 = direction;
    int dir2 = (direction+1)%2;
    int leftOrientation = orientation;
    int rightOrientation = orientation;
    if (direction==0) {
        rightOrientation *= -1;
    }
    else {
        leftOrientation *= -1;
    }
    std::vector<int> indexes(L::q);
    
    for (int iPop=0; iPop<L::q; ++iPop) {
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==leftOrientation) {
           indexes[1] = iPop;
       }
       if (L::c[iPop][dir1]==0 && L::c[iPop][dir2]==leftOrientation) {
           indexes[2] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==leftOrientation) {
           indexes[3] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==0) {
           indexes[4] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==rightOrientation) {
           indexes[5] = iPop;
       }
       if (L::c[iPop][dir1]==0 && L::c[iPop][dir2]==rightOrientation) {
           indexes[6] = iPop;
       }
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==rightOrientation) {
           indexes[7] = iPop;
       }
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==0) {
           indexes[8] = iPop;
       }
    }

    T rho = cell.computeRho();
    T uSqr = util::normSqr<T,L::d>(u);
    T fNeq[L::q];
    lbH::computeFneq(cell, fNeq, rho, u);
    
    T diff = fNeq[indexes[6]] - fNeq[indexes[2]];
    cell[indexes[4]] = lbH::equilibrium(4, rho, u, uSqr) + fNeq[indexes[8]];
    cell[indexes[3]] = lbH::equilibrium(3, rho, u, uSqr) + fNeq[indexes[7]] + 0.5*diff;
    cell[indexes[5]] = lbH::equilibrium(5, rho, u, uSqr) + fNeq[indexes[1]] - 0.5*diff;

    boundaryDynamics.staticCollide(cell, u, statistics);
    statistics.gatherStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T KeepIncomingDynamics<T,Lattice,Dynamics,direction,orientation>::getOmega() const 
{
    return boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void KeepIncomingDynamics<T,Lattice,Dynamics,direction,orientation>::setOmega(T omega_)
{
    boundaryDynamics.setOmega(omega_);
}








template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
S1Dynamics<T,Lattice,Dynamics,direction,orientation>::S1Dynamics (
        T omega_, Momenta<T,Lattice>& momenta_)
    : BasicDynamics<T,Lattice>(momenta_),
      boundaryDynamics(omega_, momenta_)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
S1Dynamics<T,Lattice,Dynamics,direction,orientation>*
    S1Dynamics<T,Lattice, Dynamics, direction, orientation>::clone() const
{
    return new S1Dynamics<T,Lattice,Dynamics,direction,orientation>(*this);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void S1Dynamics<T,Lattice,Dynamics,direction,orientation>::collide (
        Cell<T,Lattice>& cell,
        LatticeStatistics<T>& statistics )
{
    typedef lbHelpers<T,Lattice> lbH;
    typedef Lattice<T> L;

    int dir1 = direction;
    int dir2 = (direction+1)%2;
    int leftOrientation = orientation;
    int rightOrientation = orientation;
    if (direction==0) {
        rightOrientation *= -1;
    }
    else {
        leftOrientation *= -1;
    }
    std::vector<int> indexes(L::q);
    
    for (int iPop=0; iPop<L::q; ++iPop) {
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==leftOrientation) {
           indexes[1] = iPop;
       }
       if (L::c[iPop][dir1]==0 && L::c[iPop][dir2]==leftOrientation) {
           indexes[2] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==leftOrientation) {
           indexes[3] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==0) {
           indexes[4] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==rightOrientation) {
           indexes[5] = iPop;
       }
       if (L::c[iPop][dir1]==0 && L::c[iPop][dir2]==rightOrientation) {
           indexes[6] = iPop;
       }
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==rightOrientation) {
           indexes[7] = iPop;
       }
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==0) {
           indexes[8] = iPop;
       }
    }

    T rho, u[L::d], fNeq[L::q];
    cell.computeRhoU(rho,u);
    T uSqr = util::normSqr<T,L::d>(u);
    lbH::computeFneq(cell, fNeq, rho, u);
    
    cell[indexes[4]] = lbH::equilibrium(indexes[4], rho, u, uSqr) + fNeq[indexes[8]] +
                       (T)2*(fNeq[indexes[1]]+fNeq[indexes[7]]);
    cell[indexes[3]] = lbH::equilibrium(indexes[3], rho, u, uSqr) + fNeq[indexes[7]] +
                       (T)1/(T)2 * (fNeq[indexes[6]] - fNeq[indexes[2]]);
    cell[indexes[5]] = lbH::equilibrium(indexes[5], rho, u, uSqr) + fNeq[indexes[1]] +
                       (T)1/(T)2 * (fNeq[indexes[2]] - fNeq[indexes[6]]);

    boundaryDynamics.collide(cell, statistics);

    statistics.gatherStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void S1Dynamics<T,Lattice,Dynamics,direction,orientation>::staticCollide (
        Cell<T,Lattice>& cell,
        const T u[Lattice<T>::d],
        LatticeStatistics<T>& statistics )
{
    typedef lbHelpers<T,Lattice> lbH;
    typedef Lattice<T> L;

    int dir1 = direction;
    int dir2 = (direction+1)%2;
    int leftOrientation = orientation;
    int rightOrientation = orientation;
    if (direction==0) {
        rightOrientation *= -1;
    }
    else {
        leftOrientation *= -1;
    }
    std::vector<int> indexes(L::q);
    
    for (int iPop=0; iPop<L::q; ++iPop) {
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==leftOrientation) {
           indexes[1] = iPop;
       }
       if (L::c[iPop][dir1]==0 && L::c[iPop][dir2]==leftOrientation) {
           indexes[2] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==leftOrientation) {
           indexes[3] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==0) {
           indexes[4] = iPop;
       }
       if (L::c[iPop][dir1]==-orientation && L::c[iPop][dir2]==rightOrientation) {
           indexes[5] = iPop;
       }
       if (L::c[iPop][dir1]==0 && L::c[iPop][dir2]==rightOrientation) {
           indexes[6] = iPop;
       }
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==rightOrientation) {
           indexes[7] = iPop;
       }
       if (L::c[iPop][dir1]==orientation && L::c[iPop][dir2]==0) {
           indexes[8] = iPop;
       }
    }

    T rho = cell.computeRho();
    T uSqr = util::normSqr<T,L::d>(u);
    T fNeq[L::q];
    lbH::computeFneq(cell, fNeq, rho, u);
    
    T diff = fNeq[indexes[6]] - fNeq[indexes[2]];
    cell[indexes[4]] = lbH::equilibrium(4, rho, u, uSqr) + fNeq[indexes[8]];
    cell[indexes[3]] = lbH::equilibrium(3, rho, u, uSqr) + fNeq[indexes[7]] + 0.5*diff;
    cell[indexes[5]] = lbH::equilibrium(5, rho, u, uSqr) + fNeq[indexes[1]] - 0.5*diff;

    boundaryDynamics.staticCollide(cell, u, statistics);
    statistics.gatherStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T S1Dynamics<T,Lattice,Dynamics,direction,orientation>::getOmega() const 
{
    return boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void S1Dynamics<T,Lattice,Dynamics,direction,orientation>::setOmega(T omega_)
{
    boundaryDynamics.setOmega(omega_);
}


}  // namespace olb

#endif

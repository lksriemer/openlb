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

#ifndef STRAIGHT_BOUZIDI_DYNAMICS_HH
#define STRAIGHT_BOUZIDI_DYNAMICS_HH

#include "straightBouzidiDynamics.h"
#include "core/latticeDescriptors.h"
#include "core/util.h"
#include "core/lbHelpers.h"
#include <cmath>


namespace olb {

using namespace descriptors;

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
YoungDynamics<T,Lattice,Dynamics,direction,orientation>::YoungDynamics (
        T omega_, Momenta<T,Lattice>& momenta_ )
    : BasicDynamics<T,Lattice>(momenta_),
      boundaryDynamics(omega_, momenta_)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
YoungDynamics<T,Lattice,Dynamics,direction,orientation>* YoungDynamics<T,Lattice, Dynamics, direction, orientation>::clone() const
{
    return new YoungDynamics<T,Lattice,Dynamics,direction,orientation>(*this);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void YoungDynamics<T,Lattice,Dynamics,direction,orientation>::collide (
        Cell<T,Lattice>& cell,
        LatticeStatistics<T>& statistics )
{
    typedef lbHelpers<T,Lattice> lbH;
    typedef Lattice<T> L;

    std::vector<int> missingIndexes = util::subIndexOutgoing<L,direction,orientation>();
    T rho, u[L::d];
    this->momenta.computeRhoU(cell, rho, u);
    T uSqr = util::normSqr<T,L::d>(u);

    boundaryDynamics.collide(cell, statistics);
    
//     for (int iPop=1; iPop <= Lattice<T>::q/2; ++iPop) {
//         std::swap(cell[iPop], cell[iPop+Lattice<T>::q/2]);
//     }
//     
//     for (int iPop = 1; iPop < L::d; ++iPop)
//     {
// //         std::cout << missingIndexes[iPop] << ", ";
//         T c_u = T();
//         for (int iD = 0; iD < L::d; ++iD)
//         {
//             c_u += L::c[iPop][iD] * u[iD];
//         }
//         cell[iPop] += (T)2*L::invCs2*L::t[iPop]*c_u;
//     }

    for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop)
    {
        T c_u = T();
        for (int iD = 0; iD < L::d; ++iD)
        {
            c_u += L::c[missingIndexes[iPop]][iD] * u[iD];
        }
        cell[missingIndexes[iPop]] = cell[util::opposite<L>(missingIndexes[iPop])] +
                                     (T)2*L::invCs2*L::t[missingIndexes[iPop]]*c_u;
    }
    
//     T testRho, testU[L::d];
//     lbH::computeRhoU(cell,testRho,testU);
//     std::cout << rho - testRho << " " << u[0]-testU[0] << " " << u[1]-testU[1] <<"\n";
    
    
//     std::cout << "\n";

    statistics.gatherStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void YoungDynamics<T,Lattice,Dynamics,direction,orientation>::staticCollide (
        Cell<T,Lattice>& cell,
        const T u[Lattice<T>::d],
        LatticeStatistics<T>& statistics )
{
    typedef lbHelpers<T,Lattice> lbH;
    typedef Lattice<T> L;
    
    std::vector<int> missingIndexes = util::subIndexOutgoing<L,direction,orientation>();
    T rho = this->momenta.computeRho(cell);

    boundaryDynamics.staticCollide(cell, u, statistics);

    for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop)
    {
        T c_u = T();
        for (int iD = 0; iD < L::d; ++iD)
        {
            c_u += L::c[missingIndexes[iPop]][iD] * u[iD];
        }
        cell[missingIndexes[iPop]] = cell[util::opposite<L>(missingIndexes[iPop])] +
                (T)2*L::invCs2*L::t[missingIndexes[iPop]];
    }

}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T YoungDynamics<T,Lattice,Dynamics,direction,orientation>::getOmega() const 
{
    return boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void YoungDynamics<T,Lattice,Dynamics,direction,orientation>::setOmega(T omega_)
{
    boundaryDynamics.setOmega(omega_);
}


}  // namespace olb

#endif

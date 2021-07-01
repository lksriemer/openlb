/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Orestis Malaspinas, Jonas Latt
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef VISCO_OLDROYD_B_DYNAMICS_HH
#define VISCO_OLDROYD_B_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "viscoOldroydBdynamics.h"

namespace olb {

//==============================================================================//
/////////////////////////// Class ViscoOldroydBdynamics ///////////////////////////////
//==============================================================================//
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
ViscoOldroydBdynamics<T,Lattice>::ViscoOldroydBdynamics (
        T omega_, Momenta<T,Lattice>& momenta_ )
    : BasicDynamics<T,Lattice>(momenta_),
      omega(omega_)
{
	for (int iPop = Lattice<T>::q; ++iPop)
	{
		int count = 0;
		for (int iD = 0; iD < Lattice<T>::d; ++iD)
		{
			if (Lattice<T>::c[iPop][iD] != 0)
			{
				count++;
			}
		}
		if (count == 2)
			diagIndex.push_back(iPop);
		else if (count == 1)
			normIndex.push_back(iPop);
	}

}

template<typename T, template<typename U> class Lattice>
ViscoOldroydBdynamics<T,Lattice>* ViscoOldroydBdynamics<T,Lattice>::clone() const {
    return new ViscoOldroydBdynamics<T,Lattice>(*this);
}
 
template<typename T, template<typename U> class Lattice>
void ViscoOldroydBdynamics<T,Lattice>::collide (
        Cell<T,Lattice>& cell,
        LatticeStatistics<T>& statistics )
{
	using namespace util::tensorIndices2D;
    typedef Lattice<T> L;
    typedef ViscoSolventLbHelpers<T,Lattice> vLbH;
    
    T rho, u[L::d];
    this->momenta.computeRhoU(cell, rho, u);
    rho = rho - (T)1;
    	
    T* sigma = cell.getExternal(forceBeginsAt);
    
    vLbH::computeConstants(rho,sigma,A,B,C,D,E)
    
    for (int 
    vLbH::computeEquilibrium(fEq,rho,u,uSqr,A,B,C,D,E);
    
    
    
    if (cell.takesStatistics()) {
        statistics.gatherStats(rho, uSqr);
    }
}

template<typename T, template<typename U> class Lattice>
void ViscoOldroydBdynamics<T,Lattice>::staticCollide (
        Cell<T,Lattice>& cell,
        const T u[Lattice<T>::d],
        LatticeStatistics<T>& statistics )
{
    assert(false);
}

template<typename T, template<typename U> class Lattice>
T ViscoOldroydBdynamics<T,Lattice>::getOmega() const {
    return omega;
}

template<typename T, template<typename U> class Lattice>
void ViscoOldroydBdynamics<T,Lattice>::setOmega(T omega_) {
    omega = omega_;
}



}

#endif


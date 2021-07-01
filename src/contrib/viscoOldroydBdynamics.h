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
 * can be instantiated -- header file.
 */
#ifndef VISCO_OLDROYD_B_DYNAMICS_H
#define VISCO_OLDROYD_B_DYNAMICS_H

#include "../core/dynamics.h"

namespace olb {

template<typename T, template<typename U> class Lattice> class Cell;


/// Implementation of the entropic collision step
template<typename T, template<typename U> class Lattice>
class ViscoOldroydBdynamics : public BasicDynamics<T,Lattice> 
{
public:
    /// Constructor
    ViscoOldroydBdynamics(T omega_, Momenta<T,Lattice>& momenta_);
    /// Clone the object on its dynamic type.
    virtual ViscoOldroydBdynamics<T,Lattice>* clone() const;
    /// Collision step
    virtual void collide(Cell<T,Lattice>& cell,
                         LatticeStatistics<T>& statistics_);
    /// Collide with fixed velocity
    virtual void staticCollide(Cell<T,Lattice>& cell,
                               const T u[Lattice<T>::d],
                               LatticeStatistics<T>& statistics_);
    /// Get local relaxation parameter of the dynamics
    virtual T getOmega() const;
    /// Set local relaxation parameter of the dynamics
    virtual void setOmega(T omega_);
private:
    T omega;  ///< relaxation parameter
	std::vector<T> diagIndex; // indexes along the diagonal of the lattice
	std::vector<T> normIndex; // indexes along the principal axis of the lattice
};

}

#endif

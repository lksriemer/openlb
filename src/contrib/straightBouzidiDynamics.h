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

#ifndef STRAIGHT_BOUZIDI_DYNAMICS_H
#define STRAIGHT_BOUZIDI_DYNAMICS_H

#include "core/dynamics.h"

namespace olb {

/**
* Implementation of Zou-He boundary condition following
* the paper from Zou and He. This implementation is lattice independent.
*/
template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
class YoungDynamics : public BasicDynamics<T,Lattice>
{
public:
    /// Constructor
    YoungDynamics(T omega_, Momenta<T,Lattice>& momenta_);
    /// Clone the object on its dynamic type.
    virtual YoungDynamics<T, Lattice, Dynamics, direction, orientation>* clone() const;
    /// Collision step
    virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics);
    /// Collide with fixed velocity
    virtual void staticCollide(Cell<T,Lattice>& cell,
                               const T u[Lattice<T>::d],
                               LatticeStatistics<T>& statistics);
    /// Get local relaxation parameter of the dynamics
    virtual T getOmega() const;
    /// Set local relaxation parameter of the dynamics
    virtual void setOmega(T omega_);
private:
    Dynamics boundaryDynamics;
};

}

#endif

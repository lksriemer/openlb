/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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
 * Implementation of boundary cell dynamics -- header file.
 */
#ifndef HALF_WAY_BOUNCE_BACK_H
#define HALF_WAY_BOUNCE_BACK_H

#include "dynamics.h"
#include "util.h"
#include "cell.h"

namespace olb {

/// Computation of velocity momenta on a velocity boundary
template<typename T, template<typename U> class Lattice,
         int direction, int orientation>
class VelocityBM : virtual public DirichletBoundaryMomenta<T,Lattice> {
public:
    /// Default Constructor: initialization to zero
    VelocityBM();
    /// Constructor with boundary initialization
    VelocityBM(const T u_[Lattice<T>::d]);

    virtual T computeRho(Cell<T,Lattice> const& cell) const;
    virtual void computeU (
        Cell<T,Lattice> const& cell,
        T u[Lattice<T>::d] ) const;
    virtual void computeJ (
        Cell<T,Lattice> const& cell,
        T j[Lattice<T>::d] ) const;
    void computeU(T u[Lattice<T>::d]) const;
    virtual void defineRho(Cell<T,Lattice>& cell, T rho) ;
    virtual void defineU(Cell<T,Lattice>& cell,
                         const T u[Lattice<T>::d]) ;
    void defineU(const T u[Lattice<T>::d]);
    virtual void defineAllMomenta (
        Cell<T,Lattice>& cell,
        T rho, const T u[Lattice<T>::d],
        const T pi[util::TensorVal<Lattice<T> >::n] );
private:
    T u[Lattice<T>::d];   ///< value of the velocity on the boundary
};

}


#endif

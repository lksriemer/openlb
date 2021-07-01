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
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef VISCO_LATTICE_DESCRIPTORS_H
#define VISCO_LATTICE_DESCRIPTORS_H

#include <vector>

namespace olb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {
	
	struct Velocity2dDescriptor {
        static const int numScalars = 2;
        static const int numSpecies = 1;
        static const int velocityBeginsAt = 0;
        static const int sizeOfVelocity   = 2;
    };
	
	struct Conformation2dDescriptor {
        static const int numScalars = 3;
        static const int numSpecies = 1;
        static const int conformationBeginsAt = 0;
        static const int sizeOfConformation   = 3;
    };

    /// Solvent D2Q9 lattice
    template <typename T>
    struct SolventD2Q9Descriptor
    {
        typedef Conformation2dDescriptor ExternalField;
        enum { d = 2, q = 9 };      ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
    };
	
	/// Polymeric D2Q9 lattice
    template <typename T>
    struct PolymerD2Q9Descriptor
    {
        typedef Velocity2dDescriptor ExternalField;
        enum { d = 2, q = 9 };      ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
    };
}  // namespace descriptors

}  // namespace olb

#endif

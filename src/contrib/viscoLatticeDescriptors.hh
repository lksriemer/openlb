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
 *  -- generic code
 */
#ifndef VISCO_LATTICE_DESCRIPTORS_HH
#define VISCO_LATTICE_DESCRIPTORS_HH

#include "viscoLatticeDescriptors.h"

namespace olb {
namespace descriptors {

    // Solvent D2Q9 ////////////////////////////////////////////////////////////

    template<typename T>
    const int SolventD2Q9Descriptor<T>::c
        [SolventD2Q9Descriptor<T>::q][SolventD2Q9Descriptor<T>::d] =
        {
            { 0, 0},
            {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
            { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
        };
		
	// Polymer D2Q9 ////////////////////////////////////////////////////////////
	template<typename T>
    const int PolymerD2Q9Descriptor<T>::c
        [PolymerD2Q9Descriptor<T>::q][PolymerD2Q9Descriptor<T>::d] =
        {
            { 0, 0},
            {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
            { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
        };



}  // namespace descriptors

}  // namespace olb

#endif

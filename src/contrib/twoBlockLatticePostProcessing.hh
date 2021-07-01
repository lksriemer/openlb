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

#ifndef TWO_BLOCK_LATTICE_POST_PROCESSING_HH
#define TWO_BLOCK_LATTICE_POST_PROCESSING_HH

#include <cmath>
#include "core/blockLattice2D.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "core/ompManager.h"

namespace olb {

////////////////////// Class TwoBlockLatticePostProcessorGenerator2D /////////////////

template<typename T, template<typename U> class Lattice>
TwoBlockLatticePostProcessorGenerator2D<T,Lattice>::TwoBlockLatticePostProcessorGenerator2D (
            int x0One_, int x1One_, int y0One_, int y1One_,
            int x0Two_, int x1Two_, int y0Two_, int y1Two_ )
    : x0One(x0One_), x1One(x1One_), y0One(y0One_), y1One(y1One_),
	  x0Two(x0Two_), x1Two(x1Two_), y0Two(y0Two_), y1Two(y1Two_)
{ }

template<typename T, template<typename U> class Lattice>
void TwoBlockLatticePostProcessorGenerator2D<T,Lattice>::shift(int deltaXOne, int deltaYOne,
											    int deltaXTwo, int deltaYTwo) {
    x0One += deltaXOne;
    x1One += deltaXOne;
    y0One += deltaYOne;
    y1One += deltaYOne;
	
	x0Two += deltaXTwo;
    x1Two += deltaXTwo;
    y0Two += deltaYTwo;
    y1Two += deltaYTwo;
}

template<typename T, template<typename U> class Lattice>
bool TwoBlockLatticePostProcessorGenerator2D<T,Lattice>::
    extract(int x0One_, int x1One_, int y0One_, int y1One_,
		   	int x0Two_, int x1Two_, int y0Two_, int y1Two_)
{
	bool interOne, interTwo;
    int newX0One, newX1One, newY0One, newY1One;
    if ( util::intersect (
                x0One, x1One, y0One, y1One,
                x0One_, x1One_, y0One_, y1One_,
                newX0One, newX1One, newY0One, newY1One ) )
    {
        x0One = newX0One;
        x1One = newX1One;
        y0One = newY0One;
        y1One = newY1One;
        interOne = true;
    }
    else {
        interOne = false;
    }
	
	int newX0Two, newX1Two, newY0Two, newY1Two;
    if ( util::intersect (
                x0Two, x1Two, y0Two, y1Two,
                x0Two_, x1Two_, y0Two_, y1Two_,
                newX0Two, newX1Two, newY0Two, newY1Two ) )
    {
        x0Two = newX0Two;
        x1Two = newX1Two;
        y0Two = newY0Two;
        y1Two = newY1Two;
        interTwo = true;
    }
    else {
        interTwo = false;
    }
	
	if (interOne && interTwo)
	{
		return true;
	}
	else
	{
		return false;
	}
}


// ////////////////////// Class PostProcessorGenerator3D /////////////////
// 
// template<typename T, template<typename U> class Lattice>
// PostProcessorGenerator3D<T,Lattice>::PostProcessorGenerator3D (
//         int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
//     : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
// { }
// 
// template<typename T, template<typename U> class Lattice>
// void PostProcessorGenerator3D<T,Lattice>::shift (
//         int deltaX, int deltaY, int deltaZ )
// {
//     x0 += deltaX;
//     x1 += deltaX;
//     y0 += deltaY;
//     y1 += deltaY;
//     z0 += deltaZ;
//     z1 += deltaZ;
// }
// 
// template<typename T, template<typename U> class Lattice>
// bool PostProcessorGenerator3D<T,Lattice>::
//     extract(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
// {
//     int newX0, newX1, newY0, newY1, newZ0, newZ1;
//     if ( util::intersect (
//                 x0, x1, y0, y1, z0, z1,
//                 x0_, x1_, y0_, y1_, z0_, z1_,
//                 newX0, newX1, newY0, newY1, newZ0, newZ1 ) )
//     {
//         x0 = newX0;
//         x1 = newX1;
//         y0 = newY0;
//         y1 = newY1;
//         z0 = newZ0;
//         z1 = newZ1;
//         return true;
//     }
//     else {
//         return false;
//     }
// }

}  // namespace olb

#endif

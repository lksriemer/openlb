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

#include "boundaries.h"
#include "boundaries.hh"
#include "fdBoundaries2D.h"
#include "fdBoundaries2D.hh"
#include "postProcessing.h"
#include "postProcessing.hh"
#include "latticeDescriptors.h"
#include "latticeDescriptors.hh"


namespace olb {

// Instantiation with BGK dynamics

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 0, 1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 0, 1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 0, -1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 0, -1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 1, 1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 1, 1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 1, -1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 1, -1
>;


template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 0, 1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 0, 1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 0, -1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 0, -1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 1, 1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 1, 1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 1, -1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 1, -1
>;


template class ConvexVelocityCorner2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    1, 1
>;

template class ConvexVelocityCornerGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    1, 1
>;

template class ConvexVelocityCorner2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    1, -1
>;

template class ConvexVelocityCornerGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    1, -1
>;

template class ConvexVelocityCorner2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    -1, 1
>;

template class ConvexVelocityCornerGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    -1, 1
>;

template class ConvexVelocityCorner2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    -1, -1
>;

template class ConvexVelocityCornerGenerator2D <
    double, descriptors::D2Q9Descriptor,
    BGKdynamics <double, descriptors::D2Q9Descriptor>,
    -1, -1
>;



// Instantiation with regularized dynamics

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 0, 1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 0, 1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 0, -1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 0, -1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 1, 1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 1, 1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 1, -1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    VelocityBM, 1, -1
>;


template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 0, 1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 0, 1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 0, -1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 0, -1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 1, 1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 1, 1
>;

template class StraightFdBoundary2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 1, -1
>;

template class StraightFdBoundaryGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    PressureBM, 1, -1
>;


template class ConvexVelocityCorner2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    1, 1
>;

template class ConvexVelocityCornerGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    1, 1
>;

template class ConvexVelocityCorner2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    1, -1
>;

template class ConvexVelocityCornerGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    1, -1
>;

template class ConvexVelocityCorner2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    -1, 1
>;

template class ConvexVelocityCornerGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    -1, 1
>;

template class ConvexVelocityCorner2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    -1, -1
>;

template class ConvexVelocityCornerGenerator2D <
    double, descriptors::D2Q9Descriptor,
    RLBdynamics <double, descriptors::D2Q9Descriptor>,
    -1, -1
>;


}  // namespace olb

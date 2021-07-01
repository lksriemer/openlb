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

/** \file
 * Implementation of boundary cell dynamics -- template instantiation.
 */
#include "boundaries.h"
#include "boundaries.hh"
#include "boundaries3D.h"
#include "boundaries3D.hh"
#include "latticeDescriptors.h"
#include "latticeDescriptors.hh"

namespace olb {

    template class
        GenericBoundaryMomenta<double, descriptors::D3Q19Descriptor>;

    template class
        EquilibriumBM<double, descriptors::D3Q19Descriptor>;

    template class
        VelocityBM<double, descriptors::D3Q19Descriptor, 0, 1>;
    template class
        VelocityBM<double, descriptors::D3Q19Descriptor, 0, -1>;
    template class
        VelocityBM<double, descriptors::D3Q19Descriptor, 1, 1>;
    template class
        VelocityBM<double, descriptors::D3Q19Descriptor, 1, -1>;
    template class
        VelocityBM<double, descriptors::D3Q19Descriptor, 2, 1>;
    template class
        VelocityBM<double, descriptors::D3Q19Descriptor, 2, -1>;

    template class
        RegularizedVelocityBM<double, descriptors::D3Q19Descriptor, 0, 1>;
    template class
        RegularizedVelocityBM<double, descriptors::D3Q19Descriptor, 0, -1>;
    template class
        RegularizedVelocityBM<double, descriptors::D3Q19Descriptor, 1, 1>;
    template class
        RegularizedVelocityBM<double, descriptors::D3Q19Descriptor, 1, -1>;
    template class
        RegularizedVelocityBM<double, descriptors::D3Q19Descriptor, 2, 1>;
    template class
        RegularizedVelocityBM<double, descriptors::D3Q19Descriptor, 2, -1>;

    template class
        PressureBM<double, descriptors::D3Q19Descriptor, 0, 1>;
    template class
        PressureBM<double, descriptors::D3Q19Descriptor, 0, -1>;
    template class
        PressureBM<double, descriptors::D3Q19Descriptor, 1, 1>;
    template class
        PressureBM<double, descriptors::D3Q19Descriptor, 1, -1>;
    template class
        PressureBM<double, descriptors::D3Q19Descriptor, 2, 1>;
    template class
        PressureBM<double, descriptors::D3Q19Descriptor, 2, -1>;

    template class
        RegularizedPressureBM<double, descriptors::D3Q19Descriptor, 0, 1>;
    template class
        RegularizedPressureBM<double, descriptors::D3Q19Descriptor, 0, -1>;
    template class
        RegularizedPressureBM<double, descriptors::D3Q19Descriptor, 1, 1>;
    template class
        RegularizedPressureBM<double, descriptors::D3Q19Descriptor, 1, -1>;
    template class
        RegularizedPressureBM<double, descriptors::D3Q19Descriptor, 2, 1>;
    template class
        RegularizedPressureBM<double, descriptors::D3Q19Descriptor, 2, -1>;

    template class
        ConcaveCornerVelBM3D<double, descriptors::D3Q19Descriptor, 1, 1, 1>;
    template class
        ConcaveCornerVelBM3D<double, descriptors::D3Q19Descriptor, 1, 1,-1>;
    template class
        ConcaveCornerVelBM3D<double, descriptors::D3Q19Descriptor, 1,-1, 1>;
    template class
        ConcaveCornerVelBM3D<double, descriptors::D3Q19Descriptor, 1,-1,-1>;
    template class
        ConcaveCornerVelBM3D<double, descriptors::D3Q19Descriptor,-1, 1, 1>;
    template class
        ConcaveCornerVelBM3D<double, descriptors::D3Q19Descriptor,-1, 1,-1>;
    template class
        ConcaveCornerVelBM3D<double, descriptors::D3Q19Descriptor,-1,-1, 1>;
    template class
        ConcaveCornerVelBM3D<double, descriptors::D3Q19Descriptor,-1,-1,-1>;

    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 0, 1, 1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 0, 1,-1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 0,-1, 1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 0,-1,-1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 1, 1, 1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 1, 1,-1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 1,-1, 1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 1,-1,-1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 2, 1, 1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 2, 1,-1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 2,-1, 1>;
    template class
        ConcaveEdgeVelBM3D<double, descriptors::D3Q19Descriptor, 2,-1,-1>;

    template class
        StressInterpolationBM<double, descriptors::D3Q19Descriptor>;

    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, VelocityBM, 0, 1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, VelocityBM, 0, -1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, VelocityBM, 1, 1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, VelocityBM, 1, -1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, VelocityBM, 2, 1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, VelocityBM, 2, -1>;

    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, PressureBM, 0, 1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, PressureBM, 0, -1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, PressureBM, 1, 1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, PressureBM, 1, -1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, PressureBM, 2, 1>;
    template class
        StressInterpolationDirichletBM
            <double, descriptors::D3Q19Descriptor, PressureBM, 2, -1>;

}

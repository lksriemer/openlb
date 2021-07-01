/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010 Mathias J. Krause, Thomas Henn
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/> 
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
 * Input in STL format -- implementation.
 */

#include "stlReader.h"
#include <algorithm>
#include <cctype>
#include <iostream>

namespace olb {

STLreader::STLreader(const std::string& fName) {
    _innerMaterialNo = 1;
    _outerMaterialNo = 0;
    cvmlcpp::readSTL(_geometry, fName);
}

STLreader::~STLreader() {
}

void STLreader::read(BlockGeometry3D &matrix, unsigned direction,
        unsigned voxelNumber, unsigned pad, double fraction, unsigned samples) const {
    assert(direction >= 0);
    assert(direction <= 2);
    // Find dimensions
    double geometrySize = 0.;
    geometrySize = std::max(geometrySize, double(_geometry.max(direction))
            - double(_geometry.min(direction)));
    assert( geometrySize > 0.0 );
    assert( voxelNumber != 0.0 );
    double voxelSize = geometrySize / (double) voxelNumber;
    read(matrix, voxelSize, pad, fraction, samples);

}

void STLreader::read(BlockGeometry3D &matrix, double voxelSize, unsigned pad,
        double fraction, unsigned samples) const {
    if (fraction) {
        cvmlcpp::Matrix<double, 3u> voxels;
        cvmlcpp::fractionVoxelize(_geometry, voxels, voxelSize, samples, 1);
        matrix.reInit(_geometry.min(0) - (pad + 0.5) * voxelSize,
                      _geometry.min(1) - (pad + 0.5) * voxelSize,
                      _geometry.min(2) - (pad + 0.5) * voxelSize,
                      voxelSize,
                      voxels.extents()[X], voxels.extents()[Y], voxels.extents()[Z], pad);

        for (unsigned z = 0; z < voxels.extents()[Z]; ++z) {
            for (unsigned y = 0; y < voxels.extents()[Y]; ++y) {
                for (unsigned x = 0; x < voxels.extents()[X]; ++x) {
                    if (voxels[x][y][z] > fraction) {
                        matrix.setMaterial(x, y, z, _innerMaterialNo);
                    } else {
                        matrix.setMaterial(x, y, z, _outerMaterialNo);
                    }
                }
            }
        }
    } else {
        cvmlcpp::Matrix<unsigned short, 3u> voxels;
        cvmlcpp::voxelize(_geometry, voxels, voxelSize, 1);
        matrix.reInit(_geometry.min(0) - (pad + 0.5) * voxelSize,
                _geometry.min(1) - (pad + 0.5) * voxelSize, _geometry.min(2)
                        - (pad + 0.5) * voxelSize, voxelSize,
                voxels.extents()[X], voxels.extents()[Y], voxels.extents()[Z],
                pad);

        for (unsigned z = 0; z < voxels.extents()[Z]; ++z) {
            for (unsigned y = 0; y < voxels.extents()[Y]; ++y) {
                for (unsigned x = 0; x < voxels.extents()[X]; ++x) {
                    if (voxels[x][y][z] == 1)
                        matrix.setMaterial(x, y, z, _innerMaterialNo);
                    else
                        matrix.setMaterial(x, y, z, _outerMaterialNo);
                }
            }
        }
    }

}

void STLreader::setInnerMaterialNo(unsigned no) {
    _innerMaterialNo = no;
}

void STLreader::setOuterMaterialNo(unsigned no) {
    _outerMaterialNo = no;
}

} // namespace olb

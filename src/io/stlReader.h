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
 * Input in STL format -- header file.
 */

#ifndef STL_READER_H
#define STL_READER_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "../external/cvmlcpp/base/Matrix"
#include "../external/cvmlcpp/volume/Geometry"
#include "../external/cvmlcpp/volume/VolumeIO"
#include "../external/cvmlcpp/volume/Voxelizer"

#include "../core/blockGeometry3D.h"

namespace olb {

class STLreader {
public:
	/*
	 * Constructs a new STLreader from a file
	 * \param fName The STL file name
	 */
	STLreader( const std::string& fName);

	/*
	 * destructor
	 */
	~STLreader();

	/*
	 * Read a voxel mesh from the stl file
	 * \param matrix        contains the voxel mesh
	 * \param direction 
	 * \param voxelNumber   number of voxels in direction
	 * \param fraction      If fraction% of a voxel belong to the geometry the voxel belongs to the geometry
	 * \param samples       The number of samples along each dimension of the voxels.
	 */
	void read(BlockGeometry3D &matrix, unsigned direction, unsigned voxelNumber, unsigned pad=0, double fraction=0, unsigned samples=3u) const;

	/*
	 * Read a voxel mesh from the stl file
	 * \param matrix        to return the value
	 * \param voxelSize     The size of voxels.
	 * \param fraction      If fraction% of a voxel belong to the geometry the voxel belongs to the geometry
	 * \param samples       The number of samples along each dimension of the voxels.
	 */
	void read(BlockGeometry3D &matrix, double voxelSize, unsigned pad=0, double fraction=0, unsigned samples=16u) const;
    void setInnerMaterialNo(unsigned no);
    void setOuterMaterialNo(unsigned no);

private:
    cvmlcpp::Geometry<float> _geometry;
    int _innerMaterialNo, _outerMaterialNo;


};

}  // namespace olb

#endif  // STL_READER_H

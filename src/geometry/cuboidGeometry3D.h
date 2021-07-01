/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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
 * The description of a vector of 3D cuboid -- header file.
 */


#ifndef CUBOID_GEOMETRY_3D_H
#define CUBOID_GEOMETRY_3D_H

#include <vector>
#include "geometry/cuboid3D.h"
#include "io/ostreamManager.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

/// A cuboid geometry represents a voxel mesh
/** A cuboid geometry is given by a number of cuboids. To represent
 * a connected domain it is required that the distance between two
 * neighbooring cuboids is less than the smallest delta of them.
 *
 * By the class, access is provied to the cuboids. Core methods of 
 * the class are transforming lattice to physical positions in the 
 * corresponding unit systems.  
 *
 * WARNING:
 * At the moment there are only cuboids with a constant delta possible
 * and the distance between two neighbooring cuboids must be delta
 * since an interpolation operator in time and space is missing in
 * cuboidNeigbourhood and superLattice.
 *
 * This class is not intended to be derived from.
 */


template <typename T, typename S> class IndicatorF3D;

template<typename T>
class CuboidGeometry3D {

private:
  /// Vector of the cuboids
  std::vector<Cuboid3D<T> > _cuboids;
  /// Cuboid which contains all other cuboids
  Cuboid3D<T> _motherCuboid;
  /// Periodicity flag
  std::vector<bool> _periodicityOn;
  /// class specific ostream
  mutable OstreamManager clout;

public:
  /// Constructs a cuboid geometry with a cubic shape of size nX times nY times nZ with origin originR=(originX, originY, originZ) that consits of nC cuboids
  CuboidGeometry3D(T originX, T originY, T originZ, T deltaR, int nX, int nY, int nZ, int nC=1);
  /// Constructs a cuboid structure with a uniform spacing of voxelsize which consits of  nC cuboids, the cuboids not needed are removed and too big ones are shrinked
  CuboidGeometry3D(IndicatorF3D<bool,T>& indicatorF, T voxelSize, int nC=1);

  /// Read and write access to a single cuboid
  Cuboid3D<T>& get(int iC);
  /// Read access to a single cuboid
  Cuboid3D<T> const& get(int iC) const;
  /// Returns the smallest cuboid that includes all cuboids of the structure
  Cuboid3D<T> getMotherCuboid() const;
  /// Set flag to enable/disable periodicity 
  void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ);

  /// Gives for a given point (globX/globY/globZ) the related cuboidID
  /// and _p if the point is not in any of the cuboid _childrenQ
  int get_iC(T globX, T globY, T globZ, int offset = 0) const; //TODO old ones
  /// This function checks if the points (globX/globY/globZ) and
  /// (globX + orientationX/delta/globY + orientationY/delta/
  /// globZ + orientationZ/delta) is in a cuboid.
  /// It gives the related cuboidID and _p if the points are
  /// not in any of the cuboids.
  /// abs(orientationX) = abs(orientationY) = abs(orientationY) = 1
  /// must be satisfied
  int get_iC(T globX, T globY, T globZ, int orientationX, int orientationY, int orientationZ) const; //TODO old ones
  /// Returns true and the cuboid number of the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry 
  bool getC(std::vector<T> physR, int& iC) const; //TODO new one
  /// Returns true and the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry    
  bool getLatticeR(std::vector<T> physR, std::vector<int>& latticeR) const;
  /// Returns true and the floor lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getFloorLatticeR(std::vector<T> physR, std::vector<int>& latticeR) const;
  /// Returns the physical position to the given lattice position respecting periodicity for the overlap nodes which are not in the mother cuboid for the case the flag periodicityOn[iDim]=true if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  std::vector<T> getPhysR(int iCglob, int iX, int iY, int iZ) const;
  /// Returns the physical position to the given lattice position respecting periodicity for the overlap nodes which are not in the mother cuboid for the case the flag periodicityOn[iDim]=true 
  std::vector<T> getPhysR(std::vector<int> latticeR) const;

  /// Stores the iC of the neighbouring cuboids in a vector;
  void getNeighbourhood(int cuboid, std::vector<int>& neighbours, int offset = 0);
  /// Returns the number of cuboids in the structure
  int getNc() const;
  /// Returns the minimum of the ratio nX/NY in the structure
  T getMinRatio() const;
  /// Returns the maximum of the ratio nX/NY in the structure
  T getMaxRatio() const;
  /// Returns the minimum volume in the structure
  T getMinPhysVolume() const;
  /// Returns the maximum volume in the structure
  T getMaxPhysVolume() const;
  /// Returns the minimum number of nodes in the structure
  int getMinLatticeVolume() const;
  /// Returns the maximum number of nodes in the structure
  int getMaxLatticeVolume() const;
  /// Returns the minimum delta in the structure
  T getMinDeltaR() const;
  /// Returns the maximum delta in the structure
  T getMaxDeltaR() const;

  /// Sets the number of full cells of each cuboid
  void setWeights(IndicatorF3D<bool,T>& indicatorF);
  /// Adds a cuboid
  void add(Cuboid3D<T> cuboid);
  /// Splits cuboid iC, removes it and adds p cuboids
  void split(int iC, int p);
  /// Removes the cuboid iC
  void remove(int iC);
  /// Removes all cuboids where indicatorF = 0
  void remove(IndicatorF3D<bool,T>& indicatorF);
  /// Shrink all cuboids so that no empty planes are left
  void shrink(IndicatorF3D<bool,T>& indicatorF);

  /// Prints cuboid geometry details
  void print() const;
  /// Prints cuboid geometry details plus details of all cuboids
  void printExtended();
};

}  // namespace olb

#endif

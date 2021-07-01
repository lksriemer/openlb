/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2015 Thomas Henn, Mathias J. Krause
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
#include "communication/loadBalancer.h"
#include "geometry/cuboidGeometry3D.h"
#include "functors/indicatorF.h"
#include "utilities/vectorHelpers.h"
#include "octree.h"

using namespace olb::util;

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T>
class Octree;

template<typename T>
struct STLpoint {
  /// Constructor constructs
  STLpoint() : r(3, T()) {};
  /// Operator= equals
	STLpoint<T>& operator=(STLpoint<T> const& rhs) {
		r = rhs.r;
		return *this;
	};
  /// CopyConstructor copies
	STLpoint(STLpoint<T> const& rhs):r(rhs.r){};

	/// Point coordinates
  std::vector<T> r;
};

template<typename T>
struct STLtriangle {
  /** Test intersection between ray and triangle
   * \param pt Raypoint
   * \param dir Direction
   * \param q Point of intersection (if intersection occurs)
   * \param alpha Explained in http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf page 7-12
   *            q = pt + alpha * dir
   * \param rad It's complicated. Imagine you have a sphere with radius rad moving along the ray. Then q becomes the first point of the sphere to touch the triangle.
   */
  bool testRayIntersect(const std::vector<T>& pt,const std::vector<T>& dir, std::vector<T>& q, T& alpha, const T& rad = T()) const;

  /// A triangle contains 3 Points
  std::vector<STLpoint<T> > point;

  /// normal of triangle
  std::vector<T> normal;

  /// variables explained in http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf page 7-12
  std::vector<T> uBeta, uGamma;
  T d, kBeta, kGamma;

public:
  /// Constructor constructs
  STLtriangle():point(3, STLpoint<T>()), normal(3,T()), uBeta(3,T()), uGamma(3,T()), d(T()), kBeta(T()), kGamma(T()) {};
  /// CopyConstructor copies
  STLtriangle(STLtriangle<T> const& tri):point(tri.point), normal(tri.normal), uBeta(tri.uBeta), uGamma(tri.uGamma), d(tri.d), kBeta(tri.kBeta), kGamma(tri.kGamma) {};
  /// Operator= equals
  STLtriangle<T>& operator=(STLtriangle<T> const& tri) {
		point = tri.point; 
		normal = tri.normal; 
		uBeta = tri.uBeta; 
		uGamma = tri.uGamma; 
		d = tri.d; 
		kBeta = tri.kBeta; 
		kGamma = tri.kGamma;
		return *this;
	};

  /// Initializes triangle and precomputes member variables.
  void init();
  /// Returns Pt0-Pt1
  std::vector<T> getE0();
  /// Returns Pt0-Pt2
  std::vector<T> getE1();
};

template<typename T>
class STLmesh {
  /// Computes distance squared betwenn p1 and p2
  double distPoints(STLpoint<T>& p1, STLpoint<T>& p2);
  /// Filename
  const std::string& _fName;
  /// Vector of Triangles
  std::vector<STLtriangle<T> > _triangles;
  /// Min and Max points
  std::vector<T> _min, _max;
  /// I have no idea, probably max(distPoints)
  double _maxDist2;
  /// OstreamManager
  mutable OstreamManager clout;

public:
  /**
   * Constructs a new STLmesh from a file
   * \param Filename - Filename
   * \param stlSize - Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   */
  STLmesh(std::string, T stlSize = 1.);
 
  /// Returns reference to a triangle
  inline const STLtriangle<T>& getTri(size_t i) {return _triangles[i];}
  /// Returns number of triangles
  inline size_t triangleSize() const {
    return _triangles.size();
  }
  /// Returns _min
  inline std::vector<T>& getMin() {return _min;};
  /// Returns _max
  inline std::vector<T>& getMax() {return _max;};
  /// Returns maxDist squared
  inline float maxDist2() const {return _maxDist2;}
  /// Prints console output
  void print(bool full = false);
 /// Writes STL mesh in Si units
  void write(std::string fName);
};

template<typename T>
class STLreader : public IndicatorF3D<bool,T> {
private:
  /// Indicates (slow, more stable)
  void indicate1();
  /// Indicates (fast, less stable)
  void indicate2();
  /// Size of the smallest voxel
  T _voxelSize;
  /// Factor to get Si unit (m), i.e. "0.001" means mm 
  T _stlSize;
  /// The tree
  Octree<T>* _tree;
  /// The filename
  const std::string _fName;
  /// The mesh
  STLmesh<T> _mesh;
  /// Variable for output
  bool _verbose;
  /// The OstreamManager
  mutable OstreamManager clout;

public:
  /**
   * Constructs a new STLreader from a file
   * \param fName The STL file name
   * \param voxelSize Voxelsize in SI units
   * \param stlSize Conversion factor for STL (e.g. STL in mm stlSize=10^-3)
   * \param method Choose indication method
   *               0: fast, less stable
   *               1: slow, more stable (for untight STLs)
   * \param verbose Get additional information.
   */
  STLreader(const std::string fName, T voxelSize, T stlSize=1, unsigned short int method = 0, bool verbose = false);

  /// Returns whether node is inside or not.
  std::vector<bool> operator()(std::vector<T> input);

  /// Computes distance to closest triangle intersection
  bool distance(T& distance,const std::vector<T>& origin, const std::vector<T>& direction, int iC=-1);

  /// Prints console output
  void print();

  /// Writes STL mesh in Si units
  void writeSTL();

  /// Writes Octree 
  void writeOctree();

  /// Returns tree
  inline const Octree<T>* getTree() const {return _tree;};

  /// Returns mesh
  inline STLmesh<T>& getMesh() {return _mesh;};
};

}  // namespace olb

#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Thomas Henn
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
 * Octree - adapted from http://www.flipcode.com/archives/Octree_Implementation.shtml
 */

#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>

#include "core/singleton.h"

using namespace olb::util;
/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T>
class STLmesh;

template<typename T>
struct STLtriangle;

template <typename T>
class Octree
{
public:
  /*
   * Constructs Octree containing triangles of an STLmesh.
   * \param center Centerpoint
   * \param rad Radius
   * \param mesh STLmesh
   * \param maxDepth Maximal depth of tree
   * \param overlap Triangles within rad+overlap are added to this Octree
   */
  Octree(std::vector<T> center, T rad, STLmesh<T>& mesh, int maxDepth, T overlap = 0., Octree<T>* parent=NULL);
  /// Destructur destructs
  virtual ~Octree();

  /// Find the node containing the first param with remaining maxDepth
  Octree<T>* find(const std::vector<T>&,const int& maxDepth = 0);
  /// Write Octree
  void write(const std::vector<T>& pt, const std::string no);
  /// Write Octree
  void write(const int, const std::string);
  /// Write Octree
  void write(const std::string);
  /// Test intersection of ray with all triangles in Octree
  /// returns number od intersections
  int testIntersection(const std::vector<T>& pt,const std::vector<T>& dir, bool print = false);
  /// Test intersection of ray with all triangles in Octree
  /// q contains point of closest intersection to pt in direction direction
  bool closestIntersection(const std::vector<T>& pt, const std::vector<T>& direction, std::vector<T>& q, T& a);
  /// Test intersection of ray with all triangles in Octree
  /// q contains point of closest intersection to pt in direction direction
  /// tri contains triangle with closest intersection
  bool closestIntersection(const std::vector<T>& pt, const std::vector<T>& direction, std::vector<T>& q, T& a, STLtriangle<T>& tri, const T& rad = 0.);
  /// Test intersection of sphere moving along ray with radius rad
  /// q contains point of closest intersection to pt in direction direction
  /// tri contains triangle with closest intersection
  bool closestIntersectionSphere(const std::vector<T>& pt, const T& rad, const std::vector<T>& direction, std::vector<T>& q, T& a, STLtriangle<T>& tri);
  /// It's complicated. Computes intersections of a ray with triangles inside this Octree. Sets _inside depending on value of rayInside and changes rayInside depending on the number of intersections. Also takes into account if the intersections happen before or after the center.
  void checkRay(const std::vector<T>& pt,const std::vector<T>& dir, unsigned short& rayInside);
  /// Computes intersection of ray with Octree boundaries
  void intersectRayNode(const std::vector<T>& pt, const std::vector<T>& dir, std::vector<T>& s);
  /// Computes all centerpoints of Octree
  void getCenterpoints(std::vector<std::vector<T> >& pts);
  /// Collectes all leafs
  void getLeafs(std::vector<Octree<T>* >& pts);
  /// Sets Inside
  inline void setInside(bool ins) {_inside = ins;};
  /// Gets Inside
  inline bool getInside() {return _inside;};
  /// Gets Maxdepth
  inline int getMaxdepth() const {return _maxDepth;};
  /// Gets numbers of triangles contained by this Octree
  inline const std::vector<size_t>& getTriangles() const {return _triangles;};
  /// Gets centerpoint
  inline const std::vector<T>& getCenter() const {return _center;};
  /// Gets radius
  inline const T getRadius() const {return _radius;};
  /// Prints console output
  void print();

  /*
   * Atraxi: Is this world important?
   * The Doctor: Important? What's that mean, important? Six billion people live here, is that important? Here's a better question: is this world a threat to the Atraxi? Oh come on, you're monitoring the whole planet! Is this world a threat?
   * Atraxi: [after looking at a montage of world events] No.
   * The Doctor: Are the peoples of this world guilty of any crime by the laws of the Atraxi?
   * Atraxi: [after viewing another montage about earth] No.
   * The Doctor: Okay. One more, just one: is this world protected?
   */
protected:
  ///_vector _triangles contains number of triangles
  std::vector<size_t> _triangles;
  std::vector<T> _center;
  T _radius;
  STLmesh<T>& _mesh;
  int _maxDepth;
  bool _isLeaf;
  bool _inside;
  Octree<T> *_parent;
  Octree<T> *_child[8];
  mutable OstreamManager clout;
  void findTriangles(T overlap = 0.);
  bool AABBTri(const STLtriangle<T>& tri, T overlap = 0.);
};


}
#endif

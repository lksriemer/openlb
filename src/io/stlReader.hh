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

#ifndef STL_READER_HH
#define STL_READER_HH

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "core/singleton.h"
#include "communication/mpiManager.h"
#include "octree.hh"

using namespace olb::util;

/// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T>
void STLtriangle<T>::init() {
  std::vector<T> A = point[0].r;
  std::vector<T> B = point[1].r;
  std::vector<T> C = point[2].r;
  std::vector<T> b(3, 0.), c(3, 0.);
 
  T bb = 0., bc = 0., cc = 0.;
 	
  for (int i = 0; i < 3; i++) {
    b[i] = B[i] - A[i];
    c[i] = C[i] - A[i];
    bb += b[i]*b[i];
    bc += b[i]*c[i];
    cc += c[i]*c[i];
  }
  
  normal[0] = b[1]*c[2] - b[2]*c[1];
  normal[1] = b[2]*c[0] - b[0]*c[2];
  normal[2] = b[0]*c[1] - b[1]*c[0];
 	
  T norm = sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2) 
 	                 + std::pow(normal[2], 2));
  normal[0] /= norm;
  normal[1] /= norm;
  normal[2] /= norm;

  T D = 1.0 / (cc*bb - bc*bc);
  T bbD = bb*D;
  T bcD = bc*D;
  T ccD = cc*D;
	 
  kBeta = 0.;
  kGamma = 0.;
  d = 0.;
 
  for (int i = 0; i < 3; i++) {
    uBeta[i] = b[i]*ccD - c[i]*bcD;
    uGamma[i] = c[i]*bbD - b[i]*bcD;
    kBeta -= A[i]*uBeta[i];
    kGamma -= A[i]*uGamma[i];
    d += A[i]*normal[i];
  }
}

template <typename T>
std::vector<T> STLtriangle<T>::getE0() {
  std::vector<T> vec(3, T());
  vec[0] = point[0].r[0] - point[1].r[0];
  vec[1] = point[0].r[1] - point[1].r[1];
  vec[2] = point[0].r[2] - point[1].r[2];
  return vec;
}

template <typename T>
std::vector<T> STLtriangle<T>::getE1() {
  std::vector<T> vec(3, T());
  vec[0] = point[0].r[0] - point[2].r[0];
  vec[1] = point[0].r[1] - point[2].r[1];
  vec[2] = point[0].r[2] - point[2].r[2];
  return vec;
}

/* Schnitttest nach
 * http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf
 */
template <typename T>
bool STLtriangle<T>::testRayIntersect(const std::vector<T>& pt,
                                      const std::vector<T>& dir, 
                                      std::vector<T>& q, T& alpha,const T& rad) const {
  T rn = 0.;
  std::vector<T> testPt = pt + rad*normal;
  std::vector<T> help(3,0.);
  q.resize(3);
 	
  for (int i = 0; i < 3; i++) {
    rn += dir[i]*normal[i];
  }

  // Schnitttest Flugrichtung -> Ebene
  if (fabs(rn) < std::numeric_limits<T>::epsilon()) {
    return false;
  }
  alpha = d;
  for (int i = 0; i < 3; i++) {
    alpha -= testPt[i]*normal[i];
  }
  alpha /= rn;

  // Abstand Partikel Ebene
  if (alpha < 0) {
    return false;
  }
  for (int i = 0; i < 3; i++) {
    q[i] = testPt[i] + alpha*dir[i];
  }
  double beta = kBeta;
  for (int i = 0; i < 3; i++) {
    beta += uBeta[i]*q[i];
  }

  // Schnittpunkt q in der Ebene?
  if (beta < -std::numeric_limits<T>::epsilon()) {
    return false;
  }
  double gamma = kGamma;
  for (int i = 0; i < 3; i++) {
    gamma += uGamma[i]*q[i];
  }
  if (gamma < -std::numeric_limits<T>::epsilon()) {
    return false;
  }
  if (1 - beta - gamma < -std::numeric_limits<T>::epsilon()) {
   return false;
  }
  return true;
}

template <typename T>
STLmesh<T>::STLmesh(std::string fName, T stlSize): _fName(fName), _min(3,0), _max(3,0), _maxDist2(0), clout(std::cout, "STLmesh"){
  std::ifstream f(fName.c_str(), std::ios::in);
_triangles.reserve(10000);
  if (!f.good()) throw std::runtime_error("STL File not valid.");
  char buf[6];
  buf[5] = 0;
  f.read(buf, 5);
  const std::string asciiHeader = "solid";
  if (std::string(buf) == asciiHeader) {
    f.seekg(0, std::ios::beg);
    if (f.good()) {
      std::string s0, s1;
      int i = 0;
      while (!f.eof()) {
        f >> s0;
        if (s0 == "facet") {
          STLtriangle<T> tri;
          f >> s1 >> tri.normal[0] >> tri.normal[1] >> tri.normal[2];
          f >> s0 >> s1;
          f >> s0 >> tri.point[0].r[0] >> tri.point[0].r[1] >> tri.point[0].r[2];
          f >> s0 >> tri.point[1].r[0] >> tri.point[1].r[1] >> tri.point[1].r[2];
          f >> s0 >> tri.point[2].r[0] >> tri.point[2].r[1] >> tri.point[2].r[2];
          f >> s0;
          f >> s0;
          for (int k=0; k<3; k++) {
            tri.point[0].r[k] *= stlSize;
            tri.point[1].r[k] *= stlSize;
            tri.point[2].r[k] *= stlSize;
          }
          if (i==0) {
            _min.resize(3,0);
            _max.resize(3,0);

            _min[0] = tri.point[0].r[0];
            _min[1] = tri.point[0].r[1];
            _min[2] = tri.point[0].r[2];

            _max[0] = tri.point[0].r[0];
            _max[1] = tri.point[0].r[1];
            _max[2] = tri.point[0].r[2];

            _min[0] = std::min(_min[0],tri.point[1].r[0]);
            _min[1] = std::min(_min[1],tri.point[1].r[1]);
            _min[2] = std::min(_min[2],tri.point[1].r[2]);

            _max[0] = std::max(_max[0],tri.point[1].r[0]);
            _max[1] = std::max(_max[1],tri.point[1].r[1]);
            _max[2] = std::max(_max[2],tri.point[1].r[2]);

            _min[0] = std::min(_min[0],tri.point[2].r[0]);
            _min[1] = std::min(_min[1],tri.point[2].r[1]);
            _min[2] = std::min(_min[2],tri.point[2].r[2]);

            _max[0] = std::max(_max[0],tri.point[2].r[0]);
            _max[1] = std::max(_max[1],tri.point[2].r[1]);
            _max[2] = std::max(_max[2],tri.point[2].r[2]);

          } else {
            _min[0] = std::min(_min[0],tri.point[0].r[0]);
            _min[1] = std::min(_min[1],tri.point[0].r[1]);
            _min[2] = std::min(_min[2],tri.point[0].r[2]);

            _max[0] = std::max(_max[0],tri.point[0].r[0]);
            _max[1] = std::max(_max[1],tri.point[0].r[1]);
            _max[2] = std::max(_max[2],tri.point[0].r[2]);

            _min[0] = std::min(_min[0],tri.point[1].r[0]);
            _min[1] = std::min(_min[1],tri.point[1].r[1]);
            _min[2] = std::min(_min[2],tri.point[1].r[2]);

            _max[0] = std::max(_max[0],tri.point[1].r[0]);
            _max[1] = std::max(_max[1],tri.point[1].r[1]);
            _max[2] = std::max(_max[2],tri.point[1].r[2]);
 
            _min[0] = std::min(_min[0],tri.point[2].r[0]);
            _min[1] = std::min(_min[1],tri.point[2].r[1]);
            _min[2] = std::min(_min[2],tri.point[2].r[2]);

            _max[0] = std::max(_max[0],tri.point[2].r[0]);
            _max[1] = std::max(_max[1],tri.point[2].r[1]);
            _max[2] = std::max(_max[2],tri.point[2].r[2]);
          }
        
          i++;
          tri.init();
          _triangles.push_back(tri);
 
          _maxDist2 = std::max(distPoints(tri.point[0], tri.point[1]), _maxDist2);
          _maxDist2 = std::max(distPoints(tri.point[2], tri.point[1]), _maxDist2);
          _maxDist2 = std::max(distPoints(tri.point[0], tri.point[2]), _maxDist2);
        } else if (s0 == "endsolid") {
          break;
        }
      }
    }
  } else {
    f.close();
    f.open(fName.c_str(), std::ios::in | std::ios::binary);
    char comment[80];
    f.read(comment, 80);
		
    if (!f.good()) throw std::runtime_error("STL File not valid.");
		
    comment[79] = 0;
    int32_t nFacets;
    f.read(reinterpret_cast<char *>(&nFacets), sizeof(int32_t));
		
    if (!f.good()) throw std::runtime_error("STL File not valid.");

    float v[12];
    unsigned short uint16;
    for (int32_t i = 0; i < nFacets; ++i) {
      for (size_t j = 0; j < 12; ++j) {
        f.read(reinterpret_cast<char *>(&v[j]), sizeof(float));
      }
      f.read(reinterpret_cast<char *>(&uint16), sizeof(unsigned short)); 
      STLtriangle<T> tri;
      tri.normal[0] = v[0];
      tri.normal[1] = v[1];
      tri.normal[2] = v[2];
      tri.point[0].r[0] = v[3];
      tri.point[0].r[1] = v[4];
      tri.point[0].r[2] = v[5];
      tri.point[1].r[0] = v[6];
      tri.point[1].r[1] = v[7];
      tri.point[1].r[2] = v[8];
      tri.point[2].r[0] = v[9];
      tri.point[2].r[1] = v[10];
      tri.point[2].r[2] = v[11];
   			
      for (int k=0; k<3; k++) {
        tri.point[0].r[k] *= stlSize;
        tri.point[1].r[k] *= stlSize;
        tri.point[2].r[k] *= stlSize;
      }
      if (i==0) {
        _min[0] = tri.point[0].r[0];
        _min[1] = tri.point[0].r[1];
        _min[2] = tri.point[0].r[2];

        _max[0] = tri.point[0].r[0];
        _max[1] = tri.point[0].r[1];
        _max[2] = tri.point[0].r[2];

        _min[0] = std::min(_min[0],tri.point[1].r[0]);
        _min[1] = std::min(_min[1],tri.point[1].r[1]);
        _min[2] = std::min(_min[2],tri.point[1].r[2]);

        _max[0] = std::max(_max[0],tri.point[1].r[0]);
        _max[1] = std::max(_max[1],tri.point[1].r[1]);
        _max[2] = std::max(_max[2],tri.point[1].r[2]);

        _min[0] = std::min(_min[0],tri.point[2].r[0]);
        _min[1] = std::min(_min[1],tri.point[2].r[1]);
        _min[2] = std::min(_min[2],tri.point[2].r[2]);
				
        _max[0] = std::max(_max[0],tri.point[2].r[0]);
        _max[1] = std::max(_max[1],tri.point[2].r[1]);
        _max[2] = std::max(_max[2],tri.point[2].r[2]);

      } else {
        _min[0] = std::min(_min[0],tri.point[0].r[0]);
        _min[1] = std::min(_min[1],tri.point[0].r[1]);
        _min[2] = std::min(_min[2],tri.point[0].r[2]);
				
        _max[0] = std::max(_max[0],tri.point[0].r[0]);
        _max[1] = std::max(_max[1],tri.point[0].r[1]);
        _max[2] = std::max(_max[2],tri.point[0].r[2]);

        _min[0] = std::min(_min[0],tri.point[1].r[0]);
        _min[1] = std::min(_min[1],tri.point[1].r[1]);
        _min[2] = std::min(_min[2],tri.point[1].r[2]);

        _max[0] = std::max(_max[0],tri.point[1].r[0]);
        _max[1] = std::max(_max[1],tri.point[1].r[1]);
        _max[2] = std::max(_max[2],tri.point[1].r[2]);

        _min[0] = std::min(_min[0],tri.point[2].r[0]);
        _min[1] = std::min(_min[1],tri.point[2].r[1]);
        _min[2] = std::min(_min[2],tri.point[2].r[2]);
				
        _max[0] = std::max(_max[0],tri.point[2].r[0]);
        _max[1] = std::max(_max[1],tri.point[2].r[1]);
        _max[2] = std::max(_max[2],tri.point[2].r[2]);
      }
      tri.init();
      _triangles.push_back(tri);

      _maxDist2 = std::max(distPoints(tri.point[0], tri.point[1]), _maxDist2);
      _maxDist2 = std::max(distPoints(tri.point[2], tri.point[1]), _maxDist2);
      _maxDist2 = std::max(distPoints(tri.point[0], tri.point[2]), _maxDist2);
    }
  }
  f.close();
}

template <typename T>
double STLmesh<T>::distPoints(STLpoint<T>& p1, STLpoint<T>& p2) {
  return std::pow(double(p1.r[0] - p2.r[0]), 2)
       + std::pow(double(p1.r[1] - p2.r[1]), 2)
       + std::pow(double(p1.r[2] - p2.r[2]), 2);
}

template <typename T>
void STLmesh<T>::print(bool full) {
  if (full) {
    int i=0;
    clout << "Triangles: " << std::endl;
    typename std::vector<STLtriangle<T> >::iterator it = _triangles.begin();
  		
    for (; it!=_triangles.end() ; it++) {
      clout <<i++<<": "<< it->point[0].r[0] << " " << it->point[0].r[1] 
                << " " << it->point[0].r[2] << " | " << it->point[1].r[0] << " " 
                << it->point[1].r[1] << " " << it->point[1].r[2] << " | "
                << it->point[2].r[0] << " " << it->point[2].r[1] << " " 
                << it->point[2].r[2] << std::endl;
    }
  }
  clout << "nTriangles=" << _triangles.size() << "; maxDist2=" << _maxDist2 << std::endl;
  clout << "minPhysR(StlMesh)=(" << getMin()[0] << "," << getMin()[1] << "," << getMin()[2] << ")";
  clout << "; maxPhysR(StlMesh)=(" << getMax()[0] << "," << getMax()[1] << "," << getMax()[2] << ")" << std::endl;
}

template <typename T>
void STLmesh<T>::write(std::string fName) {
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if(rank==0) {
    std::string fullName = singleton::directories().getVtkOutDir() + fName + ".stl";
    std::ofstream f(fullName.c_str());
    f << "solid ascii " << fullName << "\n";
  		
    for (unsigned int i=0; i<_triangles.size(); i++) {
      f << "facet normal " << _triangles[i].normal[0] << " " 
        << _triangles[i].normal[1] << " " << _triangles[i].normal[2] << "\n";
      f << "    outer loop\n";
      f << "        vertex " << _triangles[i].point[0].r[0] << " " 
        << _triangles[i].point[0].r[1] << " " << _triangles[i].point[0].r[2] << "\n";
      f << "        vertex " << _triangles[i].point[1].r[0] << " " 
        << _triangles[i].point[1].r[1] << " " << _triangles[i].point[1].r[2] << "\n";
      f << "        vertex " << _triangles[i].point[2].r[0] << " " 
        << _triangles[i].point[2].r[1] << " " << _triangles[i].point[2].r[2] << "\n";
      f << "    endloop\n";
      f << "endfacet\n";
    }
    f.close();
  }
  /*if (_verbose)*/ clout << "Write ... OK" << std::endl;
}

/*
 * STLReader functions
 */
template <typename T>
STLreader<T>::STLreader(const std::string fName, T voxelSize, T stlSize, unsigned short int method, bool verbose): IndicatorF3D<bool,T>(1), _voxelSize(voxelSize), _stlSize(stlSize), _fName(fName), _mesh(fName, stlSize), _verbose(verbose), clout(std::cout, "STLreader") {
	if (_verbose) clout << "Voxelizing ..." << std::endl;

	std::vector<T> extension = _mesh.getMax() - _mesh.getMin();
	T max = std::max(extension[0], std::max(extension[1], extension[2]))  + _voxelSize;
	int j = 0;
	for (; _voxelSize * std::pow(2, j) < max; j++);
	std::vector<T> center(3, T());
	T radius = _voxelSize * std::pow(2, j-1);

	/// Find center of tree and move by _voxelSize/4.
	for (int i=0; i<3; i++) {
		center[i] = (_mesh.getMin()[i] + _mesh.getMax()[i]) / 2. - _voxelSize/4.;
	}

	/// Create tree
	_tree = new Octree<T>(center, radius, _mesh, j);

	/// Compute _myMin, _myMax such that they are the smallest (greatest) Voxel inside the STL.
	this->_myMin = center + _voxelSize/2.;
	this->_myMax = center - _voxelSize/2.;
	for (int i=0; i<3; i++) {
	  while (this->_myMin[i] > _mesh.getMin()[i]) {
	    this->_myMin[i] -= _voxelSize;
	  }
	  while (this->_myMax[i] < _mesh.getMax()[i]) {
	    this->_myMax[i] += _voxelSize;
	  }
	  this->_myMax[i] -= _voxelSize;
    this->_myMin[i] += _voxelSize;
	}

	/// Indicate nodes of the tree. (Inside/Outside)
	switch (method) {
	  case 1: indicate1(); break;
	  default: indicate2(); break;
	}

	if (_verbose) print();
	if (_verbose) clout << "Voxelizing ... OK" << std::endl;
}

/*
 *  Old indicate function (slower, more stable)
 *  Define three rays (X-, Y-, Z-direction) for each leaf and count intersections
 *  with STL for each ray. Odd number of intersection means inside (Majority vote).
 */
template <typename T>
void STLreader<T>::indicate1() {
	std::vector<Octree<T>* > leafs;
	_tree->getLeafs(leafs);
	typename std::vector<Octree<T>* >::iterator it = leafs.begin();
	std::vector<T> dir(3, 0), pt(3, T()), s(3, T());
	int intersections = 0;
	int inside = 0;
	Octree<T>* node = NULL;
	T step = 1./1000. * _voxelSize;
	for (; it!= leafs.end(); it++) {
		inside = 0;

		pt = (*it)->getCenter();
		intersections = 0;
		s = pt; // + step;

		/// X+ dir
		dir[0] = 1;
		dir[1] = 0;
		dir[2] = 0;
		while (s[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
			node = _tree->find(s, (*it)->getMaxdepth());
			intersections += node->testIntersection(pt, dir);
			node->intersectRayNode(pt, dir, s);
			s = s + step*dir;
		}
		inside +=(intersections%2);

		/// Y+ Test
		intersections = 0;
		s = pt; // + step;
		dir[0] = 0;
		dir[1] = 1;
		dir[2] = 0;
		while (s[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
			node = _tree->find(s, (*it)->getMaxdepth());
			intersections += node->testIntersection(pt, dir);
			node->intersectRayNode(pt, dir, s);
			s = s + step*dir;
		}
		inside +=(intersections%2);

		/// Z+ Test
		intersections = 0;
		s = pt; // + step;
		dir[0] = 0;
		dir[1] = 0;
		dir[2] = 1;
		while (s[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
			node = _tree->find(s, (*it)->getMaxdepth());
			intersections += node->testIntersection(pt, dir);
			node->intersectRayNode(pt, dir, s);
			s = s + step*dir;
		}
		inside +=(intersections%2);
		(*it)->setInside(inside>1);
	}
}

/*
 *  New indicate function (faster, less stable)
 *  Define ray in Z-direction for each Voxel in XY-layer. Indicate all nodes on the fly.
 */
template <typename T>
void STLreader<T>::indicate2() {
	T rad = _tree->getRadius();
	std::vector<T> rayPt = _tree->getCenter() - rad  + .5* _voxelSize;
	std::vector<T> pt = rayPt;
	std::vector<T> rayDir(3, T());
			rayDir[0] = 0.;
			rayDir[1] = 0.;
			rayDir[2] = 1.;
	std::vector<T> maxEdge = _tree->getCenter() +rad;

	T step = 1./1000. * _voxelSize;

	Octree<T>* node = NULL;
	unsigned short rayInside = 0;
	std::vector<T> nodeInters(3, T());
	while (pt[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
    node = _tree->find(pt);
    nodeInters = pt;
    nodeInters[2] = node->getCenter()[2] - node->getRadius();
    rayInside = 0;
    while (pt[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(pt);
      nodeInters = pt;
      nodeInters[2] = node->getCenter()[2] - node->getRadius();
      rayInside = 0;
			while (pt[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
			  node = _tree->find(pt);
				node->checkRay(nodeInters, rayDir, rayInside);
				node->intersectRayNode(pt, rayDir, nodeInters);
				pt = nodeInters + step*rayDir;
			}
			pt[2] = rayPt[2];
			pt[1] += _voxelSize;
		}
		pt[1] = rayPt[1];
		pt[0] += _voxelSize;
	}
}

template <typename T>
std::vector<bool> STLreader<T>::operator()(std::vector<T> pt) {
  std::vector<bool> ret(1, false);
  T r = _tree->getRadius();
  std::vector<T> c = _tree->getCenter();
  if (c[0]-r < pt[0] && pt[0] < c[0]+r &&
      c[1]-r < pt[1] && pt[1] < c[1]+r &&
      c[2]-r < pt[2] && pt[2] < c[2]+r) {
    ret[0] = _tree->find(pt)->getInside();
  }
  return ret;
}

template <typename T>
bool STLreader<T>::distance(T& distance, const std::vector<T>& origin, 
                            const std::vector<T>& direction, int iC) {
  Octree<T>* node = NULL;
  std::vector<T> dir = normalize(direction);
  std::vector<T> extends = _mesh.getMax() - _mesh.getMin(), pt = origin, q(3, T()), s(3, T());
  std::vector<T> center = _mesh.getMin() + 1/2. * extends;
  T step = _voxelSize/1000., a=0;
 	
  for (int i=0; i<3; i++) {
    extends[i]/=2.;
  }

  if (!(_mesh.getMin()[0] < origin[0] && origin[0] < _mesh.getMax()[0] &&
        _mesh.getMin()[1] < origin[1] && origin[1] < _mesh.getMax()[1] &&
        _mesh.getMin()[2] < origin[2] && origin[2] < _mesh.getMax()[2])) {
    T t=T(), d=T();
    bool foundQ = false;
  		
    if (dir[0] > 0) {
      d = _mesh.getMin()[0];
      t = (d - origin[0])/dir[0];
      pt = origin + (t+step)*dir;
      if(_mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1] &&
         _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
       }
    } else if (dir[0] < 0){
      d = _mesh.getMax()[0];
      t = (d - origin[0])/dir[0];
      pt = origin + (t+step)*dir;
      if(_mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1] &&
         _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }

    if (dir[1] > 0 && !foundQ) {
      d = _mesh.getMin()[1];
      t = (d - origin[1])/dir[1];
      pt = origin + (t+step)*dir;
      if(_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0] &&
         _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2] ) {
        foundQ = true;
      }
    } else if (dir[1] < 0 && !foundQ) {
      d = _mesh.getMax()[1];
      t = (d - origin[1])/dir[1];
      pt = origin + (t+step)*dir;
      if(_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0] &&
         _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }

    if (dir[2] > 0  && !foundQ) {
      d = _mesh.getMin()[2];
      t = (d - origin[2])/dir[2];
      pt = origin + (t+step)*dir;
      if(_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0] &&
         _mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]) {
        foundQ = true;
      }
    } else if (dir[2] < 0 && !foundQ) {
      d = _mesh.getMax()[2];
      t = (d - origin[2])/dir[2];
      pt = origin + (t+step)*dir;
      if(_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0] &&
         _mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]) {
        foundQ = true;
      }
    }

    if (!foundQ) {
      return false;
    }
  }

  while ((std::fabs(pt[0]-center[0]) < extends[0]) &&
         (std::fabs(pt[1]-center[1]) < extends[1]) &&
         (std::fabs(pt[2]-center[2]) < extends[2])) {
    node = _tree->find(pt);
    if (node->closestIntersection(origin, dir, q, a)) {
      std::vector<T> vek = q - origin;
      distance = norm(vek);
      return true;
    } else {
      Octree<T>* node = _tree->find(pt);
      node->intersectRayNode(pt, dir, s);
      for (int i=0; i<3; i++) {
        pt[i] = s[i] + step*dir[i];
      }
    }
  }
  clout << "Returning false" << std::endl;
  return false;
}

template <typename T>
void STLreader<T>::print() {
  _mesh.print();
  _tree->print();
  clout << "voxelSize=" << _voxelSize << "; stlSize=" << _stlSize << std::endl;
  clout << "minPhysR(VoxelMesh)=(" <<  this->_myMin[0] << "," << this->_myMin[1] << "," << this->_myMin[2] << ")";
  clout << "; maxPhysR(VoxelMesh)=(" <<  this->_myMax[0] << "," << this->_myMax[1] << "," << this->_myMax[2] << ")" << std::endl;
}

template <typename T>
void STLreader<T>::writeOctree() {
  _tree->write(_fName);
}

template <typename T>
void STLreader<T>::writeSTL() {
	_mesh.write(_fName);
}


} // namespace olb

#endif

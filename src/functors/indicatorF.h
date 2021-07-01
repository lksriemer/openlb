/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Cyril Masquelier, Mathias J. Krause
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

#ifndef INDICATOR_F_H
#define INDICATOR_F_H

#include<vector>
#include<cmath>

#include "functors/indicatorBaseF.h"
#include "io/ostreamManager.h"
#include "io/xmlReader.h"


/** \file
 * This file contains indicator functions. These return 1 if the given
 * coordinates are inside, and 0 if they are outside of the defined set.
 * Implemented are :
 - Sphere     3d
 - Cylinder   3d
 - Cone       3d
 - Pipe (not yet)
 - Cube (not yet)
 - Cuboid     2d, 3d
 - Circle     2d, 3d

 * The smoothIndicator functors return values in [0,1]. In particular there is
 * an epsilon enclosure of the set, wherein the return values are smooth and do
 * not jump from 0 to 1.

 Boolean operators allow to create unions and intersections. They can be used
 for example for initialization of a SuperGeometry.
*/

namespace olb {

template<typename S> class STLreader;


//////////////////////////////////2DIndicatorFunctors///////////////////////////

/// indicator function for a 2D-cuboid, parallel to the planes x=0, y=0;
template <typename T, typename S>
class IndicatorCuboid2D : public IndicatorF2D<T,S> {
private:
  std::vector<S>  _center;
  S               _xLength;
  S               _yLength;
public:
  /// constructs an cuboid with x axis dimension 0 to extend[0], ...
  IndicatorCuboid2D(std::vector<S> extend, std::vector<S> origin);
  /// constructs an cuboid with x axis dimension -xlength/2 to xlength/2
  IndicatorCuboid2D(S xlength, S ylength, std::vector<S> center = S());
  /// returns true if input is inside, otherwise false
  std::vector<T> operator() (std::vector<S> input);
};

/// indicator function for a 2D circle
template <typename T, typename S>
class IndicatorCircle2D : public IndicatorF2D<T,S> {
private:
  std::vector<S>  _center;
  S               _radius2;
public:
  IndicatorCircle2D(std::vector<S> center, S radius);
  std::vector<T> operator() (std::vector<S> input);
  virtual bool distance(S& distance, const std::vector<S>& origin,
                        const std::vector<S>& direction,  int iC=-1);
};



//////////////////////////////////3DIndicatorFunctors///////////////////////////

/// indicator function for a 3D circle
template <typename T, typename S>
class IndicatorCircle3D : public IndicatorF3D<T,S> {
private:
  std::vector<S>  _center;
  std::vector<S>  _normal;
  S               _radius2;
public:
  IndicatorCircle3D(std::vector<S> center, std::vector<S> normal, S radius);
  IndicatorCircle3D(S center0, S center1, S center2, S normal0, S normal1,
                    S normal2, S radius);

  std::vector<T> operator() (std::vector<S> input);
  std::vector<S> getCenter() { return _center; };
  std::vector<S> getNormal() { return _normal; };
  S getRadius() { return std::sqrt(_radius2); };
  //virtual bool distance(S& distance, std::vector<S> origin, std::vector<S> direction, int iC=-1);
};

/// indicator function for an object given by an stl file
template <typename T, typename S>
class IndicatorStl3D : public IndicatorF3D<T,S> {
private:
  STLreader<S>& _stlReader;
  S             _eps;
public:
  IndicatorStl3D(STLreader<S>& stlReader, S eps);
  std::vector<T> operator()(std::vector<S> x);
};

/// indicator function for a 3D-sphere
template <typename T, typename S>
class IndicatorSphere3D : public IndicatorF3D<T,S> {
private:
  std::vector<S>  _center;
  S               _radius2;
public:
  IndicatorSphere3D(std::vector<S> center, S radius);
  std::vector<T> operator() (std::vector<S> input);
  virtual bool distance(S& distance, const std::vector<S>& origin, const std::vector<S>& direction,
                        int iC=-1);
};

/// indicator function for a layer
template <typename T, typename S>
class IndicatorLayer3D : public IndicatorF3D<T,S> {
private:
  IndicatorIdentity3D<T,S>  _indicatorF;
  S                         _layerSize;
public:
  IndicatorLayer3D(IndicatorF3D<T,S>& indicatorF, S layerSize);
  std::vector<T> operator() (std::vector<S> input);
};

/// indicator function for a 3d-cylinder
template <typename T, typename S>
class IndicatorCylinder3D : public IndicatorF3D<T,S> {
private:
  std::vector<S>  _center1;
  std::vector<S>  _center2;
  std::vector<S>  _I;
  std::vector<S>  _J;
  std::vector<S>  _K;
  S               _length;
  S               _radius2;
  void init();
public:
  IndicatorCylinder3D(std::vector<S> center1, std::vector<S> center2, S radius);
  IndicatorCylinder3D(std::vector<S> center1, std::vector<S> normal, S radius, S eps);
  IndicatorCylinder3D(IndicatorCircle3D<T,S> circleF, S eps);
  std::vector<T> operator() (std::vector<S> input);
};

/// indicator function for a 3d frustum
template <typename T, typename S>
class IndicatorCone3D : public IndicatorF3D<T,S> {
private:
  std::vector<S>  _center1;
  std::vector<S>  _center2;
  std::vector<S>  _I;
  std::vector<S>  _J;
  std::vector<S>  _K;
  S               _length;
  S               _radius1;
  S               _radius2; // The 2nd radius is optional: if not defined, _center2 is the vertex of the cone
public:
  IndicatorCone3D(std::vector<S> center1, std::vector<S> center2,
                   S radius1, S radius2=0);
  std::vector<T> operator() (std::vector<S> input);
};

/// indicator function for a 3d-cuboid, parallel to the planes x=0, y=0, z=0;
template <typename T, typename S>
class IndicatorCuboid3D : public IndicatorF3D<T,S> {
private:
  std::vector<S>  _center;
  S               _xLength;
  S               _yLength;
  S               _zLength;
public:
  /// constructs an cuboid with x axis dimension 0 to extend[0], ...
  IndicatorCuboid3D(std::vector<S> extend, std::vector<S> origin);
  /// constructs an cuboid with x axis dimension -xlength/2 to xlength/2
  IndicatorCuboid3D(S xlength, S ylength, S zlength, std::vector<S> center);
  /// returns true if input is inside, otherwise false
  std::vector<T> operator() (std::vector<S> input);
};


template <typename T, typename S>
IndicatorCuboid3D<T,S>* createIndicatorCuboid3D(XMLreader const& params, bool verbose=false);

/// indicator function for a 3d-parallelepiped (including any cuboid)
template <typename T, typename S>
class IndicatorParallelepiped3D : public IndicatorF3D<T,S> {
private:
  std::vector<S>  _origin;   // A vertex of the parallelepiped
  std::vector<S>  _corner1;
  std::vector<S>  _corner2;
  std::vector<S>  _corner3;  // Those three vertex have a common edge with the origin.
  std::vector<S>  _I;
  std::vector<S>  _J;
  std::vector<S>  _K;
  S               _normI;
  S               _normJ;
  S               _normK;
public:
  IndicatorParallelepiped3D(std::vector<S> origin, std::vector<S> corner1,
                             std::vector<S> corner2, std::vector<S> corner3);
  std::vector<T> operator() (std::vector<S> input);
};





///////////////////////////SmoothIndicatorFunctors//////////////////////////////

/// implements a smooth circle in 2D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCircle2D : public SmoothIndicatorF2D<T,S> {
private:
  std::vector<S>  _center;
  S               _radius2;
  S               _epsilon;
public:
  SmoothIndicatorCircle2D(std::vector<S> center, S radius, S epsilon);
  std::vector<T> operator() (std::vector<S> input);
};


/// implements a smooth sphere in 3D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorSphere3D : public SmoothIndicatorF3D<T,S> {
private:
  std::vector<S>  _center;
  S               _radius2;
  S               _epsilon;
public:
  SmoothIndicatorSphere3D(std::vector<S> center, S radius, S epsilon);
  std::vector<T> operator() (std::vector<S> input);
};

/// implements a smooth cylinder in 3D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCylinder3D : public SmoothIndicatorF3D<T,S> {
private:
  std::vector<S>  _center1;
  std::vector<S>  _center2;
  std::vector<S>  _I;
  std::vector<S>  _J;
  std::vector<S>  _K;
  S               _length;
  S               _radius2;
  S               _epsilon;
public:
  SmoothIndicatorCylinder3D(std::vector<S> center1, std::vector<S> center2,
                            S radius, S epsilon);
  std::vector<T> operator() (std::vector<S> input);
};

/// implements a smooth cone in 3D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCone3D : public SmoothIndicatorF3D<T,S> {
private:
  std::vector<S>  _center1;
  std::vector<S>  _center2;
  std::vector<S>  _I;
  std::vector<S>  _J;
  std::vector<S>  _K;
  S               _length;
  S               _radius1;
  S               _radius2; // The 2nd radius is optional: if not defined, _center2 is the vertex of the cone
  S               _epsilon;
public:
  SmoothIndicatorCone3D(std::vector<S> center1, std::vector<S> center2,
                        S radius1, S radius2, S epsilon);
  std::vector<T> operator() (std::vector<S> input);
};


}

#endif


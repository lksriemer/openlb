/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Fietz, Mathias J. Krause
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


#ifndef HEURISTIC_LOAD_BALANCER_H
#define HEURISTIC_LOAD_BALANCER_H

#include "communication/mpiManager.h"
#include "geometry/cuboidGeometry2D.h"
#include "communication/cuboidNeighbourhood2D.h"
#include "geometry/cuboidGeometry3D.h"
#include "communication/cuboidNeighbourhood3D.h"
#include "communication/loadBalancer.h"
#include "io/ostreamManager.h"



namespace olb {

template <typename T> class CuboidGeometry2D;
template <typename T> class CuboidGeometry3D;

template<typename T>
class HeuristicLoadBalancer : public LoadBalancer<T> {
private:
  // Handles the MPI communication
#ifdef PARALLEL_MODE_MPI
  singleton::MpiNonBlockingHelper _mpiNbHelper;
#endif
  CuboidGeometry3D<T>* _cGeometry3d;
  CuboidGeometry2D<T>* _cGeometry2d;

public:
  /*
   * Constructs a load balancer from a given cuboid geometry
   * using a heurist
   * \param cGeometry: The cuboid geometry to base the load balance on
   * \param blockGeometry: Used to determine number of full and empty cells
   *             if given
   * \param ratioFullEmpty: Time it takes for full cells in relation to
   *              empty cells
   */
  HeuristicLoadBalancer() {};

  HeuristicLoadBalancer(CuboidGeometry3D<T>& cGeometry3d, const double ratioFullEmpty=3.7);

  HeuristicLoadBalancer(CuboidGeometry2D<T>& cGeometry2d, const double ratioFullEmpty=3.7);

  virtual void reInit(CuboidGeometry3D<T>& cGeometry3d, const double ratioFullEmpty=3.7);

  virtual void reInit(CuboidGeometry2D<T>& cGeometry2d, const double ratioFullEmpty=3.7);
};
}  // namespace olb

#endif

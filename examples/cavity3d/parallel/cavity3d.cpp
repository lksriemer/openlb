/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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

/* cavity3d.cpp:
 * This example illustrates a flow in a cuboid, lid-driven cavity.
 * This version is for parallel use. A version for sequential use
 * is also available.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor

const int N = 1; // resolution of the model
const int M = 1; // time discretization refinement
const T maxT = (T)30.; // max. simulation time in s, SI unit

void prepareGeometry(LBconverter<T> const& converter, IndicatorF3D<T>& indicator, SuperGeometry3D<T>& superGeometry) {

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  superGeometry.rename(0,2,indicator);
  superGeometry.rename(2,1,1,1,1);

  Vector<T,3> origin(T(), converter.getCharL(), T());
  Vector<T,3> extend(converter.getCharL(), converter.getLatticeL(), converter.getCharL());
  IndicatorCuboid3D<T> lid(extend,origin);

  superGeometry.rename(2,3,1,lid);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice(LBconverter<T> const& converter,
                    SuperLattice3D<T,DESCRIPTOR>& lattice,
                    Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                    SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getOmega();

  /// Material=0 -->do nothing
  lattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

  /// Material=1 -->bulk dynamics
  lattice.defineDynamics(superGeometry, 1, &bulkDynamics);

  /// Material=2,3 -->bulk dynamics, velocity boundary
  lattice.defineDynamics(superGeometry, 2, &bulkDynamics);
  lattice.defineDynamics(superGeometry, 3, &bulkDynamics);
  bc.addVelocityBoundary(superGeometry, 2, omega);
  bc.addVelocityBoundary(superGeometry, 3, omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(LBconverter<T> const&converter,
                       SuperLattice3D<T,DESCRIPTOR>& lattice, SuperGeometry3D<T>& superGeometry, int iT) {

  OstreamManager clout(std::cout,"setBoundaryValues");

  if(iT==0) {

    AnalyticalConst3D<T,T> rhoF(1);
    std::vector<T> velocity(3,T());
    AnalyticalConst3D<T,T> uF(velocity);

    lattice.iniEquilibrium(superGeometry, 1, rhoF, uF);
    lattice.iniEquilibrium(superGeometry, 2, rhoF, uF);
    lattice.iniEquilibrium(superGeometry, 3, rhoF, uF);

    lattice.defineRhoU(superGeometry, 1, rhoF, uF);
    lattice.defineRhoU(superGeometry, 2, rhoF, uF);
    lattice.defineRhoU(superGeometry, 3, rhoF, uF);

    velocity[0]=converter.getLatticeU();
    AnalyticalConst3D<T,T> u(velocity);

    lattice.defineU(superGeometry,3,u);

    /// Make the lattice ready for simulation
    lattice.initialize();
  }
}

void getResults(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                LBconverter<T> const& converter, SuperGeometry3D<T>& superGeometry, int iT, Timer<T>& timer) {

  OstreamManager clout(std::cout,"getResults");
  SuperVTKwriter3D<T> vtkWriter("cavity3d");

  const T logT     = (T)1.;
  const T vtkSave  = (T)1.;




  if (iT==0) {
    /// Writes the converter log file
    writeLogFile(converter, "cavity3d");

    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  /// Get statistics
  if (iT%converter.numTimeSteps(logT)==0 && iT>0) {
    timer.update(iT);
    timer.printStep(2);
    sLattice.getStatistics().print(iT,converter.physTime(iT));
  }

  /// Writes the VTK and GIF files
  if (iT%converter.numTimeSteps(vtkSave)==0 && iT>0) {
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( pressure );

    vtkWriter.write(iT);

    // define vector which span the gif-plane
    Vector<T,3> u(1,0,0);
    Vector<T,3> v(0,1,0);
    T tmp = T(converter.getCharL() / 2.);
    T origin[3] = {tmp,tmp,tmp};

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel(velocity);
    BlockLatticeReduction3D<T, DESCRIPTOR> planeReduction(normVel, u, v, 600, origin);
    BlockGifWriter<T> gifWriter;

    gifWriter.write(planeReduction, iT);
  }
}



int main(int argc, char **argv) {

  /// === 1st Step: Initialization ===

  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout,"main");

  LBconverter<T> converter(
    (int) 3,                               // dim
    (T)   1./30./N,                        // latticeL_
    (T)   1e-1/M,                          // latticeU_
    (T)   1./1000.,                        // charNu_
    (T)   1.,                              // charL_ = 1
    (T)   1.                               // charU_ = 1
  );

  /// === 2nd Step: Prepare Geometry ===

  /// Instantiation of a unit cube by an indicator

  std::vector<T> origin(3,T());
  std::vector<T> extend(3,converter.getCharL());
  IndicatorCuboid3D<T> cube(extend,origin);

  /// Instantiation of a cuboid geometry with weights
  int noCuboids = singleton::mpi().getSize();
  if (noCuboids<7) {
    noCuboids = 7;
  }
  CuboidGeometry3D<T> cuboidGeometry(cube, converter.getLatticeL(), noCuboids);

  /// Instantiation of a load balancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a super geometry
  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(converter, cube, superGeometry);


  /// === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

  ConstRhoBGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter.getOmega(), instances::getBulkMomenta<T,DESCRIPTOR>());

  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition(sLattice);
  createInterpBoundaryCondition3D<T,DESCRIPTOR,ConstRhoBGKdynamics<T,DESCRIPTOR> >(sBoundaryCondition);

  prepareLattice(converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry);

  /// === 4th Step: Main Loop with Timer ===

  Timer<T> timer(converter.numTimeSteps(maxT), converter.numNodes(1)*converter.numNodes(1)*converter.numNodes(1) );
  timer.start();
  int iT;

  for (iT=0; iT < converter.numTimeSteps(maxT); ++iT) {

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, sLattice, superGeometry, iT);

    /// === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, superGeometry, iT, timer);
  }

  timer.stop();
  timer.printSummary();
}


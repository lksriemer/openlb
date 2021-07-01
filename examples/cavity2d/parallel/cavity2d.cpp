/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006 - 2012 Mathias J. Krause, Jonas Fietz,
 *  Jonas Latt, Jonas Kratzke
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *••••••••••
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

/* cavity2d.cpp:
 * This example illustrates a flow in a cuboid, lid-driven cavity.
 * It also shows how to use the XML parameter files and has an
 * example description file for OpenGPI. This version is for parallel
 * use. A version for sequential use is also available.
 */


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

void prepareGeometry(LBconverter<T> const* converter,
                     SuperGeometry2D<T>& superGeometry)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  int N = converter->numNodes();
  std::vector<T> extend(2,T());
  extend[0] = (N-1)*converter->getLatticeL();
  std::vector<T> origin(2,T());
  origin[1] = (N-1)*converter->getLatticeL();
  IndicatorCuboid2D<T> lid(extend, origin);

  superGeometry.rename(0,2);

  superGeometry.rename(2,1,1,1);

  superGeometry.clean();

  /// Set material number for lid
  superGeometry.rename(2,3,1,lid);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(LBconverter<T> const* converter,
                    SuperLattice2D<T, DESCRIPTOR>& sLattice,
                    Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
                    SuperGeometry2D<T>& superGeometry) {

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter->getOmega();

  // link lattice with dynamics for collision step

  /// Material=0 -->do nothing
  sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

  /// Material=1 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 1, &bulkDynamics);

  /// Material=2,3 -->bulk dynamics, velocity boundary
  sLattice.defineDynamics(superGeometry, 2, &bulkDynamics);
  sLattice.defineDynamics(superGeometry, 3, &bulkDynamics);
  sBoundaryCondition.addVelocityBoundary(superGeometry, 2, omega);
  sBoundaryCondition.addVelocityBoundary(superGeometry, 3, omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(LBconverter<T> const* converter,
                       SuperLattice2D<T, DESCRIPTOR>& sLattice,
                       int iT, SuperGeometry2D<T>& superGeometry) {

  if(iT==0) {

    // set initial values: v = [0,0]
    AnalyticalConst2D<T,T> rhoF(1);
    std::vector<T> velocity(2,T());
    AnalyticalConst2D<T,T> uF(velocity);

    sLattice.iniEquilibrium(superGeometry, 1, rhoF, uF);
    sLattice.iniEquilibrium(superGeometry, 2, rhoF, uF);
    sLattice.iniEquilibrium(superGeometry, 3, rhoF, uF);

    sLattice.defineRhoU(superGeometry, 1, rhoF, uF);
    sLattice.defineRhoU(superGeometry, 2, rhoF, uF);
    sLattice.defineRhoU(superGeometry, 3, rhoF, uF);

    // set non-zero velocity for upper boundary cells
    velocity[0]=converter->getLatticeU();
    AnalyticalConst2D<T,T> u(velocity);

    sLattice.defineU(superGeometry,3,u);

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

void getResults(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                LBconverter<T> const* converter, int iT, Timer<T>* timer,
                const T logT, const T imSave, const T vtkSave, const T maxT,
                std::string filenameGif, std::string filenameVtk,
                const bool timerPrintSum, const int timerPrintMode,
                const int timerTimeSteps, SuperGeometry2D<T>& superGeometry) {

  OstreamManager clout(std::cout,"getResults");

  SuperVTKwriter2D<T> vtkWriter("cavity2d");
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity(sLattice, *converter);
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure(sLattice, *converter);
  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(pressure);

  if (iT==0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  /// Get statistics
  if (iT%converter->numTimeSteps(logT)==0) {
    sLattice.getStatistics().print(iT, converter->physTime(iT));
  }

  if (iT%timerTimeSteps==0) {
    timer->print(iT,timerPrintMode);
  }

  /// Writes the VTK files
  if (iT%converter->numTimeSteps(imSave)==0 && iT>0) {
    vtkWriter.write(iT);
  }

  /// Writes the Gif files
  if (iT%converter->numTimeSteps(imSave)==0 && iT>0) {
    SuperEuklidNorm2D<T,DESCRIPTOR> normVel(velocity);
    BlockLatticeReduction2D<T, DESCRIPTOR> planeReduction(normVel);
    BlockGifWriter<T> gifWriter;
//    gifWriter.write(planeReduction, 0, 1, iT, "artefacts");
    gifWriter.write(planeReduction, iT, "vel");
  }

  if (iT == converter->numTimeSteps(maxT)-1 ) {
    timer->stop();
    if (timerPrintSum==true) {
      timer->printSummary();
    }
  }

}



int main(int argc, char* argv[]) {

  /// === 1st Step: Initialization ===

  olbInit(&argc, &argv);
  OstreamManager clout(std::cout,"main");

  string fName("cavity2d.xml");
  XMLreader config(fName);

  bool multiOutput=false;
  config["Output"]["MultiOutput"].read(multiOutput);
  clout.setMultiOutput(multiOutput);

  std::string olbdir, outputdir;
  config["Application"]["OlbDir"].read(olbdir);
  config["Output"]["OutputDir"].read(outputdir);
  singleton::directories().setOlbDir(olbdir);
  singleton::directories().setOutputDir(outputdir);

  // call creator functions using xml data
  LBconverter<T>* converter = createLBconverter<T>(config);
  int N = converter->numNodes();
  Timer<T>* timer           = createTimer<T>(config, *converter, N*N);

  T logT;
  T imSave;
  T vtkSave;
  T maxT;
  int timerSkipType;
  bool timerPrintSum = true;
  int timerPrintMode = 0;
  int timerTimeSteps = 1;

  config["Application"]["PhysParam"]["MaxTime"].read(maxT);
  config["Output"]["Log"]["SaveTime"].read(logT);
  config["Output"]["VisualizationImages"]["SaveTime"].read(imSave);
  config["Output"]["ViutputsualizationVTK"]["SaveTime"].read(vtkSave);
  config["Output"]["Timer"]["SkipType"].read(timerSkipType);
  config["Output"]["Timer"]["PrintSummary"].read(timerPrintSum);
  config["Output"]["Timer"]["PrintMode"].read(timerPrintMode);

  if (timerSkipType == 0) {
    timerTimeSteps=converter->numTimeSteps(config["Output"]["Timer"]["PhysTime"].get<T>());
  } else {
    config["Output"]["Timer"]["TimeSteps"].read(timerTimeSteps);
  }

  writeLogFile(*converter, config["Output"]["Log"]["Filename"].get<std::string>());

  std::string filenameGif = config["Output"]["VisualizationImages"]["Filename"].get<std::string>();

  std::string filenameVtk = config["Output"]["VisualizationVTK"]["Filename"].get<std::string>();

  /// === 2rd Step: Prepare Geometry ===
  /// Instantiation of a cuboidGeometry with weights

  std::vector<T> extend(2,T());
  extend[0] = (N-1)*converter->getLatticeL();
  extend[1] = (N-1)*converter->getLatticeL();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter->getLatticeL(), singleton::mpi().getSize());
#else
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter->getLatticeL(), 7);
#endif

  cuboidGeometry.print();

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(converter, superGeometry);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice2D<T, DESCRIPTOR> sLattice(superGeometry);

  ConstRhoBGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter->getOmega(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sBoundaryCondition(sLattice);
  createInterpBoundaryCondition2D<T,DESCRIPTOR,ConstRhoBGKdynamics<T,DESCRIPTOR> > (sBoundaryCondition);

  prepareLattice(converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry);

  /// === 4th Step: Main Loop with Timer ===

  timer->start();

  int iT=0;
  for (iT=0; converter->physTime(iT)<maxT; ++iT) {


    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, sLattice, iT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, timer, logT, imSave, vtkSave, maxT, filenameGif, filenameVtk,
               timerPrintSum, timerPrintMode, timerTimeSteps, superGeometry);


  }

  delete converter;
  delete timer;
}

/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014 Mathias J. Krause, Thomas Henn,
 *  Cyril Masquelier
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

/* venturi3d.cpp:
 * This example examines a steady flow in a venturi tube. At the
 * main inlet, a Poiseuille profile is imposed as Dirichlet velocity
 * boundary condition, whereas at the outlet and the minor inlet 
 * a Dirichlet pressure condition is set by p=0 (i.e. rho=1).
 *
 * The example shows the usage of the Indicator functors to
 * build up a geometry and explains how to set boundary conditions
 * automatically.
 */

#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used
#include "olb3D.hh"     // Include full template code
#endif
#include <vector>
//#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor

int N = 1; // resolution of the model

SuperGeometry3D<T> prepareGeometry( LBconverter<T> const& converter) {

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  /// Definition of the geometry of the venturi
  std::vector<T> C0;
  int dummy00[] = {0,50,50};
  C0.insert(C0.begin(),dummy00, dummy00 +3);

  std::vector<T> C1;
  int dummy0[] = {5,50,50};
  C1.insert(C1.begin(),dummy0, dummy0 +3);

  std::vector<T> C2;
  int dummy1[] = {40,50,50};
  C2.insert(C2.begin(),dummy1, dummy1 +3);

  std::vector<T> C3;
  int dummy2[] = {80,50,50};
  C3.insert(C3.begin(),dummy2, dummy2 +3);

  std::vector<T> C4;
  int dummy3[] = {120,50,50};
  C4.insert(C4.begin(),dummy3, dummy3 +3);

  std::vector<T> C5;
  int dummy4[] = {160,50,50};
  C5.insert(C5.begin(),dummy4, dummy4 +3);

  std::vector<T> C6;
  int dummy5[] = {195,50,50};
  C6.insert(C6.begin(),dummy5, dummy5 +3);

  std::vector<T> C7;
  int dummy6[] = {200,50,50};
  C7.insert(C7.begin(),dummy6, dummy6 +3);

  std::vector<T> C8;
  int dummy7[] = {190,50,50};
  C8.insert(C8.begin(),dummy7, dummy7 +3);

  std::vector<T> C9;
  int dummy8[] = {115,50,50};
  C9.insert(C9.begin(),dummy8, dummy8 +3);

  std::vector<T> C10;
  int dummy9[] = {115,25,50};
  C10.insert(C10.begin(),dummy9, dummy9 +3);

  std::vector<T> C11;
  int dummy10[] = {115,5,50};
  C11.insert(C11.begin(),dummy10, dummy10 +3);

  std::vector<T> C12;
  int dummy12[] = {115,3,50};
  C12.insert(C12.begin(),dummy12, dummy12 +3);

  std::vector<T> C13;
  int dummy13[] = {115,7,50};
  C13.insert(C13.begin(),dummy13, dummy13 +3);

  T radius1 = 10 ;  // radius of the tightest part
  T radius2 = 20 ;  // radius of the widest part
  T radius3 = 4 ;   // radius of the small exit

  IndicatorCylinder3D<bool,T> inflow(C0, C1, radius2);
  IndicatorCylinder3D<bool,T> cyl1(C1, C2, radius2);
  IndicatorCone3D<bool,T> co1(C2, C3, radius2, radius1);
  IndicatorCylinder3D<bool,T> cyl2(C3, C4, radius1);
  IndicatorCone3D<bool,T> co2(C4, C5, radius1, radius2);
  IndicatorCylinder3D<bool,T> cyl3(C5, C6, radius2);
  IndicatorCylinder3D<bool,T> outflow0(C7, C8, radius2);
  IndicatorCylinder3D<bool,T> cyl4(C9, C10, radius3);
  IndicatorCone3D<bool,T> co3(C10, C11, radius3, radius1);
  IndicatorCylinder3D<bool,T> outflow1(C12, C13, radius1);

  IndicatorIdentity3D<bool,T> venturi(cyl1 + cyl2 + cyl3 + cyl4 + co1 + co2 + co3);


  /// Build CoboidGeometry from IndicatorF (weights are set, remove and shrink is done)   
  CuboidGeometry3D<T>* cuboidGeometry = new CuboidGeometry3D<T>(venturi, 1./N, 20*singleton::mpi().getSize() );

  /// Build LoadBalancer from CuboidGeometry (weights are respected) 
  HeuristicLoadBalancer<T>* loadBalancer = new HeuristicLoadBalancer<T>(*cuboidGeometry);

  /// Default instantiation of superGeometry
  SuperGeometry3D<T> superGeometry(*cuboidGeometry, *loadBalancer, 2);

  /// Set boundary voxels by rename material numbers

  superGeometry.rename(0,2,venturi);
  superGeometry.rename(2,1,1,1,1);

  superGeometry.rename(2,3,1,inflow);
  superGeometry.rename(2,4,1,outflow0);
  superGeometry.rename(2,5,1,outflow1);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();


  superGeometry.getStatistics().print();
  superGeometry.communicate();

  clout << "Prepare Geometry ... OK" << std::endl;
  return superGeometry;
}


void prepareLattice(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                    LBconverter<T> const& converter,
                    Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                    sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                    SuperGeometry3D<T>& superGeometry) {



  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getOmega();

  /// Material=0 -->do nothing
  sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

  /// Material=1 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 1, &bulkDynamics);

  /// Material=2 -->bounce back
  sLattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>());

  /// Material=3 -->bulk dynamics (inflow)
  sLattice.defineDynamics(superGeometry, 3, &bulkDynamics);

  /// Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics(superGeometry, 4, &bulkDynamics);

  /// Material=5 -->bulk dynamics (2nd outflow)
  sLattice.defineDynamics(superGeometry, 5, &bulkDynamics);

  /// Setting of the boundary conditions
  bc.addVelocityBoundary(superGeometry, 3, omega);
  bc.addPressureBoundary(superGeometry, 4, omega);
  bc.addPressureBoundary(superGeometry, 5, omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Generates a slowly increasing sinuidal inflow for the first iTMax timesteps
void setBoundaryValues(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                       LBconverter<T> const& converter, int iT, T maxPhysT,
                       SuperGeometry3D<T>& superGeometry) {

  OstreamManager clout(std::cout,"setBoundaryValues");

  // No of time steps for smooth start-up
  int iTmaxStart = converter.numTimeSteps(maxPhysT*0.8);
  int iTperiod = 50;

  if (iT==0) {
    /// Make the lattice ready for simulation
    sLattice.initialize();
  }

  else if (iT%iTperiod==0 && iT<= iTmaxStart) {
    //clout << "Set Boundary Values ..." << std::endl;

    //SinusStartScale<T,int> startScale(iTmaxStart, (T) 1);
    PolynomialStartScale<T,int> startScale(iTmaxStart, T(1));
    std::vector<int> iTvec(1,iT);
    T frac = startScale(iTvec)[0];

    // Creates and sets the Poiseuille inflow profile using functors
    std::vector<T> maxVelocity(3,0);
    maxVelocity[0] = /*2.25*/frac*converter.getLatticeU();

    RectanglePoiseuille3D<T> poiseuilleU(superGeometry, 3, maxVelocity, converter.getLatticeL(), converter.getLatticeL(), converter.getLatticeL());
    sLattice.defineU(superGeometry, 3, poiseuilleU);

    //clout << "step=" << iT << "; scalingFactor=" << frac << std::endl;
  }
  //clout << "Set Boundary Values ... ok" << std::endl;
}

void getResults(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                LBconverter<T>& converter, int iT,
                SuperGeometry3D<T>& superGeometry, Timer<double>& timer) {

  OstreamManager clout(std::cout,"getResults");
  SuperVTKwriter3D<T,DESCRIPTOR> vtkWriter("venturi3d");

  if (iT==0) {
    /// Writes the converter log file
    writeLogFile(converter, "venturi3d");
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  /// Writes the vtk files
  if (iT%converter.numTimeSteps(1.)==0) {
    // Create the data-reading functors...
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( pressure );
    vtkWriter.write(iT);
  }

  /// Writes output on the console
  if (iT%converter.numTimeSteps(1.)==0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.physTime(iT));

  }
}



int main(int argc, char* argv[]) {

  /// === 1st Step: Initialization ===

  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout,"main");
  // display messages from every single mpi process
  // clout.setMultiOutput(true);

  LBconverter<T> converter(
    (int) 3,                               // dim
    (T)   1./N,                            // latticeL_
    (T)   0.02,                            // latticeU_
    (T)   0.1,                             // charNu_
    (T)   0.1,                             // charL_ = 1
    (T)   2.                               // charU_ = 1
  );
  T maxPhysT = 200.0;                      // max. simulation time in s, SI unit

  /// === 2nd Step: Prepare Geometry ===

  SuperGeometry3D<T> superGeometry(prepareGeometry(converter) );

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

  RLBdynamics<T, DESCRIPTOR> bulkDynamics(converter.getOmega(), instances::getBulkMomenta<T, DESCRIPTOR>());

  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition(sLattice);
  createInterpBoundaryCondition3D<T, DESCRIPTOR> (sBoundaryCondition);

  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition(sLattice);
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR> (sOffBoundaryCondition);

  prepareLattice(sLattice, converter, bulkDynamics, sBoundaryCondition, sOffBoundaryCondition, superGeometry);

  Timer<double> timer(converter.numTimeSteps(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();
  getResults(sLattice, converter, 0, superGeometry, timer);

  /// === 4th Step: Main Loop with Timer ===
  for (int iT = 0; iT <= converter.numTimeSteps(maxPhysT); ++iT) {

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(sLattice, converter, iT, maxPhysT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer);
  }

  timer.stop();
  timer.printSummary();
}




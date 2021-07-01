/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Orestis Malaspinas, Jonas Latt,
 *  Mathias J. Krause, Vojtech Cvrcek, Peter Weisbrod
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

/* mrtPoiseuille2d.cpp:
 * This example examines a 2D Poseuille flow with a velocity
 * or pressure boundary at the inlet/outlet.
 * Computation of error norms via functors is shown as well as
 * the use of MRT dynamics.
 */


#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef enum {velocity, pressure} BoundaryType;

typedef double T;
#define DESCRIPTOR MRTD2Q9Descriptor


// Parameters for the simulation setup
const T lx  = 2.;
const T ly  = 1.;
const int N = 50;      // resolution of the model
const T Re = 10.;      // Reynolds number
const T maxPhysT = 5.; // max. simulation time in s, SI unit


/// Stores geometry information in form of material numbers
void prepareGeometry(LBconverter<T> const& converter,
                     SuperGeometry2D<T>& superGeometry) {

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);

  superGeometry.rename(2,1,1,1);

  std::vector<T> extend(2,T());
  std::vector<T> origin(2,T());

  /// Set material number for inflow
  extend[0] = 0; extend[1] = ly; IndicatorCuboid2D<bool,T> inflow(extend, origin);
  superGeometry.rename(2,3,1,inflow);

  /// Set material number for outflow
  origin[0] = lx; IndicatorCuboid2D<bool,T> outflow(extend, origin);
  superGeometry.rename(2,4,1,outflow);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// Set up the geometry of the simulation
void prepareLattice(LBconverter<T> const& converter,
                    SuperLattice2D<T, DESCRIPTOR>& sLattice,
                    Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
                    BoundaryType inflowBoundary, BoundaryType outflowBoundary,
                    SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  T   omega = converter.getOmega();

  /// Material=0 -->do nothing
  sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

  /// Material=1 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 1, &bulkDynamics);

  /// Material=2 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 2, &bulkDynamics);

  /// Material=3 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 3, &bulkDynamics);

  /// Material=4 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 4, &bulkDynamics);

  /// Setting of the boundary conditions
  sBoundaryCondition.addVelocityBoundary(superGeometry, 2, omega);

  if (inflowBoundary==velocity) {
    sBoundaryCondition.addVelocityBoundary(superGeometry, 3, omega);
  }
  else {
    sBoundaryCondition.addPressureBoundary(superGeometry, 3, omega);
  }

  if (outflowBoundary==velocity) {
    sBoundaryCondition.addVelocityBoundary(superGeometry, 4, omega);
  }
  else {
    sBoundaryCondition.addPressureBoundary(superGeometry, 4, omega);
  }

  /// Initial conditions
  T Lx = converter.numCells(lx);
  T Ly = converter.numCells(ly);

  T p0 =8.*converter.getLatticeNu()*converter.getLatticeU()*Lx/(Ly*Ly);
  AnalyticalLinear2D<T,T> rho(-p0/lx*DESCRIPTOR<T>::invCs2 , 0 , p0*DESCRIPTOR<T>::invCs2+1 );

  const T maxVelocity = converter.getLatticeU();
  T distance2Wall = converter.getLatticeL();
  Poiseuille2D<T> u(superGeometry, 3, maxVelocity, distance2Wall);
    
  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry, 1, rho, u);
  sLattice.iniEquilibrium(superGeometry, 1, rho, u);
  sLattice.defineRhoU(superGeometry, 2, rho, u);
  sLattice.iniEquilibrium(superGeometry, 2, rho, u);
  sLattice.defineRhoU(superGeometry, 3, rho, u);
  sLattice.iniEquilibrium(superGeometry, 3, rho, u);
  sLattice.defineRhoU(superGeometry, 4, rho, u);
  sLattice.iniEquilibrium(superGeometry, 4, rho, u);

  /// Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Compute error norms
void error(SuperGeometry2D<T>& superGeometry,
           SuperLattice2D<T, DESCRIPTOR>& sLattice,
           LBconverter<T> const& converter,
           Dynamics<T, DESCRIPTOR>& bulkDynamics) {

  OstreamManager clout(std::cout,"error");

  std::vector<int> tmp; T result; T result1;
  std::vector<double> vec(2,0);

  const T maxVelocity = converter.physVelocity(converter.getLatticeU());
  T distance2Wall = converter.getLatticeL();
  Poiseuille2D<T> uSol(superGeometry, 3, maxVelocity, distance2Wall);
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u(sLattice,converter);
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> uSolLattice(uSol,sLattice,superGeometry);

  T Lx = converter.numCells(lx);
  T Ly = converter.numCells(ly);
  T p0 =8.*converter.getLatticeNu()*converter.getLatticeU()*Lx/(Ly*Ly);
  AnalyticalLinear2D<T,T> pressureSol(-converter.physPressureFromRho(p0*DESCRIPTOR<T>::invCs2+1)/lx , 0 , converter.physPressureFromRho(p0*DESCRIPTOR<T>::invCs2+1) );
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure(sLattice,converter);
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> pressureSolLattice(pressureSol,sLattice,superGeometry);

  PoiseuilleStrainRate2D<T,T> sSol(converter, ly);
  SuperLatticePhysStrainRate2D<T,DESCRIPTOR> s(sLattice,converter);
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> sSolLattice(sSol,sLattice,superGeometry);


  // velocity error
  SuperL1Norm2D<T,DESCRIPTOR> uL1Norm(uSolLattice-u,superGeometry,1);
  SuperL1Norm2D<T,DESCRIPTOR> uSolL1Norm(uSolLattice,superGeometry,1);
  result = uL1Norm(tmp)[0]; result1=result/uSolL1Norm(tmp)[0];
  clout << "velocity-L1-error(abs)=" << result << "; velocity-L1-error(rel)=" << result1 << std::endl;

  SuperL2Norm2D<T,DESCRIPTOR> uL2Norm(uSolLattice-u,superGeometry,1);
  SuperL2Norm2D<T,DESCRIPTOR> uSolL2Norm(uSolLattice,superGeometry,1);
  result = uL2Norm(tmp)[0]; result1=result/uSolL2Norm(tmp)[0];
  clout << "velocity-L2-error(abs)=" << result << "; velocity-L2-error(rel)=" << result1 << std::endl;

  SuperLinfNorm2D<T,DESCRIPTOR> uLinfNorm(uSolLattice-u,superGeometry,1);
  result = uLinfNorm(tmp)[0];
  clout << "velocity-Linf-error(abs)=" << result << std::endl;

  // density error
  SuperL1Norm2D<T,DESCRIPTOR> pressureL1Norm(pressureSolLattice-pressure,superGeometry,1);
  SuperL1Norm2D<T,DESCRIPTOR> pressureSolL1Norm(pressureSolLattice,superGeometry,1);
  result = pressureL1Norm(tmp)[0]; result1=result/pressureSolL1Norm(tmp)[0];
  clout << "pressure-L1-error(abs)=" << result << "; pressure-L1-error(rel)=" << result1 << std::endl;

  SuperL2Norm2D<T,DESCRIPTOR> pressureL2Norm(pressureSolLattice-pressure,superGeometry,1);
  SuperL2Norm2D<T,DESCRIPTOR> pressureSolL2Norm(pressureSolLattice,superGeometry,1);
  result = pressureL2Norm(tmp)[0]; result1=result/pressureSolL2Norm(tmp)[0];
  clout << "pressure-L2-error(abs)=" << result << "; pressure-L2-error(rel)=" << result1 << std::endl;

  SuperLinfNorm2D<T,DESCRIPTOR> pressureLinfNorm(pressureSolLattice-pressure,superGeometry,1);
  result = pressureLinfNorm(tmp)[0];
  clout << "pressure-Linf-error(abs)=" << result << std::endl;

  // strain rate error
  SuperL1Norm2D<T,DESCRIPTOR> sL1Norm(sSolLattice-s,superGeometry,1);
  SuperL1Norm2D<T,DESCRIPTOR> sSolL1Norm(sSolLattice,superGeometry,1);
  result = sL1Norm(tmp)[0]; result1=result/sSolL1Norm(tmp)[0];
  clout << "strainRate-L1-error(abs)=" << result << "; strainRate-L1-error(rel)=" << result1 << std::endl;

  SuperL2Norm2D<T,DESCRIPTOR> sL2Norm(sSolLattice-s,superGeometry,1);
  SuperL2Norm2D<T,DESCRIPTOR> sSolL2Norm(sSolLattice,superGeometry,1);
  result = sL2Norm(tmp)[0]; result1=result/sSolL2Norm(tmp)[0];
  clout << "strainRate-L2-error(abs)=" << result << "; strainRate-L2-error(rel)=" << result1 << std::endl;

  SuperLinfNorm2D<T,DESCRIPTOR> sLinfNorm(sSolLattice-s,superGeometry,1);
  result = sLinfNorm(tmp)[0];
  clout << "strainRate-Linf-error(abs)=" << result << std::endl;
}

/// Output to console and files
void getResults(SuperLattice2D<T,DESCRIPTOR>& sLattice, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                LBconverter<T>& converter, int iT,
                SuperGeometry2D<T>& superGeometry, Timer<double>& timer) {

  OstreamManager clout(std::cout,"getResults");

  SuperVTKwriter2D<T,DESCRIPTOR> vtkWriter("mrtPoiseuille2d");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( pressure );

  const int vtkIter  = converter.numTimeSteps(maxPhysT/20.);
  const int statIter = converter.numTimeSteps(maxPhysT/20.);

  if (iT==0) {
    /// Writes the converter log file
    writeLogFile(converter, "mrtPoiseuille2d");

    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);  
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the vtk files and profile text file
  if (iT%vtkIter==0) {
    vtkWriter.write(iT);
    ofstream *ofile = 0;
    if (singleton::mpi().isMainProcessor()) {
      ofile = new ofstream((singleton::directories().getLogOutDir()+"centerVel.dat").c_str());
    }
    T Ly = converter.numCells(ly);
    for (int iY=0; iY<=Ly; ++iY) {
      T dx = converter.getDeltaX();
      const T maxVelocity = converter.physVelocity(converter.getLatticeU());
      std::vector<T> point(2,T());
      point[0] = lx/2.;
      point[1] = (T)iY/Ly;
      T distance2Wall = converter.getLatticeL();
      Poiseuille2D<T> uSol(superGeometry, 3, maxVelocity, distance2Wall);
      T analytical = uSol(point)[0];
      SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
      AnalyticalFfromSuperLatticeF2D<T, DESCRIPTOR> intpolateVelocity(velocity, true);
      T numerical = intpolateVelocity(point)[0];
      if (singleton::mpi().isMainProcessor()) {
        *ofile << iY*dx << " " << analytical
               << " " << numerical << "\n";
      }
    }
    delete ofile;
  }

  /// Writes output on the console
  if (iT%statIter==0) {
    /// Timer console output
    timer.update(iT);
    timer.printStep();

    /// Lattice statistics console output
    sLattice.getStatistics().print(iT,converter.physTime(iT));

    /// Error norms
    error(superGeometry, sLattice, converter, bulkDynamics);
  }
}


int main(int argc, char* argv[]) {

  /// === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout,"main");
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  const BoundaryType inflowBoundary = velocity;
  const BoundaryType outflowBoundary = pressure;

  LBconverter<T> converter(
    (int) 2,                               // dim
    1./N,                                  // latticeL_
    1./N,                                  // latticeU_
    (T)   1./Re,                           // charNu_
    (T)   1.                               // charL_ = 1,
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T()); extend[0] = lx; extend[1] = ly;
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<bool,T> cuboid(extend, origin);

  /// Instantiation of a cuboidGeometry with weights
  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
  #else
    const int noOfCuboids = 7;
  #endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getLatticeL(), noOfCuboids);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(converter, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice(superGeometry);

  MRTdynamics<T, DESCRIPTOR> bulkDynamics (
    converter.getOmega(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  // Choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sBoundaryCondition(sLattice);
  //createLocalBoundaryCondition2D<T,DESCRIPTOR,MRTdynamics<T,DESCRIPTOR> >(sBoundaryCondition);
  createInterpBoundaryCondition2D<T,DESCRIPTOR,MRTdynamics<T,DESCRIPTOR> >(sBoundaryCondition);

  prepareLattice(converter, sLattice, bulkDynamics, sBoundaryCondition, inflowBoundary, outflowBoundary, superGeometry);

  /// === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer(converter.numTimeSteps(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (int iT = 0; iT < converter.numTimeSteps(maxPhysT); ++iT) {

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    /// === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, bulkDynamics, converter, iT, superGeometry, timer);
  }

  timer.stop();
  timer.printSummary();
}

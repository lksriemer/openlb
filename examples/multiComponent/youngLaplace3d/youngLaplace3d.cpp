/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Robin Trunk
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

/* youngLaplace3d.cpp
 * In this example a Young-Laplace test is performed. A spherical domain
 * of fluid 2 is immersed in fluid 1. A diffusive interface forms and the
 * surface tension can be calculated using the Laplace pressure relation.
 * The pressure difference is calculated between a point in the middle of
 * the circular domain and a point furthest away from it in the
 * computational domain (here left bottom corner).
 *
 * This example shows the simplest case for the free-energy model with two
 * fluid components.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
typedef D3Q19<CHEM_POTENTIAL,FORCE> DESCRIPTOR;

// Parameters for the simulation setup
const int N  = 100;
const T nx   = 100.;
const T radius = 0.25 * nx;
const T alpha = 1.5;     // Interfacial width         [lattice units]
const T kappa1 = 0.0075; // For surface tensions      [lattice units]
const T kappa2 = 0.005;  // For surface tensions      [lattice units]
const T gama = 1.;       // For mobility of interface [lattice units]

const int maxIter  = 60000;
const int vtkIter  = 200;
const int statIter = 200;


void prepareGeometry( SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,1 );

  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( SuperLattice<T, DESCRIPTOR>& sLattice1,
                     SuperLattice<T, DESCRIPTOR>& sLattice2,
                     UnitConverter<T, DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  sLattice1.defineDynamics<ForcedBGKdynamics>(superGeometry, 1);
  sLattice2.defineDynamics<FreeEnergyBGKdynamics>(superGeometry, 1);

  // bulk initial conditions
  // define spherical domain for fluid 2
  std::vector<T> v( 3,T() );
  AnalyticalConst3D<T,T> zeroVelocity( v );

  AnalyticalConst3D<T,T> one ( 1. );
  IndicatorSphere3D<T> sphere( {nx/T(2), nx/T(2), nx/T(2)}, radius );
  SmoothIndicatorSphere3D<T,T> smoothSphere( sphere, 10.*alpha );

  AnalyticalIdentity3D<T,T> rho( one );
  AnalyticalIdentity3D<T,T> phi( one - smoothSphere - smoothSphere );

  sLattice1.iniEquilibrium( superGeometry, 1, rho, zeroVelocity );
  sLattice2.iniEquilibrium( superGeometry, 1, phi, zeroVelocity );

  sLattice1.setParameter<descriptors::OMEGA>( converter.getLatticeRelaxationFrequency() );
  sLattice2.setParameter<descriptors::OMEGA>( converter.getLatticeRelaxationFrequency() );
  sLattice2.setParameter<collision::FreeEnergy::GAMMA>(gama);

  sLattice1.initialize();
  sLattice2.initialize();

  {
    auto& communicator = sLattice1.getCommunicator(stage::PostPostProcess());
    communicator.requestField<POPULATION>();
    communicator.requestOverlap(sLattice1.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice2.getCommunicator(stage::PostPostProcess());
    communicator.requestField<POPULATION>();
    communicator.requestOverlap(sLattice2.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults( SuperLattice<T, DESCRIPTOR>& sLattice2,
                 SuperLattice<T, DESCRIPTOR>& sLattice1, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                 UnitConverter<T, DESCRIPTOR> converter)
{

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "youngLaplace3d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice1 );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice1 );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice1.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice2.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    AnalyticalConst3D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperLatticeDensity3D<T, DESCRIPTOR> density1( sLattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity3D<T, DESCRIPTOR> density2( sLattice2 );
    density2.getName() = "phi";

    SuperIdentity3D<T,T> c1 (half*(density1+density2));
    c1.getName() = "density-fluid-1";
    SuperIdentity3D<T,T> c2 (half*(density1-density2));
    c2.getName() = "density-fluid-2";

    vtmWriter.addFunctor( density1 );
    vtmWriter.addFunctor( density2 );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );
    vtmWriter.write( iT );

    // calculate bulk pressure, pressure difference and surface tension
    if (iT%statIter==0) {
      AnalyticalConst3D<T,T> two_( 2. );
      AnalyticalConst3D<T,T> onefive_( 1.5 );
      AnalyticalConst3D<T,T> k1_( kappa1 );
      AnalyticalConst3D<T,T> k2_( kappa2 );
      AnalyticalConst3D<T,T> cs2_( 1./descriptors::invCs2<T,DESCRIPTOR>() );
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> two(two_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> onefive(onefive_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> k1(k1_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> k2(k2_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> cs2(cs2_, sLattice1);

      // Calculation of bulk pressure:
      // c_1 = density of fluid 1; c_2 = density of fluid 2
      // p_bulk = rho*c_s^2 + kappa1 * (3/2*c_1^4 - 2*c_1^3 + 0.5*c_1^2)
      //                    + kappa2 * (3/2*c_2^4 - 2*c_2^3 + 0.5*c_2^2)
      SuperIdentity3D<T,T> bulkPressure ( density1*cs2
                                          + k1*( onefive*c1*c1*c1*c1 - two*c1*c1*c1 + half*c1*c1 )
                                          + k2*( onefive*c2*c2*c2*c2 - two*c2*c2*c2 + half*c2*c2 ) );

      AnalyticalFfromSuperF3D<T, T> interpolPressure( bulkPressure, true, 1);
      T position[3] = {nx/T(2), nx/T(2), nx/T(2)};
      T pressureIn = 0.;
      T pressureOut = 0.;
      interpolPressure(&pressureIn, position);
      position[0] = (T(N)/T(100))*converter.getPhysDeltaX();
      position[1] = (T(N)/T(100))*converter.getPhysDeltaX();
      position[2] = (T(N)/T(100))*converter.getPhysDeltaX();
      interpolPressure(&pressureOut, position);

      clout << "Pressure Difference: " << pressureIn-pressureOut << "  ;  ";
      clout << "Surface Tension: " << radius*(pressureIn-pressureOut)/2 << std::endl;
      clout << "Analytical Pressure Difference: " << alpha/(3.*radius) * (kappa1 + kappa2) << "  ;  ";
      clout << "Analytical Surface Tension: " << alpha/6. * (kappa1 + kappa2) << std::endl;
    }
  }
}


int main( int argc, char *argv[] )
{

  // === 1st Step: Initialization ===

  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    (T)   N, // resolution
    (T)   1., // lattice relaxation time (tau)
    (T)   nx, // charPhysLength: reference length of simulation geometry
    (T)   1.e-6, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1. // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  std::vector<T> extend = { nx, nx, nx };
  std::vector<T> origin = { 0, 0, 0 };
  IndicatorCuboid3D<T> cuboid(extend,origin);
#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, converter.getPhysDeltaX() );
#endif

  // set periodic boundaries to the domain
  cuboidDecomposition.setPeriodicity({ true, true, true });

  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  loadBalancer.print();

  // Instantiation of superGeometry
  SuperGeometry<T,3> superGeometry( cuboidDecomposition,loadBalancer,2 );

  prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice1( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLattice2( superGeometry );

  prepareLattice( sLattice1, sLattice2, converter, superGeometry );

  // Add the lattice couplings
  clout << "Add lattice coupling" << std::endl;
  // The chemical potential coupling must come before the force coupling
  SuperLatticeCoupling coupling1(
  ChemicalPotentialCoupling3D{},
  names::A{}, sLattice1,
  names::B{}, sLattice2);

  coupling1.template setParameter<ChemicalPotentialCoupling3D::ALPHA>(alpha);
  coupling1.template setParameter<ChemicalPotentialCoupling3D::KAPPA1>(kappa1);
  coupling1.template setParameter<ChemicalPotentialCoupling3D::KAPPA2>(kappa2);
  //coupling1.template setParameter<FreeEnergyChemicalPotentialCoupling3DM::KAPPA3>(kappa3);

  SuperLatticeCoupling coupling2(
  ForceCoupling3D{},
  names::A{}, sLattice2,
  names::B{}, sLattice1);

  {
    auto& communicator = sLattice1.getCommunicator(stage::PostCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestOverlap(sLattice1.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice2.getCommunicator(stage::PreCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestField<RhoStatistics>();
    communicator.requestOverlap(sLattice2.getOverlap());
    communicator.exchangeRequests();
  }

  sLattice1.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  sLattice2.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  clout << "Add lattice coupling ... OK!" << std::endl;

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

for ( iT=0; iT<=maxIter; ++iT ) {
    // Computation and output of the results
    getResults( sLattice2, sLattice1, iT, superGeometry, timer, converter );

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();

    // Execute coupling between the two lattices
    sLattice1.executePostProcessors(stage::PreCoupling());

    sLattice1.getCommunicator(stage::PreCoupling()).communicate();
    coupling1.execute();
    sLattice1.getCommunicator(stage::PostCoupling()).communicate();

    sLattice2.executePostProcessors(stage::PreCoupling());

    sLattice2.getCommunicator(stage::PreCoupling()).communicate();
    coupling2.execute();
    sLattice2.getCommunicator(stage::PostCoupling()).communicate();

  }

  timer.stop();
  timer.printSummary();

}

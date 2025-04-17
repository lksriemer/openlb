/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2006-2019 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod, Adrian Kummerländer
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

/* cylinder2d.h:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */

#ifndef CYLINDER_2D_H
#define CYLINDER_2D_H

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;

#define BOUZIDI

// Parameters for the simulation setup
const T Re = 20.;       // Reynolds number
const T maxPhysT = 16;  // max. simulation time in s, SI unit
const T physInterval = 0.1;  // interval for the convergence check in s
const T residuum = 1e-3;      // residuum for the convergence check

struct GeometryParameters {
  T L;
  T lengthX;
  T lengthY;
  T centerCylinderX;
  T centerCylinderY;
  T radiusCylinder;

  GeometryParameters(int N) {
    L = 0.1 / N;
    lengthX = 2.2;
    lengthY = 0.41 + L;
    centerCylinderX = 0.2;
    centerCylinderY = 0.2 + L / 2.0;
    radiusCylinder = 0.05;
  }
};


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry,
                      std::shared_ptr<IndicatorF2D<T>> circle,
                      GeometryParameters geomParams)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend( geomParams.lengthX, geomParams.lengthY );
  Vector<T,2> origin;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1} );

  // Set material number for inflow
  extend[0] = 2.*geomParams.L;
  origin[0] = -geomParams.L;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = geomParams.lengthX - geomParams.L;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );
  // Set material number for cylinder
  superGeometry.rename( 1,5, circle );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry,
                     std::shared_ptr<IndicatorF2D<T>> circle)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  // Setting of the boundary conditions
  // if boundary conditions are chosen to be local
  //boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 3);
  //boundary::set<boundary::LocalPressure>(sLattice, superGeometry, 4);

  //if boundary conditions are chosen to be interpolatedy, 3);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  // if boundary conditions are chosen to be Zou He type
  //boundary::set<boundary::ZouHeVelocity>(sLattice, superGeometry, 3);
  //boundary::set<boundary::ZouHePressure>(sLattice, superGeometry, 4);

  // Material=5 -->bouzidi / bounce back
  #ifdef BOUZIDI
  setBouzidiBoundary(sLattice, superGeometry, 5, *circle);
  #else
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 5);
  #endif

  // Initial conditions
  AnalyticalConst2D<T,T> rhoF( 1 );
  std::vector<T> velocity( 2,T( 0 ) );
  AnalyticalConst2D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,2>& superGeometry,
                        T L )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
  int iTupdate = 5;

  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = L/2.;
    Poiseuille2D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall );

    sLattice.defineU( superGeometry, 3, poiseuilleU );

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
double getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
                 Gnuplot<T>& gplot,
                 bool eoc,
                 GeometryParameters geomParams)
{
  OstreamManager clout( std::cout,"getResults" );

  // Gnuplot constructor (must be static!)
  // for real-time plotting: gplot("name", true) // experimental!
  // static Gnuplot<T> gplot( "drag" );
  T dragCoefficient = 0.;

  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );

  const int vtkIter  = converter.getLatticeTime( .3 );
  const int statIter = converter.getLatticeTime( .1 );

  if (! eoc)
  {
    SuperVTMwriter2D<T> vtmWriter( "cylinder2d" );
    VTUpointWriter<T,T,2> vtuWriter( "Points", true );
    SuperLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR> quality( sLattice, converter);
    SuperRoundingF2D<T> roundedQuality ( quality, RoundingMode::NearestInteger);
    SuperDiscretizationF2D<T> discretization ( roundedQuality, 0., 2. );

    Vector<T,2> pos1={geomParams.centerCylinderX - geomParams.radiusCylinder, geomParams.centerCylinderY};
    Vector<T,2> pos2={geomParams.centerCylinderX + geomParams.radiusCylinder, geomParams.centerCylinderY};

    vtmWriter.addFunctor( quality );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );

    vtuWriter.addPoint( pos1 );
    vtuWriter.addPoint( pos2 );
    vtuWriter.addFunctor(  pressure );
    vtuWriter.addFunctor(  velocity );

    T point[2] = {};
    point[0] = geomParams.centerCylinderX + 3*geomParams.radiusCylinder;
    point[1] = geomParams.centerCylinderY;
    AnalyticalFfromSuperF2D<T> intpolateP( pressure, true );
    T p;
    intpolateP( &p,point );

    if ( iT == 0 ) {
      // Writes the geometry, cuboid no. and rank no. as vti file for visualization
      SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
      SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
      vtmWriter.write( cuboid );
      vtmWriter.write( rank );

      vtmWriter.createMasterFile();
      vtuWriter.createMasterFile();
    }

    // Writes the vtk files
    if ( iT%vtkIter == 0 && iT > 0 ) {
      vtmWriter.write( iT );
      vtuWriter.write( iT );

      {
        SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
        BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
        // write output as JPEG
        heatmap::write(planeReduction, iT);
      }
      {
        BlockReduction2D2D<T> planeReduction( discretization, 600, BlockDataSyncMode::ReduceOnly );
        heatmap::plotParam<T> jpeg_scale;
        jpeg_scale.name = "quality";
        jpeg_scale.colour = "blackbody";
        heatmap::write( planeReduction, iT, jpeg_scale );
      }
    }
  }

  if ( iT%statIter == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T> intpolatePressure( pressure, true );
    SuperLatticePhysDrag2D<T,DESCRIPTOR> drag( sLattice, superGeometry, 5, converter );

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = geomParams.centerCylinderX - geomParams.radiusCylinder;
    point1[1] = geomParams.centerCylinderY;

    point2[0] = geomParams.centerCylinderX + geomParams.radiusCylinder;
    point2[1] = geomParams.centerCylinderY;

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    int input[3] = {};
    T _drag[drag.getTargetDim()];
    drag( _drag,input );

    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << std::endl;
    dragCoefficient = _drag[0];

    if (!eoc) {
      // set data for gnuplot: input={xValue, yValue(s), names (optional), position of key (optional)}
      gplot.setData( converter.getPhysTime( iT ), {_drag[0], 5.58}, {"drag(openLB)", "drag(schaeferTurek)"}, "bottom right", {'l','l'} );

      // every (iT%vtkIter) write an png of the plot
      if ( iT%( vtkIter ) == 0 ) {
        // writes pngs: input={name of the files (optional), x range for the plot (optional)}
        gplot.writePNG( iT, maxPhysT );
      }
    }
  }

  // write pdf at last time step
  if ( iT == converter.getLatticeTime( maxPhysT ) - 1 ) {
    // writes pdf
    gplot.writePDF();
  }

  return dragCoefficient;
}

double simulateCylinder( int N, Gnuplot<T>& gplot, bool eoc )
  {
    // === 1st Step: Initialization ===
    // initialize( &argc, &argv );
    singleton::directories().setOutputDir( "./tmp/" );
    OstreamManager clout( std::cout,"main" );

    GeometryParameters geomParams(N);

    UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      int {N},                        // resolution: number of voxels per charPhysL
      (T)   0.56,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)   2.0*geomParams.radiusCylinder,       // charPhysLength: reference length of simulation geometry
      (T)   0.2,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)   0.2*2.*geomParams.radiusCylinder/Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)   1.0                       // physDensity: physical density in __kg / m^3__
    );
    // Prints the converter log as console output
    converter.print();
    // Writes the converter log in a file
    converter.write("cylinder2d");

    // === 2rd Step: Prepare Geometry ===
    Vector<T,2> extend( geomParams.lengthX, geomParams.lengthY );
    Vector<T,2> origin;
    IndicatorCuboid2D<T> cuboid( extend, origin );

    // Instantiation of a cuboidDecomposition with weights
  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
  #else
    const int noOfCuboids = 7;
  #endif
    CuboidDecomposition2D<T> cuboidDecomposition( cuboid, geomParams.L, noOfCuboids );

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

    // Instantiation of a superGeometry
    SuperGeometry<T,2> superGeometry( cuboidDecomposition, loadBalancer );

    Vector<T,2> center( geomParams.centerCylinderX, geomParams.centerCylinderY );
    std::shared_ptr<IndicatorF2D<T>> circle = std::make_shared<IndicatorCircle2D<T>>( center, geomParams.radiusCylinder );

    prepareGeometry( converter, superGeometry, circle, geomParams );

    // === 3rd Step: Prepare Lattice ===
    SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

    //prepareLattice and set boundaryConditions
    prepareLattice( sLattice, converter, superGeometry, circle );

    // === 4th Step: Main Loop with Timer ===
    clout << "starting simulation..." << std::endl;
    util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
    // util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
    timer.start();

    T drag = 0;

    for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
      // === 5th Step: Definition of Initial and Boundary Conditions ===
      setBoundaryValues( sLattice, converter, iT, superGeometry, geomParams.L );

      // === 6th Step: Collide and Stream Execution ===
      sLattice.collideAndStream();

      // === 7th Step: Computation and Output of the Results ===
      if (iT % converter.getLatticeTime( .1 ) == 0)
      {
        drag = getResults( sLattice, converter, iT, superGeometry, timer, gplot, eoc, geomParams );
      }

    }

    timer.stop();
    timer.printSummary();

    return drag;
  }


#endif

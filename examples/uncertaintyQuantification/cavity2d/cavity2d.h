#include <olb.h>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;
using BulkDynamics = ConstRhoBGKdynamics<T,DESCRIPTOR>;

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1} );
  superGeometry.clean();

  T eps = converter.getConversionFactorLength();
  Vector<T,2> extend( T( 1 ) + 2*eps, 2*eps );
  Vector<T,2> origin( T() - eps, T( 1 ) - eps );
  IndicatorCuboid2D<T> lid( extend, origin );
  // Set material number for lid
  superGeometry.rename( 2,3,1,lid );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);

  // Material=2,3 -->bulk dynamics, velocity boundary
  boundary::set<boundary::InterpolatedVelocity<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedVelocity<T,DESCRIPTOR,BulkDynamics>>(sLattice, superGeometry, 3);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}


void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T, DESCRIPTOR>& sLattice,
                        int iT, SuperGeometry<T,2>& superGeometry )
{

  if ( iT==0 ) {
    // set initial values: v = [0,0]
    AnalyticalConst2D<T,T> rhoF( 1 );
    std::vector<T> velocity( 2,T() );
    AnalyticalConst2D<T,T> uF( velocity );

    auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3});
    sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );
    sLattice.defineRhoU( bulkIndicator, rhoF, uF );

    // set non-zero velocity for upper boundary cells
    velocity[0] = converter.getCharLatticeVelocity();
    AnalyticalConst2D<T,T> u( velocity );
    sLattice.defineU( superGeometry, 3, u );

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT, util::Timer<T> timer,
                 const T logT, const T maxPhysT, const T imSave, const T vtkSave, const T gnuplotSave,
                 std::string filenameGif, std::string filenameVtk, std::string filenameGnuplot,
                 const int timerPrintMode,
                 SuperGeometry<T,2>& superGeometry, bool converged )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( filenameVtk );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    SuperLatticeDiscreteNormal2D<T, DESCRIPTOR> discreteNormal( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );
    SuperLatticeDiscreteNormalType2D<T, DESCRIPTOR> discreteNormalType( sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );

    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.write( discreteNormal );
    vtmWriter.write( discreteNormalType );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%converter.getLatticeTime( logT )==0 || converged ) {
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
    timer.print( iT,timerPrintMode );
  }

  // Writes the VTK files
  if ( ( iT%converter.getLatticeTime( vtkSave )==0 && iT>0 ) || converged ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );

    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );
  }

  // Writes the Gif files
  if ( ( iT%converter.getLatticeTime( imSave )==0 && iT>0 ) || converged ) {
    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperEuklidNorm2D<T,DESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output of velocity as JPEG
    heatmap::write(planeReduction, iT);
  }

  // Output for x-velocity along y-position at the last time step
  //if ( iT == converter.getLatticeTime( maxPhysT ) || converged ) {
  if ( (iT%converter.getLatticeTime( gnuplotSave ) == 0 && iT>0 ) || converged ) {
    // Gives access to velocity information on lattice
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocityField( sLattice, converter );
    // Interpolation functor with velocityField information
    AnalyticalFfromSuperF2D<T> interpolation( velocityField, true, 1 );

    Vector<T,17> y_coord( {128, 125, 124, 123, 122, 109, 94, 79, 64, 58, 36, 22, 13, 9, 8, 7, 0} );
    // Ghia, Ghia and Shin, 1982: "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method";  Table 1
    Vector<T,17> vel_ghia_RE1000( { 1.0,     0.65928, 0.57492, 0.51117, 0.46604,
                                    0.33304, 0.18719, 0.05702,-0.06080,-0.10648,
                                    -0.27805,-0.38289,-0.29730,-0.22220,-0.20196,
                                    -0.18109, 0.0
                                  } );
    Vector<T,17> vel_ghia_RE100( {1.0,     0.84123, 0.78871, 0.73722, 0.68717,
                                  0.23151, 0.00332,-0.13641,-0.20581,-0.21090,
                                  -0.15662,-0.10150,-0.06434,-0.04775,-0.04192,
                                  -0.03717, 0.0
                                 } );
    Vector<T,17> vel_simulation;

    // Gnuplot interface to create plots
    Gnuplot<T> gplot( filenameGnuplot + "_iT" + std::to_string(iT) );
    // Define comparison values
    Vector<T,17> comparison = vel_ghia_RE1000;

    for ( int nY = 0; nY < 17; ++nY ) {
      // 17 data points evenly distributed between 0 and 1 (height)
      T position[2] = {0.5, y_coord[nY]/T(128)};
      T velocity[2] = {T(), T()};
      // Interpolate velocityField at "position" and save it in "velocity"
      interpolation( velocity, position );
      // Save value of velocity (in x-direction) in "vel_simulation" for every position "nY"
      vel_simulation[nY] = velocity[0];
      // Set data for plot output
      gplot.setData( position[1], {vel_simulation[nY],comparison[nY]}, {"simulated","Ghia"} );
    }
    // Create PNG file
    gplot.writePNG();
    // Console output with results
    clout << "absoluteErrorL2(line)=" << norm(vel_simulation - comparison) / 17. << "; relativeErrorL2(line)=" << norm(vel_simulation - comparison) / norm(comparison) << std::endl;
  }
}



void simulateCavity2d(double physViscosity, int N, int sample)
{

  // === 1st Step: Initialization ===
  OstreamManager clout( std::cout,"sample" );

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    N,                 // resolution
    0.5384, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    1.0,             // charPhysLength: reference length of simulation geometry
    1.0,               // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    physViscosity,                 // physViscosity: physical kinematic viscosity in __m^2 / s__
    1.0                 // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("cavity2d");

  // copy from cavity2d.xml file
  T logT = 1;
  T imSave = 5;
  T vtkSave = 1;
  T gnuplotSave = 5;
  T maxPhysT = 100;
  int timerPrintMode = 0;

  std::string filenameGif = "cavity2dimage";
  std::string filenameVtk = "cavity2dvtk";
  std::string filenameGnuplot = "centerVelocityX";

  // === 2nd Step: Prepare Geometry ===
  Vector<T,2> extend( 1,1 );
  Vector<T,2> origin( 0,0 );
  IndicatorCuboid2D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, 1.0 / N, singleton::mpi().getSize() );
#else
  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, 1.0 / N, 1 );
#endif

  cuboidDecomposition.print();

  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,2> superGeometry( cuboidDecomposition, loadBalancer );
  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  prepareLattice( converter, sLattice, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  int interval = converter.getLatticeTime( 1 /*config["Application"]["ConvergenceCheck"]["interval"].get<T>()*/ );
  T epsilon = 1e-3;
  util::ValueTracer<T> converge( interval, epsilon );

  timer.start();
  for ( std::size_t iT=0; iT <= converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice, converter, iT, timer, logT, maxPhysT, imSave, vtkSave, gnuplotSave, filenameGif, filenameVtk,
                  filenameGnuplot, timerPrintMode, superGeometry, converge.hasConverged() );
      break;
    }
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, timer, logT, maxPhysT, imSave, vtkSave, gnuplotSave, filenameGif, filenameVtk,
                filenameGnuplot, timerPrintMode, superGeometry, converge.hasConverged() );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }
  timer.stop();
  timer.printSummary();
}

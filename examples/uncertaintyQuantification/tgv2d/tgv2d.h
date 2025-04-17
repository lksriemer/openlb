#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = D2Q9<>;

#define DNS
// #define SMAGORINSKY

#ifdef DNS
  using BulkDynamics = ConstRhoBGKdynamics<T,DESCRIPTOR>;
#elif defined(SMAGORINSKY)
  using BulkDynamics = SmagorinskyBGKdynamics<T,DESCRIPTOR>;
#endif

// Parameters for the simulation setup
//const int N = 33;       // resolution of the model
const T num_vortice = 2;
const T maxPhysT = 100;  // util::max. simulation time in s, SI unit
const T startPhysT = 0;  // util::max. simulation time in s, SI unit
const T smagoConst = 0.1;

const T vtkSave = 10; //saves vtk files every vtkSave seconds

template <typename T, typename _DESCRIPTOR>
class Tgv2D : public AnalyticalF2D<T,T> {

protected:
  T u0;

// initial solution of the TGV
public:
  Tgv2D(UnitConverter<T,_DESCRIPTOR> const& converter) : AnalyticalF2D<T,T>(2)
  {
    u0 = converter.getCharLatticeVelocity();
  };

  bool operator()(T output[], const T input[]) override
  {
    const T x = input[0];
    const T y = input[1];

    output[0] = -u0 * util::cos(x) * util::sin(y);
    output[1] =  u0 * util::sin(x) * util::cos(y);

    return true;
  };
};


  void save_velocity_field(std::string dir, int nx, int ny, std::vector<std::vector<T>> u, std::vector<std::vector<T>> v, int idx) {
    std::string filename_u = dir + "/u_" + std::to_string(idx) + ".dat";
    // std::cout << "filename_u: " << filename_u << std::endl;
    std::ofstream outputFile_u(filename_u);
    if (!outputFile_u) {
      std::cerr << "Error opening the file: " << filename_u << std::endl;
      return;
    }
    outputFile_u << std::setprecision(20);
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        outputFile_u << u[i][j] << "\t";
      }
      outputFile_u << "\n";
    }
    outputFile_u.close();

    std::string filename_v = dir + "/v_" + std::to_string(idx) + ".dat";
    std::ofstream outputFile_v(filename_v);
    if (!outputFile_v) {
      std::cerr << "Error opening the file: " << filename_v << std::endl;
      return;
    }
    outputFile_v << std::setprecision(20);
    for (int i = 0; i < nx; ++i) {
      for (int j = 0; j < ny; ++j) {
        outputFile_v << v[i][j] << "\t";
      }
      outputFile_v << "\n";
    }

    outputFile_v.close();

  }

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,1 );
  superGeometry.communicate();

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

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics<T,DESCRIPTOR>>(superGeometry, 0);

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  #ifdef SMAGORINSKY
    sLattice.setParameter<collision::LES::Smagorinsky>(smagoConst);// Smagorisky Constant, for ConsistentStrainSmagorinsky smagoConst = 0.033
  #endif

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout(std::cout,"setBoundaryValues");
  AnalyticalConst2D<T,T> rhoF( 1 );
  Tgv2D<T,DESCRIPTOR> uSol(converter);

  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uSol );
  sLattice.defineRhoU( bulkIndicator, rhoF, uSol );

  // Make the lattice ready for simulation
  sLattice.initialize();
}


void saveResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                  UnitConverter<T,DESCRIPTOR> const& converter,
                  SuperGeometry<T,2>& superGeometry,
                  std::string dir,
                  int resolution,
                  int idx,
                  int iT, T td ) {

  OstreamManager clout( std::cout,"output" );
  const int UQIter = td * 0.1;

  if (iT % UQIter == 0) {
    std::vector<std::vector<double>> u(resolution, std::vector<double>(resolution));
    std::vector<std::vector<double>> v(resolution, std::vector<double>(resolution));

    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocityField( sLattice, converter );
    AnalyticalFfromSuperF2D<T> intpolateVelocity( velocityField, true, 1 );
    // SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );

    T dx = converter.getPhysDeltaX();
    // clout << "resolution: " << resolution << std::endl;

    for (int i = 0; i < resolution; ++i) {
      for (int j = 0; j < resolution; ++j) {
        T position[2] = {i * dx, j * dx};
        T velocity[2];
        intpolateVelocity(velocity, position);
        u[i][j] = velocity[0];
        v[i][j] = velocity[1];
      }
    }
    std::string outputDir = dir + "dataFiles/" + std::to_string(iT / UQIter);

    save_velocity_field( outputDir, resolution, resolution, u, v, idx);
    clout << "saved " << iT / UQIter << "\ttime step: " << iT << std::endl;
  }

}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer, T td )
{
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter("tgv2d", 1, false );

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);

    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }
  // if (iT % converter.getLatticeTime(vtkSave) == 0) {
  int step = 0.1 * td;
  if (iT % step == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write(iT);

    timer.update(iT);
    timer.printStep(2);
    sLattice.getStatistics().print(iT,converter.getPhysTime( iT ));

    // int rank = 0;
    // #ifdef PARALLEL_MODE_MPI
    //   rank = singleton::mpi().getRank();
    // #endif
    // Only rank 0 writes to the file
    // if (rank == 0) {
    //   iTList.push_back(iT);
    // }
  }
}

// void simulateTGV(double Re, std::vector<double> random, double L, double Ma, int N, std::string foldPath, int idx, std::vector<std::vector<double>>&u, std::vector<std::vector<double>>&v)
void simulateTGV(double Re, double charL, double tau, int N, std::string foldPath, int idx)
{

  // === 1st Step: Initialization ===
  //initialize( &argc, &argv );
  OstreamManager clout( std::cout,"main" );
  clout << "start sample " << idx << std::endl;

  // UnitConverter<T,DESCRIPTOR> converter(
  //   (T)   charL / N,
  //   (T)   charL / N * Ma / std::sqrt(3.0),
  //   (T)   charL,                            // charPhysLength:
  //   (T)   0.01,  // charPhysVelocity: in __m / s__
  //   (T)   0.01 * charL / Re,                    // physViscosity: in __m^2 / s__
  //   (T)   1.0                           // physDensity:  in __kg / m^3__
  // );
  double mu = 0.01 * charL / Re;
    UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      int {N},                        // resolution: number of voxels per charPhysL
      (T)   tau,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)   charL,       // charPhysLength: reference length of simulation geometry
      (T)   0.01,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)   mu, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)   1.0                       // physDensity: physical density in __kg / m^3__
    );

  //clout << u1 << "\t" << miu << std::endl;

  // Prints the converter log as console output
  // converter.print();
  // Writes the converter log in a file
  converter.write("tgv2d");

  // === 2nd Step: Prepare Geometry ===
  #ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  Vector<T,2> extend( charL,charL );
  Vector<T,2> origin( 0,0 );
  IndicatorCuboid2D<T> cuboid( extend, origin );

  CuboidDecomposition2D<T> cuboidDecomposition( cuboid, converter.getPhysDeltaX(), noOfCuboids );
  cuboidDecomposition.setPeriodicity({true, true});

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
  setBoundaryValues( converter, sLattice, superGeometry );
  timer.start();
  std::size_t iT=0;

  T td = 1.0 / (( 0.01 * charL / Re ) * (2 * converter.getPhysDeltaX() * converter.getPhysDeltaX()));
  clout << "td: " << td << std::endl;
  //T td = 6388;
  // while(iT <= converter.getLatticeTime( maxPhysT ))
  for ( iT = 0; iT <= 0.5 * td; ++iT )
  {
    if ( idx == 0 ) {
      getResults( sLattice, converter, iT, superGeometry, timer, td);
      // clout << "saved " << iT << std::endl;
    }
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    //getResults( sLattice, *converter, iT, timer, logT, maxPhysT, imSave, vtkSave, filenameGif, filenameVtk,
    //            timerPrintMode, timerTimeSteps, superGeometry, converge.hasConverged() );
    //converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
    // saveResults( sLattice, converter, superGeometry, N, idx, iT );
    saveResults( sLattice, converter, superGeometry, foldPath, N, idx, iT, td );

    // iT++;
  }
  timer.stop();
  timer.printSummary();
  // delete &converter;
  // delete &timer;
}

#include "cylinder2d.h"
#include <time.h>
#include "../uq/uq.h"
#include "../uq/postprocessing2D.h"

using namespace std;


// #define USE_GSL
int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  OstreamManager clout( std::cout,"main" );




  string foldPath = "./uq/";
  string command;
  command = "mkdir -p " + foldPath;
  int result = std::system(command.c_str());
  if (result != 0) {
      clout << "Command failed with return code: " << result << std::endl;
  }
  else {
      clout << "finish mkdir" << std::endl;
  }


  /// Simulation Parameter
  int N = 10;
  double u0 = 0.2;
  // UQ parameters
  int order = 1;
  int nq = 10;

  // Parse command-line arguments
  if (argc > 1) {
    order = std::stoi(argv[1]);
  }
  if (argc > 2) {
    nq = std::stoi(argv[2]);
  }

        // === 2rd Step: Prepare Geometry ===
    GeometryParameters geomParams(N);
    Vector<T,2> extend( geomParams.lengthX, geomParams.lengthY );
    Vector<T,2> origin;
    IndicatorCuboid2D<T> cuboid( extend, origin );

    // Instantiation of a cuboidGeometry with weights
  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
  #else
    const int noOfCuboids = 4;
  #endif
    CuboidDecomposition2D<T> cuboidDecomposition( cuboid, 0.1/N, noOfCuboids );

    // Instantiation of a loadBalancer
    HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );

    // Instantiation of a superGeometry
    SuperGeometry<T,2> superGeometry( cuboidDecomposition, loadBalancer );


    // UncertaintyQuantification uq(UQMethod::MonteCarlo);
    // unsigned int seed = 123456;
    // uq.initializeMonteCarlo(nq, 1, Distribution(DistributionType::Uniform, 0.8 * u0, 1.2 * u0), seed);
    // uq.initializeMonteCarlo(1, 1, {Distribution(DistributionType::Normal, u0, 0.1 * u0)}, seed);

    UncertaintyQuantification uq(UQMethod::GPC);
    uq.initializeGPC(order, nq, {Distribution(DistributionType::Uniform, 0.8 * u0, 1.2 * u0)});
    // uq.initializeGPC(order, nq, {Distribution(DistributionType::Normal, u0, 0.1 * u0)});
    clout << "Using order: " << order << ", nq: " << nq << std::endl;

    std::vector<std::vector<double>> samples;
    uq.getSamplingPoints(samples);

  bool eoc = false;

  clout << "Start simulation" << std::endl;

  std::vector<double> drag(samples.size(), 0.0);

  for(auto n = 0; n < samples.size(); ++n)
  {
    string subFoldPath = foldPath + std::to_string(n) + "/tmp/";
    string command;
    command = "mkdir -p " + subFoldPath;
    int result = std::system(command.c_str());
    if (result != 0) {
        clout << "Command failed with return code: " << result << std::endl;
    }
    else {
        clout << "finish mkdir" << std::endl;
    }

    singleton::directories().setOutputDir( subFoldPath );

    clout << "physvelocity = " << samples[n][0] << std::endl;

    static Gnuplot<T> gplot( "drag" );

    drag[n] = simulateCylinder(N, samples[n][0], gplot, eoc);

  }

  clout << std::setprecision(20) << std::fixed;
  clout << "Drag coefficients: ";
  clout << "mean: " << uq.mean(drag) << ", ";
  clout << "std: " << uq.std(drag) << std::endl;

  // computeMeanAndStdAndWriteVTI2D<T, DESCRIPTOR>(uq, foldPath, "cylinder2d", "physVelocity", cuboidGeometry, superGeometry);

  return 0;
}

#include "cavity2d.h"
#include <time.h>
#include "../uq/uq.h"
#include "../uq/postprocessing2D.h"

using namespace std;

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

  int N = 64;
  double physViscosity = 0.01;

  // UncertaintyQuantification uq(UQMethod::MonteCarlo);
  UncertaintyQuantification uq(UQMethod::GPC);
  // uq.initializeMonteCarlo(10, 1, Distribution(DistributionType::Uniform, 0.9, 1.1));
  uq.initializeGPC(2, 5, {Distribution(DistributionType::Uniform, 0.8 * physViscosity, 1.2 * physViscosity)});
  std::vector<std::vector<double>> samples;
  uq.getSamplingPoints(samples);

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

  clout << "Start simulation" << std::endl;

  for(int n = 0; n < samples.size(); ++n)
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
    clout << "physViscosity = " << samples[n][0] << std::endl;
    simulateCavity2d(samples[n][0], N, n);

  }

  // computeMeanAndStdAndWriteVTI2D<T, DESCRIPTOR>(uq, foldPath, "cavity2d", "physVelocity", cuboidGeometry, superGeometry);

  // end = clock();

}

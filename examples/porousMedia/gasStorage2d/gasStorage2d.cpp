/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Alice Raeli, Luiz Czelusniak, Tim Bingert
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

/* gasStorage2d.cpp
 * This example shows a liquid (water) entering inside a porous media
 * filled with gas (compressed hydrogen at 100bar) and replacing the
 * gas, only to a certain extent as some of the gas remains trapped in
 * the porous domain which gives the effective porosity for such removal
 * of stored gas processes.
 *
 * To run this example, download the appropriate vti file with the
 * the same name as the example from the OpenLB website:
 * https://openlb.net/data/gas_storage/gasStorage2d.vti
 * and execute something like 'porousCode gasStorage2d.vti "Tiff Scalars"'.
 */

#include <olb.h>

using namespace std;
using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using NSDESCRIPTOR = D2Q9<RHO, NABLARHO, FORCE, EXTERNAL_FORCE, TAU_EFF, STATISTIC, SCALAR>;
using ACDESCRIPTOR = D2Q9<CONV_POPS, FORCE, SOURCE, SOURCE_OLD, VELOCITY, OLD_PHIU, STATISTIC,
                          PSI, NORMGRADPSI, SCALAR, PSI0, THETA, BOUNDARY>;
using NSBulkDynamics = MultiPhaseIncTRTdynamics<T,NSDESCRIPTOR>;
using ACBulkDynamics = AllenCahnBGKdynamics<T, ACDESCRIPTOR>;
using Coupling = LiangPostProcessor;

// Parameters in physical units
Vector<T, 2> length;                           // domain length [m]
T inletLength = 0.1;                           // length of inlet region; % of total length
T outletLength = 0.1;                          // length of outlet region; % of total length

T pressureDrop = 50000.;                       // pressure drop [N.m-2]

const T densityGas = 7.1;                      // gas density [kg.m-3]
const T densityLiquid = 992;                   // liquid density [kg.m-3]

const T viscosityGas = 1.34e-6;                // gas kinematic viscosity [m2.s-1]
const T viscosityLiquid = 5.5e-7;              // liquid kinematic viscosity [m2.s-1]

const T surfaceTension = 72e-3;                // surface tension [N.m-1]

const T contactAngle = M_PI * 40. / 180.;      // Contact angle [radians]

const T maxPhysTime  = 0.012;                  // max simulation time [s]
const T vtkIter      = maxPhysTime/400;        // write simulation output [s]
const T statIter     = maxPhysTime/400;        // simulation statistics time [s]

// Parameters in lattice units
const T interfaceThickness = 6.;               // interface thickness [l.u.]
const T tau_mobil = 0.8;                       // relaxation time order parameter [l.u.]

void prepareGeometry(MultiPhaseUnitConverterFromRelaxationTime<T, NSDESCRIPTOR> const& converter,
                     IndicatorBlockData2Dvti<T>& indicator,
                     SuperGeometry<T, 2>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1, 1});

  // internal domain, solid rocks - not fluid
  superGeometry.rename(1, 5, indicator);

  T dx = converter.getPhysDeltaX();

  // inlet and outlet
  IndicatorCuboid2D<T> inlet( dx, length[1] - dx, { -inletLength, length[1]/2. }, 0 );
  superGeometry.rename(2, 3, inlet);

  IndicatorCuboid2D<T> outlet( dx, length[1] - dx, { outletLength + length[0] - 0.5*dx, length[1]/2. }, 0);
  superGeometry.rename(2, 4, outlet);

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.outerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template <typename STAGE>
void signedDistanceFunction( SuperLattice<T,ACDESCRIPTOR>& sLatticeAC,
                             SuperGeometry<T,2>& superGeometry, T w )
{
  T max = 1.;
  while ( max <= 1.5*w ) {
    int in[2];
    sLatticeAC.getCommunicator(STAGE{}).communicate();
    sLatticeAC.executePostProcessors(STAGE{});
    SuperLatticeExternalScalarField2D<T, ACDESCRIPTOR, PSI> psi( sLatticeAC );
    SuperMax2D<T,T> Max_psi_(psi, superGeometry, 1);
    Max_psi_(&max, in);
  }
}

T helperConvectiveU(SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                    SuperGeometry<T, 2>& superGeometry,
                    MultiPhaseUnitConverterFromRelaxationTime<T, NSDESCRIPTOR> const& converter)
{
  T dx = converter.getPhysDeltaX();
  IndicatorCuboid2D<T> beforeOutlet_( 1.1*dx, length[1], { length[0] + outletLength - 2.*dx, length[1]/2. }, 0 );
  SuperIndicatorFfromIndicatorF2D<T> beforeOutlet(beforeOutlet_, superGeometry);
  int in[2];
  T uMax[2];
  SuperLatticeVelocity2D<T, NSDESCRIPTOR> u(sLatticeNS);
  SuperMax2D<T,T> uMax_(u, beforeOutlet);
  uMax_(uMax, in);
  return uMax[0];
}

template <typename SuperLatticeCoupling>
void prepareLattice(SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                    SuperLattice<T, ACDESCRIPTOR>& sLatticeAC,
                    SuperLatticeCoupling& coupling,
                    MultiPhaseUnitConverterFromRelaxationTime<T, NSDESCRIPTOR> const& converter,
                    SuperGeometry<T,2>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // define lattice Dynamics
  sLatticeNS.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeNS.defineDynamics<NSBulkDynamics>(superGeometry, 1);

  sLatticeAC.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeAC.defineDynamics<ACBulkDynamics>(superGeometry, 1);

  // initial conditions
  AnalyticalConst2D<T, T> one(1.);
  AnalyticalConst2D<T, T> two(2.);
  AnalyticalConst2D<T, T> rhov(densityGas/converter.getConversionFactorDensity());
  AnalyticalConst2D<T, T> rhol(densityLiquid/converter.getConversionFactorDensity());
  T tau_g = converter.computeRelaxationTimefromPhysViscosity( viscosityGas );
  AnalyticalConst2D<T, T> tauv( tau_g );
  T tau_l = converter.computeRelaxationTimefromPhysViscosity( viscosityLiquid );
  AnalyticalConst2D<T, T> taul( tau_l );
  AnalyticalConst2D<T, T> angleInside( contactAngle );
  AnalyticalConst2D<T, T> angleOutside( M_PI/2. );
  AnalyticalConst2D<T, T> p0(0.);
  AnalyticalConst2D<T, T> u0(0, 0);

  auto bulk   = superGeometry.getMaterialIndicator(1);
  auto wall   = superGeometry.getMaterialIndicator({2, 5});
  auto inlet  = superGeometry.getMaterialIndicator(3);
  auto outlet = superGeometry.getMaterialIndicator(4);
  auto fluid  = superGeometry.getMaterialIndicator({1, 3});
  auto all    = superGeometry.getMaterialIndicator({0, 1, 2, 3, 4, 5});

  T dx = converter.getPhysDeltaX();

  // water enters on the left, gas in the porous rock
  SmoothIndicatorFactoredCuboid2D<T, T> phi0( {-inletLength, length[1]/2.},
                                              2.2*inletLength, 2.*length[1],
                                              interfaceThickness*dx/2., 0, {0, 0}, 0, 1. );
  AnalyticalIdentity2D<T, T> rho0(rhov + (rhol - rhov) * phi0);
  AnalyticalIdentity2D<T, T> tau0(tauv + (taul - tauv) * phi0);

  SmoothIndicatorFactoredCuboid2D<T,T> fringe( {-inletLength, length[1]/2.},
                                                2.*(inletLength+length[0]+0.5*outletLength), 0,
                                                0.25*outletLength, 0, {0,0}, 0, 1. );
  sLatticeNS.defineField<descriptors::SCALAR>(all, fringe);

  sLatticeNS.defineField<descriptors::RHO>(all, rho0);
  sLatticeNS.defineField<descriptors::TAU_EFF>(all, tau0);
  sLatticeAC.defineField<descriptors::OLD_PHIU>(all, u0);
  sLatticeAC.defineField<descriptors::BOUNDARY>(wall, two);
  sLatticeAC.defineField<descriptors::BOUNDARY>(superGeometry.getMaterialIndicator({1, 3, 4}), one);
  sLatticeAC.defineField<descriptors::THETA>(superGeometry.getMaterialIndicator(2), angleOutside);
  sLatticeAC.defineField<descriptors::THETA>(superGeometry.getMaterialIndicator(5), angleInside);

  boundary::set<boundary::BounceBackIncompressible>(sLatticeNS, wall);
  boundary::set<boundary::PhaseFieldCurvedWall>(sLatticeAC, wall);
  boundary::set<boundary::RegularizedTemperature>(sLatticeAC, inlet);
  boundary::set<boundary::IncompressibleZouHePressure>(sLatticeNS, inlet);
  boundary::set<boundary::IncompressibleZouHePressure>(sLatticeNS, outlet);
  boundary::set<boundary::PhaseFieldConvective>(sLatticeAC, outlet);

  sLatticeAC.defineRhoU(all, phi0, u0);
  sLatticeAC.iniEquilibrium(all, phi0, u0);
  sLatticeNS.defineRhoU(all, p0, u0);
  sLatticeNS.iniEquilibrium(all, p0, u0);

  sLatticeAC.addPostProcessor<stage::InitOutlet>(outlet, meta::id<SetOutletCells<1, 0>>());
  sLatticeAC.addPostProcessor<stage::PreCoupling>(fluid,meta::id<RhoStatistics>());

  sLatticeAC.addPostProcessor<stage::PreCoupling>(meta::id<initialPsi>{});
  sLatticeAC.addPostProcessor<stage::IterativePostProcess>(bulk,meta::id<normGradPsi>{});
  sLatticeAC.addPostProcessor<stage::IterativePostProcess>(meta::id<psiEvolve>{});
  sLatticeAC.setParameter<psiEvolve::DELTAT>(0.7);
  setSignedDistanceBoundary<T,ACDESCRIPTOR>(sLatticeAC, wall);
  setSignedDistanceBoundary<T,ACDESCRIPTOR>(sLatticeAC, inlet);
  setSignedDistanceBoundary<T,ACDESCRIPTOR>(sLatticeAC, outlet);

  coupling.template setParameter<Coupling::SIGMA>( surfaceTension / converter.getConversionFactorSurfaceTension() );
  coupling.template setParameter<Coupling::W>( interfaceThickness );
  coupling.template setParameter<Coupling::TAUS>({tau_g, tau_l});
  coupling.template setParameter<Coupling::RHOS>( {densityGas/converter.getConversionFactorDensity(),
                                                   densityLiquid/converter.getConversionFactorDensity()});
  coupling.template setParameter<Coupling::SWITCH>(1);

  sLatticeNS.setParameter<descriptors::OMEGA>( 1. / tau_l );
  T maxRhoGradient = (densityLiquid-densityGas)/converter.getConversionFactorDensity()/interfaceThickness;
  sLatticeNS.setParameter<collision::ITRT::TAU_MINUS>( T(1.5) );
  sLatticeNS.setParameter<collision::ITRT::MAXNABLARHO>( maxRhoGradient );
  sLatticeAC.setParameter<descriptors::OMEGA>( 1. / tau_mobil );
  sLatticeAC.addPostProcessor<stage::PhiLimiter>(fluid,meta::id<dispersionLimiter>{});
  sLatticeAC.setParameter<descriptors::EPSILON>( 1.5*interfaceThickness );
  sLatticeAC.setParameter<descriptors::INTERFACE_WIDTH>( interfaceThickness );

  T uMax = helperConvectiveU(sLatticeNS, superGeometry, converter);
  sLatticeAC.setParameter<descriptors::MAX_VELOCITY>(uMax);

  {
    auto& communicator = sLatticeAC.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }
  {
    auto& communicator =
        sLatticeAC.getCommunicator(stage::IterativePostProcess());
    communicator.requestOverlap(1);
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }

  sLatticeAC.executePostProcessors(stage::InitOutlet());
  sLatticeAC.executePostProcessors(stage::PreCoupling());
  signedDistanceFunction<stage::IterativePostProcess>(sLatticeAC,superGeometry,interfaceThickness);

  sLatticeNS.initialize();
  sLatticeAC.initialize();

  sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();
  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                       SuperLattice<T, ACDESCRIPTOR>& sLatticeAC,
                       MultiPhaseUnitConverterFromRelaxationTime<T, NSDESCRIPTOR> const& converter,
                       SuperGeometry<T, 2>& superGeometry, int iT)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  int iTmaxStart = 10000;
  int  iTupdate   = 1;
  auto inlet      = superGeometry.getMaterialIndicator(3);

  // Ramp pressure gradient of domain by increasing inlet pressure
  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    // Smooth start curve, polynomial
    PolynomialStartScale<T, T> StartScale(iTmaxStart, T(1));

    T iTvec[1] = {T(iT)};
    T frac[1]  = {};
    StartScale(frac, iTvec);

    const T imposedPressure = pressureDrop / converter.getConversionFactorPressure() * frac[0];
    AnalyticalConst2D<T, T> pressure(imposedPressure);
    sLatticeNS.defineRho(inlet, pressure);
  }
}

void getResults(SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                SuperLattice<T, ACDESCRIPTOR>& sLatticeAC, int iT,
                SuperGeometry<T, 2>& superGeometry, util::Timer<T>& timer,
                UnitConverter<T, NSDESCRIPTOR> converter)
{
  OstreamManager      clout(std::cout, "getResults");
  SuperVTMwriter2D<T> vtmWriter("gasStorage2d");

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(sLatticeNS);
    SuperLatticeRank2D<T, NSDESCRIPTOR>   rank(sLatticeNS);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }
  // Get statistics
  if (iT % converter.getLatticeTime(statIter) == 0) {
    // Timer console output
    timer.update(iT);
    timer.printStep();
    sLatticeNS.getStatistics().print(iT, converter.getPhysTime(iT));
    sLatticeAC.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  // Writes the VTK files
  if (iT % converter.getLatticeTime(vtkIter) == 0) {
    SuperLatticeDensity2D<T, NSDESCRIPTOR> p_total(sLatticeNS);
    p_total.getName() = "p_total";

    SuperLatticeExternalScalarField2D<T, NSDESCRIPTOR, RHO> rho_L(sLatticeNS);
    AnalyticalConst2D<T,T> ConversionDensity_(converter.getConversionFactorDensity());
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> ConversionDensity(ConversionDensity_, sLatticeNS);
    SuperIdentity2D<T,T> rho(rho_L * ConversionDensity);
    rho.getName() = "rho";

    SuperLatticeVelocity2D<T, NSDESCRIPTOR> velocity(sLatticeNS);
    velocity.getName() = "u";

    SuperLatticeExternalScalarField2D<T, ACDESCRIPTOR, BOUNDARY> boundary(sLatticeAC);
    boundary.getName() = "boundary";

    vtmWriter.addFunctor(p_total);
    vtmWriter.addFunctor(rho);
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(boundary); // can be used to cut the rock parts in i.e. Paraview
    vtmWriter.write(iT);
  }
}

int main(int argc, char* argv[])
{
  ///- Code init
  initialize(&argc, &argv);

  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  ///-Import vti
  std::string vtiFile;
  std::string arrayName;

  std::vector<std::string> cmdInput;
  if (argc > 1) {
    cmdInput.assign(argv + 1, argv + argc);
    vtiFile = cmdInput[0];
  }
  if (argc > 2) {
    arrayName = cmdInput[1];
  }
  else {
    clout << "Define <filename> <arrayname>(in that order)" << std::endl;
    return 1;
  }

  BlockVTIreader2D<T,T> vtiReader( vtiFile, arrayName ); // to ensure that the data persists

  T scalingFactor     = 0.0000004; // scale down
  auto cuboidSample = vtiReader.getCuboid();
  T    deltaRsample = scalingFactor;

  Vector<int, 2> extentSample     = cuboidSample.getExtent();
  Vector<T, 2>   originSamplePhys = cuboidSample.getOrigin() * scalingFactor;

  Vector<T, 2> extentSamplePhys = {deltaRsample * T(extentSample[0]),
                                   deltaRsample * T(extentSample[1])};
  for (unsigned i = 0; i < 2; ++i) {
    length[i] = extentSamplePhys[i];
  }

  // === 1st Step: Unit Converter ===
  MultiPhaseUnitConverterFromRelaxationTime<T,NSDESCRIPTOR> converter(
    (T)   1000,                      // resolution
    (T)   0.52,                      // lattice relaxation time
    (T)   densityLiquid/500.,        // lattice density
    (T)   length[1],                 // charPhysLength: reference length of simulation geometry
    (T)   viscosityLiquid,           // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   densityLiquid              // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a cuboidGeometry with weights
  BlockData<2,T,T>& blocco = vtiReader.getBlockData();
  BlockDataF2D<T,T> blockData(blocco);
  IndicatorBlockData2Dvti<T> ind(blocco, extentSamplePhys, originSamplePhys,
                                 deltaRsample, false);

  // Indicator containing porous medium, inlet and outlet zones
  inletLength *= length[0]; // computing inlet size
  outletLength *= length[0]; // computing outlet size
  IndicatorCuboid2D<T> domain( length[0] + inletLength + outletLength, length[1],
                               { ( length[0] + outletLength - inletLength )/2., length[1]/2. }, 0);

  // Creating cuboid geometry from domain indicator and voxel size
#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition2D<T> cuboidDecomposition( domain, converter.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidDecomposition2D<T> cuboidDecomposition( domain, converter.getPhysDeltaX() );
#endif
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  loadBalancer.print();
  // Instantiation of a superGeometry
  SuperGeometry<T, 2> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(converter, ind, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, NSDESCRIPTOR> sLatticeNS(superGeometry);
  SuperLattice<T, ACDESCRIPTOR> sLatticeAC(superGeometry);
  SuperLatticeCoupling coupling(LiangPostProcessor {},
                                names::NavierStokes {}, sLatticeNS,
                                names::Component1 {}, sLatticeAC);
  coupling.restrictTo(superGeometry.getMaterialIndicator({1, 4}));

  SuperLatticeCoupling velocityCoupling(VelocityCoupling{},
                                         names::NavierStokes {}, sLatticeNS,
                                         names::Component1 {}, sLatticeAC);
  velocityCoupling.restrictTo(superGeometry.getMaterialIndicator({3}));

  prepareLattice(sLatticeNS, sLatticeAC, coupling, converter, superGeometry);

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  int maxIter = converter.getLatticeTime(maxPhysTime)+1;
  util::Timer<T> timer(maxIter, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (iT = 0; iT <= maxIter; ++iT) {
    setBoundaryValues(sLatticeNS, sLatticeAC, converter, superGeometry, iT);

    // Collide and stream (and coupling) execution
    sLatticeNS.collideAndStream();
    sLatticeAC.collideAndStream();

    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeAC.executePostProcessors(stage::PreCoupling());
    if ( iT%100==0 ) {
      signedDistanceFunction<stage::IterativePostProcess>(sLatticeAC,superGeometry,interfaceThickness);
      sLatticeAC.executePostProcessors(stage::PhiLimiter());
    }
    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

    coupling.execute();

    velocityCoupling.execute();

    T uMax = helperConvectiveU(sLatticeNS, superGeometry, converter);
    sLatticeAC.setParameter<descriptors::MAX_VELOCITY>(uMax);

    // Computation and output of the results
    getResults(sLatticeNS, sLatticeAC, iT, superGeometry, timer, converter);

    if (std::isnan(sLatticeNS.getStatistics().getAverageEnergy())) {
      clout << "Code Diverged iteration: " << iT << endl;
      break;
    }
  }
  timer.stop();
  timer.printSummary();

  return 0;
}

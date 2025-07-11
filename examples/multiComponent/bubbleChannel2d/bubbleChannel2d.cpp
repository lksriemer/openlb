/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Tim Bingert
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

/* bubbleChannel2d.cpp
 * This example shows an advecting bubble of air in a channel of
 * water. The bubble has the size of 40 micrometers and has the
 * exact surface tension of an air-water mixture. This high surface
 * tension creates the necessity for a trick at the outlet, using
 * a fringe zone. This setup achieves a nice shape of the bubble
 * at the artificial outlet, which is roughly five bubble diameters
 * before the actual domain outlet.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;
using NSDESCRIPTOR =
    D2Q9<RHO, NABLARHO, FORCE, EXTERNAL_FORCE, TAU_EFF, STATISTIC, SCALAR>;
using ACDESCRIPTOR   = D2Q9<CONV_POPS, FORCE, SOURCE, VELOCITY, OLD_PHIU,
                            STATISTIC, THETA, PSI, NORMGRADPSI, SCALAR, PSI0>;
using NSBulkDynamics = MPIncTRTdynamics<T, NSDESCRIPTOR>;
using ACBulkDynamics = AllenCahnBGKdynamics<T, ACDESCRIPTOR>;
using Coupling       = LiangPostProcessor;

// Parameters for the simulation setup
const int Ny       = 160;    // domain resolution y [lattice units]
const int Nx       = 6 * Ny; // domain resolution x [lattice units]
const int diameter = Ny / 4; // lattice bubble diameter [lattice units]
const T   L_char   = 40e-6;  // physLength bubble diameter [physical units]
const T   C_rho    = 100.;   // conversion factor density [physical units]
const T   tau_l    = 0.53; // lattice relaxation time H2O liquid [lattice units]
const T   tau_v    = 13. * (tau_l - 0.5) +
                0.5; // lattice relaxation time air gas [lattice units]
const T tau_mobil =
    0.8; // lattice relaxation time for interface mobility [lattice units]
const T viscosityH2O   = 9e-7;  // physViscosity H2O liquid [physical units]
const T surfaceTension = 0.072; // physSurfaceTension H2O-air [physical units]
const std::vector<T> rhos = {
    1.2, 1000};        // physDensities air and water [physical units]
T         Re    = 100; // Reynolds number of channel flow []
const T   theta = M_PI * 40. / 180.; // contact angle (<90 is wetting) [radians]
const T   w     = 6.; // lattice interface thickness [lattice units]
int       maxIter  = 200000;
const int vtkIter  = 100;
const int statIter = 100;

void prepareGeometry(SuperGeometry<T, 2>& superGeometry, T dx)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  // bulk, MN=1
  superGeometry.rename(2, 1, {1, 1});

  // inlet, MN=3
  std::vector<T>       origin = {-0.5 * dx, dx};
  std::vector<T>       extend = {1.0 * dx, dx * (Ny - 2.)};
  IndicatorCuboid2D<T> inlet(extend, origin);
  superGeometry.rename(2, 3, inlet);

  // bottom wall, MN=4
  origin[0] = dx * (Nx - 0.5);
  IndicatorCuboid2D<T> outlet(extend, origin);
  superGeometry.rename(2, 4, outlet);

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

template <typename STAGE>
void signedDistanceFunction(SuperLattice<T, ACDESCRIPTOR>& sLatticeAC,
                            SuperGeometry<T, 2>& superGeometry, T w)
{
  T max = 1.;

  while (max <= 1.5 * w) {
    int in[2];
    sLatticeAC.getCommunicator(STAGE {}).communicate();
    sLatticeAC.executePostProcessors(STAGE {});
    SuperLatticeExternalScalarField2D<T, ACDESCRIPTOR, PSI> psi(sLatticeAC);
    SuperMax2D<T, T> Max_psi_(psi, superGeometry, 1);
    Max_psi_(&max, in);
  }
}

T helperConvectiveU(SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                    SuperGeometry<T, 2>& superGeometry, T dx)
{
  std::vector<T>                     origin = {dx * (Nx - 2.), dx};
  std::vector<T>                     extend = {1.5 * dx, dx * (Ny - 2.)};
  IndicatorCuboid2D<T>               beforeOutlet_(extend, origin);
  SuperIndicatorFfromIndicatorF2D<T> beforeOutlet(beforeOutlet_, superGeometry);
  int                                in[2];
  T                                  uMax[2];
  SuperLatticeVelocity2D<T, NSDESCRIPTOR> u(sLatticeNS);
  SuperMax2D<T, T>                        uMax_(u, beforeOutlet);
  uMax_(uMax, in);
  return uMax[0];
}

template <typename SuperLatticeCoupling>
void prepareLattice(
    SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
    SuperLattice<T, ACDESCRIPTOR>& sLatticeAC, SuperLatticeCoupling& coupling,
    MultiPhaseUnitConverterFromRelaxationTime<T, NSDESCRIPTOR> const& converter,
    SuperGeometry<T, 2>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // conversion properties
  T dx        = converter.getPhysDeltaX();
  T C_sigma   = converter.getConversionFactorSurfaceTension();
  T C_density = converter.getConversionFactorDensity();
  T C_p       = C_sigma / dx;
  T sigma     = surfaceTension / C_sigma;
  clout << "lattice surface tension: " << sigma << std::endl;

  // define lattice Dynamics
  sLatticeNS.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeNS.defineDynamics<NSBulkDynamics>(superGeometry, 1);
  sLatticeAC.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeAC.defineDynamics<ACBulkDynamics>(superGeometry, 1);

  // initial conditions
  AnalyticalConst2D<T, T> two(2.);
  AnalyticalConst2D<T, T> one(1.);
  AnalyticalConst2D<T, T> zero(0.);
  std::shared_ptr<AnalyticalF2D<T,T>> zero_ptr = std::make_shared<AnalyticalConst2D<T,T>>(0.);
  AnalyticalConst2D<T, T> rhov(rhos[0] / C_density);
  AnalyticalConst2D<T, T> rhol(rhos[1] / C_density);
  AnalyticalConst2D<T, T> tauv(tau_v);
  AnalyticalConst2D<T, T> taul(tau_l);
  AnalyticalConst2D<T, T> theta0(theta);
  AnalyticalConst2D<T, T> u0(0, 0);

  auto bulk   = superGeometry.getMaterialIndicator(1);
  auto wall   = superGeometry.getMaterialIndicator(2);
  auto inlet  = superGeometry.getMaterialIndicator(3);
  auto outlet = superGeometry.getMaterialIndicator(4);
  auto fluid  = superGeometry.getMaterialIndicator({1, 3});
  auto all    = superGeometry.getMaterialIndicator({0, 1, 2, 3, 4});

  // Start bubble part

  std::vector<T> pos0 = {dx * 2. * diameter, dx * Ny / 2. * 0.6};
  std::vector<T> pos1 = {dx * 2. * diameter, dx * Ny / 2. * 1.4};

  CircularInterface2D<T> phi0(pos0, dx * diameter / 2., dx * w, 1., true);
  CircularInterface2D<T> phi1(pos1, dx * diameter / 2., dx * w, 1., true);

  auto phi = AnalyticalIdentity2D<T, T>(one);
  // auto phi = AnalyticalIdentity2D<T, T>(one + phi0 - one + phi1 - one);

  AnalyticalIdentity2D<T, T>           rho0(rhov + (rhol - rhov) * phi);
  AnalyticalIdentity2D<T, T>           tau0(tauv + (taul - tauv) * phi);
  std::shared_ptr<AnalyticalF2D<T, T>> bubblePressure0(new LaplacePressure2D<T>(
      pos0, dx * diameter / 2., dx * w, sigma * C_sigma / C_p));
  std::shared_ptr<AnalyticalF2D<T, T>> bubblePressure1(new LaplacePressure2D<T>(
      pos1, dx * diameter / 2., dx * w, sigma * C_sigma / C_p));

  std::shared_ptr<AnalyticalF2D<T, T>> halfYLPressure(
      new AnalyticalConst2D<T, T>(2. * sigma / diameter / 2.));
  auto p0 = zero_ptr; // bubblePressure0 + bubblePressure1 + halfYLPressure + halfYLPressure;

  // End bubble part

  SmoothIndicatorFactoredCuboid2D<T, T> fringe(
      {0., dx * Ny / 2.}, dx * 2. * (Nx - 5.0 * diameter), 0,
      dx * 1.5 * diameter, 0, {0, 0}, 0, 1.);

  sLatticeNS.defineField<descriptors::RHO>(all, rho0);
  sLatticeNS.defineField<descriptors::TAU_EFF>(all, tau0);
  sLatticeNS.defineField<descriptors::SCALAR>(all, fringe);
  sLatticeNS.defineField<descriptors::SCALAR>(wall, two);
  sLatticeAC.defineField<descriptors::OLD_PHIU>(all, u0);
  sLatticeAC.defineField<descriptors::THETA>(wall, theta0);

  // walls
  std::vector<T>       origin = {-2. * dx, dx * 0.5};
  std::vector<T>       extend = {dx * (Nx + 4), dx * (Ny - 2)};
  IndicatorCuboid2D<T> wallLocation(extend, origin);
  setBouzidiBoundary(sLatticeNS, superGeometry, 2, wallLocation);
  setBouzidiPhaseField(sLatticeAC, superGeometry, 2, wallLocation);

  // inlet and outlet
  boundary::set<boundary::IncompressibleZouHeVelocity>(sLatticeNS, inlet);
  boundary::set<boundary::RegularizedTemperature>(sLatticeAC, inlet);
  boundary::set<boundary::IncompressibleZouHePressure>(sLatticeNS, outlet);
  boundary::set<boundary::PhaseFieldConvective>(sLatticeAC, outlet);

  const T        maxVelocity = Re / Ny * ((tau_l - 0.5) / 3.);
  const T        radius      = T(0.5) * dx * (Ny - 2);
  std::vector<T> axisPoint(2, T());
  axisPoint[0] = dx * Nx / 2.;
  axisPoint[1] = dx * (Ny - 1.) / 2.;
  std::vector<T> axisDirection(2, T());
  axisDirection[0] = 1;
  axisDirection[1] = 0;
  Poiseuille2D<T> poiseuille(axisPoint, axisDirection, maxVelocity, radius);

  sLatticeAC.defineRhoU(all, phi, poiseuille);
  sLatticeAC.iniEquilibrium(all, phi, poiseuille);
  sLatticeNS.defineRhoU(all, *p0, poiseuille);
  sLatticeNS.iniEquilibrium(all, *p0, poiseuille);

  sLatticeAC.addPostProcessor<stage::InitOutlet>(
      outlet, meta::id<SetOutletCells<1, 0>>());
  sLatticeAC.addPostProcessor<stage::PreCoupling>(fluid,
                                                  meta::id<RhoStatistics>());
  sLatticeAC.addPostProcessor<stage::PhiLimiter>(
      meta::id<dispersionLimiter> {});

  sLatticeAC.addPostProcessor<stage::PreCoupling>(meta::id<initialPsi> {});
  sLatticeAC.addPostProcessor<stage::IterativePostProcess>(
      bulk, meta::id<normGradPsi> {});
  sLatticeAC.addPostProcessor<stage::IterativePostProcess>(
      meta::id<psiEvolve> {});
  sLatticeAC.setParameter<psiEvolve::DELTAT>(0.7);
  setSignedDistanceBoundary<T, ACDESCRIPTOR>(sLatticeAC, wall);
  setSignedDistanceBoundary<T, ACDESCRIPTOR>(sLatticeAC, inlet);
  setSignedDistanceBoundary<T, ACDESCRIPTOR>(sLatticeAC, outlet);

  coupling.template setParameter<Coupling::SIGMA>(sigma);
  coupling.template setParameter<Coupling::W>(w);
  coupling.template setParameter<Coupling::TAUS>({tau_v, tau_l});
  coupling.template setParameter<Coupling::RHOS>(
      {rhos[0] / C_density, rhos[1] / C_density});
  coupling.template setParameter<Coupling::SWITCH>(1.);

  sLatticeNS.setParameter<descriptors::OMEGA>(1. / tau_l);
  sLatticeAC.setParameter<descriptors::OMEGA>(1. / tau_mobil);
  sLatticeAC.setParameter<descriptors::INTERFACE_WIDTH>(w);
  sLatticeAC.setParameter<descriptors::EPSILON>(1.5 * w);
  sLatticeNS.setParameter<collision::TRT::MAGIC>((tau_l - T(0.5)) *
                                                 (1.6 - 0.5));

  T uMax = helperConvectiveU(sLatticeNS, superGeometry, dx);
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
  signedDistanceFunction<stage::IterativePostProcess>(sLatticeAC, superGeometry,
                                                      w);

  sLatticeNS.initialize();
  sLatticeAC.initialize();

  sLatticeAC.iniEquilibrium(all, phi, poiseuille);
  sLatticeNS.iniEquilibrium(all, *p0, poiseuille);

  sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void insert_bubble(
    std::vector<T>& pos, double diameter,
    MultiPhaseUnitConverterFromRelaxationTime<T, NSDESCRIPTOR> const& converter,
    SuperLattice<T, NSDESCRIPTOR>&              sLatticeNS,
    SuperLattice<T, ACDESCRIPTOR>&              sLatticeAC,
    SuperGeometry<T, 2>& superGeometry)
{
  OstreamManager clout(std::cout, "main");

  clout << "Starting bubble insertion" << std::endl;
  // TODO:
  // The current approach is to set the entire fields again, 
  // changed at where the new bubble is supposed to be
  // Maybe instead only change values there?

  // TODO:
  // CUrrently all current values are passed into this fucntion
  // How to actually get them? From the sLattices?

  // TODO:
  // Check if there already is a bubble too close to the new bubble

  auto all    = superGeometry.getMaterialIndicator({1});

  // AnalyticalConst2D<T, T> tauv_anal(tau_v);
  // auto tauv = std::make_shared<SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR>>(tauv_anal);
  // AnalyticalConst2D<T, T> taul(tau_l);

  // TODO: Do we even need to do thigs with tau?
  // If yes, how to get it? This way?
  // SuperLatticeExternalScalarField2D<T, NSDESCRIPTOR, TAU_EFF> tau(sLatticeNS);

  std::shared_ptr<AnalyticalF2D<T,T>> one = std::make_shared<AnalyticalConst2D<T,T>>(1.);

  T dx        = converter.getPhysDeltaX();
  T C_sigma   = converter.getConversionFactorSurfaceTension();
  T C_density = converter.getConversionFactorDensity();
  T C_p       = C_sigma / dx;
  T sigma     = surfaceTension / C_sigma;

  SuperLatticeDensity2D<T,NSDESCRIPTOR>   p_raw(sLatticeNS);
  SuperLatticeVelocity2D<T,NSDESCRIPTOR>  u_ns_raw (sLatticeNS);
  SuperLatticeDensity2D<T,ACDESCRIPTOR>   phi_raw(sLatticeAC);
  SuperLatticeVelocity2D<T,ACDESCRIPTOR>  u_ac_raw (sLatticeAC);

  std::shared_ptr<AnalyticalF2D<T,T>> u_ns = std::make_shared<AnalyticalFfromSuperF2D<T,T>>(u_ns_raw);
  std::shared_ptr<AnalyticalF2D<T,T>> u_ac = std::make_shared<AnalyticalFfromSuperF2D<T,T>>(u_ac_raw);

  std::shared_ptr<AnalyticalF2D<T,T>> p = std::make_shared<AnalyticalFfromSuperF2D<T,T>>(p_raw);
  // auto phi_anal = AnalyticalFfromSuperF2D<T,T>(phi_raw);
  // auto phi_anal_ptr = std::make_shared<AnalyticalFfromSuperF2D<T,T>>(phi_raw);
  // auto phi = std::dynamic_pointer_cast<AnalyticalF2D<T,T>>(phi_anal_ptr);
  std::shared_ptr<AnalyticalF2D<T,T>> phi = std::make_shared<AnalyticalFfromSuperF2D<T,T>>(phi_raw);

  

  // auto rhov = std::make_shared<AnalyticalConst2D<T, T>>(rhos[0] / C_density);
  // auto tauv = std::make_shared<AnalyticalConst2D<T, T>>(tau_v);

  // Create the new bubbble phi, add it into the current phi

  std::shared_ptr<AnalyticalF2D<T,T>> phi_new_bubble = std::make_shared<CircularInterface2D<T>>(pos, dx * diameter / 2., dx * w, 1.,
                                        true);

  auto phi_new = phi + phi_new_bubble - one;

  // Update tau
  // AnalyticalIdentity2D<T, T> tau_new(tauv + (tau - tauv) * phi_new_bubble);
  // TODO: Set tau_new
  // How to do that? with defineField?
  // sLatticeNS.defineField<descriptors::TAU_EFF>(all, tau_new);


  // Now create the pressure / velocity for our new bubbble

  std::shared_ptr<AnalyticalF2D<T, T>> bubblePressure(new LaplacePressure2D<T>(
      pos, dx * diameter / 2., dx * w, sigma * C_sigma / C_p));
  std::shared_ptr<AnalyticalF2D<T, T>> halfYLPressure(
      new AnalyticalConst2D<T, T>(2. * sigma / diameter / 2.));

  // This is keeping the old pressure where the bubble is formed
  auto p_new = p + bubblePressure + halfYLPressure;

  // AnalyticalIdentity2D<T,T>  phi_new(phi_expr);
  // AnalyticalIdentity2D<T,T>  p_new  (p_expr);

  clout << "Finalizing bubble insertion" << std::endl;


  std::shared_ptr<AnalyticalF2D<T,T>> rhov = std::make_shared<AnalyticalConst2D<T, T>>(rhos[0] / C_density);
  std::shared_ptr<AnalyticalF2D<T,T>> rhol = std::make_shared<AnalyticalConst2D<T, T>>(rhos[1] / C_density);
  auto rho = rhov + (rhol - rhov) * phi_new; // rho is density

  // sLatticeNS.defineField<descriptors::RHO>(all, *rho);
  // sLatticeAC.defineField<descriptors::RHO>(all, *phi_new);

  auto pop_functor_ns = std::make_shared<IncompressibleEquilibriumPopulations2D<T, NSDESCRIPTOR>>(rho, p_new, u_ns);
  auto pop_functor_ac = std::make_shared<FirstOrderEquilibriumPopulations2D<T, ACDESCRIPTOR>>(phi_new, u_ac);

  sLatticeNS.defineField<descriptors::POPULATION>(all, *pop_functor_ns);
  sLatticeAC.defineField<descriptors::POPULATION>(all, *pop_functor_ac);

  // It says it could take a superF as u ... but it cant afaik, in practice, sadly
  // sLatticeAC.defineRho(all, *phi_new);
  // sLatticeAC.iniEquilibrium(all, *phi_new, *u_ac);
  // sLatticeNS.defineRho(all, *p_new);
  // sLatticeNS.iniEquilibrium(all, *p_new, *u_ns);

  clout << "Finished bubble insertion" << std::endl;
}

void getResults(SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                SuperLattice<T, ACDESCRIPTOR>& sLatticeAC, int iT,
                SuperGeometry<T, 2>& superGeometry, util::Timer<T>& timer,
                UnitConverter<T, NSDESCRIPTOR> converter)
{
  OstreamManager      clout(std::cout, "getResults");
  SuperVTMwriter2D<T> vtmWriter("bubbleChannel2d");
  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(sLatticeNS);
    SuperLatticeRank2D<T, NSDESCRIPTOR>   rank(sLatticeNS);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }
  // Get statistics
  if (iT % statIter == 0) {
    // Timer console output
    timer.update(iT);
    timer.printStep();
    sLatticeNS.getStatistics().print(iT, converter.getPhysTime(iT));
    sLatticeAC.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  // Writes the VTK files
  if (iT % vtkIter == 0) {
    SuperLatticePhysIncPressure2D<T, NSDESCRIPTOR> p_total(sLatticeNS,
                                                           converter);
    p_total.getName() = "p_total";
    SuperLatticeExternalScalarField2D<T, NSDESCRIPTOR, RHO> rho_L(sLatticeNS);
    AnalyticalConst2D<T, T> C_rho_(converter.getConversionFactorDensity());
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> C_rho_field(C_rho_,
                                                                sLatticeNS);
    SuperIdentity2D<T, T>                           rho(C_rho_field * rho_L);
    rho.getName() = "rho";
    SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity(sLatticeNS, converter);
    velocity.getName() = "u";
    SuperLatticeField2D<T, ACDESCRIPTOR, STATISTIC> phi(sLatticeAC);
    phi.getName() = "phi";
    SuperLatticeExternalScalarField2D<T, ACDESCRIPTOR, PSI> psi(sLatticeAC);
    psi.getName() = "psi";
    SuperLatticeExternalScalarField2D<T, NSDESCRIPTOR, SCALAR> scale(
        sLatticeNS);
    scale.getName() = "scale";

      SuperLatticeExternalScalarField2D<T, NSDESCRIPTOR, TAU_EFF> tau(sLatticeNS);
      tau.getName() = "tau_eff";

    SuperLatticeDensity2D<T,NSDESCRIPTOR>   p_raw(sLatticeNS);
  SuperLatticeVelocity2D<T,NSDESCRIPTOR>  u_ns_raw (sLatticeNS);
  SuperLatticeDensity2D<T,ACDESCRIPTOR>   phi_raw(sLatticeAC);
  SuperLatticeVelocity2D<T,ACDESCRIPTOR>  u_ac_raw (sLatticeAC);

  T dx        = converter.getPhysDeltaX();
  std::vector<T> pos = {dx * 2. * diameter, dx * Ny / 2. * 1};
  std::shared_ptr<AnalyticalF2D<T,T>> phi_new_bubble_anal = std::make_shared<CircularInterface2D<T>>(pos, dx * diameter / 2., dx * w, 1.,
                                        true);
  SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> phi_new_bubble(phi_new_bubble_anal, sLatticeNS);       

  p_raw.getName() = "p_raw";
  u_ns_raw.getName() = "u_ns_raw";
  phi_raw.getName() = "phi_raw";
  u_ac_raw.getName() = "u_ac_raw";

  phi_new_bubble.getName() = "phi_new_bubble";

  vtmWriter.addFunctor(p_raw);
  vtmWriter.addFunctor(u_ns_raw);
  vtmWriter.addFunctor(phi_raw);
  vtmWriter.addFunctor(u_ac_raw);

  vtmWriter.addFunctor(phi_new_bubble);

      vtmWriter.addFunctor(tau);
    vtmWriter.addFunctor(scale);

    vtmWriter.addFunctor(p_total);
    vtmWriter.addFunctor(rho);
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(phi);
    vtmWriter.addFunctor(psi);
    vtmWriter.write(iT);
  }
}

void simulate(T Re)
{
  OstreamManager clout(std::cout, "main");
  // === 1st Step: Initialization ===
  MultiPhaseUnitConverterFromRelaxationTime<T, NSDESCRIPTOR> converter(
      int {diameter},     // resolution
      (T)tau_l,           // lattice relaxation time
      (T)rhos[1] / C_rho, // lattice density heavier phase
      (T)L_char, // charPhysLength: reference length of simulation geometry in __m__
      (T)viscosityH2O, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)rhos[1]       // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  T u_phys = Re * viscosityH2O / (L_char / diameter * (Ny - 1.));
  T dx     = converter.getPhysDeltaX();
  clout << "Re = " << Re << std::endl;
  clout << "We = " << rhos[0] * u_phys * u_phys * L_char / surfaceTension
        << std::endl;
  clout << "u_l = " << Re / Ny * ((tau_l - 0.5) / 3.) << std::endl;
  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a cuboidGeometry with weights
  Vector<T, 2>         extend(Nx * dx, Ny * dx);
  Vector<T, 2>         origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);
#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, dx,
                                               singleton::mpi().getSize());
#else
  CuboidDecomposition2D<T> cuboidDecomposition(cuboid, dx);
#endif
  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  loadBalancer.print();
  // Instantiation of superGeometry
  SuperGeometry<T, 2> superGeometry(cuboidDecomposition, loadBalancer);
  prepareGeometry(superGeometry, dx);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, NSDESCRIPTOR> sLatticeNS(superGeometry);
  SuperLattice<T, ACDESCRIPTOR> sLatticeAC(superGeometry);
  SuperLatticeCoupling coupling(LiangPostProcessor {}, names::NavierStokes {},
                                sLatticeNS, names::Component1 {}, sLatticeAC);
  coupling.restrictTo(superGeometry.getMaterialIndicator({1, 4}));

  SuperLatticeCoupling velocityCoupling(VelocityCoupling {},
                                        names::NavierStokes {}, sLatticeNS,
                                        names::Component1 {}, sLatticeAC);
  velocityCoupling.restrictTo(superGeometry.getMaterialIndicator(3));
  prepareLattice(sLatticeNS, sLatticeAC, coupling, converter, superGeometry);


  // === 4th Step: Main Loop with Timer ===
  int iT  = 0;
  maxIter = Nx / (Re / Ny * ((tau_l - 0.5) / 3.)) * 1.3;
  clout << "starting simulation..." << std::endl;
  clout << "with " << maxIter << " timesteps " << std::endl;
  util::Timer<T> timer(maxIter, superGeometry.getStatistics().getNvoxel());
  timer.start();
  for (iT = 0; iT <= maxIter; ++iT) {

    // Collide and stream (and coupling) execution
    sLatticeNS.collideAndStream();
    sLatticeAC.collideAndStream();

    
    
    // Bubble
    if(iT == 100){
      std::vector<T> pos = {dx * 2. * diameter, dx * Ny / 2. * 1};

      insert_bubble(pos, diameter, converter, sLatticeNS, sLatticeAC, superGeometry);
    }



    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeAC.executePostProcessors(stage::PreCoupling());

    if (iT % 200 == 0) {
      signedDistanceFunction<stage::IterativePostProcess>(sLatticeAC,
                                                          superGeometry, w);
      sLatticeAC.executePostProcessors(stage::PhiLimiter());
    }
    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

    coupling.execute();
    velocityCoupling.execute();

    T uMax = helperConvectiveU(sLatticeNS, superGeometry, dx);
    sLatticeAC.setParameter<descriptors::MAX_VELOCITY>(uMax);

    // Computation and output of the results
    getResults(sLatticeNS, sLatticeAC, iT, superGeometry, timer, converter);
    if (std::isnan(sLatticeNS.getStatistics().getAverageEnergy())) {
      clout << "NS lattice has NaN average energy" << std::endl;
      break;
    }
    
    
  }
  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  if (argc > 1) {
    Re = atof(argv[1]);
  }

  singleton::directories().setOutputDir("./tmp/");
  simulate(Re);
}
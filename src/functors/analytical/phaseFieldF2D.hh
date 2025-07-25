/*  This file is part of the OpenLB library
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

#ifndef ANALYTICAL_PHASE_FIELD_F_2D_HH
#define ANALYTICAL_PHASE_FIELD_F_2D_HH

#include<vector>
#include<cmath>
#include<string>

#include "phaseFieldF2D.h"
#include "functors/genericF.h"
#include "analyticalF.h"
#include "functors/lattice/superBaseF2D.h"
#include "geometry/superGeometry.h"

#include "core/superLattice2D.h"
#include "utilities/vectorHelpers.h"  // for normalize


namespace olb {

template <typename T>
CircularInterface2D<T>::CircularInterface2D(std::vector<T> center, T radius, T interfaceWidth, T factor, bool bubble) : AnalyticalF2D<T,T>(1)
{
  this->getName() = "CircularInterface2D";
  _center.resize(2);
  for (int i = 0; i < 2; ++i) {
    _center[i] = center[i];
  }
  _interfaceWidth = interfaceWidth;
  _radius = radius;
  _factor = factor;
  _bubble = bubble;
}

template <typename T>
bool CircularInterface2D<T>::operator()(T output[], const T x[])
{
  T distToCenter2 = util::pow(_center[0]-x[0], 2) + util::pow(_center[1]-x[1], 2);
  T d = util::sqrt(distToCenter2) - _radius;
  if (_bubble) output[0] = T(_factor*(1.+tanh(2*d/_interfaceWidth))/2.);
  else output[0] = T(_factor*(1.-tanh(2*d/_interfaceWidth))/2.);
  return true;
}

template <typename T>
CircularFadingInterface2D<T>::CircularFadingInterface2D(FunctorPtr<AnalyticalF2D<T,T>> phi, std::vector<T> center, T radius, T interfaceWidth, int current_step, int fade_steps) 
: AnalyticalF2D<T,T>(1), _phi{std::move(phi)}
{
  this->getName() = "CircularFadingInterface2D";
  _center.resize(2);
  for (int i = 0; i < 2; ++i) {
    _center[i] = center[i];
  }
  _interfaceWidth = interfaceWidth;
  _radius = radius;
  _fade_steps = fade_steps;
  _current_step = current_step;
}

template <typename T>
bool CircularFadingInterface2D<T>::operator()(T output[], const T input[])
{
  constexpr double eps = 1e-12;

  T phi_val;
  _phi(&phi_val, input);

  if(phi_val > 0.999){
    output[0] = phi_val;
    return true;
  }


  int new_step = _current_step + 1; // so this starts at 2, and thats the intention
  T phi = ((phi_val * (T)_fade_steps) - (T)(_fade_steps - _current_step)) / T(_current_step);

  if(phi < 0.){
    // std::cout << typeid(*_phi).name() << std::endl;
    // printf("Phi < 0 (%lf %d %d)\n", phi_val, _fade_steps, _current_step);
    phi = 0.;
  }

  if(phi < -0.1){
    // std::cout << typeid(*_phi).name() << std::endl;
    printf("Phi < 0.1 (%lf %d %d)\n", phi_val, _fade_steps, _current_step);
  }

  if(phi > 1.){
    printf("Phi > 1 (%lf %d %d)\n", phi_val, _fade_steps, _current_step);
  }

  // Now get d from phi
  double u = phi;
  u = std::clamp(u, eps, 1.0 - eps);

    /*  log(u) - log1p(-u)  ==  log( u / (1-u) )
        log1p(x) gives full precision when x≈0,
        so the difference stays accurate even for u≈1.                */
  T d = 0.25 * _interfaceWidth * ( std::log(u) - std::log1p(-u) );

  // T d = _interfaceWidth / 2. * atanh(2 * phi - 1);
  T new_phi = T((1. + tanh(2. * d / _interfaceWidth)) / 2.);

  // Stability fix: in the interior of the bubble, the inverse is ill-conditioned
  if(phi <= 0.01){
    new_phi = eps;
  }

  T new_scaled_phi = (new_phi * T(new_step) + T(_fade_steps - new_step)) / T(_fade_steps);
  new_scaled_phi = std::clamp(new_scaled_phi, eps, 1.0 - eps);

  output[0] = new_scaled_phi;

  if(isnan(output[0])){
    printf("Got NaN in phi\n");
  }

  // if(output[0] < 0.01){
  //   printf("Writing close 0 output\n");
  // }

  return true;
}


template <typename T>
LaplacePressure2D<T>::LaplacePressure2D(std::vector<T> center, T radius, T interfaceWidth, T surfaceTension) : AnalyticalF2D<T,T>(1)
{
  this->getName() = "LaplacePressure2D";
  _center.resize(2);
  for (int i = 0; i < 2; ++i) {
    _center[i] = center[i];
  }
  _interfaceWidth = interfaceWidth;
  _radius = radius;
  _surfaceTension = surfaceTension;
}

template <typename T>
bool LaplacePressure2D<T>::operator()(T output[], const T x[])
{
  T distToCenter2 = util::pow(_center[0]-x[0], 2) + util::pow(_center[1]-x[1], 2);
  T d = util::sqrt(distToCenter2) - _radius;
  T tanh_d = tanh(2*d/_interfaceWidth);
  output[0] = _surfaceTension/(_radius*4)*tanh_d*(tanh_d*tanh_d-3.);
  return true;
}


template <typename T>
ShiftedLaplacePressure2D<T>::ShiftedLaplacePressure2D(std::vector<T> center, T radius, T interfaceWidth, T surfaceTension, T pressureShift): AnalyticalF2D<T,T>(1)
{
  this->getName() = "ShiftedLaplacePressure2D";
  _center.resize(2);
  for (int i = 0; i < 2; ++i) {
    _center[i] = center[i];
  }
  _interfaceWidth = interfaceWidth;
  _radius = radius;
  _surfaceTension = surfaceTension;
  _pressureShift = pressureShift;
}

template <typename T>
bool ShiftedLaplacePressure2D<T>::operator()(T output[], const T x[])
{
  T distToCenter2 = util::pow(_center[0]-x[0], 2) + util::pow(_center[1]-x[1], 2);
  T d = util::sqrt(distToCenter2) - _radius;
  T tanh_d = tanh(2*d/_interfaceWidth);
  output[0] = _surfaceTension/(_radius*4)*tanh_d*(tanh_d*tanh_d-3.) + _pressureShift;
  return true;
}


template <typename T, typename DESCRIPTOR>
IncompressibleEquilibriumPopulations2D<T, DESCRIPTOR>::IncompressibleEquilibriumPopulations2D(FunctorPtr<AnalyticalF2D<T,T>> density, FunctorPtr<AnalyticalF2D<T,T>> pressure, FunctorPtr<AnalyticalF2D<T,T>> velocity): AnalyticalF2D<T,T>(DESCRIPTOR::q)
  , _density{std::move(density)}, _pressure{std::move(pressure)}, _velocity{std::move(velocity)}
{
}

template <typename T, typename DESCRIPTOR>
bool IncompressibleEquilibriumPopulations2D<T, DESCRIPTOR>::operator()(T output[], const T input[])
{
  T d;
  T p;
  T v [2];
  this->_density(&d, input);
  this->_pressure(&p, input);
  this->_velocity(v, input);

  for(int iPop = 0; iPop < 9; ++iPop){
    output[iPop] = olb::equilibrium<DESCRIPTOR>::mpincompressible(iPop, d, v, p);
  }
  
  return true;
}

template <typename T, typename DESCRIPTOR>
FirstOrderEquilibriumPopulations2D<T, DESCRIPTOR>::FirstOrderEquilibriumPopulations2D(FunctorPtr<AnalyticalF2D<T,T>> density, FunctorPtr<AnalyticalF2D<T,T>> velocity): AnalyticalF2D<T,T>(DESCRIPTOR::q)
  , _density{std::move(density)}, _velocity{std::move(velocity)}
{

}

template <typename T, typename DESCRIPTOR>
bool FirstOrderEquilibriumPopulations2D<T, DESCRIPTOR>::operator()(T output[], const T input[])
{
  T d;
  T v [2];
  this->_density(&d, input);
  this->_velocity(v, input);

  for(int iPop = 0; iPop < 9; ++iPop){
    output[iPop] = olb::equilibrium<DESCRIPTOR>::firstOrder(iPop, d, v);
  }
  
  return true;
}



} // end namespace olb
#endif

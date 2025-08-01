/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Jakob Mangold, Mathias J. Krause, 2024 Julius Jeßberger
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

#ifndef SUPER_STATISTIC_F3D_H
#define SUPER_STATISTIC_F3D_H

#include "superBaseF3D.h"
#include "blockStatisticF3D.h"
#include "indicator/superIndicatorBaseF3D.h"
#include "utilities/functorPtr.h"

namespace olb {


/// SuperVarianceF3D returns the variance in each component of f on a indicated subset
template <typename T, typename W = T>
class SuperVarianceF3D : public SuperF3D<T,W> {
protected:
  FunctorPtr<SuperF3D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
  T _expectedValue;
  SuperConst3D<T> _expectedValueF;
public:
  /// Constructor for determining the standard deviation of f on a indicated subset
  /**
   * \param f          functor of which the standard deviation is to be determined
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
//  SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&&        f,
//                 FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  /// Constructor for determining the standard deviation of f on a given material
  /**
   * \param f             functor of which the average is to be determined
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
//  SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&& f,
//                 SuperGeometry<T,3>& superGeometry,
//                const int material);




  SuperVarianceF3D(FunctorPtr<SuperF3D<T,W>>&& f,
                   SuperGeometry<T,3>& superGeometry,
                   const int material,
                   T expectedValue);




  SuperVarianceF3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                   FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                   T expectedValue);


  /// Global average operator
  /**
   * Note: While this functor exposes BlockStdDeviation3D functors if possible, a call to
   * this function will not use them but calculate the global standard deviation by summing all
   * components and voxel counts.
   * Calling BlockStdDevation3D in this situation would unnecessarily complicate this as
   * we would have to weight the aggregated averages according to their share in the
   * global average.
   **/
  bool operator() (W output[], const int input[]) override;
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/// SuperStdDeviaitonF3D returns the standard deviation in each component of f on a indicated subset
template <typename T, typename W = T>
class SuperStdDeviationF3D final : public SuperVarianceF3D<T,W> {

public:
  /// Constructor for determining the standard deviation of f on a indicated subset
  /**
   * \param f          functor of which the standard deviation is to be determined
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
//  SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&&        f,
//                 FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  /// Constructor for determining the standard deviation of f on a given material
  /**
   * \param f             functor of which the average is to be determined
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
//  SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&& f,
//                 SuperGeometry<T,3>& superGeometry,
//                const int material);




  SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&& f,
                       SuperGeometry<T,3>& superGeometry,
                       const int material,
                       T expectedValue);




  SuperStdDeviationF3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                       FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                       T expectedValue);


  bool operator() (W output[], const int input[]) override;
};


}

#endif


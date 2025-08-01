/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *                          Albert Mink, Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_BASE_F_3D_H
#define SUPER_BASE_F_3D_H

#include <memory>

#include "functors/genericF.h"
#include "blockBaseF3D.h"
#include "indicator/superIndicatorBaseF3D.h"
#include "communication/superStructure.h"
#include "utilities/functorPtr.h"
#include "core/superData.h"
#include "core/superLattice.h"


/* Note: Throughout the whole source code directory functors, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

template<typename T, typename BaseType> class SuperData3D;
template<typename T, typename DESCRIPTOR> class SuperLattice;
template<typename T> class SuperStructure3D;
template<typename T> class BlockF3D;
template<typename T> class SuperIndicatorF3D;
template<typename T, typename W> class SuperIdentity3D;
template<typename T, typename DESCRIPTOR> class SuperLatticeIdentity3D;

/// represents all functors that operate on a SuperStructure<T,3> in general
template <typename T, typename W = T>
class SuperF3D : public GenericF<W,int> {
protected:
  SuperF3D(SuperStructure<T,3>& superStructure, int targetDim);

  SuperStructure<T,3>& _superStructure;
  /// Super functors may consist of several BlockF3D<W> derived functors
  /**
   * By convention: If block level functors are used at all they should
   * number exactly LoadBalancer<T>::size per process.
   **/
  std::vector<std::unique_ptr<BlockF3D<W>>> _blockF;
public:
  using identity_functor_type = SuperIdentity3D<T,W>;
  static constexpr bool isSuper = true;
  static constexpr unsigned d = 3;

  SuperF3D<T,W>& operator-(SuperF3D<T,W>& rhs);
  SuperF3D<T,W>& operator+(SuperF3D<T,W>& rhs);
  SuperF3D<T,W>& operator*(SuperF3D<T,W>& rhs);
  SuperF3D<T,W>& operator/(SuperF3D<T,W>& rhs);

  /// \return SuperF3D<T,W>::_superStructure
  SuperStructure<T,3>& getSuperStructure();
  /// \return Size of SuperF3D<T,W>::_blockF vector
  int getBlockFSize() const;
  /// \return SuperF3D<T,W>::_blockF[iCloc]
  BlockF3D<W>& getBlockF(int iCloc);

  bool operator() (W output[], const int input []);

  using GenericF<W,int>::operator();
};


/// Functor from `SuperData3D`
template<typename T, typename BaseType>
class SuperDataF3D : public SuperF3D<T,BaseType> {
protected:
  /// `SuperData3D` object this functor was created from
  SuperData<3,T,BaseType>& _superData;
public:
  /// Constructor from `SuperData3D` - stores `_superData` reference
  SuperDataF3D(SuperData<3,T,BaseType>& superData);
  /// Operator for this functor - copies data from `_superData` object into output
  bool operator() (BaseType output[], const int input[]);
  /// Getter for `_superData`
  SuperData<3,T,BaseType>& getSuperData();
};


/// identity functor for memory management
template <typename T, typename W=T>
class SuperIdentity3D : public SuperF3D<T,W> {
protected:
  FunctorPtr<SuperF3D<T,W>> _f;
public:
  SuperIdentity3D(FunctorPtr<SuperF3D<T,W>>&& f);
  bool operator() (W output[], const int input[]) override;
};


/// functor to extract one component
template <typename T, typename W=T>
class SuperExtractComponentF3D : public SuperF3D<T,W> {
protected:
  FunctorPtr<SuperF3D<T,W>> _f;
  const int _extractDim;
public:
  SuperExtractComponentF3D(FunctorPtr<SuperF3D<T,W>>&& f, int extractDim);
  int getExtractDim();
  bool operator() (W output[], const int input[]);
};


/// functor to extract one component inside an indicator
template <typename T, typename W=T>
class SuperExtractComponentIndicatorF3D : public SuperExtractComponentF3D<T,W> {
protected:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperExtractComponentIndicatorF3D(FunctorPtr<SuperF3D<T,W>>&& f, int extractDim,
                                    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  bool operator() (W output[], const int input[]) override;
};


/// functor to extract data inside an indicator
template <typename T, typename W=T>
class SuperExtractIndicatorF3D : public SuperF3D<T,W> {
protected:
  FunctorPtr<SuperF3D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperExtractIndicatorF3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                           FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  bool operator() (W output[], const int input[]);
};

/// perform explicit typecast from output type W2 to W
// user has to guarantee that cast is well-defined
template <typename T, typename W, typename W2>
class SuperTypecastF3D : public SuperF3D<T,W> {
protected:
  FunctorPtr<SuperF3D<T,W2>>        _f;
public:
  SuperTypecastF3D(FunctorPtr<SuperF3D<T,W2>>&& f);
  bool operator() (W output[], const int input[]);
};

/// identity functor for memory management
template <typename T, typename W=T>
class SuperIdentityOnSuperIndicatorF3D : public SuperF3D<T,W> {
protected:
  SuperF3D<T,W>& _f;
  SuperIndicatorF3D<T>& _indicatorF;
  W _defaultValue;
public:
  SuperIdentityOnSuperIndicatorF3D(SuperF3D<T,W>& f, SuperIndicatorF3D<T>& indicatorF, W defaultValue=0.);
  bool operator() (W output[], const int input[]);
};

/// represents all functors that operate on a SuperLattice in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, typename DESCRIPTOR>
class SuperLatticeF3D : public SuperF3D<T,T> {
protected:
  SuperLatticeF3D(SuperLattice<T,DESCRIPTOR>& superLattice, int targetDim);

  SuperLattice<T,DESCRIPTOR>& _sLattice;
public:
  using identity_functor_type = SuperLatticeIdentity3D<T,DESCRIPTOR>;

  SuperLattice<T,DESCRIPTOR>& getSuperLattice();

  bool operator() (T output [], const int input []);

  using GenericF<T,int>::operator();
};

/// identity functor for memory management
template <typename T, typename DESCRIPTOR>
class SuperLatticeIdentity3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  FunctorPtr<SuperLatticeF3D<T,DESCRIPTOR>> _f;
public:
  SuperLatticeIdentity3D(FunctorPtr<SuperLatticeF3D<T,DESCRIPTOR>>&& f);
  bool operator() (T output[], const int input[]) override;
};

/// represents all functors that operate on a DESCRIPTOR with output in Phys, e.g. physVelocity(), physForce(), physPressure()
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysF3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  SuperLatticePhysF3D(SuperLattice<T,DESCRIPTOR>& sLattice,
                      const UnitConverter<T,DESCRIPTOR>& converter, int targetDim);
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  UnitConverter<T,DESCRIPTOR> const& getConverter() const;
};

/// represents all thermal functors that operate on a DESCRIPTOR with output in Phys, e.g. physTemperature(), physHeatFlux()
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticeThermalPhysF3D : public SuperLatticeF3D<T,TDESCRIPTOR> {
protected:
  SuperLatticeThermalPhysF3D(SuperLattice<T,TDESCRIPTOR>& sLattice,
                             const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter, int targetDim);
  const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& _converter;
public:
  ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& getConverter() const;
};


template <typename T, typename DESCRIPTOR>
class ComposedSuperLatticeF3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeF3D<T,DESCRIPTOR>& _f0;
  SuperLatticeF3D<T,DESCRIPTOR>& _f1;
  SuperLatticeF3D<T,DESCRIPTOR>& _f2;
public:
  ComposedSuperLatticeF3D(SuperLatticeF3D<T,DESCRIPTOR>& f0,
                          SuperLatticeF3D<T,DESCRIPTOR>& f1,
                          SuperLatticeF3D<T,DESCRIPTOR>& f2);
  bool operator() (T output[], const int x[]) override;
};


/// Generate a SuperLatticeF from an arbitrary function
template <typename T, typename DESCRIPTOR>
class SuperLatticeFfromCallableF final : public SuperLatticeF<T,DESCRIPTOR> {
public:
  /// Generate a SuperLatticeF from an arbitrary function
  /// @tparam F function type
  /// @param f some function, e.g. lambda expression
  /// The following types of f are accepted:
  /// void(T* output, auto cell)
  /// void(T* output, auto cell, int iCuboid, const int* localCoordinates)
  /// void(T* output, int iCuboid, const int* localCoordinates)
  // iCuboid is the local cuboid index
  template <typename F>
  SuperLatticeFfromCallableF(SuperLattice<T,DESCRIPTOR>& superLattice, F&& f)
   : SuperLatticeF<T,DESCRIPTOR>(superLattice, 1)
  {
    this->getName() = "FromCallable";

    LoadBalancer<T>& load = this->getSuperStructure().getLoadBalancer();
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockLatticeFfromCallableF<T,DESCRIPTOR>(
          this->_sLattice.getBlock(iC), iC, f));
    }
  }
};


} // end namespace olb

#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2014 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef SUPER_LATTICE_LOCAL_F_3D_HH
#define SUPER_LATTICE_LOCAL_F_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm

#include "functors/superLatticeLocalF3D.h"
#include "functors/blockLatticeLocalF3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity

namespace olb {

//template <typename T, template <typename U> class DESCRIPTOR> class BlockLatticeStrainRate3D;


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeFpop3D<T,DESCRIPTOR>::SuperLatticeFpop3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,DESCRIPTOR<T>::q)
{ this->_name = "fPop"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeFpop3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;
    
    BlockLatticeFpop3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)));

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeDissipation3D<T,DESCRIPTOR>::SuperLatticeDissipation3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1), _converter(converter)
{ this->_name = "dissipation"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeDissipation3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;
    
    BlockLatticeDissipation3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), this->_converter);

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDissipation3D<T,DESCRIPTOR>::SuperLatticePhysDissipation3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1)
{ this->_name = "physDissipation"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysDissipation3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;
    
    BlockLatticePhysDissipation3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), this->_converter);

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeDensity3D<T,DESCRIPTOR>::SuperLatticeDensity3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1)
{ this->_name = "density"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeDensity3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;
    
    BlockLatticeDensity3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)) );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeVelocity3D<T,DESCRIPTOR>::SuperLatticeVelocity3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,3)
{ this->_name = "velocity"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeVelocity3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticeVelocity3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)) );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeStrainRate3D<T,DESCRIPTOR>::SuperLatticeStrainRate3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice,const LBconverter<T>& converter)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,9), _converter(converter)
{ this->_name = "strainRate"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeStrainRate3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticeStrainRate3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), this->_converter );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysStrainRate3D<T,DESCRIPTOR>::SuperLatticePhysStrainRate3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,9)
{ this->_name = "physStrainRate"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysStrainRate3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticePhysStrainRate3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), this->_converter );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeGeometry3D<T,DESCRIPTOR>::SuperLatticeGeometry3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1), _superGeometry(superGeometry),
    _material(material)
{ this->_name = "geometry"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeGeometry3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticeGeometry3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)),
                                                        this->_superGeometry.getExtendedBlockGeometry(this->_sLattice.getLoadBalancer().loc(globIC)), _material );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeRank3D<T,DESCRIPTOR>::SuperLatticeRank3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1)
{ this->_name = "rank"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeRank3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0]; // int locix= input[1]; int lociy= input[2]; int lociz= input[3];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    T rank[1];
    rank[0] = singleton::mpi().getRank() + 1;
    std::vector<T> out(rank, rank+1);  // first adress, last adress
    return out;
  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeCuboid3D<T,DESCRIPTOR>::SuperLatticeCuboid3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1)
{ this->_name = "cuboid"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeCuboid3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
  int globIC = input[0]; // int locix= input[1]; int lociy= input[2]; int lociz= input[3];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    T cuboid[1];
    cuboid[0] = globIC+1;
    std::vector<T> out(cuboid, cuboid+1);  // first adress, last adress
    return out;
  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysPressure3D<T,DESCRIPTOR>::SuperLatticePhysPressure3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1)
{ this->_name = "physPressure"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysPressure3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;
    
    BlockLatticePhysPressure3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)) , this->_converter) ;

    return blockLatticeF(inputLocal);

  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysVelocity3D<T,DESCRIPTOR>::SuperLatticePhysVelocity3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter, bool print)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,3), _print(print)
{ this->_name = "physVelocity";}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysVelocity3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticePhysVelocity3D<T,DESCRIPTOR> blockLatticeF(this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), this->_converter, _print);
    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysBoundaryForce3D<T,DESCRIPTOR>::SuperLatticePhysBoundaryForce3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,3),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "physBoundaryForce"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysBoundaryForce3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];

  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<T> force(3, T());
    if(this->_superGeometry.get(input) == _material) {
      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
        // Get direction
        const int* c = DESCRIPTOR<T>::c[iPop];
        // Get next cell located in the current direction
        // Check if the next cell is a fluid node
        if (this->_superGeometry.get(input[0], input[1] + c[0], input[2] + c[1], input[3] + c[2]) == 1) {
          int overlap = this->_sLattice.getOverlap();
          // Get f_q of next fluid cell where l = opposite(q)
          T f = this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
          // Get f_l of the boundary cell
          // Add f_q and f_opp
          f += this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
          // Update force
          force[0] -= c[0]*f;
          force[1] -= c[1]*f;
          force[2] -= c[2]*f;
        }
      }
      force[0] = this->_converter.physForce(force[0]);
      force[1] = this->_converter.physForce(force[1]);
      force[2] = this->_converter.physForce(force[2]);
      return force;
    }
    else {
      return force;
    }
  }
  else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR>::SuperLatticePhysCorrBoundaryForce3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,3),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "physCorrBoundaryForce"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];

  std::vector<T> force(3, T());
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    if(this->_superGeometry.get(input) == _material) {
      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
        // Get direction
        const int* c = DESCRIPTOR<T>::c[iPop];
        // Get next cell located in the current direction
        // Check if the next cell is a fluid node
        if (this->_superGeometry.get(input[0], input[1] + c[0], input[2] + c[1], input[3] + c[2]) == 1) {
          int overlap = this->_sLattice.getOverlap();
          // Get f_q of next fluid cell where l = opposite(q)
          T f = this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
          // Get f_l of the boundary cell
          // Add f_q and f_opp
          f += this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
          // Update force
          force[0] -= c[0]*(f-2.*DESCRIPTOR<T>::t[iPop]);
          force[1] -= c[1]*(f-2.*DESCRIPTOR<T>::t[iPop]);
          force[2] -= c[2]*(f-2.*DESCRIPTOR<T>::t[iPop]);
        }
      }
      force[0] = this->_converter.physForce(force[0]);
      force[1] = this->_converter.physForce(force[1]);
      force[2] = this->_converter.physForce(force[2]);
      return force;
    }
    else {
      return force;
    }
  }
  else {
    return force;
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePorosity3D<T,DESCRIPTOR>::SuperLatticePorosity3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1)
{ this->_name = "porosity"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePorosity3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;
    
    BlockLatticePorosity3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)) ) ;

    return blockLatticeF(inputLocal);

  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysPermeability3D<T,DESCRIPTOR>::SuperLatticePhysPermeability3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "permeability"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysPermeability3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];

  T* value = new T[1];
  int overlap = this->_sLattice.getOverlap();
  this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeExternalField(0,1,value);
  std::vector<T> result(1,this->_converter.physPermeability(value[0]));
  delete value;
  if (!(result[0]<42) &&! (result[0]>42) &&! (result[0]==42)) result[0] = 999999;
  if (std::isinf(result[0])) result[0] = 1e6;
  return result;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDarcyForce3D<T,DESCRIPTOR>::SuperLatticePhysDarcyForce3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,3),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "alphaU"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysDarcyForce3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  SuperLatticePhysPermeability3D<T,DESCRIPTOR> permeability(this->_sLattice,this->_superGeometry,this->_material,this->_converter);
  SuperLatticeVelocity3D<T,DESCRIPTOR> velocity(this->_sLattice);

  T nu = this->_converter.getCharNu();
  T K = permeability(input)[0];
  std::vector<T> u = velocity(input);

  std::vector<T> result(3,0);
  result[0] = -nu/K*u[0];
  result[1] = -nu/K*u[1];
  result[2] = -nu/K*u[2];

  return result;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeAverage3D<T,DESCRIPTOR>::SuperLatticeAverage3D
  (SuperLatticeF3D<T,DESCRIPTOR>& f, SuperGeometry3D<T>& superGeometry,
  const int material, T radius)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(), f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material), _radius(radius)
{ this->_name = "Average("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeAverage3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();

  //create boolean indicator functor isInSphere
  std::vector<T> center = cGeometry.getPhysR(input);
  IndicatorSphere3D<bool,T> isInSphere(center,_radius);

  // iterate over all cuboids & points and test for material && isInSphere
  std::vector<T> tmp( this->_n, T() );
  int numVoxels(0);
  if (this->_superGeometry.get(input) == _material) {
    for (int iC = 0; iC < load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          for (int iZ = 0; iZ < nZ; iZ++) {
            std::vector<int> testLatticeR(input); 
            testLatticeR[0] = load.glob(iC); 
            testLatticeR[1] = iX;
            testLatticeR[2] = iY;
            testLatticeR[3] = iZ;
            std::vector<T> testPhysR = cGeometry.getPhysR(testLatticeR);
            if (this->_superGeometry.get(input) == _material
                && isInSphere(testPhysR)[0] == true) {
              for (unsigned iD = 0; iD < _f(load.glob(0),0,0,0).size(); iD++) {
                tmp[iD] += _f(load.glob(iC),iX,iY,iZ)[iD];
              }
              numVoxels++;
            }
          }
        }
      }
    }

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
    for (int iD = 0; iD < _f.getTargetDim(); iD++) {
#ifdef PARALLEL_MODE_MPI
      singleton::mpi().reduceAndBcast(tmp[iD], MPI_SUM);
#endif
      if (numVoxels > 0) {
        tmp[iD] /= numVoxels;
      }
    }
  }
  return tmp;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperEuklidNorm3D<T,DESCRIPTOR>::SuperEuklidNorm3D(SuperLatticeF3D<T,DESCRIPTOR>& f)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),1), _f(f)
{
  this->_name = "EuklidNorm("+_f.getName()+")"; 
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperEuklidNorm3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  //SuperLatticeIdentity3D<T,DESCRIPTOR> ff(f);
  std::vector<T> data = _f(input);
  std::vector<T> tmp(1,0);
  for (unsigned i = 0; i < data.size(); i++) {
    tmp[0] += data[i]*data[i];
  }
  tmp[0] = sqrt(tmp[0]);
  return tmp;
}


} // end namespace olb

#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#ifndef SUPER_LATTICE_LOCAL_F_2D_HH
#define SUPER_LATTICE_LOCAL_F_2D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm

#include "functors/superLatticeLocalF2D.h"
#include "functors/blockLatticeLocalF2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity

namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeDissipation2D<T,DESCRIPTOR>::SuperLatticeDissipation2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice,1), _converter(converter)
{ this->_name = "dissipation"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeDissipation2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
////  int lociz= input[3];
//  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
//    // local coords are given, fetch local cell and compute value(s)
//    T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
//    int overlap = this->sLattice.getOverlap();
//    this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap/*, lociz+overlap*/).computeAllMomenta(rho, uTemp, pi);

//    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
//    if (util::TensorVal<DESCRIPTOR<T> >::n == 6)
//      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

//    T nuLattice = converter.getLatticeNu();
//    T omega = converter.getOmega();
//    T finalResult = PiNeqNormSqr*nuLattice*pow(omega*DESCRIPTOR<T>::invCs2,2)/rho/2.;

//    return std::vector<T>(1,finalResult);
//  } else {
//    return std::vector<T>(); // empty vector
//  }
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDissipation2D<T,DESCRIPTOR>::SuperLatticePhysDissipation2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,1)
{ this->_name = "physDissipation"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysDissipation2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
////  int lociz= input[3];
//  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
//    // local coords are given, fetch local cell and compute value(s)
//    T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
//    int overlap = this->sLattice.getOverlap();
//    this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeAllMomenta(rho, uTemp, pi);

//    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
//    if (util::TensorVal<DESCRIPTOR<T> >::n == 6)
//      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

//    T nuLattice = converter.getLatticeNu();
//    T omega = converter.getOmega();
//    T finalResult = PiNeqNormSqr*nuLattice*pow(omega*DESCRIPTOR<T>::invCs2,2)/rho/2.*converter.getCharNu()/converter.getLatticeNu()/converter.getDeltaT()/converter.getDeltaT();

//    return std::vector<T>(1,finalResult);
//  } else {
//    return std::vector<T>(); // empty vector
//  }
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeDensity2D<T,DESCRIPTOR>::SuperLatticeDensity2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice) : SuperLatticeF2D<T,DESCRIPTOR>(sLattice,1)
{ this->_name = "density"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeDensity2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(2,int());
    int overlap = this->_sLattice.getOverlap(); 
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;

    BlockLatticeDensity2D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)) );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeVelocity2D<T,DESCRIPTOR>::SuperLatticeVelocity2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice) : SuperLatticeF2D<T,DESCRIPTOR>(sLattice,2)
{ this->_name = "velocity"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeVelocity2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(2,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;

    BlockLatticeVelocity2D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)) );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}



template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysStrainRate2D<T,DESCRIPTOR>::SuperLatticePhysStrainRate2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,4)
{ this->_name = "physStrainRate"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysStrainRate2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(2,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;

    BlockLatticePhysStrainRate2D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), this->_converter );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeGeometry2D<T,DESCRIPTOR>::SuperLatticeGeometry2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice,1), _superGeometry(superGeometry),
    _material(material)
{ this->_name = "geometry"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeGeometry2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(2,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;

    BlockLatticeGeometry2D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)),
                                                        this->_superGeometry.getExtendedBlockGeometry(this->_sLattice.getLoadBalancer().loc(globIC)), _material );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeRank2D<T,DESCRIPTOR>::SuperLatticeRank2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice,1)
{ this->_name = "rank"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeRank2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
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
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeCuboid2D<T,DESCRIPTOR>::SuperLatticeCuboid2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice,1)
{ this->_name = "cuboid"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeCuboid2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0]; // int locix= input[1]; int lociy= input[2]; int lociz= input[3];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    T cuboid[1];
    cuboid[0] = globIC+1;
    std::vector<T> out(cuboid, cuboid+1);  // first adress, last adress
    return out;
  } else {
    return std::vector<T>();
  }
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysPressure2D<T,DESCRIPTOR>::SuperLatticePhysPressure2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,1)
{ this->_name = "physPressure"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysPressure2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() )
  {
    std::vector<int> inputLocal(2,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
//    inputLocal[2] = input[3] + overlap;
    
    BlockLatticePhysPressure2D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)) , this->_converter) ;

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
// // old version
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];
//  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
//    // local coords are given, fetch local cell and compute value(s)
//    int overlap = this->sLattice.getOverlap();
//    return std::vector<T>(1, this->converter.physPressureFromRho(this->sLattice
//                          .getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeRho()));
//  } else {
//    return std::vector<T>(); // empty vector
//  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysVelocity2D<T,DESCRIPTOR>::SuperLatticePhysVelocity2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2)
{ this->_name = "physVelocity"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysVelocity2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(2,T());
    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
//    inputLocal[2] = input[3] + overlap;
    BlockLatticePhysVelocity2D<T,DESCRIPTOR> blockLatticeF(this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), this->_converter);
    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysBoundaryForce2D<T,DESCRIPTOR>::SuperLatticePhysBoundaryForce2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "physBoundaryForce"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysBoundaryForce2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];

  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<T> force(2, T());
    if(this->_superGeometry.get(input) == _material) {
      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
        // Get direction
        const int* c = DESCRIPTOR<T>::c[iPop];
        // Get next cell located in the current direction
        // Check if the next cell is a fluid node
        std::vector<int> input2(input); input2[1] = input[1] + c[0]; input2[2] = input[2] + c[1];
        if (_superGeometry.get(input2) == 1) {
          int overlap = this->_sLattice.getOverlap();
          // Get f_q of next fluid cell where l = opposite(q)
          T f = this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1])[iPop];
          // Get f_l of the boundary cell
          // Add f_q and f_opp
          f += this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
          // Update force
          force[0] -= c[0]*f;
          force[1] -= c[1]*f;
        }
      }
      force[0] = this->_converter.physForce(force[0]);
      force[1] = this->_converter.physForce(force[1]);
      return force;
    }
    else {
      return force;
    }
  }
  else {
    return std::vector<T>();
  }
  //return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::SuperLatticePhysCorrBoundaryForce2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "physCorrBoundaryForce"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];

//  std::vector<T> force(3, T());
//  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
//    int globX = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosX() + locix;
//    int globY = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosY() + lociy;
//    int globZ = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosZ() + lociz;
//    if(superGeometry.getMaterial(globX, globY, globZ) == material) {
//      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
//        // Get direction
//        const int* c = DESCRIPTOR<T>::c[iPop];
//        // Get next cell located in the current direction
//        // Check if the next cell is a fluid node
//        if (superGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
//          int overlap = this->sLattice.getOverlap();
//          // Get f_q of next fluid cell where l = opposite(q)
//          T f = this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
//          // Get f_l of the boundary cell
//          // Add f_q and f_opp
//          f += this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
//          // Update force
//          force[0] -= c[0]*(f-2.*DESCRIPTOR<T>::t[iPop]);
//          force[1] -= c[1]*(f-2.*DESCRIPTOR<T>::t[iPop]);
//          force[2] -= c[2]*(f-2.*DESCRIPTOR<T>::t[iPop]);
//        }
//      }
//      force[0]=this->converter.physForce(force[0]);
//      force[1]=this->converter.physForce(force[1]);
//      force[2]=this->converter.physForce(force[2]);
//      return force;
//    }
//    else {
//      return force;
//    }
//  }
//  else {
//    return force;
//  }
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePorosity2D<T,DESCRIPTOR>::SuperLatticePorosity2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,1),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "porosity"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePorosity2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
// // new version
// // how to deal with value??

//  int globIC = input[0];
//  std::vector<int> inputLocal(3,T());
//  T* value = new T[1];
//  int overlap = this->sLattice.getOverlap();
//  inputLocal[0] = input[1] + overlap;
//  inputLocal[1] = input[2] + overlap;
//  inputLocal[2] = input[3] + overlap;
//  BlockLatticePorosity2D<T,DESCRIPTOR> blockLatticeF( this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)) );
//  std::vector<T> result(1,value[0]);
//  delete value;
//  return result;



// // old version
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];

//  T* value = new T[1];
//  int overlap = this->sLattice.getOverlap();
//  this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeExternalField(0,1,value);
//  std::vector<T> result(1,value[0]);
//  delete value;
//  return result;
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysPermeability2D<T,DESCRIPTOR>::SuperLatticePhysPermeability2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,1),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "permeability"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysPermeability2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
////  int lociz= input[3];

//  T* value = new T[1];
//  int overlap = this->sLattice.getOverlap();
//  this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap/*, lociz+overlap*/).computeExternalField(0,1,value);
//  std::vector<T> result(1,this->converter.physPermeability(value[0]));
//  delete value;
//  if (!(result[0]<42)&&!(result[0]>42)&&!(result[0]==42)) result[0] = 999999;
//  if (isinf(result[0])) result[0] = 1e6;
//  return result;
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDarcyForce2D<T,DESCRIPTOR>::SuperLatticePhysDarcyForce2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "alphaU"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysDarcyForce2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
//  SuperLatticePhysPermeability2D<T,DESCRIPTOR> permeability(this->sLattice,this->superGeometry,this->material,this->converter);
//  SuperLatticeVelocity2D<T,DESCRIPTOR> velocity(this->sLattice);

//  T nu = this->converter.getCharNu();
//  T K = permeability(input)[0];
//  std::vector<T> u = velocity(input);

//  std::vector<T> result(2,T());
//  result[0] = -nu/K*u[0];
//  result[1] = -nu/K*u[1];
////  result[2] = -nu/K*u[2];

//  return result;
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeAverage2D<T,DESCRIPTOR>::SuperLatticeAverage2D
  (SuperLatticeF2D<T,DESCRIPTOR>& f, SuperGeometry2D<T>& superGeometry,
  const int material, T radius)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(), f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material), _radius(radius)
{ this->_name = "Average("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeAverage2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
//  CuboidGeometry2D<T>& cGeometry = f.getSuperLattice2D().getCuboidGeometry();
//  LoadBalancer<T>& load = f.getSuperLattice2D().getLoadBalancer();

//  //create boolean indicator functor isInSphere
//  std::vector<T> center(3,0);
//  center[0] = (int)cGeometry.get(input[0]).get_globPosX() + input[1];
//  center[1] = (int)cGeometry.get(input[0]).get_globPosY() + input[2];
//  center[2] = (int)cGeometry.get(input[0]).get_globPosZ() + input[3];
//  SphereAnalyticalF2D<bool,T> isInSphere(center,radius);

//  // iterate over all cuboids & points and test for material && isInSphere
//  std::vector<T> tmp( this->_n, T() );
//  int numVoxels(0);
//  if (this->superGeometry.getMaterial(center[0],center[1],center[2]) == material) {
//    for (int iC=0; iC<load.size(); iC++) {
//      int nX = cGeometry.get(load.glob(iC)).getNx();
//      int nY = cGeometry.get(load.glob(iC)).getNy();
//      int nZ = cGeometry.get(load.glob(iC)).getNz();
//      for (int iX=0; iX<nX; iX++) {
//        for (int iY=0; iY<nY; iY++) {
//          for (int iZ=0; iZ<nZ; iZ++) {
//            std::vector<T> glob(3,0);
//            glob[0] = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
//            glob[1] = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
//            glob[2] = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
//            if (this->superGeometry.getMaterial(glob[0],glob[1],glob[2]) == material
//                && isInSphere(glob)[0]==true) {
//              for (unsigned iD=0; iD<f(load.glob(0),0,0,0).size(); iD++) {
//                tmp[iD]+=f(load.glob(iC),iX,iY,iZ)[iD];
//              }
//              numVoxels++;
//            }
//          }
//        }
//      }
//    }

//#ifdef PARALLEL_MODE_MPI
//    singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
//#endif
//    for (int iD=0; iD<f.getTargetDim(); iD++) {
//#ifdef PARALLEL_MODE_MPI
//      singleton::mpi().reduceAndBcast(tmp[iD], MPI_SUM);
//#endif
//      if (numVoxels>0) {
//        tmp[iD] /= numVoxels;
//      }
//    }
//  }
//  return tmp;
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperEuklidNorm2D<T,DESCRIPTOR>::SuperEuklidNorm2D(SuperLatticeF2D<T,DESCRIPTOR>& f)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),1), _f(f)
{ this->_name = "l2("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperEuklidNorm2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
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

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_LOCAL_F_3D_HH
#define BLOCK_LATTICE_LOCAL_F_3D_HH

#include<vector>
#include<cmath>

#include "functors/blockLatticeLocalF3D.h"
#include "functors/blockLatticeBaseF3D.h"
#include "functors/genericF.h"
#include "functors/analyticalF.h"
#include "functors/indicatorF.h"
#include "core/blockLattice3D.h"
#include "core/blockLatticeStructure3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "communication/mpiManager.h"

namespace olb {

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeFpop3D<T,DESCRIPTOR>::BlockLatticeFpop3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,DESCRIPTOR<T>::q)
{ this->_name = "fPop"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeFpop3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  std::vector<T> output;
  for (int iPop=0; iPop<DESCRIPTOR<T>::q; iPop++) {
    output.push_back(this->_blockLattice.get( input[0], input[1], input[2] )[iPop] + DESCRIPTOR<T>::t[iPop]);
  }
  return output;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeDissipation3D<T,DESCRIPTOR>::BlockLatticeDissipation3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,1), _converter(converter)
{ this->_name = "dissipation"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeDissipation3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1], input[2] ).computeAllMomenta(rho, uTemp, pi);

  T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6)
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

  T nuLattice = _converter.getLatticeNu();
  T omega = _converter.getOmega();
  T finalResult = PiNeqNormSqr*nuLattice*pow(omega*DESCRIPTOR<T>::invCs2,2)/rho/2.;

  return std::vector<T>(1,finalResult);
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysDissipation3D<T,DESCRIPTOR>::BlockLatticePhysDissipation3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,1), _converter(converter)
{ this->_name = "physDissipation"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysDissipation3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1], input[2] ).computeAllMomenta(rho, uTemp, pi);

  T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6)
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

  T nuLattice = _converter.getLatticeNu();
  T omega = _converter.getOmega();
  T dt = _converter.physTime();
  T finalResult = PiNeqNormSqr*nuLattice*pow(omega*DESCRIPTOR<T>::invCs2/rho,2)/2.*_converter.getCharNu()/_converter.getLatticeNu()/dt/dt;

  return std::vector<T>(1,finalResult);
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeDensity3D<T,DESCRIPTOR>::BlockLatticeDensity3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,1)
{ this->_name = "density"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeDensity3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T rho = this->_blockLattice.get( input[0] , input[1] , input[2] ).computeRho();
  std::vector<T> output( 1 , rho );
  return output;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeVelocity3D<T,DESCRIPTOR>::BlockLatticeVelocity3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,3)
{ this->_name = "velocity"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeVelocity3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
  T rho, u[3];
  this->_blockLattice.get(input[0], input[1], input[2]).computeRhoU(rho,u);
  std::vector<T> output(u, u+3); // first adress, last adress
  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeStrainRate3D<T,DESCRIPTOR>::BlockLatticeStrainRate3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,9)
{ this->_name = "strainRate"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeStrainRate3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T strainRate[9];
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1], input[2] ).computeAllMomenta(rho, uTemp, pi);

  T omega = this->_converter.getOmega();

  strainRate[0] = -pi[0]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[1] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[2] = -pi[2]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[3] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[4] = -pi[3]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[5] = -pi[4]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[6] = -pi[2]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[7] = -pi[4]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[8] = -pi[5]*omega*DESCRIPTOR<T>::invCs2/rho/2.;

  std::vector<T> output(strainRate, strainRate+9); // first adress, last adress
  return output;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysStrainRate3D<T,DESCRIPTOR>::BlockLatticePhysStrainRate3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,9)
{ this->_name = "strainRate"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysStrainRate3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T strainRate[9];
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1], input[2] ).computeAllMomenta(rho, uTemp, pi);

  T omega = this->_converter.getOmega();
  T dt = this->_converter.physTime();

  strainRate[0] = -pi[0]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[1] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[2] = -pi[2]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[3] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[4] = -pi[3]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[5] = -pi[4]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[6] = -pi[2]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[7] = -pi[4]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[8] = -pi[5]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;

  std::vector<T> output(strainRate, strainRate+9); // first adress, last adress
  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeGeometry3D<T,DESCRIPTOR>::BlockLatticeGeometry3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
  BlockGeometryStructure3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,1), _blockGeometry(blockGeometry),
    _material(material)
{ this->_name = "geometry"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeGeometry3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
  
  int materialTmp = _blockGeometry.getMaterial( input[0], input[1], input[2] );

  if (_material != -1) {
    if (_material == materialTmp) {
         return std::vector<T>(1, T(1));
    }
    else return std::vector<T>(1, T());
  }
  return std::vector<T>(1, materialTmp);
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysPressure3D<T,DESCRIPTOR>::BlockLatticePhysPressure3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,1)
{ this->_name = "physPressure"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysPressure3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T rho = this->_blockLattice.get( input[0] , input[1] , input[2] ).computeRho();
  std::vector<T> output( 1 , this->_converter.physPressureFromRho(rho) );

  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysVelocity3D<T,DESCRIPTOR>::BlockLatticePhysVelocity3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter, bool print)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3) , _print(print)
{ this->_name = "physVelocity";}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysVelocity3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T rho, u[3];
  if (_print) {
	  std::cout << input[0] << " "  << input[1] << " " << input[2] << " | " << singleton::mpi().getRank() << std::endl;
  }
  this->_blockLattice.get( input[0], input[1], input[2] ).computeRhoU(rho,u);
  u[0]=this->_converter.physVelocity(u[0]);
  u[1]=this->_converter.physVelocity(u[1]);
  u[2]=this->_converter.physVelocity(u[2]);
  std::vector<T> output(u, u+3); // first adress, last adress
  return output;
}





template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysBoundaryForce3D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material) 
{ this->_name = "physBoundaryForce"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysBoundaryForce3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];

//  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
//    int globX = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosX() + locix;
//    int globY = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosY() + lociy;
//    int globZ = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosZ() + lociz;
//    std::vector<T> force(3, T());
//    if(BlockGeometry.getMaterial(globX, globY, globZ) == material) {
//      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
//        // Get direction
//        const int* c = DESCRIPTOR<T>::c[iPop];
//        // Get next cell located in the current direction
//        // Check if the next cell is a fluid node
//        if (BlockGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
//          int overlap = this->_blockLattice.getOverlap();
//          // Get f_q of next fluid cell where l = opposite(q)
//          T f = this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
//          // Get f_l of the boundary cell
//          // Add f_q and f_opp
//          f += this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
//          // Update force
//          force[0] -= c[0]*f;
//          force[1] -= c[1]*f;
//          force[2] -= c[2]*f;
//        }
//      }
//      force[0]=this->_converter.physForce(force[0]);
//      force[1]=this->_converter.physForce(force[1]);
//      force[2]=this->_converter.physForce(force[2]);
//      return force;
//    }
//    else {
//      return force;
//    }
//  }
//  else {
//    return std::vector<T>();
//  }
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR>::BlockLatticePhysCorrBoundaryForce3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material)
{ this->_name = "physCorrBoundaryForce"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];

//  std::vector<T> force(3, T());
//  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
//    int globX = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosX() + locix;
//    int globY = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosY() + lociy;
//    int globZ = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosZ() + lociz;
//    if(BlockGeometry.getMaterial(globX, globY, globZ) == material) {
//      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
//        // Get direction
//        const int* c = DESCRIPTOR<T>::c[iPop];
//        // Get next cell located in the current direction
//        // Check if the next cell is a fluid node
//        if (BlockGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
//          int overlap = this->_blockLattice.getOverlap();
//          // Get f_q of next fluid cell where l = opposite(q)
//          T f = this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
//          // Get f_l of the boundary cell
//          // Add f_q and f_opp
//          f += this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
//          // Update force
//          force[0] -= c[0]*(f-2.*DESCRIPTOR<T>::t[iPop]);
//          force[1] -= c[1]*(f-2.*DESCRIPTOR<T>::t[iPop]);
//          force[2] -= c[2]*(f-2.*DESCRIPTOR<T>::t[iPop]);
//        }
//      }
//      force[0]=this->_converter.physForce(force[0]);
//      force[1]=this->_converter.physForce(force[1]);
//      force[2]=this->_converter.physForce(force[2]);
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
BlockLatticePorosity3D<T,DESCRIPTOR>::BlockLatticePorosity3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,1)
{ this->_name = "porosity"; }

// under construction
template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePorosity3D<T,DESCRIPTOR>::operator()(std::vector<int> input) {

  std::vector<T> output( 1 , T() );
  this->_blockLattice.get( input[0] , input[1] , input[2] ).computeExternalField(0,1,&output[0]);
  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysPermeability3D<T,DESCRIPTOR>::BlockLatticePhysPermeability3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,1),
    _blockGeometry(blockGeometry),  _material(material)
{ this->_name = "permeability"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysPermeability3D<T,DESCRIPTOR>::operator()(std::vector<int> input) {
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];

//  T* value = new T[1];
//  int overlap = this->_blockLattice.getOverlap();
//  this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeExternalField(0,1,value);
//  std::vector<T> result(1,this->_converter.physPermeability(value[0]));
//  delete value;
//  if (!(result[0]<42)&&!(result[0]>42)&&!(result[0]==42)) result[0] = 999999;
//  if (isinf(result[0])) result[0] = 1e6;
//  return result;
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysDarcyForce3D<T,DESCRIPTOR>::BlockLatticePhysDarcyForce3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material)
{ this->_name = "alphaU"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysDarcyForce3D<T,DESCRIPTOR>::operator()(std::vector<int> input) {
  BlockLatticePhysPermeability3D<T,DESCRIPTOR> permeability(this->_blockLattice,this->_blockGeometry,this->_material,this->_converter);
  BlockLatticeVelocity3D<T,DESCRIPTOR> velocity(this->_blockLattice);

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
BlockLatticeAverage3D<T,DESCRIPTOR>::BlockLatticeAverage3D
  (BlockLatticeF3D<T,DESCRIPTOR>& f, BlockGeometry3D<T>& blockGeometry,
  int material, T radius)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice3D(), f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material), _radius(radius)
{ this->_name = "Average("+f.getName()+")"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeAverage3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
//  CuboidGeometry3D<T>& cGeometry = f.getBlockLattice3D().get_cGeometry();
//  loadBalancer& load = f.getBlockLattice3D().get_load();

//  //create boolean indicator functor isInSphere
//  std::vector<T> center(3,0);
//  center[0] = (int)cGeometry.get(load.glob(input[0])).get_globPosX() + input[1];
//  center[1] = (int)cGeometry.get(load.glob(input[0])).get_globPosY() + input[2];
//  center[2] = (int)cGeometry.get(load.glob(input[0])).get_globPosZ() + input[3];
//  SphereAnalyticalF3D<bool,T> isInSphere(center,radius);

  // iterate over all cuboids & points and test for material && isInSphere
  std::vector<T> tmp( this->_n, T() );
//  int numVoxels(0);
//  if (this->blockGeometry.getMaterial(center[0],center[1],center[2]) == material) {
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
//            if (this->blockGeometry.getMaterial(glob[0],glob[1],glob[2]) == material
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
  return tmp;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockL2Norm3D<T,DESCRIPTOR>::BlockL2Norm3D(BlockLatticeF3D<T,DESCRIPTOR>& f)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice3D(),1), _f(f)
{ this->_name = "l2("+f.getName()+")"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockL2Norm3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
  std::vector<T> data = _f(input);
  std::vector<T> tmp(1,T());
  for (unsigned i=0; i<data.size(); i++) {
    tmp[0] += data[i]*data[i];
  }
  tmp[0] = sqrt(tmp[0]);
  return tmp;
}


} // end namespace olb

#endif

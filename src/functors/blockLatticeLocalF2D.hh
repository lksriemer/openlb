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

#ifndef BLOCK_LATTICE_LOCAL_F_2D_HH
#define BLOCK_LATTICE_LOCAL_F_2D_HH

#include<vector>
#include<cmath>

#include "functors/blockLatticeLocalF2D.h"
#include "functors/blockLatticeBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analyticalF.h"
#include "functors/indicatorF.h"
#include "core/blockLattice2D.h"
#include "core/blockLatticeStructure2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity

namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeDissipation2D<T,DESCRIPTOR>::BlockLatticeDissipation2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _converter(converter)
{ this->_name = "dissipation"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeDissipation2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];
//  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
//    // local coords are given, fetch local cell and compute value(s)
//    T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
//    int overlap = this->_blockLattice.getOverlap();
//    this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeAllMomenta(rho, uTemp, pi);

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
BlockLatticeDensity2D<T,DESCRIPTOR>::BlockLatticeDensity2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{ this->_name = "density"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeDensity2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T rho = this->_blockLattice.get( input[0] , input[1] ).computeRho();
  std::vector<T> output( 1 , rho );
  return output;

//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];
//  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
//    // local coords are given, fetch local cell and compute value(s)

//    int overlap = this->_blockLattice.getOverlap();
//    return std::vector<T>(1, this->_blockLattice
//                          .getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)
//                          .computeRho() );
//  } else {
//    return std::vector<T>(); // empty vector
//  }
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeVelocity2D<T,DESCRIPTOR>::BlockLatticeVelocity2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,2)
{ this->_name = "velocity"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeVelocity2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
  T rho, u[2];
  this->_blockLattice.get(input[0], input[1]).computeRhoU(rho,u);
  std::vector<T> output(u, u+2); // first adress, last adress
  return output;
}

/*template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeStrainRate2D<T,DESCRIPTOR>::BlockLatticeStrainRate2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,4)
{ this->_name = "strainRate"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeStrainRate2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T strainRate[4];
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T omega = this->_converter.getOmega();

  strainRate[0] = -pi[0]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[1] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[2] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2.;
  strainRate[3] = -pi[2]*omega*DESCRIPTOR<T>::invCs2/rho/2.;

  //cout << pi[0] << " " << pi[1] << " " << pi[2] << " " << DESCRIPTOR<T>::invCs2 << endl;

  std::vector<T> output(strainRate, strainRate+4); // first adress, last adress
  return output;
}*/

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysStrainRate2D<T,DESCRIPTOR>::BlockLatticePhysStrainRate2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,4)
{ this->_name = "strainRate"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysStrainRate2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T strainRate[4];
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T omega = this->_converter.getOmega();
  T dt = this->_converter.physTime();

  strainRate[0] = -pi[0]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[1] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[2] = -pi[1]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;
  strainRate[3] = -pi[2]*omega*DESCRIPTOR<T>::invCs2/rho/2./dt;

  std::vector<T> output(strainRate, strainRate+4); // first adress, last adress
  return output;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeGeometry2D<T,DESCRIPTOR>::BlockLatticeGeometry2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
  int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _blockGeometry(blockGeometry),
    _material(material)
{ this->_name = "geometry"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeGeometry2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  int materialTmp = _blockGeometry.getMaterial( input[0], input[1]);

  if (_material != -1) {
    if (_material == materialTmp) {
         return std::vector<T>(1, T(1));
    }
    else return std::vector<T>(1, T());
  }
  return std::vector<T>(1, materialTmp);
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeRank2D<T,DESCRIPTOR>::BlockLatticeRank2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{ this->_name = "rank"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeRank2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
//  int globIC = input[0]; // int locix= input[1]; int lociy= input[2]; int lociz= input[3];
//  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
//    T rank[1];
//    rank[0] = singleton::mpi().getRank() + 1;
//    std::vector<T> out(rank, rank+1);  // first adress, last adress
//    return out;
//  } else {
//    return std::vector<T>();
//  }
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeCuboid2D<T,DESCRIPTOR>::BlockLatticeCuboid2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{ this->_name = "cuboid"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeCuboid2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
//  int globIC = input[0]; // int locix= input[1]; int lociy= input[2]; int lociz= input[3];
//  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
//    T cuboid[1];
//    cuboid[0] = globIC+1;
//    std::vector<T> out(cuboid, cuboid+1);  // first adress, last adress
//    return out;
//  } else {
//    return std::vector<T>();
//  }
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysPressure2D<T,DESCRIPTOR>::BlockLatticePhysPressure2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1) 
{ this->_name = "physPressure"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysPressure2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T rho = this->_blockLattice.get( input[0] , input[1] ).computeRho();
  std::vector<T> output( 1 , this->_converter.physPressureFromRho(rho) );

  return output;



//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];
//  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
//    // local coords are given, fetch local cell and compute value(s)
//    int overlap = this->_blockLattice.getOverlap();
//    return std::vector<T>(1, this->_converter.physPressureFromRho(this->_blockLattice
//                          .getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeRho()));
//  } else {
//    return std::vector<T>(); // empty vector
//  }
//  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysVelocity2D<T,DESCRIPTOR>::BlockLatticePhysVelocity2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const LBconverter<T>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2)
{ this->_name = "physVelocity"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysVelocity2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T rho, u[2];
  this->_blockLattice.get( input[0], input[1] ).computeRhoU(rho,u);
  u[0]=this->_converter.physVelocity(u[0]);
  u[1]=this->_converter.physVelocity(u[1]);
  std::vector<T> output(u, u+2); // first adress, last adress
  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce2D  
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometry2D<T>& blockGeometry,
   int material, const LBconverter<T>& converter)
   : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
     _blockGeometry(blockGeometry), _material(material)
{ this->_name = "physBoundaryForce"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
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
BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::BlockLatticePhysCorrBoundaryForce2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,BlockGeometry2D<T>& blockGeometry,
   int material, const LBconverter<T>& converter)
   : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
     _blockGeometry(blockGeometry), _material(material)
{ this->_name = "physCorrBoundaryForce"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
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
BlockLatticePorosity2D<T,DESCRIPTOR>::BlockLatticePorosity2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometry2D<T>& blockGeometry,
   int material, const LBconverter<T>& converter)
   : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1),
     _blockGeometry(blockGeometry), _material(material)
{ this->_name = "porosity"; }

// under construction
template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePorosity2D<T,DESCRIPTOR>::operator()(std::vector<int> input) {
//  int globIC = input[0];
//  int locix= input[1];
//  int lociy= input[2];
//  int lociz= input[3];

//  T* value = new T[1];
//  int overlap = this->_blockLattice.getOverlap();
//  this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeExternalField(0,1,value);
//  std::vector<T> result(1,value[0]);
//  delete value;
//  return result;
  return std::vector<T>();
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysPermeability2D<T,DESCRIPTOR>::BlockLatticePhysPermeability2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometry2D<T>& blockGeometry,
   int material, const LBconverter<T>& converter)
   : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1),
     _blockGeometry(blockGeometry),  _material(material)
{ this->_name = "permeability"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysPermeability2D<T,DESCRIPTOR>::operator()(std::vector<int> input) {
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
BlockLatticePhysDarcyForce2D<T,DESCRIPTOR>::BlockLatticePhysDarcyForce2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometry2D<T>& blockGeometry,
   int material, const LBconverter<T>& converter)
   : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
     _blockGeometry(blockGeometry), _material(material)
{ this->_name = "alphaU"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysDarcyForce2D<T,DESCRIPTOR>::operator()(std::vector<int> input) {
  BlockLatticePhysPermeability2D<T,DESCRIPTOR> permeability(this->_blockLattice,this->_blockGeometry,this->_material,this->_converter);
  BlockLatticeVelocity2D<T,DESCRIPTOR> velocity(this->_blockLattice);

  T nu = this->_converter.getCharNu();
  T K = permeability(input)[0];
  std::vector<T> u = velocity(input);

  std::vector<T> result(2,0);
  result[0] = -nu/K*u[0];
  result[1] = -nu/K*u[1];

  return result;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeAverage2D<T,DESCRIPTOR>::BlockLatticeAverage2D
  (BlockLatticeF2D<T,DESCRIPTOR>& f, BlockGeometry2D<T>& blockGeometry,
   int material, T radius)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice2D(), f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material), _radius(radius)
{ this->_name = "Average("+f.getName()+")"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeAverage2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
//  CuboidGeometry2D<T>& cGeometry = f.getBlockLattice2D().get_cGeometry();
//  loadBalancer& load = f.getBlockLattice2D().get_load();

//  //create boolean indicator functor isInSphere
//  std::vector<T> center(3,0);
//  center[0] = (int)cGeometry.get(load.glob(input[0])).get_globPosX() + input[1];
//  center[1] = (int)cGeometry.get(load.glob(input[0])).get_globPosY() + input[2];
//  center[2] = (int)cGeometry.get(load.glob(input[0])).get_globPosZ() + input[3];
//  SphereAnalyticalF2D<bool,T> isInSphere(center,radius);

  // iterate over all cuboids & points and test for material && isInSphere
//  std::vector<T> tmp( this->_n, T() );
//  int numVoxels(0);
//  if (this->_blockGeometry.getMaterial(center[0],center[1],center[2]) == material) {
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
//            if (this->_blockGeometry.getMaterial(glob[0],glob[1],glob[2]) == material
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
BlockL2Norm2D<T,DESCRIPTOR>::BlockL2Norm2D(BlockLatticeF2D<T,DESCRIPTOR>& f)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice2D(),1), _f(f)
{ this->_name = "l2("+f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockL2Norm2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
  std::vector<T> data = _f(input);
  std::vector<T> tmp(1,0);
  for (unsigned i=0; i<data.size(); i++) {
    tmp[0] += data[i]*data[i];
  }
  tmp[0] = sqrt(tmp[0]);
  return tmp;
}


} // end namespace olb

#endif

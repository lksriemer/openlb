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

#ifndef BLOCK_LATTICE_INTEGRAL_F_3D_HH
#define BLOCK_LATTICE_INTEGRAL_F_3D_HH

#include<vector>
#include<cmath>

#include "functors/blockLatticeIntegralF3D.h"
#include "functors/blockLatticeLocalF3D.h"
#include "functors/blockLatticeCalcF3D.h" // for IdentityF
#include "geometry/cuboidGeometry3D.h"
#include "communication/loadBalancer.h"

namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
BlockMax3D<T,DESCRIPTOR>::BlockMax3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
  BlockGeometry3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice3D(), f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{ this->_name = "Max("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockMax3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

//  f.getBlockLattice3D().communicate();
//  CuboidGeometry3D<T>& cGeometry = f.getBlockLattice3D().get_cGeometry();
//  loadBalancer& load = f.getBlockLattice3D().get_load();

  std::vector<T> tmp(this->getTargetDim(), T() );
//  for (int i=0; i<this->getTargetDim(); i++) {

//    for (int iC=0; iC<load.size(); iC++) {
//      int nX = cGeometry.get(load.glob(iC)).getNx();
//      int nY = cGeometry.get(load.glob(iC)).getNy();
//      int nZ = cGeometry.get(load.glob(iC)).getNz();
//      for (int iX=0; iX<nX; iX++) {
//        for (int iY=0; iY<nY; iY++) {
//          for (int iZ=0; iZ<nZ; iZ++) {
//            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
//            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
//            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
//            if (blockGeometry.get_material(globX, globY, globZ) == material) {
//              if (fabs(f(load.glob(iC),iX,iY,iZ)[i]) > tmp[i]) {
//                tmp[i]=fabs(f(load.glob(iC),iX,iY,iZ)[i]);
//              }
//            }
//          }
//        }
//      }
//    }
//#ifdef PARALLEL_MODE_MPI
//    singleton::mpi().reduceAndBcast(tmp[i], MPI_MAX);
//#endif
//  }
  return tmp;

}

template <typename T, template <typename U> class DESCRIPTOR>
BlockSum3D<T,DESCRIPTOR>::BlockSum3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
  BlockGeometry3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice3D(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{ this->_name = "Sum("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockSum3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
  BlockLatticeIdentity3D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted
  std::vector<T> tmp(this->getTargetDim(), T() );
//  for (int i=0; i<this->getTargetDim(); i++) {
//    int nX = f.getBlockLattice3D().getNx();
//    int nY = f.getBlockLattice3D().getNy();
//    int nZ = f.getBlockLattice3D().getNz();
//    for (int iX=0; iX<nX; iX++) {
//      for (int iY=0; iY<nY; iY++) {
//        for (int iZ=0; iZ<nZ; iZ++) {
//          if (this->blockGeometry.get_material(iX, iY, iZ) == material) {
//            tmp[i]+=f(iX,iY,iZ)[i];
//          }
//        }
//      }
//    }
//  }
  return tmp;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockIntegral3D<T,DESCRIPTOR>::BlockIntegral3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
  BlockGeometry3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice3D(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{ this->_name = "Integral("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockIntegral3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

//  f.getBlockLattice3D().communicate();
//  CuboidGeometry3D<T>& cGeometry = f.getBlockLattice3D().get_cGeometry();
//  loadBalancer& load = f.getBlockLattice3D().get_load();

  std::vector<T> tmp(this->_n, T() );
//  for (int i=0; i<this->n; i++) {
//    for (int iC=0; iC<load.size(); iC++) {
//      int nX = cGeometry.get(load.glob(iC)).getNx();
//      int nY = cGeometry.get(load.glob(iC)).getNy();
//      int nZ = cGeometry.get(load.glob(iC)).getNz();
//      T weight = pow(this->blockGeometry.getDeltaR(),3);
//      for (int iX=0; iX<nX; iX++) {
//        for (int iY=0; iY<nY; iY++) {
//          for (int iZ=0; iZ<nZ; iZ++) {
//            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
//            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
//            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
//            if (this->blockGeometry.get_material(globX, globY, globZ) == material) {
//              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*weight;
//            }
//          }
//        }
//      }
//    }
//#ifdef PARALLEL_MODE_MPI
//    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
//#endif
//  }
  return tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockL1Norm3D<T,DESCRIPTOR>::BlockL1Norm3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
  BlockGeometry3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice3D(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{ this->_name = "L1("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockL1Norm3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  BlockLatticeIdentity3D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted

  std::vector<T> tmp(this->getTargetDim(), T() );
  for (int i=0; i<this->getTargetDim(); i++) {
    for (int iX=0; iX<_f.getBlockLattice3D().getNx(); iX++) {
      for (int iY=0; iY<_f.getBlockLattice3D().getNy(); iY++) {
        for (int iZ=0; iZ<_f.getBlockLattice3D().getNz(); iZ++) {
          if (this->_blockGeometry.getMaterial(iX, iY, iZ) == _material) {
            if (fabs(_f(iX,iY,iZ)[i]) > tmp[i]) {
              tmp[i]=fabs(_f(iX,iY,iZ)[i]);
            }
          }
        }
      }
    }
  }
  return tmp;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockL223D<T,DESCRIPTOR>::BlockL223D(BlockLatticeF3D<T,DESCRIPTOR>& f,
  BlockGeometry3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice3D(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{ this->_name = "L22("+f.getName()+")"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockL223D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

//  f.getBlockLattice3D().communicate();
//  CuboidGeometry3D<T>& cGeometry = f.getBlockLattice3D().get_cGeometry();
//  loadBalancer& load = f.getBlockLattice3D().get_load();

  std::vector<T> tmp(this->_n, T() );
//  for (int i=0; i<this->n; i++) {

//    for (int iC=0; iC<load.size(); iC++) {
//      int nX = cGeometry.get(load.glob(iC)).getNx();
//      int nY = cGeometry.get(load.glob(iC)).getNy();
//      int nZ = cGeometry.get(load.glob(iC)).getNz();
//      T weight = pow(this->blockGeometry.getDeltaR(),3);
//      for (int iX=0; iX<nX; iX++) {
//        for (int iY=0; iY<nY; iY++) {
//          for (int iZ=0; iZ<nZ; iZ++) {
//            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
//            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
//            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
//            if (this->blockGeometry.get_material(globX, globY, globZ) == material) {
//              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight;
//            }
//          }
//        }
//      }
//    }
//#ifdef PARALLEL_MODE_MPI
//    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
//#endif
//  }
  return tmp;
}


template <typename T>
BlockGeometryFaces3D<T>::BlockGeometryFaces3D(BlockGeometryStructure3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : GenericF<T,int>(7,4), _blockGeometry(blockGeometry), _material(material),
    _converter(converter)
{ this->_name = "blockGeometryFaces"; }

template <typename T>
std::vector<T> BlockGeometryFaces3D<T>::operator() (std::vector<int> input) {

  int counter[] = {0,0,0,0,0,0,0};
  if (_blockGeometry.getStatistics().getNvoxel(_material)!=0) {
    const int x0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[0];
    const int y0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[1];
    const int z0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[2];
    const int x1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[0];
    const int y1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[1];
    const int z1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[2];

    // Iterate over all cells and count the cells of the face
    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        for (int iZ = z0; iZ <= z1; iZ++) {
          // Lock at solid nodes only
          if (_blockGeometry.getMaterial(iX, iY, iZ) == _material) {
            if (_blockGeometry.getMaterial(iX-1, iY, iZ) == 1) {
              counter[0]++;
            }
            if (_blockGeometry.getMaterial(iX, iY-1, iZ) == 1) {
              counter[1]++;
            }
            if (_blockGeometry.getMaterial(iX, iY, iZ-1) == 1) {
              counter[2]++;
            }
            if (_blockGeometry.getMaterial(iX+1, iY, iZ) == 1) {
              counter[3]++;
            }
            if (_blockGeometry.getMaterial(iX, iY+1, iZ) == 1) {
              counter[4]++;
            }
            if (_blockGeometry.getMaterial(iX, iY, iZ+1) == 1) {
              counter[5]++;
            }
          }
        }
      }
    }

    T dx2 = _converter.getLatticeL()*_converter.getLatticeL();
    std::vector<T> output;
    T total = T();
    for (int i=0; i<6; i++) {
      output.push_back((T) counter[i] * dx2);
      total+=(T) counter[i] * dx2;
    }
    output.push_back(total);
    return output;
  }
  else {
    return std::vector<T>(7,T());
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysDrag3D<T,DESCRIPTOR>::BlockLatticePhysDrag3D
  (BlockLattice3D<T,DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material)
{ this->_name = "physDrag"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysDrag3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  BlockGeometryFaces3D<T> faces(_blockGeometry, _material, this->_converter);
  BlockLatticePhysBoundaryForce3D<T,DESCRIPTOR> fTemp(this->_blockLattice, _blockGeometry, _material, this->_converter);
  BlockSum3D<T,DESCRIPTOR> sumF(fTemp, _blockGeometry, _material);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());

  std::vector<T> drag(3,T());
  drag[0] = factor * sumF(input)[0] / faces(input)[0];
  drag[1] = factor * sumF(input)[1] / faces(input)[1];
  drag[2] = factor * sumF(input)[2] / faces(input)[2];
  return drag;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysCorrDrag3D<T,DESCRIPTOR>::BlockLatticePhysCorrDrag3D
  (BlockLattice3D<T,DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const LBconverter<T>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material)
{ this->_name = "physCorrDrag"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePhysCorrDrag3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  BlockGeometryFaces3D<T> faces(_blockGeometry, _material, this->_converter);

  BlockLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR> fTemp(this->_blockLattice, _blockGeometry, _material, this->_converter);
  BlockSum3D<T,DESCRIPTOR> sumF(fTemp, _blockGeometry, _material);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());

  std::vector<T> drag(3,T());
  drag[0] = factor * sumF(input)[0] / faces(input)[0];
  drag[1] = factor * sumF(input)[1] / faces(input)[1];
  drag[2] = factor * sumF(input)[2] / faces(input)[2];
  return drag;


}


} // end namespace olb

#endif

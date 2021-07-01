/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause,
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

#ifndef SUPER_LATTICE_INTEGRAL_F_2D_HH
#define SUPER_LATTICE_INTEGRAL_F_2D_HH

#include<vector>
#include<cmath>     

#include "functors/superLatticeIntegralF2D.h"
#include "functors/superLatticeCalcF2D.h" // for IdentityF
#include "functors/blockLatticeIntegralF2D.h"
#include "utilities/vectorHelpers.h"

namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
SuperMax2D<T,DESCRIPTOR>::SuperMax2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
  SuperGeometry2D<T>& superGeometry, const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),f.getTargetDim()), 
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "Max("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperMax2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity2D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted
  _f.getSuperLattice2D().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperLattice2D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice2D().getLoadBalancer();

  std::vector<T> tmp(this->getTargetDim(), T() );
  for (int i = 0; i < this->getTargetDim(); i++) {

    for (int iC = 0; iC < load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (_superGeometry.get(load.glob(iC), iX, iY) == _material) {
            if (fabs(_f(load.glob(iC),iX,iY)[i]) > tmp[i]) {
              tmp[i] = fabs(_f(load.glob(iC),iX,iY)[i]);
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(tmp[i], MPI_MAX);
#endif
  }
  return tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperSum2D<T,DESCRIPTOR>::SuperSum2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
  SuperGeometry2D<T>& superGeometry, const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "Sum("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperSum2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity2D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted
  _f.getSuperLattice2D().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperLattice2D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice2D().getLoadBalancer();

  int numVoxels(0);
  std::vector<T> tmp(this->getTargetDim(), T() );
    for (int iC=0; iC<load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          	for (int i = 0; i < this->getTargetDim(); i++) {
              tmp[i] += _f(load.glob(iC),iX,iY)[i];
           	}
            numVoxels++;
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    for (int i=0; i<this->getTargetDim(); i++) {
    	singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
    }
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  tmp.push_back(numVoxels);
  return tmp;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperIntegral2D<T,DESCRIPTOR>::SuperIntegral2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
  SuperGeometry2D<T>& superGeometry, const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "Integral("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperIntegral2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
//  SuperLatticeIdentity2D<T,DESCRIPTOR> ff(f); // exists only to prevent f from being deleted
//  f.getSuperLattice2D().communicate();
//  CuboidGeometry2D<T>& cGeometry = f.getSuperLattice2D().getCuboidGeometry();
//  LoadBalancer<T>& load = f.getSuperLattice2D().getLoadBalancer();

//  std::vector<T> tmp(this->_n, T() );
//  for (int i=0; i<this->_n; i++) {
//    for (int iC=0; iC<load.size(); iC++) {
//      int nX = cGeometry.get(load.glob(iC)).getNx();
//      int nY = cGeometry.get(load.glob(iC)).getNy();
////      int nZ = cGeometry.get(load.glob(iC)).getNz();
//      T weight = pow(this->superGeometry.getDeltaR(),3);
//      for (int iX=0; iX<nX; iX++) {
//        for (int iY=0; iY<nY; iY++) {
////          for (int iZ=0; iZ<nZ; iZ++) {
//            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
//            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
////            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
////            if (this->superGeometry.getMaterial(globX, globY) == material) {
////              tmp[i]+=f(load.glob(iC),iX,iY)[i]*weight;
////            }
////            if (this->superGeometry.getMaterial(globX, globY, globZ) == material) {
////              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*weight;
////            }

////          }
//        }
//      }
//    }
//#ifdef PARALLEL_MODE_MPI
//    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
//#endif
//  }
//  return tmp;
  return std::vector<T>();
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperL1Norm2D<T,DESCRIPTOR>::SuperL1Norm2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
  SuperGeometry2D<T>& superGeometry, const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "L1("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperL1Norm2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity2D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted
  _f.getSuperLattice2D().communicate();
  SuperEuklidNorm2D<T,DESCRIPTOR> euklidF(_f);
  CuboidGeometry2D<T>& cGeometry = _f.getSuperLattice2D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice2D().getLoadBalancer();

  int numVoxels(0);
  std::vector<T> tmp(1, T() );
  for (int iC=0; iC<load.size(); iC++) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(),3);
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          std::vector<int> inputTmp; inputTmp.push_back(load.glob(iC));
          inputTmp.push_back(iX); inputTmp.push_back(iY);
          tmp[0] += euklidF(inputTmp)[0]*weight;
          numVoxels++;
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(tmp[0], MPI_SUM);
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  tmp.push_back(numVoxels);
  return tmp;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperL2Norm2D<T,DESCRIPTOR>::SuperL2Norm2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
  SuperGeometry2D<T>& superGeometry, const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),1),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "L2Norm("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperL2Norm2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity2D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted
  _f.getSuperLattice2D().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperLattice2D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice2D().getLoadBalancer();
  SuperEuklidNorm2D<T,DESCRIPTOR> euklidF(_f);

  std::vector<T> tmp(1, T() );
  for (int i = 0; i < _f.getTargetDim(); i++) {
    for (int iC = 0; iC < load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(),3);
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
            //std::cout << _f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight << std::endl;
            tmp[0] += _f(load.glob(iC),iX,iY)[i]*_f(load.glob(iC),iX,iY)[i]*weight;
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(tmp[0], MPI_SUM);
#endif
  tmp[0] = sqrt(tmp[0]);
  return tmp;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLinfNorm2D<T,DESCRIPTOR>::SuperLinfNorm2D(SuperLatticeF2D<T,DESCRIPTOR>& f,
  SuperGeometry2D<T>& superGeometry, const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),1),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "LinfNorm("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLinfNorm2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity2D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted
  _f.getSuperLattice2D().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperLattice2D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice2D().getLoadBalancer();

  std::vector<T> tmp(4, T() );
  for (int iC = 0; iC < load.size(); iC++) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          T tmpMax = T();
          for (int i = 0; i < _f.getTargetDim(); i++) {
            tmpMax += _f(load.glob(iC),iX,iY)[i]*_f(load.glob(iC),iX,iY)[i];
          }
          if (sqrt(tmpMax) > tmp[0]) {
            tmp[0] = sqrt(tmpMax);
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(tmp[0], MPI_MAX);
#endif
  return tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperL222D<T,DESCRIPTOR>::SuperL222D(SuperLatticeF2D<T,DESCRIPTOR>& f,
  SuperGeometry2D<T>& superGeometry, const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice2D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "L22("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperL222D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
//  SuperLatticeIdentity2D<T,DESCRIPTOR> ff(f); // exists only to prevent f from being deleted
//  f.getSuperLattice2D().communicate();
//  CuboidGeometry2D<T>& cGeometry = f.getSuperLattice2D().getCuboidGeometry();
//  LoadBalancer<T>& load = f.getSuperLattice2D().getLoadBalancer();

//  std::vector<T> tmp(this->_n, T() );
//  for (int i=0; i<this->_n; i++) {

//    for (int iC=0; iC<load.size(); iC++) {
//      int nX = cGeometry.get(load.glob(iC)).getNx();
//      int nY = cGeometry.get(load.glob(iC)).getNy();
////      int nZ = cGeometry.get(load.glob(iC)).getNz();
//      T weight = pow(this->superGeometry.getDeltaR(),3);
//      for (int iX=0; iX<nX; iX++) {
//        for (int iY=0; iY<nY; iY++) {
////          for (int iZ=0; iZ<nZ; iZ++) {
//            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
//            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
////            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
////            if (this->superGeometry.getMaterial(globX, globY) == material) {
////              //std::cout << f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight << std::endl;
////              tmp[i]+=f(load.glob(iC),iX,iY)[i]*f(load.glob(iC),iX,iY)[i]*weight;
////            }
////            if (this->superGeometry.getMaterial(globX, globY, globZ) == material) {
////              //std::cout << f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight << std::endl;
////              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight;
////            }

////          }
//        }
//      }
//    }
//#ifdef PARALLEL_MODE_MPI
//    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
//#endif
//  }
//  return tmp;
  return std::vector<T>();
}



template <typename T>
SuperGeometryFaces2D<T>::SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : GenericF<T,int>(7,3), _superGeometry(superGeometry), _material(material),
    _converter(converter) 
{ this->_name = "superGeometryFaces"; }

template <typename T>
std::vector<T> SuperGeometryFaces2D<T>::operator() (std::vector<int> input) {

  _superGeometry.communicate();
  std::vector<T> counter(7,T());
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); iC++) {
    BlockGeometryFaces2D<T> f(_superGeometry.getBlockGeometry(iC), _material, _converter);
    std::vector<T> tmp = f(input);
    for (int iDim = 0; iDim < 7; iDim++) {
      counter[iDim] += tmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; iDim++) {
    singleton::mpi().reduceAndBcast( counter[iDim], MPI_SUM);
  }
#endif
  return counter;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDrag2D<T,DESCRIPTOR>::SuperLatticePhysDrag2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "physDrag"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysDrag2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{

  SuperGeometryFaces2D<T> faces(_superGeometry, _material, this->_converter);
  SuperLatticePhysBoundaryForce2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _material, this->_converter);
  SuperSum2D<T,DESCRIPTOR> sumF(f, _superGeometry, _material);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());

  std::vector<T> drag(2,T());
  drag[0] = factor * sumF(input)[0] / faces(input)[0];
  drag[1] = factor * sumF(input)[1] / faces(input)[1];

  return drag;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysCorrDrag2D<T,DESCRIPTOR>::SuperLatticePhysCorrDrag2D
  (SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "physCorrDrag"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysCorrDrag2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
//  SuperGeometryFaces2D<T> faces(superGeometry, material, this->converter);

//  SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR> f(this->sLattice, superGeometry, material, this->converter);
//  SuperSum2D<T,DESCRIPTOR> sumF(f, superGeometry, material);

//  T factor = 2. / (this->converter.getCharRho() * this->converter.getCharU() * this->converter.getCharU());

//  std::vector<T> drag(2,T());
//  drag[0] = factor * sumF(input)[0] / faces(input)[0];
//  drag[1] = factor * sumF(input)[1] / faces(input)[1];
////  drag[2] = factor * sumF(input)[2] / faces(input)[2];
//  return drag;
  return std::vector<T>();
}



template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux2D<T, DESCRIPTOR>::SuperLatticeFlux2D(SuperLatticeF2D<T, DESCRIPTOR>& f,
  SuperGeometry2D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(f.getSuperLattice2D(),4), _sg(sg), _u(2, 0.), _A(A), _n(n), _rad(radius), _h(h), _vox(0), _analyticalF(f)
{
  init(f);
  _mat.push_back(1);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux2D<T, DESCRIPTOR>::SuperLatticeFlux2D(SuperLatticeF2D<T, DESCRIPTOR>& f,
  SuperGeometry2D<T>& sg, std::vector<T>& n, std::vector<T> A, std::list<int> materials,
  T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(f.getSuperLattice2D(),4), _sg(sg), _u(2, 0.), _A(A), _n(n), _rad(radius), _h(h), _vox(0), _mat(materials), _analyticalF(f)
{
  init(f);
}


//initialization of member variables
template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux2D<T, DESCRIPTOR>::init(SuperLatticeF2D<T, DESCRIPTOR>& f) {

  if (util::nearZero(_h) ){
    _h = f.getSuperLattice2D().getCuboidGeometry().getMinDeltaR();
  }
  //normalize normal
  _n = util::normalize(_n);

  //rotation of normal 90° clockwise
  _u[0] = _n[1]; _u[1] = -_n[0];
  _u = sumVector(_h/util::norm(_u), _u, 0 , _u);

  //checking radius
  //maxPhysDist is the diameter of the geometry
  T maxPhysDist = pow(_sg.getStatistics().getMaxPhysR(1)[0]-_sg.getStatistics().getMinPhysR(1)[0],2)
                 +pow(_sg.getStatistics().getMaxPhysR(1)[1]-_sg.getStatistics().getMinPhysR(1)[1],2);
  maxPhysDist = sqrt(maxPhysDist);
  if(_rad < 0 || util::nearZero(_rad) || _rad > maxPhysDist) {
    _rad = maxPhysDist + _h;
    if (singleton::mpi().getRank() == 0) {
      std::cout << "WARNING: bad radius! Setting radius to maxPhysDist=" << _rad << std::endl;
    }
  }
}


//sums and scales three vectors component-wise
template<typename T, template<typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeFlux2D<T, DESCRIPTOR>::sumVector(T a, std::vector<T> A,
	T b, std::vector<T> u)
{
  std::vector<T> w(3, T());
  w[0] = a*A[0] + b*u[0];
  w[1] = a*A[1] + b*u[1];
  return w;
}

//check if point is inside
template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFlux2D<T, DESCRIPTOR>::checkInside(std::vector<T> physR, int iC) {
  std::vector<int> dPos(4,0);
  _sg.getCuboidGeometry().getFloorLatticeR(physR, dPos);
  int iX = dPos[1], iY = dPos[2];

  //list of material numbers of the four neighbours
  std::list<int> neighbourCellMaterial;
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY));
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY+1));
  neighbourCellMaterial.push_back(_sg.get(iC, iX+1, iY));
  neighbourCellMaterial.push_back(_sg.get(iC, iX+1, iY+1));

  //if a neighbour has none of the right material numbers of _mat it returns false
  bool interpolationPossible = false;
  std::list<int>::iterator i;
  std::list<int>::iterator j;
  for (i = neighbourCellMaterial.begin(); i != neighbourCellMaterial.end(); ++i) {
    for (j = _mat.begin(); j != _mat.end(); ++j) {
      if(*i == *j) {
        interpolationPossible = true;
      }
    }
    if(interpolationPossible!=true) {
      return false;
    }
    else interpolationPossible=false;
  }
  return true;
}

//summation(integration) of interpolated values
template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux2D<T, DESCRIPTOR>::calculate(std::vector<T>& flow, int xDir, int xStart) {
  //xDir is direction of the line
  //if positive, start at 0 (default)
  //if negative, start at -1, to prevent double summation at point 0
  if (xDir == -1) xStart = -1;

  int i = xStart;

  //while distance from origin < radius
  while(pow(i*_u[0],2)+pow(i*_u[1],2) < pow(_rad,2)) {
    std::vector<T> pos = sumVector(1, _A, i, _u);
    int iC = _sg.getCuboidGeometry().get_iC(pos[0],pos[1], 0);
    if(this->_sg.getLoadBalancer().rank(iC) == singleton::mpi().getRank()) {
      //check if point has the right material number
      if(checkInside(pos,iC)) {
        for (int k = 0; k < 2; k++) {
          if (_analyticalF.getTargetDim() == 2)
            flow[k] += _analyticalF(pos)[k];          //interpolation
          else
            flow[k] += _analyticalF(pos)[0];
        }
      _vox++;
      }
    }
    i += xDir;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeFlux2D<T, DESCRIPTOR>::operator()(std::vector<int> input) {
  this->_sLattice.communicate();

  _vox = 0;
  std::vector<T> flow(4, T()), tmpFlow(4, T());
  int iXp, iXn;
  iXp = 1; iXn = -1;

  //calculate for both direction of the line
  calculate(flow, iXp);
  calculate(flow, iXn);

  //communicate
#ifdef PARALLEL_MODE_MPI
  for (int j = 0; j < 2; j++)	{
    singleton::mpi().reduceAndBcast(flow[j], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(_vox, MPI_SUM);
#endif

  //quadrature
  for (int l = 0; l < 2; l++) {
    tmpFlow[l+2] = flow[l];	//summation of quantity
    flow[l] *= _h;          //integration
  }

  if (_analyticalF.getTargetDim() == 2)
    tmpFlow[0] = flow[0]*_n[0] + flow[1]*_n[1];	//projection of flow vector on plane normal (= flux)
  else
	tmpFlow[0] = flow[1];   //summation of quantity, if quantity is a scalar (= force)

  //number of voxel * grid length (= length)
  tmpFlow[1] = _vox*_h;

  return tmpFlow;
}


template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux2D<T, DESCRIPTOR>::print(std::string regionName,
  std::string fluxSiScaleName, std::string meanSiScaleName ) { }



template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux2D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux2D
  (SuperLattice2D<T, DESCRIPTOR>& sLattice, LBconverter<T> const& converter,
  SuperGeometry2D<T>& sg, std::vector<T>& u, std::vector<T> A, T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 4), _p(sLattice, converter),
    _fluxF(_p, sg, u, A, radius, h), clout(std::cout,"SuperLatticePhysPressureFlux2D")
{
  //this->name = "SuperLatticePhysPressureFlux2D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux2D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux2D
  (SuperLattice2D<T, DESCRIPTOR>& sLattice, LBconverter<T> const& converter,
  SuperGeometry2D<T>& sg, std::vector<T>& u, std::vector<T> A, std::list<int> materials,
  T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 4), _p(sLattice, converter),
    _fluxF(_p, sg, u, A, materials, radius, h), clout(std::cout,"SuperLatticePhysPressureFlux2D")
{
  //this->name = "SuperLatticePhysPressureFlux2D";
}


template<typename T, template<typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysPressureFlux2D<T, DESCRIPTOR>::operator()(std::vector<int> input)
{
  return _fluxF(input);
}


template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticePhysPressureFlux2D<T, DESCRIPTOR>::print(std::string regionName,
  std::string fluxSiScaleName, std::string meanSiScaleName)
{
  std::vector<int> input;
  std::vector<T> output = this->operator()(input);
  if ( regionName != "")
    clout << "regionName=" << regionName << "; regionSize[m]=" << output[1] << std::flush;
  else
    clout << "regionSize[m]=" << output[1] << std::flush;
  if (singleton::mpi().isMainProcessor() ) {
    if ( fluxSiScaleName == "MN" )
      std::cout << "; force[MN]=" << output[0]/T(1.e6) << std::flush;
    else if ( fluxSiScaleName == "kN")
      std::cout << "; force[kN]=" << output[0]/T(1.e3) << std::flush;
    else
      std::cout << "; force[N]=" << output[0] << std::flush;
    if ( meanSiScaleName == "mmHg" )
      std::cout << "; meanPressure[mmHg]=" << fabs(output[0])/output[1]/T(133.322) << std::endl;
    else
      std::cout << "; meanPressure[Pa]=" << fabs(output[0])/output[1] << std::endl;
  }
}



template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux2D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux2D
  (SuperLattice2D<T, DESCRIPTOR>& sLattice, LBconverter<T> const& converter,
  SuperGeometry2D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 4), _vel(sLattice, converter),
    _fluxF(_vel, sg, n, A, radius, h), clout(std::cout,"SuperLatticePhysVelocityFlux2D")
{
  //this->name = "SuperLatticePhysVelocityFlux2D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux2D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux2D
  (SuperLattice2D<T, DESCRIPTOR>& sLattice, LBconverter<T> const& converter,
  SuperGeometry2D<T>& sg, std::vector<T>& n, std::vector<T> A, std::list<int> materials,
  T radius, T h)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 4), _vel(sLattice, converter),
    _fluxF(_vel, sg, n, A, materials, radius, h), clout(std::cout,"SuperLatticePhysVelocityFlux2D")
{
  //this->name = "SuperLatticePhysVelocityFlux2D";
}


template<typename T, template<typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysVelocityFlux2D<T, DESCRIPTOR>::operator()(std::vector<int> input)
{
  return _fluxF(input);
}


template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticePhysVelocityFlux2D<T, DESCRIPTOR>::print(std::string regionName,
  std::string fluxSiScaleName, std::string meanSiScaleName )
{
  std::vector<int> input;
  std::vector<T> output = this->operator()(input);
  if ( regionName != "")
    clout << "regionName=" << regionName << "; regionSize[m]=" << output[1] << std::flush;
  else
    clout << "regionSize[m]=" << output[1] << std::flush;
  if (singleton::mpi().isMainProcessor() ) {
    std::cout << "; flowRate[m^2/s]=" << output[0] << std::flush;
    if ( meanSiScaleName == "mm/s" )
      std::cout << "; meanVelocity[mm/s]=" << output[0]/output[1]*T(1.e3) << std::endl;
    else
      std::cout << "; meanVelocity[m/s]=" << output[0]/output[1] << std::endl;
  }
}



#endif

} //end namespace olb 

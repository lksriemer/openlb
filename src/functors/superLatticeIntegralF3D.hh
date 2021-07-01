/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Fabian Klemens
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

#ifndef SUPER_LATTICE_INTEGRAL_F_3D_HH
#define SUPER_LATTICE_INTEGRAL_F_3D_HH

#include "functors/superLatticeIntegralF3D.h"
#include "functors/superLatticeCalcF3D.h" // for IdentityF
#include "functors/blockLatticeIntegralF3D.h"
#include "utilities/vectorHelpers.h"
#include "io/ostreamManager.h"

namespace olb {

template <typename T, template <typename U> class DESCRIPTOR>
SuperMin3D<T,DESCRIPTOR>::SuperMin3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "Min("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperMin3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();

  std::vector<T> tmp(this->getTargetDim(), std::numeric_limits<T>::max() );
  for (int i = 0; i < this->getTargetDim(); i++) {

    for (int iC=0; iC<load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          for (int iZ = 0; iZ < nZ; iZ++) {
            if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
              if (_f(load.glob(iC),iX,iY,iZ)[i] < tmp[i]) {
                tmp[i] = _f(load.glob(iC),iX,iY,iZ)[i];
              }
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(tmp[i], MPI_MIN);
#endif
  }
  return tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperMax3D<T,DESCRIPTOR>::SuperMax3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "Max("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperMax3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f); // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();

  std::vector<T> tmp(this->getTargetDim(), std::numeric_limits<T>::min() );
  for (int i = 0; i < this->getTargetDim(); i++) {

    for (int iC = 0; iC < load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          for (int iZ = 0; iZ < nZ; iZ++) {
            if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
              if (_f(load.glob(iC),iX,iY,iZ)[i] > tmp[i]) {
                tmp[i] = _f(load.glob(iC),iX,iY,iZ)[i];
              }
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
SuperSum3D<T,DESCRIPTOR>::SuperSum3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "Sum("+_f.getName()+")"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperSum3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f);    // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();

  int numVoxels(0);
  std::vector<T> tmp(this->getTargetDim(), T() );
  for (int iC = 0; iC < load.size(); iC++) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        for (int iZ = 0; iZ < nZ; iZ++) {
          if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
             for (int i = 0; i < this->getTargetDim(); i++) {
               tmp[i] += _f(load.glob(iC),iX,iY,iZ)[i];
             }
            numVoxels++;
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim(); i++) {
    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  tmp.push_back(numVoxels);
  return tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperAverage3D<T,DESCRIPTOR>::SuperAverage3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "Average("+_f.getName()+")"; }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperAverage3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f);    // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();

  int numVoxels(0);
  std::vector<T> tmp(this->getTargetDim(), T() );
  for (int iC = 0; iC < load.size(); iC++) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        for (int iZ = 0; iZ < nZ; iZ++) {
          if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
            //if (f(load.glob(iC),iX,iY,iZ)[i]!=0) std::cout<< _f(load.glob(iC),iX,iY,iZ)[i] <<std::endl;
             for (int i = 0; i < this->getTargetDim(); i++) {
               tmp[i] += _f(load.glob(iC),iX,iY,iZ)[i];
             }
            numVoxels++;
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
  for (int i = 0; i < this->getTargetDim(); i++) {
    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
    tmp[i] /= numVoxels;
  }
#endif
  tmp.push_back(numVoxels);
  return tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperIntegral3D<T,DESCRIPTOR>::SuperIntegral3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "Integral("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperIntegral3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f);    // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();

  std::vector<T> tmp(this->_n, T() );
  for (int i = 0; i < this->_n; i++) {
    for (int iC = 0; iC < load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(),3);
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          for (int iZ = 0; iZ < nZ; iZ++) {
            if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
              tmp[i] += _f(load.glob(iC),iX,iY,iZ)[i]*weight;
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
#endif
  }
  return tmp;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperL1Norm3D<T,DESCRIPTOR>::SuperL1Norm3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),2),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "L1("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperL1Norm3D<T,DESCRIPTOR>::operator() (std::vector<int> input) 
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f);    // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  SuperEuklidNorm3D<T,DESCRIPTOR> euklidF(_f);   
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();

  int numVoxels(0);
  std::vector<T> tmp(1, T() );
  for (int iC = 0; iC < load.size(); iC++) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();
    T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(),3);
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        for (int iZ = 0; iZ < nZ; iZ++) {
          if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
            std::vector<int> inputTmp; inputTmp.push_back(load.glob(iC));
            inputTmp.push_back(iX); inputTmp.push_back(iY); inputTmp.push_back(iZ);
            tmp[0] += euklidF(inputTmp)[0]*weight;
            numVoxels++;
          }
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
SuperL2Norm3D<T,DESCRIPTOR>::SuperL2Norm3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),1),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "L2Norm("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperL2Norm3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f);    // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();
  SuperEuklidNorm3D<T,DESCRIPTOR> euklidF(_f); 

  std::vector<T> tmp(1, T() );
  for (int i = 0; i < _f.getTargetDim(); i++) {
    for (int iC = 0; iC < load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(),3);
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          for (int iZ = 0; iZ < nZ; iZ++) {
            if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
              //std::cout << _f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight << std::endl;
              tmp[0] += _f(load.glob(iC),iX,iY,iZ)[i]*_f(load.glob(iC),iX,iY,iZ)[i]*weight;
            }
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
SuperLpNorm3D<T,DESCRIPTOR>::SuperLpNorm3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material, int p)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),2),
    _f(f), _superGeometry(superGeometry), _material(material), _p(p)
{ this->_name = "LpNorm("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLpNorm3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f);    // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();
  SuperEuklidNorm3D<T,DESCRIPTOR> euklidF(_f); 

  std::vector<T> tmp(2, T() );
  for (int i = 0; i < _f.getTargetDim(); i++) {
    for (int iC = 0; iC < load.size(); iC++) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      int nZ = cGeometry.get(load.glob(iC)).getNz();
      T weight = pow(cGeometry.get(load.glob(iC)).getDeltaR(),3);
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          for (int iZ = 0; iZ < nZ; iZ++) {
            if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
              tmp[0] += pow(_f(load.glob(iC),iX,iY,iZ)[i],_p)*weight;
            }
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(tmp[0], MPI_SUM);
#endif
  tmp[1] = tmp[0];
  tmp[0] = pow(tmp[0],1./_p);
  return tmp;
}
template <typename T, template <typename U> class DESCRIPTOR>
SuperLinfNorm3D<T,DESCRIPTOR>::SuperLinfNorm3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
  SuperGeometry3D<T>& superGeometry, const int material)
  : SuperLatticeF3D<T,DESCRIPTOR>(f.getSuperLattice3D(),1),
    _f(f), _superGeometry(superGeometry), _material(material)
{ this->_name = "LinfNorm("+_f.getName()+")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLinfNorm3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperLatticeIdentity3D<T,DESCRIPTOR> ff(_f);    // exists only to prevent f from being deleted
  _f.getSuperLattice3D().communicate();
  CuboidGeometry3D<T>& cGeometry = _f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperLattice3D().getLoadBalancer();

  std::vector<T> tmp(4, T() );
  for (int iC = 0; iC < load.size(); iC++) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        for (int iZ = 0; iZ < nZ; iZ++) {
          if (this->_superGeometry.get(load.glob(iC), iX, iY, iZ) == _material) {
            T tmpMax = T();
            for (int i = 0; i < _f.getTargetDim(); i++) {              
              tmpMax += _f(load.glob(iC),iX,iY,iZ)[i]*_f(load.glob(iC),iX,iY,iZ)[i];
            }
            if (sqrt(tmpMax) > tmp[0]) {
              tmp[0] = sqrt(tmpMax); 
            }
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

template <typename T>
SuperGeometryFaces3D<T>::SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : GenericF<T,int>(7,4), _superGeometry(superGeometry), _material(material),
    _converter(converter) 
{ this->_name = "superGeometryFaces"; }

template <typename T>
std::vector<T> SuperGeometryFaces3D<T>::operator() (std::vector<int> input)
{
  _superGeometry.communicate();
  std::vector<T> counter(7,T());
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); iC++) {
    BlockGeometryFaces3D<T> f(_superGeometry.getBlockGeometry(iC), _material, _converter);
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
SuperLatticePhysDrag3D<T,DESCRIPTOR>::SuperLatticePhysDrag3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,3),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "physDrag"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysDrag3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperGeometryFaces3D<T> faces(_superGeometry, _material, this->_converter);
  SuperLatticePhysBoundaryForce3D<T,DESCRIPTOR> pBoundForce(this->_sLattice, _superGeometry, _material, this->_converter);
  SuperSum3D<T,DESCRIPTOR> sumF(pBoundForce, _superGeometry, _material);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());

  std::vector<T> drag(3,T());
  drag[0] = factor * sumF(input)[0] / faces(input)[0];
  drag[1] = factor * sumF(input)[1] / faces(input)[1];
  drag[2] = factor * sumF(input)[2] / faces(input)[2];
  return drag;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysCorrDrag3D<T,DESCRIPTOR>::SuperLatticePhysCorrDrag3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const LBconverter<T>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,3),
    _superGeometry(superGeometry), _material(material)
{ this->_name = "physCorrDrag"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysCorrDrag3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  SuperGeometryFaces3D<T> faces(_superGeometry, _material, this->_converter);
  SuperLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR>  pBoundForce(this->_sLattice, _superGeometry,
                                                                 _material, this->_converter);
  SuperSum3D<T,DESCRIPTOR> sumF(pBoundForce, _superGeometry, _material);

  T factor = 2. / (this->_converter.getCharRho() * this->_converter.getCharU() * this->_converter.getCharU());

  std::vector<T> drag(3,T());
  drag[0] = factor * sumF(input)[0] / faces(input)[0];
  drag[1] = factor * sumF(input)[1] / faces(input)[1];
  drag[2] = factor * sumF(input)[2] / faces(input)[2];
  return drag;
}




template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v, std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice3D(),5), _sg(sg), _u(u), _v(v),
    _A(A), _rad(radius), _h(h), _analyticalF(f, false)
{
  _mat.push_back(1);
  // n perpendicular to u and v
  _n = sumVector( T(1)/util::norm(util::crossProduct3D(_u, _v)), util::crossProduct3D(_u, _v),
                  T(), _A, T(), _A );
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
  std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice3D(),5), _sg(sg),  _u(u), _v(v),
    _A(A), _rad(radius), _h(h), _mat(materials), _analyticalF(f, false)
{
  _n = sumVector( T(1)/util::norm(util::crossProduct3D(_u, _v)), util::crossProduct3D(_u, _v),
                  T(), _A, T(), _A);
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice3D(),5), _sg(sg), _u(3, 0.),
    _v(3, 0.), _A(A), _n(n), _rad(radius), _h(h), _analyticalF(f, false)
{
  _mat.push_back(1);
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A, std::list<int> materials,
  T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice3D(),5), _sg(sg), _u(3, 0.),
    _v(3, 0.), _A(A), _n(n), _rad(radius), _h(h), _mat(materials), _analyticalF(f, false)
{
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<bool,T>& circle, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice3D(),5), _sg(sg), _u(3, 0.),
    _v(3, 0.), _A(circle.getCenter()), _n(circle.getNormal()), _rad(circle.getRadius()),
    _h(h), _vox(0), _analyticalF(f, false)
{
  _mat.push_back(1);
  init(f);
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(SuperLatticeF3D<T, DESCRIPTOR>& f,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<bool,T>& circle, std::list<int> materials, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice3D(),5), _sg(sg), _u(3, 0.),
    _v(3, 0.), _A(circle.getCenter()), _n(circle.getNormal()), _rad(circle.getRadius()),
    _h(h), _vox(0), _mat(materials), _analyticalF(f, false)
{
  init(f);
}


//initialization of member variables
template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux3D<T, DESCRIPTOR>::init(SuperLatticeF3D<T, DESCRIPTOR>& f)
{
  this->_name = "SuperLatticeFlux3D";

  //set grid length h to lattice length
  if (util::nearZero(_h)) {
    _h = f.getSuperLattice3D().getCuboidGeometry().getMinDeltaR();
  }

  //define vector u in plane (perpendicular to n)
  if (util::nearZero(_n[2]) ) {
    _u[2] = T(1);
  } else {
    _u[0] = _n[2];
    _u[2] = -_n[0];
  }

  //define vector v in plane perpendicular to n and u
  _v = util::crossProduct3D(_n, _u);

  _n = util::normalize(_n); // normalize n
  _u = sumVector(_h/util::norm(_u), _u, T(), _v, T(), _A); // normalize u and set length h
  _v = sumVector(_h/util::norm(_v), _v, T(), _A, T(), _A); // normalize v and set length h

  //checking radius
  //maxPhysDist is the diameter of the geometry
  T maxPhysDist = pow(_sg.getStatistics().getMaxPhysR(1)[0]-_sg.getStatistics().getMinPhysR(1)[0],2)
       +pow(_sg.getStatistics().getMaxPhysR(1)[1] - _sg.getStatistics().getMinPhysR(1)[1],2)
       +pow(_sg.getStatistics().getMaxPhysR(1)[2] - _sg.getStatistics().getMinPhysR(1)[2],2);
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
std::vector<T> SuperLatticeFlux3D<T, DESCRIPTOR>::sumVector(T a, std::vector<T> A,
    T b, std::vector<T> u, T c, std::vector<T> v) {
  std::vector<T> w(3, T());
  w[0] = a*A[0] + b*u[0] + c*v[0];
  w[1] = a*A[1] + b*u[1] + c*v[1];
  w[2] = a*A[2] + b*u[2] + c*v[2];
  return w;
}

//check if point is inside
template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFlux3D<T, DESCRIPTOR>::checkInside(std::vector<T> physR, int iC)
{
  std::vector<int> dPos(4,0);
  //get nearest lattice point
  _sg.getCuboidGeometry().getFloorLatticeR(physR, dPos);
  int iX = dPos[1], iY = dPos[2], iZ = dPos[3];

  //list of material numbers of the eight neighbours of the lattice point
  std::list<int> neighbourCellMaterial;
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY, iZ));
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY, iZ+1));
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY+1, iZ));
  neighbourCellMaterial.push_back(_sg.get(iC, iX, iY+1, iZ+1));
  neighbourCellMaterial.push_back(_sg.get(iC, iX+1, iY, iZ));
  neighbourCellMaterial.push_back(_sg.get(iC, iX+1, iY, iZ+1));
  neighbourCellMaterial.push_back(_sg.get(iC, iX+1, iY+1, iZ));
  neighbourCellMaterial.push_back(_sg.get(iC, iX+1, iY+1, iZ+1));

  //if a neighbour has none of the right material numbers it returns false
  bool interpolationPossible=false;
  std::list<int>::iterator i;
  std::list<int>::iterator j;
  for (i = neighbourCellMaterial.begin(); i != neighbourCellMaterial.end(); ++i) {
    for (j = _mat.begin(); j != _mat.end(); ++j) {
      if(*i == *j) {
        interpolationPossible = true;
      }
    }
    if(interpolationPossible != true) {
      return false;
    }
    else interpolationPossible = false;
  }
  return true;
}

//summation(integration) of interpolated values
template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux3D<T, DESCRIPTOR>::calculate(std::vector<T>& flow, int xDir,
  int yDir, int xStart, int yStart)
{
  //x and y direction of the plane
  //if positive, start at 0 (default)
  //if negative, start at -1, to prevent double summation at point 0
  if (xDir == -1) xStart = -1;
  if (yDir == -1) yStart = -1;

  //i,j are the affine coordinates of the plane
  int i = xStart;
  int j = yStart;

  //while distance from origin < radius
  while(pow(i*_u[0]+j*_v[0],2)+pow(i*_u[1]+j*_v[1],2)+pow(i*_u[2]+j*_v[2],2) < pow(_rad,2)) {
    while(pow(i*_u[0]+j*_v[0],2)+pow(i*_u[1]+j*_v[1],2)+pow(i*_u[2]+j*_v[2],2) < pow(_rad,2)) {
      std::vector<T> pos = sumVector(1, _A, i, _u, j, _v);
      int iC = _sg.getCuboidGeometry().get_iC(pos[0],pos[1],pos[2], 0);
      if(iC != _sg.getCuboidGeometry().getNc() ) {
        if(this->_sg.getLoadBalancer().rank(iC) == singleton::mpi().getRank()) {
          //check if point has the right material number
          if(checkInside(pos,iC)) {
            for (int k = 0; k < 3; k++) {
              if (_analyticalF.getTargetDim() == 3)		//interpolation
                flow[k] += _analyticalF(pos)[k];		//if quantity is three dimensional
              else
                flow[k] += _analyticalF(pos)[0];		//if quantity is one dimensional
            }
            _vox++;
          }
        }
      }
      j += yDir;
    }
    i += xDir;
    j = yStart;
  }
}

//returns values
template<typename T, template<typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeFlux3D<T, DESCRIPTOR>::operator()(std::vector<int> input)
{
  this->_sLattice.communicate();

  _vox = 0;
  std::vector<T> flow(3,T()), tmpFlow(5,T());
  int iXp, iXn, iYp, iYn;
  iXp=iYp=1; iXn=iYn=-1;

  //calculate for each of the four quadrants in a plane
  calculate(flow, iXp, iYp);
  calculate(flow, iXp, iYn);
  calculate(flow, iXn, iYp);
  calculate(flow, iXn, iYn);

  //communicate
#ifdef PARALLEL_MODE_MPI
  for (int j = 0; j < 3; j++) {
    singleton::mpi().reduceAndBcast(flow[j], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(_vox, MPI_SUM);
#endif

  for (int l = 0; l < 3; l++) {
    tmpFlow[l+2] = flow[l];     //summation of quantity
    flow[l] *= _h*_h;           //integration
  }

  if (_analyticalF.getTargetDim() == 3)
    tmpFlow[0] = util::dotProduct3D(flow, _n);  //projection of flow vector on plane normal (= flux)
  else
    tmpFlow[0] = flow[2];  //summation of quantity, if quantity is one dimensional (= force)

  //number of voxel * grid length^2 (= area)
  tmpFlow[1] = _vox*_h*_h;

  return tmpFlow;
}

template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticeFlux3D<T, DESCRIPTOR>::print(std::string regionName,
  std::string fluxSiScaleName, std::string meanSiScaleName ) {
}


template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
  T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _p(sLattice, converter),
    _fluxF(_p, sg, u, v, A, radius, h), clout(std::cout,"SuperLatticePhysPressureFlux3D")
{
  this->_name = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
  std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _p(sLattice, converter),
    _fluxF(_p, sg, u, v, A, materials, radius, h), clout(std::cout,"SuperLatticePhysPressureFlux3D")
{
  this->_name = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _p(sLattice, converter),
    _fluxF(_p, sg, n, A, radius, h), clout(std::cout,"SuperLatticePhysPressureFlux3D")
{
  this->_name = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A,
  std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _p(sLattice, converter),
    _fluxF(_p, sg, n, A, materials, radius, h), clout(std::cout,"SuperLatticePhysPressureFlux3D")
{
  this->_name = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<bool,T>& circle, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _p(sLattice, converter),
    _fluxF(_p, sg, circle, h), clout(std::cout,"SuperLatticePhysPressureFlux3D")
{
 this->_name = "SuperLatticePhysPressureFlux3D";
}
template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::SuperLatticePhysPressureFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<bool,T>& circle, std::list<int> materials, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _p(sLattice, converter),
    _fluxF(_p, sg, circle, materials, h), clout(std::cout,"SuperLatticePhysPressureFlux3D")
{
  this->_name = "SuperLatticePhysPressureFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::operator()(std::vector<int> input) {
  return _fluxF(input);
}


template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticePhysPressureFlux3D<T, DESCRIPTOR>::print(std::string regionName,
  std::string fluxSiScaleName, std::string meanSiScaleName )
{
  std::vector<int> input; 
  std::vector<T> output = this->operator()(input);
  if ( regionName != "")
    clout << "regionName=" << regionName << "; regionSize[m^2]=" << output[1] << std::flush;
  else 
    clout << "regionSize[m^2]=" << output[1] << std::flush;
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
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
  T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _vel(sLattice, converter),
    _fluxF(_vel, sg, u, v, A, radius, h), clout(std::cout,"SuperLatticePhysVelocityFlux3D")
{
  this->_name = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& u, std::vector<T>& v, std::vector<T> A,
  std::list<int> materials, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _vel(sLattice, converter),
    _fluxF(_vel, sg, u, v, A, materials, radius, h), clout(std::cout,"SuperLatticePhysVelocityFlux3D")
{
  this->_name = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A, T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _vel(sLattice, converter),
    _fluxF(_vel, sg, n, A, radius, h), clout(std::cout,"SuperLatticePhysVelocityFlux3D")
{
  this->_name = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, std::vector<T>& n, std::vector<T> A, std::list<int> materials,
  T radius, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _vel(sLattice, converter),
    _fluxF(_vel, sg, n, A, materials, radius, h), clout(std::cout,"SuperLatticePhysVelocityFlux3D")
{
  this->_name = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<bool,T>& circle, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _vel(sLattice, converter),
    _fluxF(_vel, sg, circle, h), clout(std::cout,"SuperLatticePhysVelocityFlux3D")
{
  this->_name = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::SuperLatticePhysVelocityFlux3D
  (SuperLattice3D<T, DESCRIPTOR>& sLattice, LBconverter<T>& converter,
  SuperGeometry3D<T>& sg, IndicatorCircle3D<bool,T>& circle, std::list<int> materials, T h)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 5), _vel(sLattice, converter),
    _fluxF(_vel, sg, circle, materials, h), clout(std::cout,"SuperLatticePhysVelocityFlux3D")
{
  this->_name = "SuperLatticePhysVelocityFlux3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::operator()(std::vector<int> input) {
  return _fluxF(input);
}


template<typename T, template<typename U> class DESCRIPTOR>
void SuperLatticePhysVelocityFlux3D<T, DESCRIPTOR>::print(std::string regionName,
  std::string fluxSiScaleName, std::string meanSiScaleName )
{
  std::vector<int> input; 
  std::vector<T> output = this->operator()(input);
  if ( regionName != "")
    clout << "regionName=" << regionName << "; regionSize[m^2]=" << output[1] << std::flush;
  else 
    clout << "regionSize[m^2]=" << output[1] << std::flush;
  if (singleton::mpi().isMainProcessor() ) {
    if ( fluxSiScaleName == "ml/s" )
      std::cout << "; volumetricFlowRate[ml/s]=" << output[0]*T(1.e6) << std::flush;
    else if ( fluxSiScaleName == "l/s") 
      std::cout << "; volumetricFlowRate[l/s]=" << output[0]*T(1.e3) << std::flush;
    else
      std::cout << "; volumetricFlowRate[m^3/s]=" << output[0] << std::flush;
    if ( meanSiScaleName == "mm/s" )
      std::cout << "; meanVelocity[mm/s]=" << output[0]/output[1]*T(1.e3) << std::endl;  
    else
      std::cout << "; meanVelocity[m/s]=" << output[0]/output[1] << std::endl;    
  }
}

} // end namespace olb

#endif

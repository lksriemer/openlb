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

#ifndef INTERPOLATION_F_3D_HH
#define INTERPOLATION_F_3D_HH


#include "functors/interpolationF3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity

namespace olb {

/// a class used to convert a block functor to an analytical functor
template<typename T>
AnalyticalFfromBlockF3D<T>::AnalyticalFfromBlockF3D(BlockF3D<T>& f, Cuboid3D<T>& cuboid)
  : AnalyticalF3D<T, T>(f.getTargetDim()), _f(f), _cuboid(cuboid)
{
  this->getName() = "fromBlockF";
}

template<typename T>
bool AnalyticalFfromBlockF3D<T>::operator()(
  T output[], const T physC[])
{
  int latticeR[3];
  _cuboid.getFloorLatticeR(latticeR, physC);


  T physRiC[3], d[3], e[3], output_tmp[3];

  _cuboid.getPhysR(physRiC, latticeR[0], latticeR[1], latticeR[1]);
  T dr = 1/_cuboid.getDeltaR();

  // compute weights
  d[0] = (physC[0] - physRiC[0])*dr;
  d[1] = (physC[1] - physRiC[1])*dr;
  d[2] = (physC[2] - physRiC[2])*dr;

  e[0] = 1. - d[0];
  e[1] = 1. - d[1];
  e[2] = 1. - d[2];


  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = T();
    output_tmp[iD] = T();
  }

  //if (d[0]*d[0]>1||d[1]*d[1]>1||d[2]*d[2]>1)
  //cout << d[0] << " " << d[1] << " " << d[2] << endl;
  //if (locX<0||locY<0||locZ<0)
  //cout << locX << " " << locY << " " << locZ << "_overlap=" <<_overlap<< endl;
  //      std::vector<T> output(f(latticeC).size(), T());
  //std::vector<T> output_tmp(_f.getTargetDim(), T());

  //0=1=2=
  _f(output_tmp,latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*e[1]*e[2];
  }

  latticeR[1]++;
  //0=1+2=
  _f(output_tmp,latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*d[1]*e[2];
  }

  latticeR[0]++;
  latticeR[1]--;
  //0+1=2=
  _f(output_tmp,latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*e[1]*e[2];
  }

  latticeR[1]++;
  //0+1+2=
  _f(output_tmp,latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*d[1]*e[2];
  }

  latticeR[0]--;
  latticeR[1]--;
  latticeR[2]++;
  //0=1=2+
  _f(output_tmp,latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*e[1]*d[2];
  }

  latticeR[1]++;
  //0=1+2+
  _f(output_tmp,latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*d[1]*d[2];
  }

  latticeR[0]++;
  latticeR[1]--;
  //0+1=2+
  _f(output_tmp,latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*e[1]*d[2];
  }

  latticeR[1]++;
  //0+1+2+
  _f(output_tmp,latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*d[1]*d[2];
  }

  return true;
}


/// a class used to convert lattice functions to analytical functions
template<typename T>
AnalyticalFfromSuperF3D<T>::AnalyticalFfromSuperF3D(
  SuperF3D<T>& f, bool communicateToAll, int overlap)
  : AnalyticalF3D<T, T>(f.getTargetDim()), _f(f),
    _cuboidGeometry(f.getSuperStructure().getCuboidGeometry()),
    _communicateToAll(communicateToAll), _overlap(overlap)
{
  this->getName() = "fromSuperF";
  //std::cout << _f.getName() << std::endl;
  if (overlap == -1) {
    _overlap = f.getSuperStructure().getOverlap();
  }
  /*if (&(_f.getBlockF(0))!=NULL ) {
     for (int iC = 0; iC < _f.getSuperStructure().getLoadBalancer().size(); iC++ ) {
       int iCglob = _f.getSuperStructure().getLoadBalancer().glob(iC);
       this->_analyticalFfromBlockF.push_back(new AnalyticalFfromBlockF3D<T>(_f.getBlockF(iC), _f.getSuperStructure().getCuboidGeometry().get(iCglob) ) );
     }
   }*/
}


template<typename T>
AnalyticalFfromSuperF3D<T>::~AnalyticalFfromSuperF3D()
{
  if (_analyticalFfromBlockF.size() != 0) {
    for (unsigned int iC = 0; iC < _analyticalFfromBlockF.size(); ++iC) {
      delete _analyticalFfromBlockF[iC];
    }
    _analyticalFfromBlockF.resize(0);
  }
}

template<typename T>
bool AnalyticalFfromSuperF3D<T>::operator()(
  T output[], const T physC[])
{
  int latticeR[4];
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = T();
  }
  if (!(_cuboidGeometry.getLatticeR(latticeR, physC))) {
    return false;
  }
  // convert to lattice coordinates
  T d[3];

  _f.getSuperStructure().communicate();

  int locX, locY, locZ;

  int dataSize = 0;
  int dataFound = 0;

  int latticeC[4] = {};

  for (int iC = 0; iC < _f.getSuperStructure().getLoadBalancer().size(); ++iC) {
    latticeC[0] = _f.getSuperStructure().getLoadBalancer().glob(iC);
    _cuboidGeometry.get(latticeC[0]).getFloorLatticeR(latticeR, physC);
    if (latticeR[0] >= -_overlap && latticeR[0] + 1 < _cuboidGeometry.get(latticeC[0]).getNx() + _overlap &&
        latticeR[1] >= -_overlap && latticeR[1] + 1 < _cuboidGeometry.get(latticeC[0]).getNy() + _overlap &&
        latticeR[2] >= -_overlap && latticeR[2] + 1 < _cuboidGeometry.get(latticeC[0]).getNz() + _overlap ) {
      if (_analyticalFfromBlockF.size() != 0 ) {
        (*(_analyticalFfromBlockF[iC]))(output, physC);
      } else {
        locX = latticeR[0];
        locY = latticeR[1];
        locZ = latticeR[2];

        T physRiC[3];
        _cuboidGeometry.get(latticeC[0]).getPhysR(physRiC, locX, locY, locZ);

        d[0] = (physC[0] - physRiC[0]);
        d[1] = (physC[1] - physRiC[1]);
        d[2] = (physC[2] - physRiC[2]);

        d[0] /= _cuboidGeometry.get(latticeC[0]).getDeltaR();
        d[1] /= _cuboidGeometry.get(latticeC[0]).getDeltaR();
        d[2] /= _cuboidGeometry.get(latticeC[0]).getDeltaR();

        //if (d[0]*d[0]>1||d[1]*d[1]>1||d[2]*d[2]>1)
        //cout << d[0] << " " << d[1] << " " << d[2] << endl;
        //if (locX<0||locY<0||locZ<0)
        //cout << locX << " " << locY << " " << locZ << "_overlap=" <<_overlap<< endl;

        //      std::vector<T> output(f(latticeC).size(), T());
        //std::vector<T> output_tmp(_f.getTargetDim(), T());
        T output_tmp[_f.getTargetDim()];
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] = T();
          output_tmp[iD] = T();
        }

        latticeC[1] = locX;
        latticeC[2] = locY;
        latticeC[3] = locZ;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (1 - d[0]) * (1 - d[1]) * (1 - d[2]));
          output_tmp[iD] = T();
        }


        latticeC[1] = locX;
        latticeC[2] = locY + 1;
        latticeC[3] = locZ;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (1 - d[0]) * (d[1]) * (1 - d[2]));
          output_tmp[iD] = T();
        }

        latticeC[1] = locX + 1;
        latticeC[2] = locY;
        latticeC[3] = locZ;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (d[0]) * (1 - d[1]) * (1 - d[2]));
          output_tmp[iD] = T();
        }

        latticeC[1] = locX + 1;
        latticeC[2] = locY + 1;
        latticeC[3] = locZ;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (d[0]) * (d[1]) * (1 - d[2]));
          output_tmp[iD] = T();
        }

        latticeC[1] = locX;
        latticeC[2] = locY;
        latticeC[3] = locZ + 1;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (1 - d[0]) * (1 - d[1]) * (d[2]));
          output_tmp[iD] = T();
        }

        latticeC[1] = locX;
        latticeC[2] = locY + 1;
        latticeC[3] = locZ + 1;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (1 - d[0]) * (d[1]) * (d[2]));
          output_tmp[iD] = T();
        }

        latticeC[1] = locX + 1;
        latticeC[2] = locY;
        latticeC[3] = locZ + 1;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (d[0]) * (1 - d[1]) * (d[2]));
          output_tmp[iD] = T();
        }

        latticeC[1] = locX + 1;
        latticeC[2] = locY + 1;
        latticeC[3] = locZ + 1;
        _f(output_tmp,latticeC);
        for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
          output[iD] += (output_tmp[iD] * (d[0]) * (d[1]) * (d[2]));
          output_tmp[iD] = T();
        }
      }
      dataSize += _f.getTargetDim();
      dataFound ++;
      //      for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      //        output[iD] = output_tmp[iD];
      //      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  if (_communicateToAll) {
    singleton::mpi().reduceAndBcast(dataFound, MPI_SUM);
    singleton::mpi().reduceAndBcast(dataSize, MPI_SUM);
    dataSize /= dataFound;
    for (int iD = 0; iD < dataSize; ++iD) {
      singleton::mpi().reduceAndBcast(output[iD], MPI_SUM);
    }
    for (int iD = 0; iD < dataSize; ++iD) {
      output[iD]/=dataFound;
    }
  }
#endif
  if (dataFound>0) {
    return true;
  }
  return false;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::SuperLatticeFfromAnalyticalF3D(
  AnalyticalF3D<T, T>& f, SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, f.getTargetDim()),
    _f(f)
{
  this->getName() = "fromAnalyticalF";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  this->_sLattice.getCuboidGeometry().getPhysR(physR,input);
  _f(output,physR);
  return true;
}

//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::BlockLatticeFfromAnalyticalF3D(
  AnalyticalF3D<T, T>& f, BlockLattice3D<T, DESCRIPTOR>& sLattice,
  BlockGeometry3D<T>& superGeometry, CuboidGeometry3D<T>& cuboidGeometry)
  : BlockLatticeF3D<T, DESCRIPTOR>(sLattice, f.getTargetDim()),
    _f(f),
    _superGeometry(superGeometry),
    _cuboidGeometry(cuboidGeometry)
{
  this->getName() = "fromAnalyticalF";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeFfromAnalyticalF3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[3] = {};
  _superGeometry.getPhysR(physR,input[0],input[1],input[2] );
  _f(output,physR);
  return true;
}

}  // end namespace olb

#endif

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

#include<vector>    // for generic i/o

#include "functors/interpolationF3D.h"
#include "functors/genericF.h"
#include "functors/analyticalF.h"
#include "functors/indicatorF.h"
#include "core/superLattice3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity


namespace olb {


/// a class used to convert lattice functions to analytical functions
template <typename T, template <typename U> class DESCRIPTOR>
AnalyticalFfromSuperLatticeF3D<T,DESCRIPTOR>::AnalyticalFfromSuperLatticeF3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, bool communicateToAll, int overlap )
  : AnalyticalF3D<T,T>(f.getTargetDim()), _f(f),
    _cuboidGeometry(f.getSuperLattice3D().getCuboidGeometry()),
    _communicateToAll(communicateToAll), _overlap(overlap)
{
  this->_name = "fromSuperLatticeF";
}


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> AnalyticalFfromSuperLatticeF3D<T,DESCRIPTOR>::operator()
(std::vector<T> physC)
{ 
  // convert to lattice coordinates
  T d[3];

  if(_communicateToAll) {
	  _f.getSuperLattice3D().communicate();
  }

  int locX, locY, locZ;
#ifdef PARALLEL_MODE_MPI
  int dataSize = 0;
  int dataFound = 0;
#endif
  std::vector<T> pOutput(_f.getTargetDim(),T());
  std::vector<int> latticeC(4,0);

  for (int iC = 0; iC < _f.getSuperLattice3D().getLoadBalancer().size(); iC++) {
	  latticeC[0] = _f.getSuperLattice3D().getLoadBalancer().glob(iC);
    if (_cuboidGeometry.get(latticeC[0]).checkPoint(
          physC[0], physC[1], physC[2], locX, locY, locZ, _overlap-1)) {
      locX-=(_overlap-1);
      locY-=(_overlap-1);
      locZ-=(_overlap-1);

      std::vector<T> physRiC = _cuboidGeometry.get(latticeC[0]).getPhysR(locX,locY,locZ);

      d[0] = (physC[0] - physRiC[0]);
      d[1] = (physC[1] - physRiC[1]);
      d[2] = (physC[2] - physRiC[2]);

      d[0]/=_cuboidGeometry.get(latticeC[0]).getDeltaR();
      d[1]/=_cuboidGeometry.get(latticeC[0]).getDeltaR();
      d[2]/=_cuboidGeometry.get(latticeC[0]).getDeltaR();

//      cout << d[0] << " " << d[1] << " " << d[2] << endl;

//      std::vector<T> output(f(latticeC).size(), T());
      std::vector<T> output(_f.getTargetDim(), T());

      for (unsigned int iD=0; iD<output.size(); iD++) {
        latticeC[1] = locX;
        latticeC[2] = locY;
        latticeC[3] = locZ;
        output[iD] += (_f(latticeC)[iD]*(1-d[0])*(1-d[1])*(1-d[2]));

        latticeC[1] = locX;
        latticeC[2] = locY+1;
        latticeC[3] = locZ;
        output[iD] += (_f(latticeC)[iD]*(1-d[0])*(d[1])*(1-d[2]));

        latticeC[1] = locX+1;
        latticeC[2] = locY;
        latticeC[3] = locZ;
        output[iD] += (_f(latticeC)[iD]*(d[0])*(1-d[1])*(1-d[2]));

        latticeC[1] = locX+1;
        latticeC[2] = locY+1;
        latticeC[3] = locZ;
        output[iD] += (_f(latticeC)[iD]*(d[0])*(d[1])*(1-d[2]));

        latticeC[1] = locX;
        latticeC[2] = locY;
        latticeC[3] = locZ+1;
        output[iD] += (_f(latticeC)[iD]*(1-d[0])*(1-d[1])*(d[2]));

        latticeC[1] = locX;
        latticeC[2] = locY+1;
        latticeC[3] = locZ+1;
        output[iD] += (_f(latticeC)[iD]*(1-d[0])*(d[1])*(d[2]));

        latticeC[1] = locX+1;
        latticeC[2] = locY;
        latticeC[3] = locZ+1;
        output[iD] += (_f(latticeC)[iD]*(d[0])*(1-d[1])*(d[2]));

        latticeC[1] = locX+1;
        latticeC[2] = locY+1;
        latticeC[3] = locZ+1;
        output[iD] += (_f(latticeC)[iD]*(d[0])*(d[1])*(d[2]));
      }
#ifdef PARALLEL_MODE_MPI
      dataSize = _f(latticeC).size();
      dataFound = 1;
#endif
      pOutput = output;
    }
  }

#ifdef PARALLEL_MODE_MPI
  if(_communicateToAll) {
    singleton::mpi().reduceAndBcast(dataFound, MPI_SUM);
    singleton::mpi().reduceAndBcast(dataSize, MPI_SUM);
    dataSize/=dataFound;
    if (pOutput.size()==0) {
      for (int iD=0; iD<dataSize; iD++) {
        pOutput.push_back(T());
      }
    }
    for (int iD=0; iD<dataSize; iD++) {
      singleton::mpi().reduceAndBcast(pOutput[iD], MPI_SUM);
    }
    for (int iD=0; iD<dataSize; iD++) {
      pOutput[iD]/=dataFound;
    }
  }
#endif

  return pOutput;
}



template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR>::SuperLatticeFfromAnalyticalF3D(
  AnalyticalF3D<T,T>& f, SuperLattice3D<T,DESCRIPTOR>& sLattice )
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,f.getTargetDim()), _f(f)
{
  this->_name = "fromAnalyticalF";
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR>::operator()
  (std::vector<int> input)
{
  return _f(this->_sLattice.getCuboidGeometry().getPhysR(input));
}



//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeFfromAnalyticalF3D<T,DESCRIPTOR>::
BlockLatticeFfromAnalyticalF3D(AnalyticalF3D<T,T>& f,
  BlockLattice3D<T,DESCRIPTOR>& sLattice, BlockGeometry3D<T>& superGeometry,
  CuboidGeometry3D<T>& cuboidGeometry)
  : BlockLatticeF3D<T,DESCRIPTOR>(sLattice,f.getTargetDim()), _f(f),
    _superGeometry(superGeometry), _cuboidGeometry(cuboidGeometry)
{
  this->_name = "fromAnalyticalF";
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeFfromAnalyticalF3D<T,DESCRIPTOR>::operator()
  (std::vector<int> input)
{
  std::vector<T> physR = _superGeometry.getPhysR(input[0],input[1],input[2]);
  return _f(physR);
}



} // end namespace olb

#endif

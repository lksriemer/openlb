/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
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

#ifndef INTERPOLATION_F_2D_HH
#define INTERPOLATION_F_2D_HH

#include<vector>

#include "functors/interpolationF2D.h"
#include "functors/genericF.h"
#include "functors/analyticalF.h"
#include "functors/indicatorF.h"
#include "core/superLattice2D.h"
#include "dynamics/lbHelpers.h"


namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
AnalyticalFfromSuperLatticeF2D<T,DESCRIPTOR>::
AnalyticalFfromSuperLatticeF2D(SuperLatticeF2D<T,DESCRIPTOR>& f, bool communicateToAll,
  int overlap)
  : AnalyticalF2D<T,T>(f.getTargetDim()), _f(f),
    _cg(_f.getSuperLattice2D().getCuboidGeometry()), _communicateToAll(communicateToAll),
    _overlap(overlap)
{
  this->_name = "fromSuperLatticeF";
}


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> AnalyticalFfromSuperLatticeF2D<T,DESCRIPTOR>::operator() (std::vector<T> physC)
{
  // convert to lattice coordinates
  T d[2];

  if(_communicateToAll) {
	  _f.getSuperLattice2D().communicate();
  }

  int locX, locY;
#ifdef PARALLEL_MODE_MPI
  int dataSize = 0;
  int dataFound = 0;
#endif
  std::vector<T> pOutput;
  std::vector<int> latticeC(3,0);

  for (int iC = 0; iC < _f.getSuperLattice2D().getLoadBalancer().size(); iC++) {
    latticeC[0] = _f.getSuperLattice2D().getLoadBalancer().glob(iC);
    if (_cg.get(latticeC[0]).checkPoint(physC[0], physC[1], locX, locY, _overlap-1)) {
      locX -= (_overlap-1);
      locY -= (_overlap-1);

      std::vector<T> physRiC = _cg.get(latticeC[0]).getPhysR(locX,locY);

      d[0] = (physC[0] - physRiC[0]);
      d[1] = (physC[1] - physRiC[1]);

      d[0] /= _cg.get(latticeC[0]).getDeltaR();
      d[1] /= _cg.get(latticeC[0]).getDeltaR();

      std::vector<T> output(_f.getTargetDim(), T());

      for (unsigned int iD = 0; iD < output.size(); iD++) {
        latticeC[1] = locX;
        latticeC[2] = locY;
        output[iD] += (_f(latticeC)[iD]*(1-d[0])*(1-d[1]));

        latticeC[1] = locX;
        latticeC[2] = locY+1;
        output[iD] += (_f(latticeC)[iD]*(1-d[0])*(d[1]));

        latticeC[1] = locX+1;
        latticeC[2] = locY;
        output[iD] += (_f(latticeC)[iD]*(d[0])*(1-d[1]));

        latticeC[1] = locX+1;
        latticeC[2] = locY+1;
        output[iD] += (_f(latticeC)[iD]*(d[0])*(d[1]));
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
    dataSize /= dataFound;
    if (pOutput.size() == 0) {
      for (int iD = 0; iD < dataSize; iD++) {
        pOutput.push_back(T());
      }
    }
    for (int iD = 0; iD < dataSize; iD++) {
      singleton::mpi().reduceAndBcast(pOutput[iD], MPI_SUM);
    }
    for (int iD = 0; iD < dataSize; iD++) {
      pOutput[iD] /= dataFound;
    }
  }
#endif

  return pOutput;
}



template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR>::
SuperLatticeFfromAnalyticalF2D(AnalyticalF2D<T,T>& f, SuperLattice2D<T,DESCRIPTOR>& sLattice,
  SuperGeometry2D<T>& sg)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, f.getTargetDim()), _f(f), _sg(sg)
{
  this->_name = "fromAnalyticalF";
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
	// convert to physical coordinates
	std::vector<T> physCoordinate;
	physCoordinate.push_back(this->_sg.getPhysR(input)[0]);
	physCoordinate.push_back(this->_sg.getPhysR(input)[1]);

	return _f(physCoordinate);
}



} // end namespace olb

#endif

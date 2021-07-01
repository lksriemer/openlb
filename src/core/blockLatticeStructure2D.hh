/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
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

/** \file
 * Dynamics for a generic 2D block structure -- header file.
 */
#ifndef BLOCK_LATTICE_STRUCTURE_2D_HH
#define BLOCK_LATTICE_STRUCTURE_2D_HH

#include <vector>
#include "blockLatticeStructure2D.h"


namespace olb {


template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineDynamics(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  Dynamics<T,Lattice>* dynamics)
{
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        get(iX,iY).defineDynamics(dynamics);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineRho(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& rho)
{
  T rhoTmp;
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        std::vector<T> physCoordinate = blockGeometry.getPhysR(iX,iY);
        rhoTmp = rho(physCoordinate)[0];
        get(iX,iY).defineRho(rhoTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineU(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& u)
{
  T uTmp[2];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        std::vector<T> physCoordinate = blockGeometry.getPhysR(iX,iY);
        for (int i=0; i<2; i++) {
          uTmp[i] = u(physCoordinate)[i];
        }
        get(iX,iY).defineU(uTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineRhoU(BlockGeometryStructure2D<T>& blockGeometry,
  int material, AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  T rhoTmp;
  T uTmp[2];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        std::vector<T> physCoordinate = blockGeometry.getPhysR(iX,iY);
        rhoTmp = rho(physCoordinate)[0];
        for (int i=0; i<2; i++) {
          uTmp[i] = u(physCoordinate)[i];
        }
        get(iX,iY).defineRhoU(rhoTmp,uTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::definePopulations(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& Pop)
{
  T PopTmp[Lattice<T>::q];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        std::vector<T> physCoordinate = blockGeometry.getPhysR(iX,iY);
        for (int i=0; i<Lattice<T>::q; i++) {
          PopTmp[i] = Pop(physCoordinate)[i];
        }
        get(iX,iY).definePopulations(PopTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineExternalField(
  BlockGeometryStructure2D<T>& blockGeometry, int material, int fieldBeginsAt,
  int sizeOfField, AnalyticalF2D<T,T>& field)
{
  T* fieldTmp = new T [sizeOfField];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        std::vector<T> physCoordinate = blockGeometry.getPhysR(iX,iY);
        for (int i=0; i<sizeOfField; i++) {
          fieldTmp[i] = field(physCoordinate)[i];
        }
        get(iX,iY).defineExternalField(fieldBeginsAt,sizeOfField, fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}


template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::iniEquilibrium(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  AnalyticalF2D<T,T>& rho , AnalyticalF2D<T,T>& u)
{
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        std::vector<T> physCoordinate = blockGeometry.getPhysR(iX,iY);
        T uTmp[] = {u(physCoordinate)[0],u(physCoordinate)[1]};
        get(iX,iY).iniEquilibrium(rho(physCoordinate)[0],uTmp);
      }
    }
  }
}

}  // namespace olb

#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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
 * The description of a vector of 3D cuboid  -- generic implementation.
 */


#ifndef CUBOID_GEOMETRY_3D_HH
#define CUBOID_GEOMETRY_3D_HH

#include <vector>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "geometry/cuboid3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "functors/indicatorF.h"

namespace olb {

////////////////////// Class CuboidGeometry3D /////////////////////////

template <typename T, typename S> class IndicatorF3D;

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(T originX, T originY, T originZ, T deltaR, int nX, int nY, int nZ, int nC)
  : _motherCuboid(originX, originY, originZ, deltaR, nX, nY, nZ), 
    _periodicityOn(3,bool(false)), 
    clout(std::cout,"CuboidGeometry3D")
{
  add(_motherCuboid);
  split(0, nC);
}

template<typename T>
CuboidGeometry3D<T>::CuboidGeometry3D(IndicatorF3D<bool,T>& indicator, T voxelSize, int nC) 
  : _motherCuboid( indicator.getMin()[0],  indicator.getMin()[1],  indicator.getMin()[2], voxelSize, (int)((indicator.getMax()[0]-indicator.getMin()[0])/voxelSize+1.5), (int)((indicator.getMax()[1]-indicator.getMin()[1])/voxelSize+1.5), (int)((indicator.getMax()[2]-indicator.getMin()[2])/voxelSize+1.5)), 
    _periodicityOn(3,bool(false)),
    clout(std::cout,"CuboidGeometry3D") 
{
  // Identity prevents tmp from being deleted!
  IndicatorIdentity3D<bool,T> indicatorF(indicator);

  add(_motherCuboid);
  split(0, nC);
  shrink(indicatorF);
}


template<typename T>
Cuboid3D<T>& CuboidGeometry3D<T>::get(int iC) {
  return _cuboids[iC];
}

template<typename T>
Cuboid3D<T> const& CuboidGeometry3D<T>::get(int iC) const {
  return _cuboids[iC];
}

template<typename T>
Cuboid3D<T> CuboidGeometry3D<T>::getMotherCuboid() const {
  return _motherCuboid;
}

template<typename T>
void CuboidGeometry3D<T>::setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ) {
  _periodicityOn.resize(3);
  _periodicityOn[0] = periodicityX;
  _periodicityOn[1] = periodicityY;
  _periodicityOn[2] = periodicityZ;
}


template<typename T>
int CuboidGeometry3D<T>::get_iC(T x, T y, T z, int offset) const {
  unsigned i;
  for (i=0; i<_cuboids.size(); i++) {
    if (_cuboids[i].checkPoint(x, y, z, offset)) return (int)i;
  }
  return (int)i;
}

template<typename T>
int CuboidGeometry3D<T>::get_iC(T x, T y, T z,
                                int orientationX, int orientationY, int orientationZ) const {
  unsigned i;
  for (i=0; i<_cuboids.size(); i++) {
    if (_cuboids[i].checkPoint(x, y, z) &&
        _cuboids[i].checkPoint(x + orientationX/_cuboids[i].getDeltaR(),
                               y + orientationY/_cuboids[i].getDeltaR(),
                               z + orientationZ/_cuboids[i].getDeltaR())) {
      return (int)i;
    }
  }
  return (int)i;
}

template<typename T>
bool CuboidGeometry3D<T>::getC(std::vector<T> physR, int& iC) const {
  int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp<getNc()) {
    iC=iCtmp; 
    return true;
  }
  else return false;
}

template<typename T>
bool CuboidGeometry3D<T>::getLatticeR(std::vector<T> physR, std::vector<int>& latticeR) const {
 int iCtmp = get_iC(physR[0], physR[1], physR[2]);
 if (iCtmp<getNc()) {
   latticeR[0] = iCtmp; 
   latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] )/_cuboids[latticeR[0]].getDeltaR() +.5);
   latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] )/_cuboids[latticeR[0]].getDeltaR() +.5);
   latticeR[3] = (int)floor( (physR[2] - _cuboids[latticeR[0]].getOrigin()[2] )/_cuboids[latticeR[0]].getDeltaR() +.5);
   return true;
  }
  else return false;
}

template<typename T>
bool CuboidGeometry3D<T>::getFloorLatticeR(std::vector<T> physR, std::vector<int>& latticeR) const {
  int iCtmp = get_iC(physR[0], physR[1], physR[2]);
  if (iCtmp<getNc()) {
    latticeR[0]=iCtmp; 
    latticeR[1] = (int)floor( (physR[0] - _cuboids[latticeR[0]].getOrigin()[0] )/_cuboids[latticeR[0]].getDeltaR() );
    latticeR[2] = (int)floor( (physR[1] - _cuboids[latticeR[0]].getOrigin()[1] )/_cuboids[latticeR[0]].getDeltaR() );
    latticeR[3] = (int)floor( (physR[2] - _cuboids[latticeR[0]].getOrigin()[2] )/_cuboids[latticeR[0]].getDeltaR() );
    return true;}
  else return false;
}

template<typename T>
std::vector<T> CuboidGeometry3D<T>::getPhysR(int iCglob, int iX, int iY, int iZ) const {
  std::vector<T> physR = _cuboids[iCglob].getPhysR(iX,iY,iZ);
  for (int iDim=0; iDim<3; iDim++) {
    if (_periodicityOn[iDim]) {
      //std::cout << iDim << _periodicityOn[iDim] <<":"<< _motherCuboid.getDeltaR()*(_motherCuboid.getExtend()[iDim]) << std::endl;
      physR[iDim] = remainder( physR[iDim] - _motherCuboid.getOrigin()[iDim]
                                       + _motherCuboid.getDeltaR()*(_motherCuboid.getExtend()[iDim]) ,
                           _motherCuboid.getDeltaR()*(_motherCuboid.getExtend()[iDim]));
      // solving the rounding error problem for double 
      if( physR[iDim]*physR[iDim] < 0.001*_motherCuboid.getDeltaR()*_motherCuboid.getDeltaR() ) {
        if( physR[iDim] > 0 ) physR[iDim] = _motherCuboid.getDeltaR()*_motherCuboid.getExtend()[iDim];
        else  physR[iDim] = T();
      }
      // make it to mod instead remainer
      if( physR[iDim] < 0 ) physR[iDim] += _motherCuboid.getDeltaR()*_motherCuboid.getExtend()[iDim];
      // add origin
      physR[iDim] += _motherCuboid.getOrigin()[iDim];
    }
  }  
  return physR;
}

template<typename T>
std::vector<T> CuboidGeometry3D<T>::getPhysR(std::vector<int> latticeR) const {
 return getPhysR(latticeR[0], latticeR[1], latticeR[2], latticeR[3]);
}


template<typename T>
void CuboidGeometry3D<T>::getNeighbourhood(int cuboid, std::vector<int>& neighbours, int overlap) {
  neighbours.clear();
  for (int iC=0; iC<getNc(); iC++) {
    if(cuboid == iC) continue;
    T globX = get(iC).getOrigin()[0];
    T globY = get(iC).getOrigin()[1];
    T globZ = get(iC).getOrigin()[2];
    T nX = get(iC).getNx();
    T nY = get(iC).getNy();
    T nZ = get(iC).getNz();
    T deltaR = get(iC).getDeltaR();
    if(get(cuboid).checkInters(globX-overlap*deltaR,
    		globX+(nX+overlap-1)*deltaR,
			globY-overlap*deltaR,
			globY+(nY+overlap-1)*deltaR,
			globZ-overlap*deltaR,
			globZ+(nZ+overlap-1)*deltaR, overlap)) {
      neighbours.push_back(iC);
    }
  }
}

template<typename T>
int CuboidGeometry3D<T>::getNc() const {
  return _cuboids.size();
}

template<typename T>
T CuboidGeometry3D<T>::getMinRatio() const {
  T minRatio = 1.;
  for (unsigned i=0; i<_cuboids.size(); i++) {
    if((T)_cuboids[i].getNx()/(T)_cuboids[i].getNy() < minRatio) {
      minRatio = (T)_cuboids[i].getNx()/(T)_cuboids[i].getNy();
    }
    if((T)_cuboids[i].getNy()/(T)_cuboids[i].getNz() < minRatio) {
      minRatio = (T)_cuboids[i].getNy()/(T)_cuboids[i].getNz();
    }
    if((T)_cuboids[i].getNz()/(T)_cuboids[i].getNx() < minRatio) {
      minRatio = (T)_cuboids[i].getNz()/(T)_cuboids[i].getNx();
    }
  }
  return minRatio;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxRatio() const {
  T maxRatio = 1.;
  for (unsigned i=0; i<_cuboids.size(); i++) {
    if((T)_cuboids[i].getNx()/(T)_cuboids[i].getNy() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNx()/(T)_cuboids[i].getNy();
    }
    if((T)_cuboids[i].getNy()/(T)_cuboids[i].getNz() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNy()/(T)_cuboids[i].getNz();
    }
    if((T)_cuboids[i].getNz()/(T)_cuboids[i].getNx() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNz()/(T)_cuboids[i].getNx();
    }
  }
  return maxRatio;
}

template<typename T>
T CuboidGeometry3D<T>::getMinPhysVolume() const {
  T minVolume = _cuboids[0].getPhysVolume();
  for (unsigned i=0; i<_cuboids.size(); i++) {
    if(_cuboids[i].getPhysVolume() < minVolume) {
      minVolume = _cuboids[i].getPhysVolume();
    }
  }
  return minVolume;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxPhysVolume() const {
  T maxVolume = _cuboids[0].getPhysVolume();
  for (unsigned i=0; i<_cuboids.size(); i++) {
    if(_cuboids[i].getPhysVolume() > maxVolume) {
      maxVolume = _cuboids[i].getPhysVolume();
    }
  }
  return maxVolume;
}

template<typename T>
int CuboidGeometry3D<T>::getMinLatticeVolume() const {
  int minNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i=0; i<_cuboids.size(); i++) {
    if(_cuboids[i].getLatticeVolume() < minNodes) {
      minNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return minNodes;
}

template<typename T>
int CuboidGeometry3D<T>::getMaxLatticeVolume() const {
  int maxNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i=0; i<_cuboids.size(); i++) {
    if(_cuboids[i].getLatticeVolume() > maxNodes) {
      maxNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return maxNodes;
}

template<typename T>
T CuboidGeometry3D<T>::getMinDeltaR() const {
  T minDelta = _cuboids[0].getDeltaR();
  for (unsigned i=0; i<_cuboids.size(); i++) {
    if(_cuboids[i].getDeltaR() < minDelta) {
      minDelta = _cuboids[i].getDeltaR();
    }
  }
  return minDelta;
}

template<typename T>
T CuboidGeometry3D<T>::getMaxDeltaR() const {
  T maxDelta = _cuboids[0].getDeltaR();
  for (unsigned i=0; i<_cuboids.size(); i++) {
    if(_cuboids[i].getDeltaR() > maxDelta) {
      maxDelta = _cuboids[i].getDeltaR();
    }
  }
  return maxDelta;
}

template<typename T>
void CuboidGeometry3D<T>::add(Cuboid3D<T> cuboid) {

  _cuboids.push_back(cuboid);
}

template<typename T>
void CuboidGeometry3D<T>::remove(int iC) {

  _cuboids.erase(_cuboids.begin() + iC);
}


template<typename T>
void CuboidGeometry3D<T>::remove(IndicatorF3D<bool,T>& indicatorF) {

  //IndicatorIdentity3D<bool,T> tmpIndicatorF(indicatorF);

  std::vector<bool> allZero;
  std::vector<int> latticeR(4,0);
  for (unsigned iC=0; iC < _cuboids.size(); iC++) {
    latticeR[0] = iC;
    allZero.push_back(true);
    for (int iX=0; iX<_cuboids[iC].getNx(); iX++) {
      for (int iY=0; iY<_cuboids[iC].getNy(); iY++) {
        for (int iZ=0; iZ<_cuboids[iC].getNz(); iZ++) {
          latticeR[1] = iX; latticeR[2] = iY; latticeR[3] = iZ;
          if (indicatorF(getPhysR(latticeR))[0])
            allZero[iC] = 0;
        }
      }
    }
  }
  for (int iC=_cuboids.size()-1; iC>=0; iC--) {
    if (allZero[iC] ) {
      remove(iC);
    }
  }
}

template<typename T>
void CuboidGeometry3D<T>::shrink(IndicatorF3D<bool,T>& indicatorF) {

  //IndicatorIdentity3D<bool,T> tmpIndicatorF(indicatorF);
  int newX,newY,newZ,maxX,maxY,maxZ;
  int nC = getNc();
  std::vector<int> latticeR(4,0);
  for (int iC = nC-1; iC >= 0; iC--) {
    latticeR[0] = iC;
    int fullCells = 0;
    int xN = get(iC).getNx();
    int yN = get(iC).getNy();
    int zN = get(iC).getNz();
    maxX = 0; maxY = 0; maxZ = 0;
    newX = xN-1; newY = yN-1; newZ = zN-1;
    for (int iX=0; iX<xN; iX++) {
      for (int iY=0; iY<yN; iY++) {
        for (int iZ=0; iZ<zN; iZ++) {
          latticeR[1] = iX; latticeR[2] = iY; latticeR[3] = iZ;
          if (indicatorF(getPhysR(latticeR))[0]) {
            fullCells++;
            maxX=std::max(maxX,iX); maxY=std::max(maxY,iY); maxZ=std::max(maxZ,iZ);
            newX=std::min(newX,iX); newY=std::min(newY,iY); newZ=std::min(newZ,iZ);
          }
        }
      }
    }
//    if (maxX+2 < xN) maxX+=2; else if (maxX+1 < xN) maxX+=1;
//    if (maxY+2 < yN) maxY+=2; else if (maxY+1 < yN) maxY+=1;
//    if (maxZ+2 < zN) maxZ+=2; else if (maxZ+1 < zN) maxZ+=1;
//
//    if (newX-2 >= 0) newX-=2; else if (newX-1 >= 0) newX-=1;
//    if (newY-2 >= 0) newY-=2; else if (newY-1 >= 0) newY-=1;
//    if (newZ-2 >= 0) newZ-=2; else if (newZ-1 >= 0) newZ-=1;

    if (fullCells > 0) {
      get(iC).setWeight(fullCells);
      _cuboids[iC].resize(newX, newY, newZ, maxX-newX+1, maxY-newY+1, maxZ-newZ+1);
    }
    else {
      remove(iC);
    }
  }
}


template<typename T>
void CuboidGeometry3D<T>::split(int iC, int p) {

  Cuboid3D<T> temp(_cuboids[iC].getOrigin()[0],_cuboids[iC].getOrigin()[1],
                   _cuboids[iC].getOrigin()[2],  _cuboids[iC].getDeltaR(),
                   _cuboids[iC].getNx(), _cuboids[iC].getNy(), _cuboids[iC].getNz());
  temp.divide(p, _cuboids);
  remove(iC);
}

template<typename T>
void CuboidGeometry3D<T>::setWeights(IndicatorF3D<bool,T>& indicatorF) {

  //IndicatorIdentity3D<bool,T> tmpIndicatorF(indicatorF);
  int xN, yN, zN;
  int nC = getNc();
  std::vector<int> latticeR(4,0);
  for( int iC = nC-1; iC >= 0; iC--) { // assemble neighbourhood information
    latticeR[0] = iC;
    xN  = get(iC).getNx();
    yN  = get(iC).getNy();
    zN  = get(iC).getNz();
    int fullCells = 0;
    for (int iX=0; iX<xN; iX++) {
      for (int iY=0; iY<yN; iY++) {
        for (int iZ=0; iZ<zN; iZ++) {
          latticeR[1] = iX; latticeR[2] = iY; latticeR[3] = iZ;
          if (indicatorF(getPhysR(latticeR))[0]) {
            fullCells++;
          }
        }
      }
    }
    if (fullCells > 0) {
      get(iC).setWeight(fullCells);
    }
    else {
      remove(iC);
    }
  }
}


template<typename T>
void CuboidGeometry3D<T>::print() const {
  clout << "---Cuboid Stucture Statistics---" << std::endl;
  clout << " Number of Cuboids: " << "\t" << getNc() << std::endl;
  clout << " Delta (min): " << "\t" << "\t" << getMinDeltaR() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxDeltaR() << std::endl;
  clout << " Ratio (min): " << "\t" << "\t" << getMinRatio() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxRatio() << std::endl;
  clout << " Nodes (min): " << "\t" << "\t" << getMinLatticeVolume() << std::endl;
  clout << "       (max): " << "\t" << "\t" << getMaxLatticeVolume() << std::endl;
  clout << "--------------------------------" << std::endl;
}

template<typename T>
void CuboidGeometry3D<T>::printExtended() {
  clout << "Mothercuboid :" << std::endl;
  getMotherCuboid().print();

  for (int iC=0; iC<getNc(); iC++) {
    clout << "Cuboid #" << iC << ": " << std::endl;
    get(iC).print();
  }
}


}  // namespace olb

#endif

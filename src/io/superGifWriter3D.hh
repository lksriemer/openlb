/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias Krause
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

#ifndef SUPER_GIF_WRITER_3D_HH
#define SUPER_GIF_WRITER_3D_HH

//#include <cstdlib>
//#include <cmath>
#include <iostream>
//#include <sstream>
//#include <fstream>
#include <string>
#include <vector>
#include "colormaps.h"
//#include "core/singleton.h"
//#include "communication/mpiManager.h"
#include "geometry/cuboidGeometry3D.h"
#include "utilities/vectorHelpers.h"
#include "functors/interpolationF3D.h"

namespace olb {

////////// class BlockGifWriter3D ////////////////////////////////////////

template<typename T, template <typename U> class DESCRIPTOR>
SuperGifWriter3D<T,DESCRIPTOR>::SuperGifWriter3D(std::string const& map)
  : clout(std::cout,"ImageWriter"), _colorRange(1024), _numColors(1024),
    _colorMap( graphics::mapGenerators::generateMap<T>(map) )
{ }

template<typename T, template <typename U> class DESCRIPTOR>
void SuperGifWriter3D<T,DESCRIPTOR>::write( SuperLatticeF3D<T,DESCRIPTOR>& f,
  std::vector<int> normal)
{
//  util::print<int>( vecU(normal) );
//  util::print<int>( vecV(normal) );

//  CuboidGeometry3D<T> const& cGeometry = f.getSuperLattice3D().getCuboidGeometry();
//  std::cout << "motherC getNx: " <<  cGeometry.getMotherCuboid().getNx() << std::endl;
//  std::cout << "motherC getNy: " <<  cGeometry.getMotherCuboid().getNy() << std::endl;
//  std::cout << "motherC getNz: " <<  cGeometry.getMotherCuboid().getNz() << std::endl;
  setFrameExt( f, normal );
  writeXYplane( f );

}


// determines from the given normal, which plane is demanded
// stores extension of the slice, therefore one extension is allways 0
// only normals: e1,e2,e3 are allowed.
template<typename T, template <typename U> class DESCRIPTOR>
void SuperGifWriter3D<T,DESCRIPTOR>::setFrameExt( SuperLatticeF3D<T,DESCRIPTOR>& f,
  std::vector<int> normal )
{

  CuboidGeometry3D<T> const& cGeometry = f.getSuperLattice3D().getCuboidGeometry();
//  LoadBalancer<T>& load = f.getSuperLattice3D().getLoadBalancer();
//  for (int iCloc=0; iCloc<load.size(); iCloc++)
//  { // every thread gets his cuboids via load.glob(iC)
//    // hidden parallelization !!
//    int nx = cGeometry.get(load.glob(iCloc)).getNx();
//    int ny = cGeometry.get(load.glob(iCloc)).getNy();
//    int nz = cGeometry.get(load.glob(iCloc)).getNz();
//  }

  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if( rank == 0 )
  {
    _extNx = cGeometry.getMotherCuboid().getNx();
    _extNy = cGeometry.getMotherCuboid().getNy();
    _extNz = cGeometry.getMotherCuboid().getNz();
    if(normal[0] == 1 && normal[1] == 0 && normal[2] == 0)
    {
  //      std::cout << "Y-Z Plane" << std::endl;
      _extNx = T(0);
    }
    else if(normal[0] == 0 && normal[1] == 1 && normal[2] == 0)
    {
  //      std::cout << "X-Z Plane" << std::endl;
      _extNy = T(0);
    }
    else if(normal[0] == 0 && normal[1] == 0 && normal[2] == 1)
    {
  //      std::cout << "X-Y Plane" << std::endl;
      _extNz = T(0);
    }
  }

}


template<typename T, template <typename U> class DESCRIPTOR>
void SuperGifWriter3D<T,DESCRIPTOR>::writeXYplane( SuperLatticeF3D<T,DESCRIPTOR>& f )
{

  std::string fullName = singleton::directories().getImageOutDir() + f.getName()+".ppm";
  std::ofstream fout(fullName.c_str());
  // write header
  fout << "P3\n";
  // dimension of image
  fout << _extNx << " " << _extNy << "\n";
  // dynamic range
  fout << (_colorRange-1) << "\n";

  AnalyticalFfromSuperLatticeF3D<T,DESCRIPTOR> ana(f,true);
  
  int iC = 0;

  int dim = f.getTargetDim();
  std::vector<double> tmpVec(4,T(0));
  for( int iY = _extNy ; iY >= 0 ; iY -- )
  {
    for( int iX = 0; iX <= _extNx; iX ++ )
    {
//      coordinate: _extNx/2 + iX, _extNy/2 + iY, _extNz/2 + iZ
      for( int iDim = 0; iDim < dim; iDim ++ )
      {
        tmpVec[0]=iC; tmpVec[1]=iX; tmpVec[2]=iY; tmpVec[3]=0;
        T evaluated = ana(tmpVec)[iDim];

        if (evaluated <   T()) evaluated = T();
        if (evaluated >= T(1)) evaluated = (T)(_numColors-1)/(T)_numColors;
        graphics::rgb<T> color = _colorMap.get(evaluated);
        fout << (int) (color.r*(_colorRange-1)) << " "
             << (int) (color.g*(_colorRange-1)) << " "
             << (int) (color.b*(_colorRange-1)) << "\n";
      }
    }
  }

}


// functor evalutation
//          for (int iDim=0; iDim<dim; ++iDim) {
//            tmpVec[0] = iX; tmpVec[1] = iY; tmpVec[2] = iZ;
//            T evaluated = f(tmpVec)[iDim]; 
//            if(singleton::mpi().getRank()==0) fout << evaluated << " ";



//// vecOut spans plane
//// problem of the sign for planes
//// geometry extens from 0,10 on x axis, but what parameters of s,t ???
//template<typename T, template <typename U> class DESCRIPTOR>
//std::vector<int> BlockGifWriter3D<T,DESCRIPTOR>::vecU( std::vector<int> normal )
//{
////   c++11
////  std::vector<T> vecOut = {-normal[1],normal[0],T(0)};
//  std::vector<int> vecOut( 3,int(0) );
//  vecOut[0] = -normal[1];
//  vecOut[1] = normal[0];
////  vecOut[2] = T(0); not necessary due to instanciation

//  // if vecOut equals zero, than vecOut = (-n3, 0, n1)
//  if( vecOut[0] == 0 && vecOut[1] == 0 && vecOut[2] == 0 )
//  {
//    vecOut[0] = -normal[2];
//    vecOut[1] = T(0);
//    vecOut[2] = normal[1];
//  }
//  return vecOut;
//}

//// vecOut spans plane
//template<typename T, template <typename U> class DESCRIPTOR>
//std::vector<int> BlockGifWriter3D<T,DESCRIPTOR>::vecV( std::vector<int> normal )
//{
////   c++11
////  std::vector<T> vecOut = {T(0),-normal[2],normal[1]};
//  std::vector<int> vecOut( 3,int(0) );
////  vecOut[0] = T(0); not necessary due to instanciation
//  vecOut[1] = -normal[2];
//  vecOut[2] = normal[1];

//  // if vecOut equals zero, than vecOut = (-n3, 0, n1)
//  if( vecOut[0] == 0 && vecOut[1] == 0 && vecOut[2] == 0 )
//  {
//    vecOut[0] = -normal[2];
//    vecOut[1] = T(0);
//    vecOut[2] = normal[1];
//  }
//  return vecOut;
//}


}  // namespace olb

#endif

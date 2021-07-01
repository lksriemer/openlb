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

#ifndef BLOCK_GIF_WRITER_3D_HH
#define BLOCK_GIF_WRITER_3D_HH

//#include <cstdlib>
//#include <cmath>
#include <iostream>
//#include <sstream>
//#include <fstream>
#include <string>
#include <vector>
#include "colormaps.h"
#include "core/singleton.h"
#include "communication/mpiManager.h"
#include "geometry/cuboidGeometry3D.h"
#include "utilities/vectorHelpers.h"

namespace olb {

////////// class BlockGifWriter3D ////////////////////////////////////////

template<typename T, template <typename U> class DESCRIPTOR>
BlockGifWriter3D<T,DESCRIPTOR>::BlockGifWriter3D(std::string const& map)
  : clout(std::cout,"ImageWriter"), _colorRange(1024), _numColors(1024),
    _colorMap( graphics::mapGenerators::generateMap<T>(map) )
{ }

template<typename T, template <typename U> class DESCRIPTOR>
void BlockGifWriter3D<T,DESCRIPTOR>::write( BlockLatticeF3D<T,DESCRIPTOR>& f,
  std::vector<int> normal)
{
//  util::print<int>( vecU(normal) );
//  util::print<int>( vecV(normal) );

  writeGif( f, normal);

}


// determines from the given normal, which plane is demanded
// stores extension of the slice, therefore one extension is allways 0
// only normals: e1,e2,e3 are allowed.
template<typename T, template <typename U> class DESCRIPTOR>
void BlockGifWriter3D<T,DESCRIPTOR>::setFrame( BlockLatticeF3D<T,DESCRIPTOR>& f,
  std::vector<int> normal )
{

  if(singleton::mpi().getRank()==0)
  {
    if(normal[0] == 1 && normal[1] == 0 && normal[2] == 0){
//      std::cout << "Y-Z Plane" << std::endl;
      _extNx = T(0);
      _extNy = f.getBlockLattice3D().getNy() -1;
      _extNz = f.getBlockLattice3D().getNz() -1;
    }
    else if(normal[0] == 0 && normal[1] == 1 && normal[2] == 0){
//      std::cout << "X-Z Plane" << std::endl;
      _extNx = f.getBlockLattice3D().getNx() -1;
      _extNy = T(0);
      _extNz = f.getBlockLattice3D().getNz() -1;
    }
    else if(normal[0] == 0 && normal[1] == 0 && normal[2] == 1){
//      std::cout << "X-Y Plane" << std::endl;
      _extNx = f.getBlockLattice3D().getNx() -1;
      _extNy = f.getBlockLattice3D().getNy() -1;
      _extNz = T(0);
    }
  }
}


template<typename T, template <typename U> class DESCRIPTOR>
void BlockGifWriter3D<T,DESCRIPTOR>::writeGif( BlockLatticeF3D<T,DESCRIPTOR>& f,
   std::vector<int> normal )
{

/*
origin in terms of functors is given by:
  _extNx/2
  _extNy/2
  _extNz/2
*/

//   erbt von AnalyticalF3D Funktor.
//   AnalyticalFfromSuperLatticeF3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
//                                  SuperGeometry3D<T>& superGeometry,
//                                  bool communicateToAll=false, int overlap=2);

//  AnalyticalFfromSuperLatticeF3D<T,DESCRIPTOR> ana();


//  for(int iX = -_imageNx/2; iX <= _imageNx/2; iX ++)
//  {
//    for(int iY = -_imageNy/2; iY <= _imageNy/2; iY ++)
//    {
//      for(int iZ = -_imageNz/2; iZ <= _imageNz/2; iZ ++)
//      {
////        coordinate: _extNx/2 + iX, _extNy/2 + iY, _extNz/2 + iZ
//        
//      }
//    }
//  }

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

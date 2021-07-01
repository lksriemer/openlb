/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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
 * A method to write vtk data for block geometries
 * (only for uniform grids) -- header file.
 *
 * One can add functors via addFunctors. To write added functors
 * call write(int iT=0).
 * 
 * To write a functor without adding him, call 
 * call write(BlockLatticeF3D<T,DESCRIPTOR>& f, int iT).
 *
 * default output type is binary, to change it have a look on the
 * constructor.
 * 
 */

#ifndef BLOCK_VTK_WRITER_3D_H
#define BLOCK_VTK_WRITER_3D_H

#include "io/ostreamManager.h"
//#include "functors/blockLatticeBaseF3D.h"
#include "functors/blockLatticeCalcF3D.h"     // BlockLatticeIdentity

namespace olb {

template<typename T, template <typename U> class DESCRIPTOR>
class BlockVTKwriter3D {
public:
  BlockVTKwriter3D( std::string name, bool binary = true );
  ~BlockVTKwriter3D();
  ///  method calls preamble(), pointData(), data() and coresponding
  ///  closing methods.
  ///  writes given functor
  void write( BlockLatticeF3D<T,DESCRIPTOR>& f, int iT = 0 );
  ///  writes functors stored at pointerVec
  void write(int iT=0);
  ///  put functor to _pointerVec
  ///  to simplify writing process of several functors
  void addFunctor( BlockLatticeF3D<T,DESCRIPTOR>& f );
  ///  to clear stored functors
  void clearAddedFunctors();
private:
  ///  writes <VTKFile .... >, <ImageData ... > and <Piece ... >
  void preamble( std::string& fullName, int nx,int ny,int nz, 
                 T originX=0, T originY=0, T originZ=0);
  ///  writes </Piece>, </ImageData> and  </VTKFile>
  void closePreamble( std::string& fullName );
  ///  writes <PointData Scalar="..." >
  void pointData( std::string& fullName );
  ///  writes </PointData>
  void closePointData( std::string& fullName  );
  ///  writes all functors stored at pointerVec
  void dataArray(std::string& fullNameVti, int nx, int ny, int nz);
  void dataArrayBinary(std::string& fullNameVti, int nx, int ny, int nz);
  ///  writes instantaniously given functor, without adding to _pointerVec
  void dataFunctorToOneFile( std::string& fullNameVti, 
                            BlockLatticeF3D<T,DESCRIPTOR>& f,
                             int nx, int ny, int nz );
  void dataFunctorToOneFileBinary( std::string& fullNameVti, 
                             BlockLatticeF3D<T,DESCRIPTOR>& f,
                             int nx, int ny, int nz );
private:
  mutable OstreamManager clout;
  ///  output files are called "_name + iT + .vti"
  std::string _name;
  std::vector< BlockLatticeIdentity3D<T,DESCRIPTOR>* > _pointerVec;
  int _offset;
  bool _binary;
};

}  // namespace olb


#endif

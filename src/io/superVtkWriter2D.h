/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014  Albert Mink, Mathias J. Krause
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
 * A method to write vtk data for cuboid geometries
 * (only for uniform grids) -- header file.
 *
 *
 * In .pvd files, there are only links/references to the files with the real 
 * data - the .vti files 
 * 
 * .pvd file structur
 * there exists one so called ...master.pvd
 * and for every timestep a ...iT.pvd file
 *
 * for every timestep iT and every cuboid iC
 * there is a .vti file which stores the real
 * data.
 * 
 *
 */

#ifndef SUPER_VTK_WRITER_2D_H
#define SUPER_VTK_WRITER_2D_H

#include <sstream>
#include <iomanip>
#include <vector>
#include "io/ostreamManager.h"
#include "functors/superLatticeBaseF2D.h"
#include "geometry/cuboid2D.h"

namespace olb {

template<typename T, template <typename U> class DESCRIPTOR>
class SuperVTKwriter2D {
public:
  SuperVTKwriter2D( std::string name, bool binary = true );
  ///  writes functors stored in pointerVec
  ///  every thread writes a vti file with data from his cuboids
  ///  the vti files are linked in a pvd file
  void write(int iT=0);
  ///  writes functor instantaneously, same vti-pvd file structure as above
  void write(SuperLatticeF2D<T,DESCRIPTOR>& f, int iT = 0);
  ///  have to be called before writting added functors, in order to create
  ///  master.pvd file
  void createMasterFile();
  ///  put functor to _pointerVec
  ///  to simplify writing process of several functors
  void addFunctor(SuperLatticeF2D<T,DESCRIPTOR>& f);
  ///  to clear stored functors
  void clearAddedFunctors();
private:
  ///  performes <VTKFile ...>, <ImageData ...> and <PieceExtent ...>
  void preambleVTI(const std::string& fullName, int x0, int y0, int x1, int y1, 
                   T originX,T originY, T originZ, T delta);
  ///  performes </ImageData> and </VTKFile>
  void closeVTI(const std::string& fullNamePiece);
  ///  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  ///  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);
  ///  performes <PointData ...>
  void pointData(const std::string& fullName);
  ///  performes </PointData>
  void closePointData(const std::string& fullName);
  ///  performes <DataSet timestep= ... file=namePiece />
  ///  *** nasty function ***
  void dataPVD(int iT, int iC, const std::string& fullNamePVD,
               const std::string& namePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVDmaster(int iT, int iC, const std::string& fullNamePVDMaster,
                     const std::string& namePiece);
  ///  writes given functor f
  ///  !!functor f have to be saved with Identity!!
  void dataArray(const std::string& fullName, SuperLatticeF2D<T,DESCRIPTOR>& f,
                 int iC, int nx, int ny);
  ///  writes given functor f, base64
  ///  !!functor f have to be saved with Identity!!
  void dataArrayBinary(const std::string& fullName, SuperLatticeF2D<T,DESCRIPTOR>& f,
                       int iC, int nx, int ny);
  ///  performes </Piece>
  void closePiece(const std::string& fullNamePiece);
private:
  mutable OstreamManager clout;
  ///  default is false, call createMasterFile() and it will be true
  bool _createFile;
  ///  determines the name of .vti and .pvd per iT
  std::string _name;
  std::vector< SuperLatticeF2D<T,DESCRIPTOR>* > _pointerVec;
  bool _binary;
};

}  // namespace olb


#endif

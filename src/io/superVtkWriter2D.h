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
 */

#ifndef SUPER_VTK_WRITER_2D_H
#define SUPER_VTK_WRITER_2D_H

#include <sstream>
#include <iomanip>
#include <vector>
#include "io/ostreamManager.h"
#include "functors/superBaseF2D.h"


namespace olb {

/** SuperVTKwriter2D writes any SuperF2D to vtk-based output files.
 *
 * In .pvd files, there are only links/references to the files with the real
 * data - the .vti files
 *
 * .pvd file structur
 * there exists one file named after the name given in the ctor
 * and for every timestep a xxx.pvd file named after the correspondig time step.
 *
 * For every timestep iT and every cuboid iC
 * there is a .vti file which stores the real
 * data.
 */
template<typename T>
class SuperVTKwriter2D {
public:
  /** construtor
   *  \param name determines beginning of output file name
   *  \param binary corresponds to the vtk data format, ASCII or binary
   */
  SuperVTKwriter2D( std::string name, bool binary=true );
  ///  writes functors stored in pointerVec
  ///  every thread writes a vti file with data from his cuboids
  ///  the vti files are linked in a pvd file
  void write(int iT=0);
  ///  writes functor instantaneously, same vti-pvd file structure as above
  void write(SuperF2D<T>& f, int iT=0);
  ///  have to be called before writting added functors, in order to create
  ///  master.pvd file
  void createMasterFile();
  ///  put functor to _pointerVec
  ///  to simplify writing process of several functors
  void addFunctor(SuperF2D<T>& f);
  ///  to clear stored functors
  void clearAddedFunctors();
private:
  ///  performes <VTKFile ...>, <ImageData ...>, <PieceExtent ...> and <PointData ...>
  void preambleVTI(const std::string& fullName, int x0, int y0, int x1, int y1,
                   T originX,T originY, T delta);
  ///  performes </ImageData> and </VTKFile>
  void closeVTI(const std::string& fullNamePiece);
  ///  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  ///  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);
  ///  performes <DataSet timestep= ... file=namePiece />
  ///  *** nasty function ***
  void dataPVD(int iT, int iC, const std::string& fullNamePVD, const std::string& namePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVDmaster(int iT, int iC, const std::string& fullNamePVDMaster,
                     const std::string& namePiece);
  ///  writes given functor f
  void dataArray(const std::string& fullName, SuperF2D<T>& f,
                 int iC, int nx, int ny);
  ///  writes given functor f, base64
  void dataArrayBinary(const std::string& fullName, SuperF2D<T>& f,
                       int iC, int nx, int ny);
  ///  performes </PointData> and </Piece>
  void closePiece(const std::string& fullNamePiece);
private:
  mutable OstreamManager clout;
  ///  default is false, call createMasterFile() and it will be true
  bool _createFile;
  ///  determines the name of .vti and .pvd per iT
  std::string _name;
  ///  default is true, may be changed at constructor
  bool _binary;
  ///  holds added functor, to simplify the use of write function
  std::vector< SuperF2D<T>* > _pointerVec;
};

}  // namespace olb


#endif

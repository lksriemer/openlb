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
 * (only for uniform grids) -- generic implementation.
 */

#ifndef SUPER_VTK_WRITER_2D_HH
#define SUPER_VTK_WRITER_2D_HH

#include <fstream>
#include <iostream>
#include <string>
#include "core/singleton.h"
#include "communication/loadBalancer.h"
#include "geometry/cuboidGeometry2D.h"
#include "communication/mpiManager.h"
#include "io/superVtkWriter2D.h"
#include "io/base64.h"
#include "io/fileName.h"

namespace olb {

template<typename T>
SuperVTKwriter2D<T>::SuperVTKwriter2D ( std::string name, bool binary )
  : clout( std::cout,"SuperVTKwriter2D" ), _createFile(false), _name(name), _binary(binary)
{ }

template<typename T>
void SuperVTKwriter2D<T>::write(int iT)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  //  !!!!!!!!!!! check whether _pointerVec is empty
  if ( _pointerVec.empty() ) {
    clout << "Error: Did you add a Functor ?";
  } else {
    // DO WORK
    //  auto it = _pointerVec.begin();
    // to get first element _pointerVec
    // problem if functors with different SuperLattice are stored
    // since till now, there is only one origin
    auto it = _pointerVec.cbegin();
    CuboidGeometry2D<T> const& cGeometry = (**it).getSuperStructure().getCuboidGeometry();
    // no gaps between vti files (cuboids)
    (**it).getSuperStructure().communicate();
    LoadBalancer<T>& load = (**it).getSuperStructure().getLoadBalancer();

    //-------------------------------------
    // PVD
    // master.pvd owns all
    //-------------------------------------
    if ( rank == 0 ) {
      // master only
      std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                      + createFileName( _name ) + ".pvd";
      std::string fullNamePVD = singleton::directories().getVtkOutDir()
                                + "data/" + createFileName( _name, iT ) + ".pvd";

      preamblePVD(fullNamePVD);         // timestep
      for (int iC = 0; iC < cGeometry.getNc(); iC++) {
        // cuboid
        std::string namePieceData = "data/" + createFileName( _name, iT, iC) + ".vti";
        std::string namePiece = createFileName( _name, iT, iC) + ".vti";
        // puts name of .vti piece to a .pvd file [fullNamePVD]
        dataPVD( iT, iC, fullNamePVD, namePiece );
        // adds a namePiece to master.pvd file.
        // To do so we overwrite closePVD() and add new entry.
        dataPVDmaster( iT, iC, fullNamePVDmaster, namePieceData );
      } // cuboid
      closePVD(fullNamePVD);            // timestep
    } // master only
    //-------------------------------------

    //-------------------------------------
    // VTI
    // each process writes his cuboids
    //-------------------------------------
    int originLatticeR[3] = {int()};
    for (int iCloc = 0; iCloc < load.size(); iCloc++) {
      // cuboid
      int nx = cGeometry.get(load.glob(iCloc)).getNx();
      int ny = cGeometry.get(load.glob(iCloc)).getNy();
      // to be changed into the following line once local refinement has been implemented
      // double deltaX = cGeometry.get(load.glob(iCloc)).getDeltaR();
      T delta = cGeometry.getMotherCuboid().getDeltaR();

      std::string fullNameVTI = singleton::directories().getVtkOutDir() + "data/"
                                + createFileName( _name, iT, load.glob(iCloc) ) + ".vti";

      originLatticeR[0] = load.glob(iCloc);
      T originPhysR[2] = {T()};
      cGeometry.getPhysR(originPhysR,originLatticeR);

      // blow up by 1, to ensure overlap at output.
      // no more gaps in paraview!
      preambleVTI(fullNameVTI, -1,-1,nx,ny, originPhysR[0],originPhysR[1], delta);
      if ( _binary ) {
        // iterate the stored functors and write them
        for ( auto it = _pointerVec.cbegin(); it != _pointerVec.cend(); ++it) {
          // write <DataArray ..../> for functor, e.g. velocity, pressure, ...
          dataArrayBinary(fullNameVTI, (**it), load.glob(iCloc), nx,ny);
        }
      } else {

        // iterate the stored functors and write them
        for ( auto it = _pointerVec.cbegin(); it != _pointerVec.cend(); ++it) {
          // write <DataArray ..../> for functor, e.g. velocity, pressure, ...
          dataArray(fullNameVTI, (**it), load.glob(iCloc), nx,ny);
        }
      }
      closePiece(fullNameVTI);
      closeVTI(fullNameVTI);
    } // cuboid
  } // END WORK
  ////-------------------------------------
}

template<typename T>
void SuperVTKwriter2D<T>::write(SuperF2D<T>& f, int iT)
{
  CuboidGeometry2D<T> const& cGeometry = f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = f.getSuperStructure().getLoadBalancer();
  f.getSuperStructure().communicate();
  T delta   = cGeometry.getMotherCuboid().getDeltaR();

  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  //  T originX = 0; //cGeometry.getMotherCuboid().get_globPosX();
  //  T originY = 0; //cGeometry.getMotherCuboid().get_globPosY();

  // write a pvd file, which links all vti files
  // each vti file is written by one thread, which may own severals cuboids
  if ( rank == 0 ) {
    // master only
    std::string namePVD = singleton::directories().getVtkOutDir() +
                          createFileName( f.getName(), iT )  + ".pvd";

    preamblePVD(namePVD);
    for (int iC = 0; iC < cGeometry.getNc(); iC++) {
      std::string nameVTI =  "data/" + createFileName( f.getName(), iT, iC)
                             + ".vti";
      // puts name of .vti piece to a .pvd file [fullNamePVD]
      dataPVD( iT, iC, namePVD, nameVTI );
    }
    closePVD(namePVD);
  } // master only

  int originLatticeR[3] = {int()};
  for (int iCloc = 0; iCloc < load.size(); iCloc++) {
    int nx = cGeometry.get(load.glob(iCloc)).getNx();
    int ny = cGeometry.get(load.glob(iCloc)).getNy();
    // get dimension/extent for each cuboid
    originLatticeR[0] = load.glob(iCloc);
    T originPhysR[2] = {T()};
    cGeometry.getPhysR(originPhysR,originLatticeR);

    std::string fullNameVTI = singleton::directories().getVtkOutDir() + "data/"
                              + createFileName( f.getName(), iT, load.glob(iCloc) ) + ".vti";

    preambleVTI(fullNameVTI, -1,-1,nx,ny,originPhysR[0],originPhysR[1], delta);
    if ( _binary ) {
      dataArrayBinary(fullNameVTI, f, load.glob(iCloc), nx,ny);
    } else {
      dataArray(fullNameVTI, f, load.glob(iCloc), nx,ny);
    }
    closePiece(fullNameVTI);
    closeVTI(fullNameVTI);
  }
}

template<typename T>
void SuperVTKwriter2D<T>::createMasterFile()
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if ( rank == 0 ) {
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir()
                                    + createFileName( _name ) + ".pvd";
    preamblePVD(fullNamePVDmaster);
    closePVD(fullNamePVDmaster);
    _createFile = true;
  }
}

template<typename T>
void SuperVTKwriter2D<T>::addFunctor(SuperF2D<T>& f)
{
  _pointerVec.push_back(&f);
}

template<typename T>
void SuperVTKwriter2D<T>::clearAddedFunctors()
{
  _pointerVec.clear();
}




////////////////////private member functions///////////////////////////////////
template<typename T>
void SuperVTKwriter2D<T>::preambleVTI (const std::string& fullName,
                                       int x0,int y0,int x1,int y1, T originX,T originY, T delta)
{
  std::ofstream fout(fullName.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n";
  fout << "<ImageData WholeExtent=\""
       << x0 <<" "<< x1 <<" "
       << y0 <<" "<< y1 <<" "
       << 0 <<" "<< 0
       << "\" Origin=\"" << originX << " " << originY << " " << "0"
       << "\" Spacing=\"" << delta << " " << delta << " " << delta << "\">\n";
  fout << "<Piece Extent=\""
       << x0 <<" "<< x1 <<" "
       << y0 <<" "<< y1 <<" "
       << 0 <<" "<< 0 <<"\">\n";
  fout << "<PointData>\n";
  fout.close();
}

template<typename T>
void SuperVTKwriter2D<T>::closeVTI(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</ImageData>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T>
void SuperVTKwriter2D<T>::preamblePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }

  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n"
       << "<Collection>\n";
  fout.close();
}

template<typename T>
void SuperVTKwriter2D<T>::closePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T>
void SuperVTKwriter2D<T>::dataPVD(int iT, int iC,
                                  const std::string& fullNamePVD, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNamePVD << std::endl;
  }

  fout << "<DataSet timestep=\"" << iT << "\" "
       << "group=\"\" part=\" " <<  iC << "\" "
       << "file=\"" << namePiece << "\"/>\n";
  fout.close();
}

template<typename T>
void SuperVTKwriter2D<T>::dataPVDmaster(int iT, int iC,
                                        const std::string& fullNamePVDMaster, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVDMaster.c_str(), std::ios::in | std::ios::out | std::ios::ate);
  if (fout) {
    fout.seekp(-25,std::ios::end);    // jump -25 form the end of file to overwrite closePVD

    fout << "<DataSet timestep=\"" << iT << "\" "
         << "group=\"\" part=\" " <<  iC << "\" "
         << "file=\"" << namePiece << "\"/>\n";
    fout.close();
    closePVD(fullNamePVDMaster);
  } else {
    clout << "Error: could not open " << fullNamePVDMaster << std::endl;
  }
}

template<typename T>
void SuperVTKwriter2D<T>::dataArray(const std::string& fullName,
                                    SuperF2D<T>& f, int iC, int nx, int ny)
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  fout << "<DataArray " ;
  if (f.getTargetDim() == 1) {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\">\n";
  } else {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "NumberOfComponents=\"" << f.getTargetDim() << "\">\n";
  }

  // since cuboid has been blowed up by 1 [every dimension]
  // looping from -1 to ny (= ny+1, as passed)
  int i[3];
  i[0] = iC;
  T evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = T();
  }
  for (i[2] = -1; i[2] < ny+1; ++i[2]) {
    for (i[1] = -1; i[1] < nx+1; ++i[1]) {
      f(evaluated,i);
      for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
        fout << evaluated[iDim] << " ";
      }
    }
  }
  fout << "\n</DataArray>\n";

  fout.close();
}

template<typename T>
void SuperVTKwriter2D<T>::dataArrayBinary(const std::string& fullName,
    SuperF2D<T>& f, int iC, int nx, int ny)
{
  const char* fileName = fullName.c_str();

  std::ofstream fout( fileName, std::ios::out | std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fileName << std::endl;
  }

  fout << "<DataArray ";
  if (f.getTargetDim() == 1) {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "format=\"binary\" encoding=\"base64\">";
  } else {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "format=\"binary\" encoding=\"base64\" "
         << "NumberOfComponents=\"" << f.getTargetDim() << "\">\n";
  }
  fout.close();

  std::ofstream* ofstr;   // only used for binary output // passed to Base64Encoder
  ofstr = new std::ofstream( fileName, std::ios::out | std::ios::app | std::ios::binary );
  if (!ofstr) {
    clout << "Error: could not open " << fileName << std::endl;
  }

  size_t fullSize = f.getTargetDim() * (nx+2) * (ny+2);    //  how many numbers to write
  size_t binarySize = size_t( fullSize*sizeof(float) );
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(*ofstr, 1);
  unsigned int uintBinarySize = (unsigned int)binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<float>* dataEncoder = new Base64Encoder<float>( *ofstr, fullSize );

  int i[3];
  i[0] = iC;
  T evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = T();
  }
  for (i[2] = -1; i[2] < ny+1; ++i[2]) {
    for (i[1] = -1; i[1] < nx+1; ++i[1]) {
      f(evaluated,i);
      for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
        const float tmpC = float( evaluated[iDim] );
        dataEncoder->encode( &tmpC, 1 );
      }
    }
  }
  ofstr->close();
  delete ofstr;

  std::ofstream ffout( fileName,  std::ios::out | std::ios::app );
  if (!ffout) {
    clout << "Error: could not open " << fileName << std::endl;
  }
  ffout << "\n</DataArray>\n";
  ffout.close();
  delete dataEncoder;
}

template<typename T>
void SuperVTKwriter2D<T>::closePiece(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullNamePiece << std::endl;
  }
  fout << "</PointData>\n";
  fout << "</Piece>\n";
  fout.close();
}

}  // namespace olb

#endif

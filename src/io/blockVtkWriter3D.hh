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
 * A method to write vtk data for cuboid geometries
 * (only for uniform grids) -- generic implementation.
 */

#ifndef BLOCK_VTK_WRITER_3D_HH
#define BLOCK_VTK_WRITER_3D_HH

#include <fstream>
#include <iostream>
#include "communication/mpiManager.h"
#include "core/singleton.h"
#include "io/imageWriter.h"
#include "io/blockVtkWriter3D.h"
#include "io/base64.h"

namespace olb {

template<typename T, template<typename U> class DESCRIPTOR>
class BlockLatticeF3D;


////////// class VTKwriter3D ////////////////////////////////////////

template<typename T, template <typename U> class DESCRIPTOR>
BlockVTKwriter3D<T,DESCRIPTOR>::BlockVTKwriter3D
  ( std::string name, bool binary ) : clout( std::cout,"BlockVTKwriter3D" )
{ 
  _name = name;
  _binary = binary;
}

template<typename T, template <typename U> class DESCRIPTOR>
BlockVTKwriter3D<T,DESCRIPTOR>::~BlockVTKwriter3D ()
{ clearAddedFunctors(); }

//  iteration on _pointerVec is realized by function
//  dataArray() respective dataArrayBinary()
template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::write(int iT)
{  
  if( _pointerVec.empty() ) { 
    // secure code. doesn't appear on console ??
    clout << "Error: Please add functor via addFunctor()";
  } else {
    //  auto it = _pointerVec.begin();
    typename std::vector<BlockLatticeIdentity3D<T,DESCRIPTOR>* >::iterator it;
    it = _pointerVec.begin();
    BlockLatticeIdentity3D<T,DESCRIPTOR> f( **it );

    T originX = 0;
    T originY = 0;
    T originZ = 0;
    int nx = f.getBlockLattice3D().getNx() -1;
    int ny = f.getBlockLattice3D().getNy() -1;
    int nz = f.getBlockLattice3D().getNz() -1;

    std::string fullNameVti = singleton::directories().getVtkOutDir()
                              + graphics::createFileName( f.getName(), iT ) + ".vti";

    preamble( fullNameVti, nx,ny,nz, originX,originY,originZ );
    pointData( fullNameVti );
    if( _binary ) {
      dataArrayBinary( fullNameVti, nx,ny,nz);
    } else { 
      dataArray( fullNameVti, nx,ny,nz);
    }
    closePointData( fullNameVti );
    closePreamble( fullNameVti );
  }  
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::write(BlockLatticeF3D<T,DESCRIPTOR>& f, int iT)
{
  T originX = 0;
  T originY = 0;
  T originZ = 0;
  int nx = f.getBlockLattice3D().getNx() -1;
  int ny = f.getBlockLattice3D().getNy() -1;
  int nz = f.getBlockLattice3D().getNz() -1;

  std::string fullNameVti = singleton::directories().getVtkOutDir()
                            + graphics::createFileName( f.getName(), iT ) + ".vti";

  preamble( fullNameVti, nx,ny,nz, originX,originY,originZ );
  pointData( fullNameVti );
  if( _binary ) {
    dataFunctorToOneFileBinary( fullNameVti, f, nx,ny,nz );
  } else { 
    dataFunctorToOneFile( fullNameVti, f, nx,ny,nz );
  }

  closePointData( fullNameVti );
  closePreamble( fullNameVti );
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::addFunctor(BlockLatticeF3D<T,DESCRIPTOR>& f)
{  
  BlockLatticeIdentity3D<T,DESCRIPTOR>* tmp = new BlockLatticeIdentity3D<T,DESCRIPTOR>(f);
  _pointerVec.push_back(tmp); 
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::clearAddedFunctors()
{ 
  typename std::vector< BlockLatticeIdentity3D<T,DESCRIPTOR>* >::iterator it;
  for( it = _pointerVec.begin(); it != _pointerVec.end(); it++)
  {
    delete *it;
  }
  _pointerVec.clear();
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::preamble(std::string& fullName,
  int nx,int ny,int nz, T originX, T originY, T originZ)
{
  if(singleton::mpi().getRank()==0)
  {
    std::ofstream fout(fullName.c_str());
    if (!fout) clout << "Error: could not open " << fullName << std::endl;

    T spacing = double(1.0/nx);

    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
         << "byte_order=\"LittleEndian\">\n";
    fout << "<ImageData WholeExtent=\"" 
            << "0" <<" "<< nx <<" "
            << "0" <<" "<< ny <<" "
            << "0" <<" "<< nz
         << "\" Origin=\"" << originX << " " << originY << " " << originZ
         << "\" Spacing=\"" << spacing << " " << spacing << " " << spacing << "\">\n";

    fout << "<Piece Extent=\""
         << 0 <<" "<< nx <<" "
         << 0 <<" "<< ny <<" "
         << 0 <<" "<< nz <<"\">\n";

    fout.close();
  }
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::closePreamble(std::string& fullNamePiece)
{
  if(singleton::mpi().getRank()==0)
  {
    std::ofstream fout(fullNamePiece.c_str(), std::ios::app );
    if (!fout) clout << "Error: could not open " << fullNamePiece << std::endl;
    fout << "</Piece>\n";
    fout << "</ImageData>\n";
    fout << "</VTKFile>\n";
    fout.close();
  }
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::pointData(std::string& fullName)
{ 
  if(singleton::mpi().getRank()==0)
  {
    std::ofstream fout(fullName.c_str(), std::ios::app);
    if (!fout) clout << "Error: could not open " << fullName << std::endl;
  
    fout << "<PointData>\n";
    fout.close();
  }
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::closePointData(std::string& fullName)
{
  if(singleton::mpi().getRank()==0)
  {
    std::ofstream fout(fullName.c_str(), std::ios::app);
    if (!fout) clout << "Error: could not open " << fullName << std::endl;
    fout << "</PointData>\n";
    fout.close();
  }
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::dataArray(std::string& fullNameVti,
  int nx, int ny, int nz)
{
  std::ofstream fout(fullNameVti.c_str(), std::ios::app);
  if (!fout) clout << "Error: could not open " << fullNameVti << std::endl;

//  for( auto it = _pointerVec.begin(); it != _pointerVec.end(); it++)
  typename std::vector< BlockLatticeIdentity3D<T,DESCRIPTOR>* >::iterator it;
  for( it = _pointerVec.begin(); it != _pointerVec.end(); it++)
  { //  functors
    BlockLatticeIdentity3D<T,DESCRIPTOR> f( **it );
    int dim = f.getTargetDim();

    if(singleton::mpi().getRank()==0) {
      fout << "<DataArray " ;
      if (dim == 1) {
        fout << "type=\"Float32\" Name=\"" << f.getName() <<"\">\n";
      } else {
        fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
             << "NumberOfComponents=\"" << dim <<"\">\n";
      }
    }

    std::vector<int> tmpVec( 3,int(0) );
    for (int iZ=0; iZ<nz+1; ++iZ) {
      for (int iY=0; iY<ny+1; ++iY) {
        for (int iX=0; iX<nx+1; ++iX) {
          for (int iDim=0; iDim<dim; ++iDim) {
            tmpVec[0] = iX; tmpVec[1] = iY; tmpVec[2] = iZ;
            T evaluated = f(tmpVec)[iDim]; 
            if(singleton::mpi().getRank()==0) fout << evaluated << " ";
          }
        }
      }
    }

    if(singleton::mpi().getRank()==0) {
      fout << "\n";
      fout << "</DataArray>\n";
    }
  } //  functors
  fout.close();
}

//  uses base64 encoder to write binary output
//  first number is written by a seperate sizeEncoder
//  this number indicates how many numbers will be stored.
//  then dataEncoder will be called to write output.
//  !!code is fixed to float functor values!!
template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::dataArrayBinary(std::string& fullNameVti,
  int nx, int ny, int nz)
{
  const char* fileName = fullNameVti.c_str();

  typename std::vector< BlockLatticeIdentity3D<T,DESCRIPTOR>* >::iterator it;
  for( it = _pointerVec.begin(); it != _pointerVec.end(); it++)
  { //  functors
    BlockLatticeIdentity3D<T,DESCRIPTOR> f( **it );
    int dim = f.getTargetDim();

    if(singleton::mpi().getRank()==0) {
      std::ofstream fout(fileName,  std::ios::out | std::ios::app);
      if (!fout) clout << "Error: could not open " << fileName << std::endl;

      fout << "<DataArray ";
      if (dim == 1) {
        fout << "type=\"Float32\" Name=\"" << f.getName() << "\" " 
             << "format=\"binary\" encoding=\"base64\">\n";
      } else {
        fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
             << "format=\"binary\" encoding=\"base64\" "
             << "NumberOfComponents=\"" << dim <<"\">\n";
      }
      fout.close();
    }

    std::ofstream* ofstr;   // only used for binary output // passed to Base64Encoder
    ofstr = new std::ofstream( fileName, std::ios::out | std::ios::app | std::ios::binary );
    if (!ofstr) clout << "Error: could not open " << fileName << std::endl;

    size_t fullSize = dim * (1 + nx) * (1 + ny) * (1 + nz); //  take care about DIMENSION
    size_t binarySize = size_t( fullSize * sizeof(float) );
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(*ofstr, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float>* dataEncoder = 0;
    dataEncoder = new Base64Encoder<float>( *ofstr, fullSize );
    std::vector<int> tmpVec( 3,int(0) );
    for (int iZ=0; iZ<nz+1; ++iZ) {
      for (int iY=0; iY<ny+1; ++iY) {
        for (int iX=0; iX<nx+1; ++iX) {
          for (int iDim=0; iDim<dim; ++iDim) {
            tmpVec[0] = iX; tmpVec[1] = iY; tmpVec[2] = iZ;
            const float evaluated = float( f(tmpVec)[iDim] );
            if(singleton::mpi().getRank()==0) dataEncoder->encode( &evaluated, 1 );
          }
        }
      }
    }
    ofstr->close();
    
    if(singleton::mpi().getRank()==0) {
      std::ofstream fout(fileName,  std::ios::out | std::ios::app);
      if (!fout) clout << "Error: could not open " << fileName << std::endl;
      fout << "</DataArray>\n";
      fout.close();
    }
  }//  functors
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::dataFunctorToOneFile
    (std::string& fullNameVti, BlockLatticeF3D<T,DESCRIPTOR>& f,
     int nx, int ny, int nz)
{
  std::ofstream fout(fullNameVti.c_str(), std::ios::app);
  if (!fout) clout << "Error: could not open " << fullNameVti << std::endl;

  int dim = f.getTargetDim();

  if(singleton::mpi().getRank()==0) {
    fout << "<DataArray " ;
    if (dim == 1) {
      fout << "type=\"Float32\" Name=\"" << f.getName() <<"\">\n";
    } else {
      fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
           << "NumberOfComponents=\"" << dim <<"\">\n";
    }
  }

  std::vector<int> tmpVec( 3,int(0) );
  for (int iZ=0; iZ<nz+1; ++iZ) {
    for (int iY=0; iY<ny+1; ++iY) {
      for (int iX=0; iX<nx+1; ++iX) {
        for (int iDim=0; iDim<dim; ++iDim) {
          tmpVec[0] = iX; tmpVec[1] = iY; tmpVec[2] = iZ;
          T evaluated = f(tmpVec)[iDim]; 
          if(singleton::mpi().getRank()==0) fout << evaluated << " ";
        }
      }
    }
  }
  if(singleton::mpi().getRank()==0) {
    fout << "\n";
    fout << "</DataArray>\n";
  }

  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter3D<T,DESCRIPTOR>::dataFunctorToOneFileBinary
    (std::string& fullNameVti, BlockLatticeF3D<T,DESCRIPTOR>& f,
     int nx, int ny, int nz)
{
  const char* fileName = fullNameVti.c_str();
  std::ofstream fout(fileName, std::ios::app);
  if (!fout) clout << "Error: could not open " << fileName << std::endl;

  int dim = f.getTargetDim();

  if(singleton::mpi().getRank()==0) {
    fout << "<DataArray " ;
    if (dim == 1) {
      fout << "type=\"Float32\" Name=\"" << f.getName() << "\" " 
           << "format=\"binary\" encoding=\"base64\">\n";
    }
    if (dim != 1) {
      fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
           << "format=\"binary\" encoding=\"base64\" "
           << "NumberOfComponents=\"" << dim <<"\">\n";
    }
  }
  fout.close();

  std::ofstream* ofstr;   // only used for binary output // passed to Base64Encoder
  ofstr = new std::ofstream( fileName, std::ios::out | std::ios::app | std::ios::binary );
  if (!ofstr) clout << "Error: could not open " << fileName << std::endl;

  size_t fullSize = dim * (1 + nx) * (1 + ny) * (1 + nz);
  size_t binarySize = size_t( fullSize * sizeof(float) );
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(*ofstr, 1);
  unsigned int uintBinarySize = (unsigned int)binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<float>* dataEncoder = 0;
  dataEncoder = new Base64Encoder<float>( *ofstr, fullSize );
  std::vector<int> tmpVec( 3,int(0) );
  for (int iZ=0; iZ<nz+1; ++iZ) {
    for (int iY=0; iY<ny+1; ++iY) {
      for (int iX=0; iX<nx+1; ++iX) {
        for (int iDim=0; iDim<dim; ++iDim) {
          tmpVec[0] = iX; tmpVec[1] = iY; tmpVec[2] = iZ;
          const float evaluated = float( f(tmpVec)[iDim] );
          if(singleton::mpi().getRank()==0) dataEncoder->encode( &evaluated, 1 );
        }
      }
    }
  }
  ofstr->close();

  if(singleton::mpi().getRank()==0) {
    std::ofstream fout(fileName,  std::ios::out | std::ios::app);
    if (!fout) clout << "Error: could not open " << fileName << std::endl;
    fout << "\n";
    fout << "</DataArray>\n";
    fout.close();
  }
}



}  // namespace olb

#endif

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

#ifndef BLOCK_VTK_WRITER_2D_HH
#define BLOCK_VTK_WRITER_2D_HH

#include <fstream>
#include <iostream>
#include "communication/mpiManager.h"
#include "core/singleton.h"
#include "io/imageWriter.h"
#include "io/blockVtkWriter2D.h"
#include "io/base64.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, template<typename U> class DESCRIPTOR>
class BlockLatticeF2D;


////////// class VTKwriter2D ////////////////////////////////////////

template<typename T, template <typename U> class DESCRIPTOR>
BlockVTKwriter2D<T,DESCRIPTOR>::BlockVTKwriter2D
  ( std::string name, bool binary ) : clout( std::cout,"BlockVTKwriter2D" )
{ 
  _name = name;
  _binary = binary;
}

template<typename T, template <typename U> class DESCRIPTOR>
BlockVTKwriter2D<T,DESCRIPTOR>::~BlockVTKwriter2D ()
{ clearAddedFunctors(); }

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter2D<T,DESCRIPTOR>::addFunctor(BlockLatticeF2D<T,DESCRIPTOR>& f)
{  
  BlockLatticeIdentity2D<T,DESCRIPTOR>* tmp = new BlockLatticeIdentity2D<T,DESCRIPTOR>(f);
  _pointerVec.push_back(tmp); 
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter2D<T,DESCRIPTOR>::clearAddedFunctors()
{ 
  typename std::vector< BlockLatticeIdentity2D<T,DESCRIPTOR>* >::iterator it;
  for( it = _pointerVec.begin(); it != _pointerVec.end(); it++)
  {
    delete *it;
  }
  _pointerVec.clear();
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter2D<T,DESCRIPTOR>::preamble(std::string& fullName,
  int nx,int ny, T originX, T originY, T originZ)
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
            << 0 <<" "<< nx <<" "
            << 0 <<" "<< ny <<" "
            << 0 <<" "<< 0
         << "\" Origin=\"" << originX << " " << originY << " " << originZ
         << "\" Spacing=\"" << spacing << " " << spacing << " " << spacing << "\">\n";

    fout << "<Piece Extent=\""
         << 0 <<" "<< nx <<" "
         << 0 <<" "<< ny <<" "
         << 0 <<" "<< 0 <<"\">\n";

    fout.close();
  }
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter2D<T,DESCRIPTOR>::closePreamble(std::string& fullNamePiece)
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
void BlockVTKwriter2D<T,DESCRIPTOR>::pointData(std::string& fullName)
{ 
  if(singleton::mpi().getRank()==0)
  {
    std::ofstream fout(fullName.c_str(), std::ios::app);
    if (!fout) clout << "Error: could not open " << fullName << std::endl;
  
    fout << "<PointData>\n";
  }
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter2D<T,DESCRIPTOR>::closePointData(std::string& fullName)
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
void BlockVTKwriter2D<T,DESCRIPTOR>::dataArray(std::string& fullNameVti,
  int nx, int ny)
{
  std::ofstream fout(fullNameVti.c_str(), std::ios::app);
  if (!fout) clout << "Error: could not open " << fullNameVti << std::endl;

//  for( auto it = _pointerVec.begin(); it != _pointerVec.end(); it++)
  typename std::vector< BlockLatticeIdentity2D<T,DESCRIPTOR>* >::iterator it;
  for( it = _pointerVec.begin(); it != _pointerVec.end(); it++)
  { //  functors
    BlockLatticeIdentity2D<T,DESCRIPTOR> f( **it );
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

    std::vector<int> tmpVec( 2,int(0) );
    for (int iY=0; iY<ny+1; ++iY) {
      for (int iX=0; iX<nx+1; ++iX) {
        for (int iDim=0; iDim<dim; ++iDim) {
          tmpVec[0] = iX; tmpVec[1] = iY;
          T evaluated = f(tmpVec)[iDim]; 
          if(singleton::mpi().getRank()==0) fout << evaluated << " ";
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
void BlockVTKwriter2D<T,DESCRIPTOR>::dataArrayBinary(std::string& fullNameVti,
  int nx, int ny)
{
  const char* fileName = fullNameVti.c_str();

  typename std::vector< BlockLatticeIdentity2D<T,DESCRIPTOR>* >::iterator it;
  for( it = _pointerVec.begin(); it != _pointerVec.end(); it++)
  { //  functors
    BlockLatticeIdentity2D<T,DESCRIPTOR> f( **it );
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

    size_t fullSize = dim * (1 + nx) * (1 + ny) * (1); //  take care about DIMENSION
    size_t binarySize = size_t( fullSize * sizeof(float) );
    // writes first number, which have to be the size(byte) of the following data
    Base64Encoder<unsigned int> sizeEncoder(*ofstr, 1);
    unsigned int uintBinarySize = (unsigned int)binarySize;
    sizeEncoder.encode(&uintBinarySize, 1);
    //  write numbers from functor
    Base64Encoder<float>* dataEncoder = 0;
    dataEncoder = new Base64Encoder<float>( *ofstr, fullSize );
    
    std::vector<int> tmpVec( 2,int(0) );
    for (int iY=0; iY<ny+1; ++iY) {
      for (int iX=0; iX<nx+1; ++iX) {
        for (int iDim=0; iDim<dim; ++iDim) {
          tmpVec[0] = iX; tmpVec[1] = iY;
          const float evaluated = float( f(tmpVec)[iDim] );
          if(singleton::mpi().getRank()==0) dataEncoder->encode( &evaluated, 1 );
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
void BlockVTKwriter2D<T,DESCRIPTOR>::dataFunctorToOneFile
    (std::string& fullNameVti, BlockLatticeF2D<T,DESCRIPTOR>& f,
     int nx, int ny)
{
  BlockLatticeIdentity2D<T,DESCRIPTOR> tmp(f);
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

  std::vector<int> tmpVec( 2,int(0) );
  for (int iY=0; iY<ny+1; ++iY) {
    for (int iX=0; iX<nx+1; ++iX) {
      for (int iDim=0; iDim<dim; ++iDim) {
        tmpVec[0] = iX; tmpVec[1] = iY;
        T evaluated = f(tmpVec)[iDim];
          if(singleton::mpi().getRank()==0) fout << evaluated << " ";
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
void BlockVTKwriter2D<T,DESCRIPTOR>::dataFunctorToOneFileBinary
    (std::string& fullNameVti, BlockLatticeF2D<T,DESCRIPTOR>& f,
     int nx, int ny)
{
  BlockLatticeIdentity2D<T,DESCRIPTOR> tmp(f);
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

  size_t fullSize = dim * (1 + nx) * (1 + ny) * (1);
//  size_t fullSize = 3 * (1 + nx) * (1 + ny) ;
  size_t binarySize = size_t( fullSize * sizeof(float) );
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(*ofstr, 1);
  unsigned int uintBinarySize = (unsigned int)binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<float>* dataEncoder = 0;
  dataEncoder = new Base64Encoder<float>( *ofstr, fullSize );
  
  std::vector<int> tmpVec( 2,int(0) );
  for (int iY=0; iY<ny+1; ++iY) {
    for (int iX=0; iX<nx+1; ++iX) {
      for (int iDim=0; iDim<dim; ++iDim) {
        tmpVec[0] = iX; tmpVec[1] = iY;
        const float evaluated = float( f(tmpVec)[iDim] );
        if(singleton::mpi().getRank()==0) dataEncoder->encode( &evaluated, 1 );
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

//  iteration on _pointerVec is realized by function
//  dataArray() respective dataArrayBinary()
template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter2D<T,DESCRIPTOR>::write(int iT)
{  
  if( _pointerVec.empty() ) { 
    // secure code. doesn't appear on console ??
    clout << "Error: Please add functor via addFunctor()";
  } else {
    //  auto it = _pointerVec.begin();
    typename std::vector<BlockLatticeIdentity2D<T,DESCRIPTOR>* >::iterator it;
    it = _pointerVec.begin();
    BlockLatticeIdentity2D<T,DESCRIPTOR> f( **it );

    T originX = 0;
    T originY = 0;
    T originZ = 0;
    int nx = f.getBlockLattice2D().getNx() -1;
    int ny = f.getBlockLattice2D().getNy() -1;

    std::string fullNameVti = singleton::directories().getVtkOutDir()
                              + graphics::createFileName( _name, iT ) + ".vti";

    preamble( fullNameVti, nx,ny, originX,originY,originZ );
    pointData( fullNameVti );
    if( _binary )  dataArrayBinary( fullNameVti, nx,ny);
    if( !_binary ) dataArray( fullNameVti, nx,ny);
    closePointData( fullNameVti );
    closePreamble( fullNameVti );
  }  
}

template<typename T, template <typename U> class DESCRIPTOR>
void BlockVTKwriter2D<T,DESCRIPTOR>::write(BlockLatticeF2D<T,DESCRIPTOR>& f, int iT)
{
  T originX = 0;
  T originY = 0;
  T originZ = 0;
  int nx = f.getBlockLattice2D().getNx() -1;
  int ny = f.getBlockLattice2D().getNy() -1;

  std::string fullNameVti = singleton::directories().getVtkOutDir()
                              + graphics::createFileName( f.getName(), iT ) + ".vti";
//                              + graphics::createFileName( _name, iT ) + ".vti";

  preamble( fullNameVti, nx,ny, originX,originY,originZ );
  pointData( fullNameVti );
  if( _binary )  dataFunctorToOneFileBinary( fullNameVti, f, nx,ny );
  if( !_binary ) dataFunctorToOneFile( fullNameVti, f, nx,ny );
  closePointData( fullNameVti );
  closePreamble( fullNameVti );
}

}  // namespace olb

#endif

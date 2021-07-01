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

#ifndef SUPER_VTK_WRITER_3D_HH
#define SUPER_VTK_WRITER_3D_HH

#include <fstream>
//#include <cstdlib>    // for creating folders
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "core/singleton.h"
#include "communication/loadBalancer.h"
#include "geometry/cuboidGeometry3D.h"
#include "communication/mpiManager.h"
#include "functors/superLatticeCalcF3D.h" // for IdentityF
#include "io/imageWriter.h"
#include "io/superVtkWriter3D.h"
#include "io/base64.h"

using namespace std;
namespace olb {

template<typename T, template<typename U> class DESCRIPTOR>
class SuperLatticeF3D;


//////////////// class SuperVTKwriter3D ////////////////////////////////////////

//////////////////////public member functions///////////////////////////////////
template<typename T, template <typename U> class DESCRIPTOR>
SuperVTKwriter3D<T,DESCRIPTOR>::SuperVTKwriter3D
  ( std::string name, bool binary ) : clout( std::cout,"SuperVTKwriter3D" )
{ 
  _name = name;
  _binary = binary;
  _createFile = false;  // true if createMasterFile() has been called

//  std::string consoleCommand = std::string("mkdir -p ") + std::string("tmp/data");
//  if( !system(consoleCommand.c_str()) )
//  {
//    std::cout << "smth went wrong, exists tmp/data ??" << std::endl;
//  }
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::write(int iT)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  //  !!!!!!!!!!! check whether _pointerVec is empty
  if( _pointerVec.empty() ) {
    clout << "Error: Did you add a Functor ?";
  } else
  { // DO WORK 
    //  auto it = _pointerVec.begin();
    // to get first element _pointerVec
    // problem if functors with different SuperLattice are stored
    // since till now, there is only one origin
    typename std::vector<SuperLatticeF3D<T,DESCRIPTOR>* >::iterator it;
    it = _pointerVec.begin();
    SuperLatticeIdentity3D<T,DESCRIPTOR> f( **it );
    CuboidGeometry3D<T> const& cGeometry = f.getSuperLattice3D().getCuboidGeometry();
    // no gaps between vti files (cuboids)
    f.getSuperLattice3D().communicate();
    LoadBalancer<T>& load = f.getSuperLattice3D().getLoadBalancer();

    //-------------------------------------
    // PVD
    // master.pvd owns all
    //-------------------------------------
    if( rank == 0 )
    { // master only
      std::string fullNamePVDmaster = singleton::directories().getVtkOutDir() 
                                      + graphics::createFileName( _name ) + ".pvd";    
      std::string fullNamePVD = singleton::directories().getVtkOutDir()
                                + "data/" + graphics::createFileName( _name, iT ) + ".pvd";
      preamblePVD(fullNamePVD);
      for (int iC=0; iC<cGeometry.getNc(); iC++) 
      {
        std::string namePieceData = "data/" + graphics::createFileName( _name, iT, iC) + ".vti";
        std::string namePiece = graphics::createFileName( _name, iT, iC) + ".vti";

        // puts name of .vti piece to a .pvd file [fullNamePVD]
        dataPVD( iT, iC, fullNamePVD, namePiece );
        // puts name of .vti piece to the master .pvd file [fullNamePVDmaster]
        // adds vti piece name, via deleting the ending of the master pvd file
        // and performs the closing as well.
        dataPVDmaster( iT, iC, fullNamePVDmaster, namePieceData );
      }
      closePVD(fullNamePVD);
    } // master only
    //-------------------------------------

    //-------------------------------------
    // VTI
    // each process writes his cuboids
    //-------------------------------------
    for (int iCloc=0; iCloc<load.size(); iCloc++)
    { // cuboid
      int nx = cGeometry.get(load.glob(iCloc)).getNx();
      int ny = cGeometry.get(load.glob(iCloc)).getNy();
      int nz = cGeometry.get(load.glob(iCloc)).getNz();
      // to be changed into the following line once local refinement has been implemented
      // double deltaX = cGeometry.get(load.glob(iCloc)).getDeltaR();
      T delta = cGeometry.getMotherCuboid().getDeltaR();

      std::string fullNameVTI = singleton::directories().getVtkOutDir() + "data/"
                  + graphics::createFileName( _name, iT, load.glob(iCloc) ) + ".vti";

      // get dimension/extent for each cuboid
      std::vector<int> originLatticeR(4,int());
      originLatticeR[0]=load.glob(iCloc);
      std::vector<T> originPhysR = cGeometry.getPhysR(originLatticeR);

      preambleVTI(fullNameVTI, -1,-1,-1,nx,ny,nz,
                  originPhysR[0],originPhysR[1],originPhysR[2], delta);
      pointData(fullNameVTI);

      typename std::vector<SuperLatticeF3D<T,DESCRIPTOR>* >::iterator it;
      for( it = _pointerVec.begin(); it != _pointerVec.end(); it++)
      {
        // get functor and prevent deleting with use of Identity
        SuperLatticeIdentity3D<T,DESCRIPTOR> f( **it );
        // write <DataArray ..../> for functor, e.g. velocity, pressure, ...
        if( _binary ) {
          dataArrayBinary(fullNameVTI, f, load.glob(iCloc), nx+1,ny+1,nz+1);
        } else {
          dataArray(fullNameVTI, f, load.glob(iCloc), nx+1,ny+1,nz+1);
        }
      }

      closePointData(fullNameVTI);
      closePiece(fullNameVTI);
      closeVTI(fullNameVTI);
    } // cuboid
  } // END WORK
  ////-------------------------------------
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::write(SuperLatticeF3D<T,DESCRIPTOR>& f, int iT)
{
  CuboidGeometry3D<T> const& cGeometry = f.getSuperLattice3D().getCuboidGeometry();
  LoadBalancer<T>& load = f.getSuperLattice3D().getLoadBalancer();
  // no gaps between vti files (cuboids)
  f.getSuperLattice3D().communicate();

  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif

  // write a pvd file, which links all vti files
  // each vti file is written by one thread, which may own severals cuboids
  if( rank == 0 )
  { // master only
    std::string namePVD = singleton::directories().getVtkOutDir() + 
                          graphics::createFileName( _name, f.getName(), iT )  + ".pvd";

    preamblePVD(namePVD);
    for (int iC=0; iC<cGeometry.getNc(); iC++) 
    {
      std::string nameVTI = "data/" + graphics::createFileName( _name, f.getName(), iT, iC)
                            + ".vti";
      // puts name of .vti piece to a .pvd file [fullNamePVD]
      dataPVD( iT, iC, namePVD, nameVTI );
    }
    closePVD(namePVD);
  } // master only

  for (int iCloc=0; iCloc<load.size(); iCloc++)
  { // cuboid
    int nx = cGeometry.get(load.glob(iCloc)).getNx();
    int ny = cGeometry.get(load.glob(iCloc)).getNy();
    int nz = cGeometry.get(load.glob(iCloc)).getNz();
    // to be changed into the following line once local refinement has been implemented
    // double deltaX = cGeometry.get(load.glob(iCloc)).getDeltaR();
    T delta = cGeometry.getMotherCuboid().getDeltaR();

    std::string fullNameVTI = singleton::directories().getVtkOutDir() + "data/"
                + graphics::createFileName( _name, f.getName(), iT, load.glob(iCloc) ) + ".vti";

    // get dimension/extent for each cuboid
    std::vector<int> originLatticeR(4,int());
    originLatticeR[0]=load.glob(iCloc);
    std::vector<T> originPhysR = cGeometry.getPhysR(originLatticeR);

    preambleVTI(fullNameVTI, -1,-1,-1, nx,ny,nz,
                originPhysR[0],originPhysR[1],originPhysR[2], delta);
    pointData(fullNameVTI);
    if( _binary ) {
      dataArrayBinary(fullNameVTI, f,  load.glob(iCloc), nx+1,ny+1,nz+1);
    } else {
      dataArray(fullNameVTI, f, load.glob(iCloc), nx+1,ny+1,nz+1);
    }

    closePointData(fullNameVTI);
    closePiece(fullNameVTI);
    closeVTI(fullNameVTI);
  } // cuboid
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::createMasterFile()
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if( rank == 0 ) {
    std::string fullNamePVDmaster = singleton::directories().getVtkOutDir() 
                + graphics::createFileName( _name ) + ".pvd"; 
    preamblePVD(fullNamePVDmaster);
    closePVD(fullNamePVDmaster);
    _createFile = true;
  }
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::addFunctor(SuperLatticeF3D<T,DESCRIPTOR>& f)
{ _pointerVec.push_back(&f); }

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::clearAddedFunctors()
{ _pointerVec.clear(); }




////////////////////private member functions///////////////////////////////////
template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::preambleVTI (const std::string& fullName,
  int x0,int y0,int z0,int x1,int y1,int z1, T originX,T originY,T originZ, T delta)
{
  std::ofstream fout(fullName.c_str(), std::ios::trunc);
  if (!fout) clout << "Error: could not open " << fullName << std::endl;

  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
          << "byte_order=\"LittleEndian\">\n";
  fout << "<ImageData WholeExtent=\""
          << x0 <<" "<< x1 <<" "
          << y0 <<" "<< y1 <<" "
          << z0 <<" "<< z1 
       << "\" Origin=\"" << originX << " " << originY << " " << originZ
       << "\" Spacing=\"" << delta << " " << delta << " " << delta << "\">\n";
  fout << "<Piece Extent=\"" 
          << x0 <<" "<< x1 <<" "
          << y0 <<" "<< y1 <<" "
          << z0 <<" "<< z1 <<"\">\n";
  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::closeVTI(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app );
  if (!fout) clout << "Error: could not open " << fullNamePiece << std::endl;
  fout << "</ImageData>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::preamblePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::trunc);
  if (!fout) clout << "Error: could not open " << fullNamePVD << std::endl;

  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"Collection\" version=\"0.1\" " 
          << "byte_order=\"LittleEndian\">\n"
       << "<Collection>\n";
  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::closePVD(const std::string& fullNamePVD)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app );
  if (!fout) clout << "Error: could not open " << fullNamePVD << std::endl;
  fout << "</Collection>\n";
  fout << "</VTKFile>\n";
  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::pointData(const std::string& fullName)
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) clout << "Error: could not open " << fullName << std::endl;

  fout << "<PointData>\n";
  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::closePointData(const std::string& fullName)
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) clout << "Error: could not open " << fullName << std::endl;
  fout << "</PointData>\n";
  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::dataPVD(int iT, int iC,
  const std::string& fullNamePVD, const std::string& namePiece)
{
  std::ofstream fout(fullNamePVD.c_str(), std::ios::app);
  if (!fout) clout << "Error: could not open " << fullNamePVD << std::endl;

  fout << "<DataSet timestep=\"" << iT << "\" "
          << "group=\"\" part=\" " <<  iC << "\" "
          << "file=\"" << namePiece << "\"/>\n";
  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::dataPVDmaster(int iT, int iC,
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
  }
  else { clout << "Error: could not open " << fullNamePVDMaster << std::endl; }
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::dataArray( const std::string& fullName,
  SuperLatticeF3D<T,DESCRIPTOR>& f, int iC, int nx, int ny, int nz)
{
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) clout << "Error: could not open " << fullName << std::endl;

  int dim = f.getTargetDim();
  fout << "<DataArray " ;
  if (dim == 1) {
    fout << "type=\"Float32\" Name=\"" << f.getName() <<"\">\n";
  }
  if (dim != 1) {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "NumberOfComponents=\"" << dim <<"\">\n";
  }

  // since cuboid has been blowed up by 1 [every dimension]
  // looping from -1 to nz (= nz+1, as passed)
  std::vector<int> tmpVec( 4,int(0) );
  for (int iZ=-1; iZ<nz; ++iZ) {
    for (int iY=-1; iY<ny; ++iY) {
      for (int iX=-1; iX<nx; ++iX) {
        for (int iDim=0; iDim<dim; ++iDim) {
          tmpVec[0]=iC; tmpVec[1]=iX; tmpVec[2]=iY; tmpVec[3]=iZ;
          //  tmpVec = {iC,iX,iY,iZ};  // std=c++11 
          fout <<  f(tmpVec)[iDim] << " ";
        }
      }
    }
  }
  fout << "</DataArray>\n";

  fout.close();
}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::dataArrayBinary(const std::string& fullName,
  SuperLatticeF3D<T,DESCRIPTOR>& f, int iC, int nx, int ny, int nz)
{
  const char* fileName = fullName.c_str();

  std::ofstream fout( fileName, std::ios::out | std::ios::app );
  if (!fout) clout << "Error: could not open " << fileName << std::endl;

  int dim = f.getTargetDim();
  fout << "<DataArray ";
  if (dim == 1) {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" " 
         << "format=\"binary\" encoding=\"base64\">";
  }
  if (dim != 1) {
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "format=\"binary\" encoding=\"base64\" "
         << "NumberOfComponents=\"" << dim 
         <<"\">\n";
  }
  fout.close();

  std::ofstream* ofstr;   // only used for binary output // passed to Base64Encoder
  ofstr = new std::ofstream( fileName, std::ios::out | std::ios::app | std::ios::binary );
  if (!ofstr) clout << "Error: could not open " << fileName << std::endl;

  size_t fullSize = dim * (nx+1) * (ny+1) * (nz+1);    //  how many numbers to write
  size_t binarySize = size_t( fullSize*sizeof(float) );
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(*ofstr, 1);
  unsigned int uintBinarySize = (unsigned int)binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<float>* dataEncoder = new Base64Encoder<float>( *ofstr, fullSize );
  std::vector<int> tmpVec( 4,int(0) );
  for (int iZ=-1; iZ<nz; ++iZ) {
    for (int iY=-1; iY<ny; ++iY) {
      for (int iX=-1; iX<nx; ++iX) {
        for (int iDim=0; iDim<dim; ++iDim) {
          tmpVec[0] = iC; tmpVec[1] = iX; tmpVec[2] = iY; tmpVec[3] = iZ;
          const float tmp = float( f(tmpVec)[iDim] );
          dataEncoder->encode( &tmp, 1 );
        }
      }
    }
  }
  ofstr->close();
  
  std::ofstream ffout( fileName,  std::ios::out | std::ios::app );
  if (!ffout) clout << "Error: could not open " << fileName << std::endl;
  ffout << "</DataArray>\n";
  ffout.close();

}

template<typename T, template <typename U> class DESCRIPTOR>
void SuperVTKwriter3D<T,DESCRIPTOR>::closePiece(const std::string& fullNamePiece)
{
  std::ofstream fout(fullNamePiece.c_str(), std::ios::app );
  if (!fout) clout << "Error: could not open " << fullNamePiece << std::endl;
  fout << "</Piece>\n";
  fout.close();
}

}  // namespace olb

#endif

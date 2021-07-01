/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Thomas Henn, Mathias J. Krause
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

#ifndef SUPERPARTICLESYSVTUOUT_H
#define SUPERPARTICLESYSVTUOUT_H

#include <stddef.h>
#include <sys/types.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

#include "core/singleton.h"
#include "core/units.h"
#include "io/base64.h"
#include "io/base64.hh"
#include "io/fileName.h"
#include "io/ostreamManager.h"
#include "superParticleSystem3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class SuperParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE>
class SuperParticleSysVtuWriter {
public:
  //  SuperParticleSysVtuWriter() = default;
  SuperParticleSysVtuWriter(SuperParticleSystem3D<T, PARTICLETYPE>&,
                            std::string const,
                            unsigned short properties,
                            bool binary = true);
  SuperParticleSysVtuWriter(const SuperParticleSysVtuWriter<T, PARTICLETYPE>& rhs);
  SuperParticleSysVtuWriter(const SuperParticleSysVtuWriter<T, PARTICLETYPE>&& rhs);

  void write(int iT = 0);

  int numofpsys();

  void set(unsigned short);

  enum particleProperties
  : unsigned short {velocity = 1, radius = 2, mass = 4, force = 8, cuboid = 16, active = 32};

private:
  ///  performes <VTKFile ...>, <ImageData ...> and <PieceExtent ...>
  void preambleVTU(const std::string& fullName);
  ///  performes </ImageData> and </VTKFile>
  void closeVTU(const std::string& fullNamePiece);
  ///  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  ///  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVD(int iT, int iR, const std::string& fullNamePVD,
               const std::string& namePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  void dataPVDmaster(int iT, int iR, const std::string& fullNamePVDMaster,
                     const std::string& namePiece);
  ///  writes functors stored at pointerVec
  void dataArray(const std::string& fullName);
  ///  writes functors stored at pointerVec
  void dataArrayBinary(const std::string& fullName);
  ///  performes <VTKFile...> and <ImageData ...>
  void preambleOneFile(const std::string& fullName);
  ///  writes instantaniously given functor, without adding to _pointerVec
  void writePieceToOneFile(const std::string& fullName);

  void createMasterFile();

private:
  unsigned short _properties;
  bool _haveMaster;
  SuperParticleSystem3D<T, PARTICLETYPE>& _psys;
  std::string _name;
  bool _binary;
  mutable OstreamManager clout;
};

}  // namespace OLB

#endif /* SUPERPARTICLESYSVTUOUT_H */

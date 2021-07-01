/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathen
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

#ifndef ADM_LATTICE_DESCRIPTOR_H
#define ADM_LATTICE_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"


namespace olb {

namespace descriptors {

// 2D Descriptors for ADM


struct ADM2Ddescriptor {
  static const int numScalars = 3;
  static const int numSpecies = 2;
  static const int rhoIsAt     = 0;
  static const int sizeOfRho   = 1;
  static const int velocityBeginsAt = 1;
  static const int sizeOfVelocity   = 2;
};

struct ADM2DdescriptorBase {
  typedef ADM2Ddescriptor ExternalField;
};

template <typename T> struct ADMD2Q9Descriptor
    : public D2Q9DescriptorBase<T>, public ADM2DdescriptorBase {
};




////////////////////////////////////////////////////////////////////////////////
// extended descriptor for ADM

struct ADM3dDescriptor {
  static const int numScalars = 4;
  static const int numSpecies = 2;
  static const int rhoIsAt     = 0;
  static const int sizeOfRho   = 1;
  static const int velocityBeginsAt = 1;
  static const int sizeOfVelocity   = 3;
};

struct ADM3dDescriptorBase {
  typedef ADM3dDescriptor ExternalField;
};

template <typename T> struct ADMD3Q19Descriptor
    : public D3Q19DescriptorBase<T>, public ADM3dDescriptorBase {
};


////////////////////////////////////////
/// ADM filter descriptors for forced fields

//// Forced 2D ADM filter scheme
struct ForcedADM2Ddescriptor {
  static const int numScalars = 5;
  static const int numSpecies = 3;
  static const int rhoIsAt     = 0;
  static const int sizeOfRho   = 1;
  static const int velocityBeginsAt = 1;
  static const int sizeOfVelocity   = 2;
  static const int forceBeginsAt    = 3;
  static const int sizeOfForce      = 2;
};

struct ForcedADM2DdescriptorBase {
  typedef ForcedADM2Ddescriptor ExternalField;
};

template <typename T> struct ForcedADMD2Q9Descriptor
    : public D2Q9DescriptorBase<T>, public ForcedADM2DdescriptorBase {
};


//// Forced 3D ADM filter scheme

struct ForcedADM3dDescriptor {
  static const int numScalars = 7;
  static const int numSpecies = 2;
  static const int rhoIsAt     = 0;
  static const int sizeOfRho   = 1;
  static const int velocityBeginsAt = 1;
  static const int sizeOfVelocity   = 3;
  static const int forceBeginsAt    = 4;
  static const int sizeOfForce      = 3;
};

struct ForcedADM3dDescriptorBase {
  typedef ForcedADM3dDescriptor ExternalField;
};

template <typename T> struct ForcedADMD3Q19Descriptor
    : public D3Q19DescriptorBase<T>, public ForcedADM3dDescriptorBase {
};


} // namespace descriptors

} // namespace olb

#endif

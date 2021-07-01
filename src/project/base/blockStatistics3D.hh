/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
 *  E-mail: jonas.latt@gmail.com
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

#ifndef BLOCK_STATISTICS_3D_HH
#define BLOCK_STATISTICS_3D_HH

#include <cmath>
#include "cell.h"
#include "blockStatistics3D.h"
#include "finiteDifference.h"
#include "util.h"

namespace olb {


/////// Class BlockStatistics3D  /////////////////////////////

template<typename T, template<typename U> class Lattice>
BlockStatistics3D<T,Lattice>::BlockStatistics3D (
    BlockStructure3D<T,Lattice> const& block_ )
    : block(block_),
      velField       (block.getNx(), block.getNy(), block.getNz()),
      momentumField  (block.getNx(), block.getNy(), block.getNz()),
      pressureField  (block.getNx(), block.getNy(), block.getNz()),
      velNormField   (block.getNx(), block.getNy(), block.getNz()),
      vortNormField  (block.getNx(), block.getNy(), block.getNz()),
      vortField      (block.getNx(), block.getNy(), block.getNz()),
      strainRateField(block.getNx(), block.getNy(), block.getNz()),
      stressField    (block.getNx(), block.getNy(), block.getNz()),
      divRhoUField   (block.getNx(), block.getNy(), block.getNz()),
      poissonField   (block.getNx(), block.getNy(), block.getNz()),
      populationField(block.getNx(), block.getNy(), block.getNz())
{
    iniFlags();
}

template<typename T, template<typename U> class Lattice>
BlockStatistics3D<T,Lattice>::~BlockStatistics3D() {
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::reset() const {
    iniFlags();
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::iniFlags() const {
    velFieldComputed = false;
    momentumFieldComputed = false;
    pressureFieldComputed = false;
    velNormFieldComputed = false;
    vortNormFieldComputed = false;
    vortFieldComputed = false;
    strainRateFieldComputed = false;
    stressFieldComputed = false;
    divRhoUFieldComputed = false;
    poissonFieldComputed = false;
    populationFieldComputed = false;
}

template<typename T, template<typename U> class Lattice>
TensorField3D<T,3> const&
    BlockStatistics3D<T,Lattice>::getVelocity() const
{
    computeVelocityField();
    return velField;
}

template<typename T, template<typename U> class Lattice>
TensorField3D<T,3> const&
    BlockStatistics3D<T,Lattice>::getMomentum() const
{
    computeMomentumField();
    return momentumField;
}

template<typename T, template<typename U> class Lattice>
ScalarField3D<T> const&
    BlockStatistics3D<T,Lattice>::getPressure() const
{
    computePressureField();
    return pressureField;
}

template<typename T, template<typename U> class Lattice>
TensorField3D<T,3> const&
    BlockStatistics3D<T,Lattice>::getVorticity() const
{
    computeVorticityField();
    return vortField;
}

template<typename T, template<typename U> class Lattice>
ScalarField3D<T> const&
    BlockStatistics3D<T,Lattice>::getVelocityNorm() const
{
    computeVelocityNormField();
    return velNormField;
}

template<typename T, template<typename U> class Lattice>
ScalarField3D<T> const&
    BlockStatistics3D<T,Lattice>::getVorticityNorm() const
{
    computeVorticityNormField();
    return vortNormField;
}

template<typename T, template<typename U> class Lattice>
TensorField3D<T,6> const&
    BlockStatistics3D<T,Lattice>::getStrainRate() const
{
    computeStrainRateField();
    return strainRateField;
}

template<typename T, template<typename U> class Lattice>
TensorField3D<T,6> const&
    BlockStatistics3D<T,Lattice>::getStress() const
{
    computeStressField();
    return stressField;
}

template<typename T, template<typename U> class Lattice>
ScalarField3D<T> const&
    BlockStatistics3D<T,Lattice>::getDivRhoU() const
{
    computeDivRhoUField();
    return divRhoUField;
}

template<typename T, template<typename U> class Lattice>
ScalarField3D<T> const&
    BlockStatistics3D<T,Lattice>::getPoissonTerm() const
{
    computePoissonTerm();
    return poissonField;
}

template<typename T, template<typename U> class Lattice>
TensorField3D<T, Lattice<T>::q > const&
    BlockStatistics3D<T,Lattice>::getPopulations() const
{
    computePopulations();
    return populationField;
}


template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::computeMeanEnstrophy() const {
    computeVelocityField();
    int nx = velField.getNx()-1;
    int ny = velField.getNy()-1;
    int nz = velField.getNz()-1;

    T enstrophy = T();
    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            for (int iZ=0; iZ<nz; ++iZ) {
                T dxuy = (
                      velField.get(iX+1,iY+1,iZ  )[1]
                    + velField.get(iX+1,iY  ,iZ  )[1]
                    + velField.get(iX+1,iY+1,iZ+1)[1]
                    + velField.get(iX+1,iY  ,iZ+1)[1]
                    - velField.get(iX  ,iY+1,iZ  )[1]
                    - velField.get(iX  ,iY  ,iZ  )[1]
                    - velField.get(iX  ,iY+1,iZ+1)[1]
                    - velField.get(iX  ,iY  ,iZ+1)[1] ) / (T)4;
                T dxuz = (
                      velField.get(iX+1,iY+1,iZ  )[2]
                    + velField.get(iX+1,iY  ,iZ  )[2]
                    + velField.get(iX+1,iY+1,iZ+1)[2]
                    + velField.get(iX+1,iY  ,iZ+1)[2]
                    - velField.get(iX  ,iY+1,iZ  )[2]
                    - velField.get(iX  ,iY  ,iZ  )[2]
                    - velField.get(iX  ,iY+1,iZ+1)[2]
                    - velField.get(iX  ,iY  ,iZ+1)[2] ) / (T)4;
                T dyux = (
                      velField.get(iX  ,iY+1,iZ  )[0]
                    + velField.get(iX+1,iY+1,iZ  )[0]
                    + velField.get(iX  ,iY+1,iZ+1)[0]
                    + velField.get(iX+1,iY+1,iZ+1)[0]
                    - velField.get(iX  ,iY  ,iZ  )[0]
                    - velField.get(iX+1,iY  ,iZ  )[0]
                    - velField.get(iX  ,iY  ,iZ+1)[0]
                    - velField.get(iX+1,iY  ,iZ+1)[0] ) / (T)4;
                T dyuz = (
                      velField.get(iX  ,iY+1,iZ  )[2]
                    + velField.get(iX+1,iY+1,iZ  )[2]
                    + velField.get(iX  ,iY+1,iZ+1)[2]
                    + velField.get(iX+1,iY+1,iZ+1)[2]
                    - velField.get(iX  ,iY  ,iZ  )[2]
                    - velField.get(iX+1,iY  ,iZ  )[2]
                    - velField.get(iX  ,iY  ,iZ+1)[2]
                    - velField.get(iX+1,iY  ,iZ+1)[2] ) / (T)4;
                T dzux = (
                      velField.get(iX  ,iY  ,iZ+1)[0]
                    + velField.get(iX+1,iY+1,iZ+1)[0]
                    + velField.get(iX  ,iY+1,iZ+1)[0]
                    + velField.get(iX+1,iY  ,iZ+1)[0]
                    - velField.get(iX  ,iY  ,iZ  )[0]
                    - velField.get(iX+1,iY  ,iZ  )[0]
                    - velField.get(iX  ,iY+1,iZ  )[0]
                    - velField.get(iX+1,iY+1,iZ  )[0] ) / (T)4;
                T dzuy = (
                      velField.get(iX  ,iY  ,iZ+1)[1]
                    + velField.get(iX+1,iY+1,iZ+1)[1]
                    + velField.get(iX  ,iY+1,iZ+1)[1]
                    + velField.get(iX+1,iY  ,iZ+1)[1]
                    - velField.get(iX  ,iY  ,iZ  )[1]
                    - velField.get(iX+1,iY  ,iZ  )[1]
                    - velField.get(iX  ,iY+1,iZ  )[1]
                    - velField.get(iX+1,iY+1,iZ  )[1] ) / (T)4;

                T omegaX = dyuz - dzuy;
                T omegaY = dzux - dxuz;
                T omegaZ = dxuy - dyux;

                enstrophy += omegaX*omegaX+omegaY*omegaY+omegaZ*omegaZ;
            }
        }
    }
    enstrophy /= (2*nx*ny*nz);
    return enstrophy;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computeVelocityField() const {
    if (velFieldComputed) return;
    velField.construct();
    for (int iX=0; iX<velField.getNx(); ++iX) {
        for (int iY=0; iY<velField.getNy(); ++iY) {
            for (int iZ=0; iZ<velField.getNz(); ++iZ) {
                block.get(iX,iY,iZ).computeU(velField.get(iX,iY,iZ));
            }
        }
    }
    velFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computeMomentumField() const {
    if (momentumFieldComputed) return;
    momentumField.construct();
    for (int iX=0; iX<momentumField.getNx(); ++iX) {
        for (int iY=0; iY<momentumField.getNy(); ++iY) {
            for (int iZ=0; iZ<momentumField.getNz(); ++iZ) {
                T rho;
                block.get(iX,iY,iZ).computeRhoU (
                        rho, momentumField.get(iX,iY,iZ) );
                for (int iD=0; iD<Lattice<T>::d; ++iD) {
                    momentumField.get(iX,iY,iZ)[iD] *= rho;
                }
            }
        }
    }
    momentumFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computePressureField() const {
    if (pressureFieldComputed) return;
    pressureField.construct();
    for (int iX=0; iX<pressureField.getNx(); ++iX) {
        for (int iY=0; iY<pressureField.getNy(); ++iY) {
            for (int iZ=0; iZ<momentumField.getNz(); ++iZ) {
                pressureField.get(iX,iY,iZ) =
                    block.get(iX,iY,iZ).computeRho();
                pressureField.get(iX,iY,iZ) -= (T)1;
                pressureField.get(iX,iY,iZ) /= Lattice<T>::invCs2;
            }
        }
    }
    pressureFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computeVelocityNormField() const {
    if (velNormFieldComputed) return;
    velNormField.construct();
    computeVelocityField();
    for (int iEl=0; iEl<velNormField.getSize(); ++iEl) {
        velNormField[iEl] = sqrt(util::normSqr<T,3>(velField[iEl]));
    }
    velNormFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computeVorticityNormField() const {
    if (vortNormFieldComputed) return;
    vortNormField.construct();
    computeVorticityField();
    for (int iEl=0; iEl<vortNormField.getSize(); ++iEl) {
        vortNormField[iEl] = sqrt(util::normSqr<T,3>(vortField[iEl]));
    }
    vortNormFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkVorticityX(int iX, int iY, int iZ) const {
    OLB_PRECONDITION(velFieldComputed);
    OLB_PRECONDITION(iX>=1 && iX<=vortField.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=vortField.getNy()-2);
    OLB_PRECONDITION(iZ>=1 && iZ<=vortField.getNz()-2);

    T dyuz = bulkYderiv(iX,iY,iZ, 2, velField);
    T dzuy = bulkZderiv(iX,iY,iZ, 1, velField);

    return dyuz - dzuy;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkVorticityY(int iX, int iY, int iZ) const {
    OLB_PRECONDITION(velFieldComputed);
    OLB_PRECONDITION(iX>=1 && iX<=vortField.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=vortField.getNy()-2);
    OLB_PRECONDITION(iZ>=1 && iZ<=vortField.getNz()-2);

    T dzux = bulkZderiv(iX,iY,iZ, 0, velField);
    T dxuz = bulkXderiv(iX,iY,iZ, 2, velField);

    return dzux - dxuz;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkVorticityZ(int iX, int iY, int iZ) const {
    OLB_PRECONDITION(iX>=1 && iX<=vortField.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=vortField.getNy()-2);
    OLB_PRECONDITION(iZ>=1 && iZ<=vortField.getNz()-2);

    T dxuy = bulkXderiv(iX,iY,iZ, 1, velField);
    T dyux = bulkYderiv(iX,iY,iZ, 0, velField);

    return dxuy - dyux;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkXderiv (
        int iX, int iY, int iZ, int iD,
        TensorField3D<T,3> const& field) const
{
    OLB_PRECONDITION(iX>=1 && iX<=field.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=field.getNy()-2);
    OLB_PRECONDITION(iZ>=1 && iZ<=field.getNz()-2);

    T dxu = fd::centralGradient(field.get(iX+1,iY,iZ)[iD],
                                field.get(iX-1,iY,iZ)[iD]);
    return dxu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkYderiv (
        int iX, int iY, int iZ, int iD,
        TensorField3D<T,3> const& field) const
{
    OLB_PRECONDITION(iX>=1 && iX<=field.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=field.getNy()-2);
    OLB_PRECONDITION(iZ>=1 && iZ<=field.getNz()-2);

    T dyu = fd::centralGradient(field.get(iX,iY+1,iZ)[iD],
                                field.get(iX,iY-1,iZ)[iD]);
    return dyu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkZderiv (
        int iX, int iY, int iZ, int iD,
        TensorField3D<T,3> const& field) const
{
    OLB_PRECONDITION(iX>=1 && iX<=field.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=field.getNy()-2);
    OLB_PRECONDITION(iZ>=1 && iZ<=field.getNz()-2);

    T dzu = fd::centralGradient(field.get(iX,iY,iZ+1)[iD],
                                field.get(iX,iY,iZ-1)[iD]);
    return dzu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkDeriv (
        int iX, int iY, int iZ, int iAlpha, int iBeta,
        TensorField3D<T,3> const& field) const
{
    switch(iAlpha) {
        case 0:
            return bulkXderiv(iX,iY,iZ, iBeta, field);
        case 1:
            return bulkYderiv(iX,iY,iZ, iBeta, field);
        case 2:
            return bulkZderiv(iX,iY,iZ, iBeta, field);
        default:
            OLB_ASSERT( false, "iAlpha>2!");
            return T();
    }
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkStrain (
    int iX, int iY, int iZ, int iAlpha, int iBeta) const
{
    OLB_PRECONDITION( momentumFieldComputed );
    return ( bulkDeriv(iX,iY,iZ, iAlpha,iBeta, momentumField) +
             bulkDeriv(iX,iY,iZ, iBeta,iAlpha, momentumField) ) / (T)2;

}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkDivRhoU(int iX, int iY, int iZ) const {
    OLB_PRECONDITION( momentumFieldComputed );
    return bulkDeriv(iX,iY,iZ, 0,0, momentumField) +
           bulkDeriv(iX,iY,iZ, 1,1, momentumField) +
           bulkDeriv(iX,iY,iZ, 2,2, momentumField);
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::bulkPoisson(int iX, int iY, int iZ) const {
    OLB_PRECONDITION( velFieldComputed );

    T dxux = bulkDeriv(iX,iY,iZ, 0,0, velField);
    T dxuy = bulkDeriv(iX,iY,iZ, 0,1, velField);
    T dxuz = bulkDeriv(iX,iY,iZ, 0,2, velField);
    T dyux = bulkDeriv(iX,iY,iZ, 1,0, velField);
    T dyuy = bulkDeriv(iX,iY,iZ, 1,1, velField);
    T dyuz = bulkDeriv(iX,iY,iZ, 1,2, velField);
    T dzux = bulkDeriv(iX,iY,iZ, 2,0, velField);
    T dzuy = bulkDeriv(iX,iY,iZ, 2,1, velField);
    T dzuz = bulkDeriv(iX,iY,iZ, 2,2, velField);

    return dxux*dxux + dyuy*dyuy + dzuz*dzuz +
           (T)2*( dxuy*dyux + dxuz*dzux + dyuz*dzuy );
}


template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryXderiv (
        int iX, int iY, int iZ, int iD,
        TensorField3D<T,3> const& field) const
{
    OLB_PRECONDITION(iX>=0 && iX<=field.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=field.getNy()-1);
    OLB_PRECONDITION(iZ>=0 && iZ<=field.getNz()-1);

    T dxu;
    
    if (iX==0) {
        dxu = fd::boundaryGradient(field.get(iX,iY,iZ)[iD],
                                   field.get(iX+1,iY,iZ)[iD],
                                   field.get(iX+2,iY,iZ)[iD]);
    }
    else if (iX==field.getNx()-1) {
        dxu = -fd::boundaryGradient(field.get(iX,iY,iZ)[iD],
                                    field.get(iX-1,iY,iZ)[iD],
                                    field.get(iX-2,iY,iZ)[iD]);
    }
    else {
        dxu = fd::centralGradient(field.get(iX+1,iY,iZ)[iD],
                                  field.get(iX-1,iY,iZ)[iD]);
    }

    return dxu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryYderiv (
        int iX, int iY, int iZ, int iD,
        TensorField3D<T,3> const& field) const
{
    OLB_PRECONDITION(iX>=0 && iX<=field.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=field.getNy()-1);
    OLB_PRECONDITION(iZ>=0 && iZ<=field.getNz()-1);

    T dyu;

    if (iY==0) {
        dyu = fd::boundaryGradient(field.get(iX,iY,iZ)[iD],
                                   field.get(iX,iY+1,iZ)[iD],
                                   field.get(iX,iY+2,iZ)[iD]);
    }
    else if (iY==field.getNy()-1) {
        dyu = -fd::boundaryGradient(field.get(iX,iY,iZ)[iD],
                                    field.get(iX,iY-1,iZ)[iD],
                                    field.get(iX,iY-2,iZ)[iD]);
    }
    else {
        dyu = fd::centralGradient(field.get(iX,iY+1,iZ)[iD],
                                  field.get(iX,iY-1,iZ)[iD]);
    }
   
    return dyu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryZderiv (
        int iX, int iY, int iZ, int iD,
        TensorField3D<T,3> const& field) const
{
    OLB_PRECONDITION(iX>=0 && iX<=field.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=field.getNy()-1);
    OLB_PRECONDITION(iZ>=0 && iZ<=field.getNz()-1);

    T dzu;

    if (iZ==0) {
        dzu = fd::boundaryGradient(field.get(iX,iY,iZ)[iD],
                                   field.get(iX,iY,iZ+1)[iD],
                                   field.get(iX,iY,iZ+2)[iD]);
    }
    else if (iZ==field.getNz()-1) {
        dzu = -fd::boundaryGradient(field.get(iX,iY,iZ)[iD],
                                    field.get(iX,iY,iZ-1)[iD],
                                    field.get(iX,iY,iZ-2)[iD]);
    }
    else {
        dzu = fd::centralGradient(field.get(iX,iY,iZ+1)[iD],
                                  field.get(iX,iY,iZ-1)[iD]);
    }
   
    return dzu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryDeriv (
    int iX, int iY, int iZ, int iAlpha, int iBeta,
    TensorField3D<T,3> const& field) const
{
    switch(iAlpha) {
        case 0:
            return boundaryXderiv(iX,iY,iZ, iBeta, field);
        case 1:
            return boundaryYderiv(iX,iY,iZ, iBeta, field);
        case 2:
            return boundaryZderiv(iX,iY,iZ, iBeta, field);
        default:
            OLB_ASSERT( false, "iAlpha>2!");
            return T();
    }
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryStrain (
    int iX, int iY, int iZ, int iAlpha, int iBeta) const
{
    OLB_PRECONDITION( momentumFieldComputed );
    return ( boundaryDeriv(iX,iY,iZ,iAlpha,iBeta, momentumField) +
             boundaryDeriv(iX,iY,iZ,iBeta,iAlpha, momentumField) ) / (T)2;

}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryDivRhoU (
        int iX, int iY, int iZ) const
{
    OLB_PRECONDITION( momentumFieldComputed );
    return boundaryDeriv(iX,iY,iZ, 0,0, momentumField) +
           boundaryDeriv(iX,iY,iZ, 1,1, momentumField) +
           boundaryDeriv(iX,iY,iZ, 2,2, momentumField);
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryPoisson (
        int iX, int iY, int iZ) const
{
    OLB_PRECONDITION( velFieldComputed );

    T dxux = boundaryDeriv(iX,iY,iZ, 0,0, velField);
    T dxuy = boundaryDeriv(iX,iY,iZ, 0,1, velField);
    T dxuz = boundaryDeriv(iX,iY,iZ, 0,2, velField);
    T dyux = boundaryDeriv(iX,iY,iZ, 1,0, velField);
    T dyuy = boundaryDeriv(iX,iY,iZ, 1,1, velField);
    T dyuz = boundaryDeriv(iX,iY,iZ, 1,2, velField);
    T dzux = boundaryDeriv(iX,iY,iZ, 2,0, velField);
    T dzuy = boundaryDeriv(iX,iY,iZ, 2,1, velField);
    T dzuz = boundaryDeriv(iX,iY,iZ, 2,2, velField);

    return dxux*dxux + dyuy*dyuy + dzuz*dzuz +
           (T)2*( dxuy*dyux + dxuz*dzux + dyuz*dzuy );
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryVorticityX (
        int iX, int iY, int iZ) const
{
    OLB_PRECONDITION(velFieldComputed);
    OLB_PRECONDITION(iX>=0 && iX<=vortField.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=vortField.getNy()-1);
    OLB_PRECONDITION(iZ>=0 && iZ<=vortField.getNz()-1);

    T dyuz = boundaryYderiv(iX,iY,iZ, 2, velField);
    T dzuy = boundaryZderiv(iX,iY,iZ, 1, velField);

    return dyuz - dzuy;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryVorticityY (
        int iX, int iY, int iZ) const
{
    OLB_PRECONDITION(velFieldComputed);
    OLB_PRECONDITION(iX>=0 && iX<=vortField.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=vortField.getNy()-1);
    OLB_PRECONDITION(iZ>=0 && iZ<=vortField.getNz()-1);

    T dzux = boundaryZderiv(iX,iY,iZ, 0, velField);
    T dxuz = boundaryXderiv(iX,iY,iZ, 2, velField);

    return dzux - dxuz;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics3D<T,Lattice>::boundaryVorticityZ (
        int iX, int iY, int iZ) const
{
    OLB_PRECONDITION(velFieldComputed);
    OLB_PRECONDITION(iX>=0 && iX<=vortField.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=vortField.getNy()-1);
    OLB_PRECONDITION(iZ>=0 && iZ<=vortField.getNz()-1);

    T dxuy = boundaryXderiv(iX,iY,iZ, 1, velField);
    T dyux = boundaryYderiv(iX,iY,iZ, 0, velField);

    return dxuy - dyux;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computeVorticityField() const {
    if (vortFieldComputed) return;
    vortField.construct();
    computeVelocityField();

    int nx = vortField.getNx();
    int ny = vortField.getNy();
    int nz = vortField.getNz();

    for (int iX=1; iX<nx-1; ++iX) {
        for (int iY=1; iY<ny-1; ++iY) {
            for (int iZ=1; iZ<nz-1; ++iZ) {
                vortField.get(iX,iY,iZ)[0] = bulkVorticityX(iX,iY,iZ);
                vortField.get(iX,iY,iZ)[1] = bulkVorticityY(iX,iY,iZ);
                vortField.get(iX,iY,iZ)[2] = bulkVorticityZ(iX,iY,iZ);
            }
        }
    }

    for (int iX=0; iX<nx-1; ++iX) {
        for (int iY=0; iY<ny-1; ++iY) {
            vortField.get(iX,iY,0)[0] = boundaryVorticityX(iX,iY,0);
            vortField.get(iX,iY,0)[1] = boundaryVorticityY(iX,iY,0);
            vortField.get(iX,iY,0)[2] = boundaryVorticityZ(iX,iY,0);
            vortField.get(iX,iY,nz-1)[0] = boundaryVorticityX(iX,iY,nz-1);
            vortField.get(iX,iY,nz-1)[1] = boundaryVorticityY(iX,iY,nz-1);
            vortField.get(iX,iY,nz-1)[2] = boundaryVorticityZ(iX,iY,nz-1);
        }
    }

    for (int iX=0; iX<nx-1; ++iX) {
        for (int iZ=0; iZ<nz-1; ++iZ) {
            vortField.get(iX,0,iZ)[0] = boundaryVorticityX(iX,0,iZ);
            vortField.get(iX,0,iZ)[1] = boundaryVorticityY(iX,0,iZ);
            vortField.get(iX,0,iZ)[2] = boundaryVorticityZ(iX,0,iZ);
            vortField.get(iX,ny-1,iZ)[0] = boundaryVorticityX(iX,ny-1,iZ);
            vortField.get(iX,ny-1,iZ)[1] = boundaryVorticityY(iX,ny-1,iZ);
            vortField.get(iX,ny-1,iZ)[2] = boundaryVorticityZ(iX,ny-1,iZ);
        }
    }

    for (int iY=0; iY<ny-1; ++iY) {
        for (int iZ=0; iZ<nz-1; ++iZ) {
            vortField.get(0,iY,iZ)[0] = boundaryVorticityX(0,iY,iZ);
            vortField.get(0,iY,iZ)[1] = boundaryVorticityY(0,iY,iZ);
            vortField.get(0,iY,iZ)[2] = boundaryVorticityZ(0,iY,iZ);
            vortField.get(nx-1,iY,iZ)[0] = boundaryVorticityX(nx-1,iY,iZ);
            vortField.get(nx-1,iY,iZ)[1] = boundaryVorticityY(nx-1,iY,iZ);
            vortField.get(nx-1,iY,iZ)[2] = boundaryVorticityZ(nx-1,iY,iZ);
        }
    }

    vortFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computeStrainRateField() const {
    if (strainRateFieldComputed) return;
    strainRateField.construct();
    computeMomentumField();

    int nx = vortField.getNx();
    int ny = vortField.getNy();
    int nz = vortField.getNz();

    int iPi = 0;
    for (int iAlpha=0; iAlpha<3; ++iAlpha) {
        for (int iBeta=iAlpha; iBeta<3; ++iBeta) {

            for (int iX=1; iX<nx-1; ++iX) {
                for (int iY=1; iY<ny-1; ++iY) {
                    for (int iZ=1; iZ<nz-1; ++iZ) {
                        strainRateField.get(iX,iY,iZ)[iPi] =
                            bulkStrain(iX,iY,iZ, iAlpha,iBeta);
                    }
                }
            }

            for (int iX=0; iX<nx; ++iX) {
                for (int iY=0; iY<ny; ++iY) {
                    strainRateField.get(iX,iY,0)[iPi] =
                        boundaryStrain(iX,iY,0, iAlpha,iBeta);
                    strainRateField.get(iX,iY,nz-1)[iPi] =
                        boundaryStrain(iX,iY,nz-1, iAlpha,iBeta);
                }
            }

            for (int iX=0; iX<nx; ++iX) {
                for (int iZ=0; iZ<nz; ++iZ) {
                    strainRateField.get(iX,0,iZ)[iPi] =
                        boundaryStrain(iX,0,iZ, iAlpha,iBeta);
                    strainRateField.get(iX,ny-1,iZ)[iPi] =
                        boundaryStrain(iX,ny-1,iZ, iAlpha,iBeta);
                }
            }

            for (int iY=0; iY<ny; ++iY) {
                for (int iZ=0; iZ<nz; ++iZ) {
                    strainRateField.get(0,iY,iZ)[iPi] =
                        boundaryStrain(0,iY,iZ, iAlpha,iBeta);
                    strainRateField.get(nx-1,iY,iZ)[iPi] =
                        boundaryStrain(nx-1,iY,iZ, iAlpha,iBeta);
                }
            }

            ++iPi;
        }
    }

    strainRateFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computeStressField() const {
    if (stressFieldComputed) return;
    stressField.construct();

    for (int iX=0; iX<velField.getNx(); ++iX) {
        for (int iY=0; iY<velField.getNy(); ++iY) {
            for (int iZ=0; iZ<velField.getNz(); ++iZ) {
                block.get(iX,iY,iZ).computeStress (
                        stressField.get(iX,iY,iZ) );
                T omega = block.get(iX,iY,iZ).getDynamics()->getOmega();
                for (int iPi=0; iPi<6; ++iPi) {
                    stressField.get(iX,iY,iZ)[iPi] *=
                        -omega / (T)2 * Lattice<T>::invCs2;
                }
            }
        }
    }

    stressFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computeDivRhoUField() const {
    if (divRhoUFieldComputed) return;
    divRhoUField.construct();
    computeMomentumField();

    int nx = divRhoUField.getNx();
    int ny = divRhoUField.getNy();
    int nz = divRhoUField.getNz();

    for (int iX=1; iX<nx-1; ++iX) {
        for (int iY=1; iY<ny-1; ++iY) {
            for (int iZ=1; iZ<nz-1; ++iZ) {
                divRhoUField.get(iX,iY,iZ) = bulkDivRhoU(iX,iY,iZ);
            }
        }
    }

    for (int iX=0; iX<nx-1; ++iX) {
        for (int iY=0; iY<ny-1; ++iY) {
            divRhoUField.get(iX,iY,0) = boundaryDivRhoU(iX,iY,0);
            divRhoUField.get(iX,iY,nz-1) = boundaryDivRhoU(iX,iY,nz-1);
        }
    }

    for (int iX=0; iX<nx-1; ++iX) {
        for (int iZ=0; iZ<nz-1; ++iZ) {
            divRhoUField.get(iX,0,iZ) = boundaryDivRhoU(iX,0,iZ);
            divRhoUField.get(iX,ny-1,iZ) = boundaryDivRhoU(iX,ny-1,iZ);
        }
    }

    for (int iY=0; iY<ny-1; ++iY) {
        for (int iZ=0; iZ<nz-1; ++iZ) {
            divRhoUField.get(0,iY,iZ) = boundaryDivRhoU(0,iY,iZ);
            divRhoUField.get(nx-1,iY,iZ) = boundaryDivRhoU(nx-1,iY,iZ);
        }
    }

    divRhoUFieldComputed = true;
}


template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computePoissonTerm() const {
    if (poissonFieldComputed) return;
    poissonField.construct();
    computeVelocityField();

    int nx = poissonField.getNx();
    int ny = poissonField.getNy();
    int nz = poissonField.getNz();

    for (int iX=1; iX<nx-1; ++iX) {
        for (int iY=1; iY<ny-1; ++iY) {
            for (int iZ=1; iZ<nz-1; ++iZ) {
                poissonField.get(iX,iY,iZ) = bulkPoisson(iX,iY,iZ);
            }
        }
    }

    for (int iX=0; iX<nx-1; ++iX) {
        for (int iY=0; iY<ny-1; ++iY) {
            poissonField.get(iX,iY,0) = boundaryPoisson(iX,iY,0);
            poissonField.get(iX,iY,nz-1) = boundaryPoisson(iX,iY,nz-1);
        }
    }

    for (int iX=0; iX<nx-1; ++iX) {
        for (int iZ=0; iZ<nz-1; ++iZ) {
            poissonField.get(iX,0,iZ) = boundaryPoisson(iX,0,iZ);
            poissonField.get(iX,ny-1,iZ) = boundaryPoisson(iX,ny-1,iZ);
        }
    }

    for (int iY=0; iY<ny-1; ++iY) {
        for (int iZ=0; iZ<nz-1; ++iZ) {
            poissonField.get(0,iY,iZ) = boundaryPoisson(0,iY,iZ);
            poissonField.get(nx-1,iY,iZ) = boundaryPoisson(nx-1,iY,iZ);
        }
    }

    poissonFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics3D<T,Lattice>::computePopulations() const {
    if (populationFieldComputed) return;
    populationField.construct();

    int nx = populationField.getNx();
    int ny = populationField.getNy();
    int nz = populationField.getNz();

    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            for (int iZ=0; iZ<nz; ++iZ) {
                for (int iPop=0; iPop<Lattice<T>::q; ++iPop) {
                    populationField.get(iX,iY,iZ)[iPop] = block.get(iX,iY,iZ)[iPop];
                }
            }
        }
    }

    populationFieldComputed = true;
}



}  // namespace olb


#endif

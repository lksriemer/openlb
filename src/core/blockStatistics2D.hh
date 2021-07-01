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

#ifndef BLOCK_STATISTICS_2D_HH
#define BLOCK_STATISTICS_2D_HH

#include <cmath>
#include "cell.h"
#include "blockStatistics2D.h"
#include "finiteDifference.h"
#include "util.h"

namespace olb {


/////// Class BlockStatistics2D  /////////////////////////////

template<typename T, template<typename U> class Lattice>
BlockStatistics2D<T,Lattice>::BlockStatistics2D (
    BlockStructure2D<T,Lattice> const& block_ )
    : block(block_),
      velField(block.getNx(), block.getNy()),
      momentumField(block.getNx(), block.getNy()),
      pressureField(block.getNx(), block.getNy()),
      velNormField(block.getNx(), block.getNy()),
      vortField(block.getNx(), block.getNy()),
      strainRateField(block.getNx(), block.getNy()),
      stressField(block.getNx(), block.getNy()),
      divRhoUField(block.getNx(), block.getNy()),
      poissonField(block.getNx(), block.getNy()),
      populationField(block.getNx(), block.getNy())
{
    iniFlags();
}

template<typename T, template<typename U> class Lattice>
BlockStatistics2D<T,Lattice>::~BlockStatistics2D() {
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::reset() const {
    iniFlags();
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::iniFlags() const {
    velFieldComputed = false;
    momentumFieldComputed = false;
    pressureFieldComputed = false;
    velNormFieldComputed = false;
    vortFieldComputed = false;
    strainRateFieldComputed = false;
    stressFieldComputed = false;
    divRhoUFieldComputed = false;
    poissonFieldComputed = false;
    populationFieldComputed = false;
}

template<typename T, template<typename U> class Lattice>
TensorField2D<T,2> const&
    BlockStatistics2D<T,Lattice>::getVelocity() const
{
    computeVelocityField();
    return velField;
}

template<typename T, template<typename U> class Lattice>
TensorField2D<T,2> const&
    BlockStatistics2D<T,Lattice>::getMomentum() const
{
    computeMomentumField();
    return momentumField;
}

template<typename T, template<typename U> class Lattice>
ScalarField2D<T> const&
    BlockStatistics2D<T,Lattice>::getPressure() const
{
    computePressureField();
    return pressureField;
}

template<typename T, template<typename U> class Lattice>
ScalarField2D<T> const&
    BlockStatistics2D<T,Lattice>::getVelocityNorm() const
{
    computeVelocityNormField();
    return velNormField;
}

template<typename T, template<typename U> class Lattice>
ScalarField2D<T> const&
    BlockStatistics2D<T,Lattice>::getVorticity() const
{
    computeVorticityField();
    return vortField;
}

template<typename T, template<typename U> class Lattice>
TensorField2D<T,3> const&
    BlockStatistics2D<T,Lattice>::getStrainRate() const
{
    computeStrainRateField();
    return strainRateField;
}

template<typename T, template<typename U> class Lattice>
TensorField2D<T,3> const&
    BlockStatistics2D<T,Lattice>::getStrainRateFromStress() const
{
    computeStrainRateFieldFromStress();
    return stressField;
}

template<typename T, template<typename U> class Lattice>
ScalarField2D<T> const&
    BlockStatistics2D<T,Lattice>::getDivRhoU() const
{
    computeDivRhoUField();
    return divRhoUField;
}

template<typename T, template<typename U> class Lattice>
ScalarField2D<T> const&
    BlockStatistics2D<T,Lattice>::getPoissonTerm() const
{
    computePoissonTerm();
    return poissonField;
}

template<typename T, template<typename U> class Lattice>
TensorField2D<T, Lattice<T>::q > const&
    BlockStatistics2D<T,Lattice>::getPopulations() const
{
    computePopulationField();
    return populationField;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::computeMeanEnstrophy() const {
    computeVelocityField();
    int nx = velField.getNx()-1;
    int ny = velField.getNy()-1;

    T enstrophy = T();
    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
           T dyux = (
                 velField.get(iX  ,iY+1)[0]
               + velField.get(iX+1,iY+1)[0]
               - velField.get(iX  ,iY  )[0]
               - velField.get(iX+1,iY  )[0] ) / (T)2;
           T dxuy = (
                 velField.get(iX+1,iY  )[1]
               + velField.get(iX+1,iY+1)[1]
               - velField.get(iX  ,iY  )[1]
               - velField.get(iX  ,iY+1)[1] ) / (T)2;
           T omega = dxuy - dyux;
           enstrophy += omega*omega;
        }
    }
    enstrophy /= (2*nx*ny);
    return enstrophy;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::computeMeanEnstrophy2() const {
    computeVorticityField();
    int nx = vortField.getNx();
    int ny = vortField.getNy();

    T enstrophy = T();
    for (int iX=1; iX<nx-1; ++iX) {
        for (int iY=1; iY<ny-1; ++iY) {
           enstrophy += util::sqr(vortField.get(iX,iY));
        }
    }

    for (int iX=1; iX<nx-1; ++iX) {
        enstrophy += 0.5* (
                util::sqr(vortField.get(iX,0)) +
                util::sqr(vortField.get(iX,ny-1)) );
    }
    for (int iY=1; iY<ny-1; ++iY) {
        enstrophy += 0.5* (
                util::sqr(vortField.get(0,iY)) +
                util::sqr(vortField.get(nx-1,iY)) );
    }
    enstrophy += 0.25 * util::sqr(vortField.get(0,0));
    enstrophy += 0.25 * util::sqr(vortField.get(0,ny-1));
    enstrophy += 0.25 * util::sqr(vortField.get(nx-1,0));
    enstrophy += 0.25 * util::sqr(vortField.get(nx-1,ny-1));

    enstrophy /= 2*(nx-1)*(ny-1);

    return enstrophy;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computeVelocityField() const {
    if (velFieldComputed) return;
    velField.construct();
    for (int iX=0; iX<velField.getNx(); ++iX) {
        for (int iY=0; iY<velField.getNy(); ++iY) {
            block.get(iX,iY).computeU(velField.get(iX,iY));
        }
    }
    velFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computeMomentumField() const {
    if (momentumFieldComputed) return;
    momentumField.construct();
    for (int iX=0; iX<momentumField.getNx(); ++iX) {
        for (int iY=0; iY<momentumField.getNy(); ++iY) {
            T rho;
            block.get(iX,iY).computeRhoU(rho, momentumField.get(iX,iY));
            for (int iD=0; iD<Lattice<T>::d; ++iD) {
                momentumField.get(iX,iY)[iD] *= rho;
            }
        }
    }
    momentumFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computePressureField() const {
    if (pressureFieldComputed) return;
    pressureField.construct();
    for (int iX=0; iX<pressureField.getNx(); ++iX) {
        for (int iY=0; iY<pressureField.getNy(); ++iY) {
            pressureField.get(iX,iY) = block.get(iX,iY).computeRho();
            pressureField.get(iX,iY) -= (T)1;
            pressureField.get(iX,iY) /= Lattice<T>::invCs2;
        }
    }
    pressureFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computeVelocityNormField() const {
    if (velNormFieldComputed) return;
    velNormField.construct();
    computeVelocityField();
    for (int iEl=0; iEl<velNormField.getNx()*velNormField.getNy(); ++iEl) {
        velNormField[iEl] = sqrt(util::normSqr<T,2>(velField[iEl]));
    }
    velNormFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::bulkVorticity(int iX, int iY) const {
    OLB_PRECONDITION(velFieldComputed);
    OLB_PRECONDITION(iX>=1 && iX<=vortField.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=vortField.getNy()-2);

    T dyux = fd::centralGradient(velField.get(iX,iY+1)[0],
                                 velField.get(iX,iY-1)[0]);
    T dxuy = fd::centralGradient(velField.get(iX+1,iY)[1],
                                 velField.get(iX-1,iY)[1]);
    return dxuy - dyux;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::bulkXderiv (
        int iX, int iY, int iD, TensorField2D<T,2> const& field) const
{
    OLB_PRECONDITION(iX>=1 && iX<=field.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=field.getNy()-2);

    T dxu = fd::centralGradient(field.get(iX+1,iY)[iD],
                                field.get(iX-1,iY)[iD]);
    return dxu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::bulkYderiv (
        int iX, int iY, int iD, TensorField2D<T,2> const& field) const
{
    OLB_PRECONDITION(iX>=1 && iX<=field.getNx()-2);
    OLB_PRECONDITION(iY>=1 && iY<=field.getNy()-2);

    T dyu = fd::centralGradient(field.get(iX,iY+1)[iD],
                                field.get(iX,iY-1)[iD]);
    return dyu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::bulkDeriv (
    int iX, int iY, int iAlpha, int iBeta,
    TensorField2D<T,2> const& field) const
{
    switch(iAlpha) {
        case 0:
            return bulkXderiv(iX,iY,iBeta,field);
        case 1:
            return bulkYderiv(iX,iY,iBeta,field);
        default:
            OLB_ASSERT( false, "iAlpha>1!");
            return T();
    }
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::bulkStrain (
    int iX, int iY, int iAlpha, int iBeta) const
{
    OLB_PRECONDITION( momentumFieldComputed );
    return ( bulkDeriv(iX,iY,iAlpha,iBeta, momentumField) +
             bulkDeriv(iX,iY,iBeta,iAlpha, momentumField) ) / (T)2;

}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::bulkDivRhoU(int iX, int iY) const {
    OLB_PRECONDITION( momentumFieldComputed );
    return bulkDeriv(iX,iY, 0,0, momentumField) +
           bulkDeriv(iX,iY, 1,1, momentumField);
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::bulkPoisson(int iX, int iY) const {
    OLB_PRECONDITION( velFieldComputed );
    T dxux = bulkDeriv(iX,iY,0,0, velField);
    T dxuy = bulkDeriv(iX,iY,0,1, velField);
    T dyux = bulkDeriv(iX,iY,1,0, velField);
    T dyuy = bulkDeriv(iX,iY,1,1, velField);
    return dxux*dxux + (T)2*dxuy*dyux + dyuy*dyuy;
}


template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::boundaryXderiv (
        int iX, int iY, int iD, TensorField2D<T,2> const& field ) const
{
    OLB_PRECONDITION(iX>=0 && iX<=field.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=field.getNy()-1);

    T dxu;
    
    if (iX==0) {
        dxu = fd::boundaryGradient(field.get(iX,iY)[iD],
                                   field.get(iX+1,iY)[iD],
                                   field.get(iX+2,iY)[iD]);
    }
    else if (iX==vortField.getNx()-1) {
        dxu = -fd::boundaryGradient(field.get(iX,iY)[iD],
                                    field.get(iX-1,iY)[iD],
                                    field.get(iX-2,iY)[iD]);
    }
    else {
        dxu = fd::centralGradient(field.get(iX+1,iY)[iD],
                                  field.get(iX-1,iY)[iD]);
    }

    return dxu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::boundaryYderiv (
        int iX, int iY, int iD, TensorField2D<T,2> const& field ) const
{
    OLB_PRECONDITION(iX>=0 && iX<=field.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=field.getNy()-1);

    T dyu;

    if (iY==0) {
        dyu = fd::boundaryGradient(field.get(iX,iY)[iD],
                                   field.get(iX,iY+1)[iD],
                                   field.get(iX,iY+2)[iD]);
    }
    else if (iY==vortField.getNy()-1) {
        dyu = -fd::boundaryGradient(field.get(iX,iY)[iD],
                                    field.get(iX,iY-1)[iD],
                                    field.get(iX,iY-2)[iD]);
    }
    else {
        dyu = fd::centralGradient(field.get(iX,iY+1)[iD],
                                  field.get(iX,iY-1)[iD]);
    }
   
    return dyu;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::boundaryDeriv (
    int iX, int iY, int iAlpha, int iBeta,
    TensorField2D<T,2> const& field ) const
{
    switch(iAlpha) {
        case 0:
            return boundaryXderiv(iX,iY,iBeta, field);
        case 1:
            return boundaryYderiv(iX,iY,iBeta, field);
        default:
            OLB_ASSERT( false, "iAlpha>1!");
            return T();
    }
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::boundaryStrain (
    int iX, int iY, int iAlpha, int iBeta) const
{
    OLB_PRECONDITION( momentumFieldComputed );
    return ( boundaryDeriv(iX,iY,iAlpha,iBeta, momentumField) +
             boundaryDeriv(iX,iY,iBeta,iAlpha, momentumField) ) / (T)2;

}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::boundaryDivRhoU(int iX, int iY) const {
    OLB_PRECONDITION( momentumFieldComputed );
    return boundaryDeriv(iX,iY,0,0, momentumField) +
           boundaryDeriv(iX,iY,1,1, momentumField);
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::boundaryPoisson(int iX, int iY) const {
    OLB_PRECONDITION( velFieldComputed );
    T dxux = boundaryDeriv(iX,iY,0,0, velField);
    T dxuy = boundaryDeriv(iX,iY,0,1, velField);
    T dyux = boundaryDeriv(iX,iY,1,0, velField);
    T dyuy = boundaryDeriv(iX,iY,1,1, velField);
    return dxux*dxux + (T)2*dxuy*dyux + dyuy*dyuy;
}

template<typename T, template<typename U> class Lattice>
T BlockStatistics2D<T,Lattice>::boundaryVorticity(int iX, int iY) const {
    OLB_PRECONDITION(velFieldComputed);
    OLB_PRECONDITION(iX>=0 && iX<=vortField.getNx()-1);
    OLB_PRECONDITION(iY>=0 && iY<=vortField.getNy()-1);

    T dyux = boundaryYderiv(iX,iY, 0, velField);
    T dxuy = boundaryXderiv(iX,iY, 1, velField);

    return dxuy - dyux;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computeVorticityField() const {
    if (vortFieldComputed) return;
    vortField.construct();
    computeVelocityField();

    int nx = vortField.getNx();
    int ny = vortField.getNy();

    for (int iX=1; iX<nx-1; ++iX) {
        for (int iY=1; iY<ny-1; ++iY) {
            vortField.get(iX,iY) = bulkVorticity(iX,iY);
        }
    }

    for (int iX=1; iX<nx-1; ++iX) {
        vortField.get(iX,0) = boundaryVorticity(iX,0);
        vortField.get(iX,ny-1) = boundaryVorticity(iX,ny-1);
    }

    for (int iY=0; iY<ny; ++iY) {
        vortField.get(0,iY) = boundaryVorticity(0,iY);
        vortField.get(nx-1,iY) = boundaryVorticity(nx-1,iY);
    }

    vortFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computeStrainRateField() const {
    if (strainRateFieldComputed) return;
    strainRateField.construct();
    computeMomentumField();

    int nx = vortField.getNx();
    int ny = vortField.getNy();

    int iPi = 0;
    for (int iAlpha=0; iAlpha<2; ++iAlpha) {
        for (int iBeta=iAlpha; iBeta<2; ++iBeta) {

            for (int iX=1; iX<nx-1; ++iX) {
                for (int iY=1; iY<ny-1; ++iY) {
                    strainRateField.get(iX,iY)[iPi] =
                        bulkStrain(iX,iY, iAlpha,iBeta);
                }
            }

            for (int iX=1; iX<nx-1; ++iX) {
                strainRateField.get(iX,0)[iPi] =
                    boundaryStrain(iX,0, iAlpha,iBeta);
                strainRateField.get(iX,ny-1)[iPi] =
                    boundaryStrain(iX,ny-1, iAlpha,iBeta);
            }

            for (int iY=0; iY<ny; ++iY) {
                strainRateField.get(0,iY)[iPi] =
                    boundaryStrain(0,iY, iAlpha,iBeta);
                strainRateField.get(nx-1,iY)[iPi] =
                    boundaryStrain(nx-1,iY, iAlpha,iBeta);
            }

            ++iPi;
        }
    }

    strainRateFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computeStrainRateFieldFromStress() const {
    if (stressFieldComputed) return;
    stressField.construct();

    for (int iX=0; iX<velField.getNx(); ++iX) {
        for (int iY=0; iY<velField.getNy(); ++iY) {
            block.get(iX,iY).computeStress(stressField.get(iX,iY));
            T omega = block.get(iX,iY).getDynamics()->getOmega();
            for (int iPi=0; iPi<3; ++iPi) {
                stressField.get(iX,iY)[iPi] *=
                    -omega / (T)2 * Lattice<T>::invCs2;
            }
        }
    }

    stressFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computeDivRhoUField() const {
    if (divRhoUFieldComputed) return;
    divRhoUField.construct();
    computeMomentumField();

    int nx = divRhoUField.getNx();
    int ny = divRhoUField.getNy();

    for (int iX=1; iX<nx-1; ++iX) {
        for (int iY=1; iY<ny-1; ++iY) {
            divRhoUField.get(iX,iY) = bulkDivRhoU(iX,iY);
        }
    }

    for (int iX=1; iX<nx-1; ++iX) {
        divRhoUField.get(iX,0) = boundaryDivRhoU(iX,0);
        divRhoUField.get(iX,ny-1) = boundaryDivRhoU(iX,ny-1);
    }

    for (int iY=0; iY<ny; ++iY) {
        divRhoUField.get(0,iY) = boundaryDivRhoU(0,iY);
        divRhoUField.get(nx-1,iY) = boundaryDivRhoU(nx-1,iY);
    }

    divRhoUFieldComputed = true;
}


template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computePoissonTerm() const {
    if (poissonFieldComputed) return;
    poissonField.construct();
    computeVelocityField();

    int nx = poissonField.getNx();
    int ny = poissonField.getNy();

    for (int iX=1; iX<nx-1; ++iX) {
        for (int iY=1; iY<ny-1; ++iY) {
            poissonField.get(iX,iY) = bulkPoisson(iX,iY);
        }
    }

    for (int iX=1; iX<nx-1; ++iX) {
        poissonField.get(iX,0) = boundaryPoisson(iX,0);
        poissonField.get(iX,ny-1) = boundaryPoisson(iX,ny-1);
    }

    for (int iY=0; iY<ny; ++iY) {
        poissonField.get(0,iY) = boundaryPoisson(0,iY);
        poissonField.get(nx-1,iY) = boundaryPoisson(nx-1,iY);
    }

    poissonFieldComputed = true;
}

template<typename T, template<typename U> class Lattice>
void BlockStatistics2D<T,Lattice>::computePopulationField() const {
    if (populationFieldComputed) return;
    populationField.construct();

    int nx = populationField.getNx();
    int ny = populationField.getNy();
    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            for (int iPop=0; iPop<Lattice<T>::q; ++iPop) {
                populationField.get(iX,iY)[iPop] = block.get(iX,iY)[iPop];
            }
        }
    }

    populationFieldComputed = true;
}

}  // namespace olb


#endif

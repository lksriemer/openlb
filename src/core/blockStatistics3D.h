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

#ifndef BLOCK_STATISTICS_3D_H
#define BLOCK_STATISTICS_3D_H

#include "blockLattice3D.h"
#include "blockLatticeView3D.h"
#include "dataFields3D.h"

namespace olb {

template<typename T, template<typename U> class Lattice>
class BlockStatistics3D {
public:
    BlockStatistics3D(BlockStructure3D<T,Lattice> const& block_);
    ~BlockStatistics3D();

    void reset() const;

    TensorField3D<T,3> const& getVelocity() const;
    TensorField3D<T,3> const& getMomentum() const;
    ScalarField3D<T>   const& getPressure() const;
    TensorField3D<T,3> const& getVorticity() const;
    ScalarField3D<T>   const& getVelocityNorm() const;
    ScalarField3D<T>   const& getVorticityNorm() const;
    TensorField3D<T,6> const& getStrainRate() const;
    TensorField3D<T,6> const& getStrainRateFromStress() const;
    ScalarField3D<T>   const& getDivRhoU() const;
    ScalarField3D<T>   const& getPoissonTerm() const;
    TensorField3D<T,Lattice<T>::q > const& getPopulations() const;

    T computeMeanEnstrophy() const;

    int getNx() const { return block.getNx(); }
    int getNy() const { return block.getNy(); }
    int getNz() const { return block.getNz(); }
private:
    void iniFlags() const;
    void computeVelocityField() const;
    void computeMomentumField() const;
    void computePressureField() const;
    void computeVelocityNormField() const;
    void computeVorticityNormField() const;
    void computeVorticityField() const;
    void computeStrainRateField() const;
    void computeStrainRateFieldFromStress() const;
    void computeDivRhoUField() const;
    void computePoissonTerm() const;
    void computePopulations() const;
    T bulkVorticityX(int iX, int iY, int iZ) const;
    T bulkVorticityY(int iX, int iY, int iZ) const;
    T bulkVorticityZ(int iX, int iY, int iZ) const;
    T boundaryVorticityX(int iX, int iY, int iZ) const;
    T boundaryVorticityY(int iX, int iY, int iZ) const;
    T boundaryVorticityZ(int iX, int iY, int iZ) const;
    T bulkXderiv(int iX, int iY, int iZ, int iD, TensorField3D<T,3> const& field) const;
    T bulkYderiv(int iX, int iY, int iZ, int iD, TensorField3D<T,3> const& field) const;
    T bulkZderiv(int iX, int iY, int iZ, int iD, TensorField3D<T,3> const& field) const;
    T bulkDeriv(int iX, int iY, int iZ, int iAlpha, int iBeta, TensorField3D<T,3> const& field) const;
    T bulkStrain(int iX, int iY, int iZ, int iAlpha, int iBeta) const;
    T bulkDivRhoU(int iX, int iY, int iZ) const;
    T bulkPoisson(int iX, int iY, int iZ) const;
    T boundaryXderiv(int iX, int iY, int iZ, int iD, TensorField3D<T,3> const& field) const;
    T boundaryYderiv(int iX, int iY, int iZ, int iD, TensorField3D<T,3> const& field) const;
    T boundaryZderiv(int iX, int iY, int iZ, int iD, TensorField3D<T,3> const& field) const;
    T boundaryDeriv(int iX, int iY, int iZ, int iAlpha, int iBeta, TensorField3D<T,3> const& field) const;
    T boundaryStrain(int iX, int iY, int iZ, int iAlpha, int iBeta) const;
    T boundaryDivRhoU(int iX, int iY, int iZ) const;
    T boundaryPoisson(int iX, int iY, int iZ) const;
private:
    BlockStructure3D<T,Lattice> const& block;
    mutable TensorField3D<T,3> velField;
    mutable TensorField3D<T,3> momentumField;
    mutable ScalarField3D<T>   pressureField;
    mutable ScalarField3D<T>   velNormField;
    mutable ScalarField3D<T>   vortNormField;
    mutable TensorField3D<T,3> vortField;
    mutable TensorField3D<T,6> strainRateField;
    mutable TensorField3D<T,6> stressField;
    mutable ScalarField3D<T>   divRhoUField;
    mutable ScalarField3D<T>   poissonField;
    mutable TensorField3D<T, Lattice<T>::q > populationField;
    mutable bool velFieldComputed;
    mutable bool momentumFieldComputed;
    mutable bool pressureFieldComputed;
    mutable bool velNormFieldComputed;
    mutable bool vortNormFieldComputed;
    mutable bool vortFieldComputed;
    mutable bool strainRateFieldComputed;
    mutable bool stressFieldComputed;
    mutable bool divRhoUFieldComputed;
    mutable bool poissonFieldComputed;
    mutable bool populationFieldComputed;
};


}  // namespace olb;

#endif

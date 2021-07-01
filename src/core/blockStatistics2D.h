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

#ifndef BLOCK_STATISTICS_2D_H
#define BLOCK_STATISTICS_2D_H

#include "blockLattice2D.h"
#include "blockLatticeView2D.h"
#include "dataFields2D.h"

namespace olb {

template<typename T, template<typename U> class Lattice>
class BlockStatistics2D {
public:
    BlockStatistics2D(BlockStructure2D<T,Lattice> const& block_);
    ~BlockStatistics2D();

    void reset() const;

    TensorField2D<T,2> const& getVelocity() const;
    TensorField2D<T,2> const& getMomentum() const;
    ScalarField2D<T> const& getPressure() const;
    ScalarField2D<T> const& getVelocityNorm() const;
    ScalarField2D<T> const& getVorticity() const;
    TensorField2D<T,3> const& getStrainRate() const;
    TensorField2D<T,3> const& getStress() const;
    ScalarField2D<T> const& getDivRhoU() const;
    ScalarField2D<T> const& getPoissonTerm() const;
    TensorField2D<T,Lattice<T>::q > const& getPopulations() const;

    T computeMeanEnstrophy() const;
    T computeMeanEnstrophy2() const;

    int getNx() const { return block.getNx(); }
    int getNy() const { return block.getNy(); }
private:
    void iniFlags() const;
    void computeVelocityField() const;
    void computeMomentumField() const;
    void computePressureField() const;
    void computeVelocityNormField() const;
    void computeVorticityField() const;
    void computeStrainRateField() const;
    void computeStrainRateFieldFromStress() const;
    void computeDivRhoUField() const;
    void computePoissonTerm() const;
    void computePopulationField() const;
    T bulkVorticity(int iX, int iY) const;
    T boundaryVorticity(int iX, int iY) const;
    T bulkXderiv(int iX, int iY, int iD, TensorField2D<T,2> const& field) const;
    T bulkYderiv(int iX, int iY, int iD, TensorField2D<T,2> const& field) const;
    T bulkDeriv(int iX, int iY, int iAlpha, int iBeta, TensorField2D<T,2> const& field) const;
    T bulkStrain(int iX, int iY, int iAlpha, int iBeta) const;
    T bulkDivRhoU(int iX, int iY) const;
    T bulkPoisson(int iX, int iY) const;
    T boundaryXderiv(int iX, int iY, int iD, TensorField2D<T,2> const& field) const;
    T boundaryYderiv(int iX, int iY, int iD, TensorField2D<T,2> const& field) const;
    T boundaryDeriv(int iX, int iY, int iAlpha, int iBeta, TensorField2D<T,2> const& field) const;
    T boundaryStrain(int iX, int iY, int iAlpha, int iBeta) const;
    T boundaryDivRhoU(int iX, int iY) const;
    T boundaryPoisson(int iX, int iY) const;
private:
    BlockStructure2D<T,Lattice> const& block;
    mutable TensorField2D<T,2> velField;
    mutable TensorField2D<T,2> momentumField;
    mutable ScalarField2D<T>   pressureField;
    mutable ScalarField2D<T>   velNormField;
    mutable ScalarField2D<T>   vortField;
    mutable TensorField2D<T,3> strainRateField;
    mutable TensorField2D<T,3> stressField;
    mutable ScalarField2D<T>   divRhoUField;
    mutable ScalarField2D<T>   poissonField;
    mutable TensorField2D<T,Lattice<T>::q > populationField;
    mutable bool velFieldComputed;
    mutable bool momentumFieldComputed;
    mutable bool pressureFieldComputed;
    mutable bool velNormFieldComputed;
    mutable bool vortFieldComputed;
    mutable bool strainRateFieldComputed;
    mutable bool stressFieldComputed;
    mutable bool divRhoUFieldComputed;
    mutable bool poissonFieldComputed;
    mutable bool populationFieldComputed;
};


}  // namespace olb;

#endif

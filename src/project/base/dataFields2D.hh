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

/** \file
 * Scalar, vector and tensor fields for 2D data analysis -- generic implementation.
 */

#ifndef DATA_FIELDS_2D_HH
#define DATA_FIELDS_2D_HH

#include <cmath>
#include <algorithm>
#include <limits>
#include "dataFields2D.h"

namespace olb {

/////// Class ScalarField2D //////////////////////////////////

template<typename T>
ScalarField2D<T>::ScalarField2D(int nx_, int ny_)
    : nx(nx_), ny(ny_), rawData(0), field(0)
{ }

template<typename T>
ScalarField2D<T>::~ScalarField2D() {
    deConstruct();
}

template<typename T>
ScalarField2D<T>::ScalarField2D(ScalarField2D<T> const& rhs) {
    nx = rhs.nx;
    ny = rhs.ny;
    rawData = 0;
    field   = 0;
    if (rhs.isConstructed()) {
        construct();
        for (int iData=0; iData<getSize(); ++iData) {
            (*this)[iData] = rhs[iData];
        }
    }
}

template<typename T>
ScalarField2D<T>& ScalarField2D<T>::operator=(ScalarField2D<T> const& rhs) {
    ScalarField2D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
bool ScalarField2D<T>::isConstructed() const {
    return rawData;
}

template<typename T>
void ScalarField2D<T>::construct() {
    if (!isConstructed()) {
        allocateMemory();
    }
}

template<typename T>
void ScalarField2D<T>::deConstruct() {
    if (isConstructed()) {
        releaseMemory();
    }
}

template<typename T>
void ScalarField2D<T>::reset() {
    OLB_PRECONDITION(isConstructed());
    for (int index=0; index<nx*ny; ++index) {
        (*this)[index] = T();
    }
}


template<typename T>
void ScalarField2D<T>::swap(ScalarField2D<T>& rhs) {
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}

template<typename T>
void ScalarField2D<T>::allocateMemory() {
    rawData = new T[nx*ny];
    field   = new T* [nx];
    for (int iX=0; iX<nx; ++iX) {
        field[iX] = rawData + iX*ny;
    }
}

template<typename T>
void ScalarField2D<T>::releaseMemory() {
    delete [] rawData; rawData = 0;
    delete [] field;
}


//////// Class TensorField2D //////////////////////////////////

template<typename T, int nDim>
TensorField2D<T,nDim>::TensorField2D(int nx_, int ny_)
    : nx(nx_), ny(ny_), rawData(0), field(0)
{
    for (int iDim=0; iDim<nDim; ++iDim) {
        components[iDim] = new ScalarField2D<T>(nx, ny);
    }
}

template<typename T, int nDim>
TensorField2D<T,nDim>::~TensorField2D() {
    deConstruct();
    for (int iDim=0; iDim<nDim; ++iDim) {
        delete components[iDim];
    }
}

template<typename T, int nDim>
TensorField2D<T,nDim>::TensorField2D(TensorField2D<T,nDim> const& rhs) {
    nx = rhs.nx;
    ny = rhs.ny;
    rawData = 0;
    field   = 0;
    for (int iDim=0; iDim<nDim; ++iDim) {
        components[iDim] = new ScalarField2D<T>(nx, ny);
    }
    if (rhs.isConstructed()) {
        construct();
        for (int iData=0; iData<nx*ny; ++iData) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                (*this)[iData][iDim] = rhs[iData][iDim];
            }
        }
    }
}

template<typename T, int nDim>
TensorField2D<T,nDim>& TensorField2D<T,nDim>::operator=(TensorField2D<T,nDim> const& rhs) {
    TensorField2D<T,nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T, int nDim>
bool TensorField2D<T,nDim>::isConstructed() const {
    return rawData;
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::construct() {
    if (!isConstructed()) {
        allocateMemory();
    }
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::deConstruct() {
    if (isConstructed()) {
        releaseMemory();
    }
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::reset() {
    OLB_PRECONDITION(isConstructed());
    for (int index=0; index<nx*ny; ++index) {
        for (int iDim=0; iDim<nDim; ++iDim) {
            (*this)[index][iDim] = T();
        }
    }
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::swap(TensorField2D<T,nDim>& rhs) {
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
    for (int iDim=0; iDim<nDim; ++iDim) {
        std::swap(components[iDim], rhs.components[iDim]);
    }
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::allocateMemory() {
    rawData = new Tensor[nx*ny];
    field   = new Tensor* [nx];
    for (int iX=0; iX<nx; ++iX) {
        field[iX] = rawData + iX*ny;
    }
}

template<typename T, int nDim>
void TensorField2D<T,nDim>::releaseMemory() {
    delete [] rawData; rawData = 0;
    delete [] field;
}

template<typename T, int nDim>
ScalarField2D<T> const& TensorField2D<T,nDim>::extractComponent(int whichDim)
    const
{
    components[whichDim]->construct();
    for (int iEl=0; iEl<nx*ny; ++iEl) {
        (*components[whichDim])[iEl] = (*this)[iEl][whichDim];
    };
    return *components[whichDim];
}


//////// Free functions //////////////////////////////////

template<typename T, template<typename U> class ScalarField>
T computeMin(ScalarField<T> const& field) {
    T minimum = std::numeric_limits<T>::max();
    for (int iData=0; iData<field.getSize(); ++iData) {
        if (field[iData] < minimum) {
            minimum = field[iData];
        }
    }
    return minimum;
}

template<typename T, template<typename U> class ScalarField>
T computeMax(ScalarField<T> const& field) {
    T maximum = std::numeric_limits<T>::min();
    for (int iData=0; iData<field.getSize(); ++iData) {
        if (field[iData] > maximum) {
            maximum = field[iData];
        }
    }
    return maximum;
}

template<typename T, template<typename U> class ScalarField>
T computeAverage(ScalarField<T> const& field) {
    T average = T();
    for (int iData=0; iData<field.getSize(); ++iData) {
        average += field[iData];
    }
    average /= (T)(field.getSize());
    return average;
}

template<typename T, template<typename U> class ScalarField>
T computeRMS(ScalarField<T> const& field) {
    return sqrt(computeNormSqr(field));
}

template<typename T, template<typename U> class ScalarField>
T computeNormSqr(ScalarField<T> const& field) {
    T normSqr = T();
    for (int iData=0; iData<field.getSize(); ++iData) {
        normSqr += field[iData]*field[iData];
    }
    normSqr /= (T)(field.getSize());
    return normSqr;
}

}  // namespace olb

#endif  

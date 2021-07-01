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
 * Scalar, vector and tensor fields for 3D data analysis -- generic implementation.
 */
#ifndef DATA_FIELDS_3D_HH
#define DATA_FIELDS_3D_HH

#include <cmath>
#include <algorithm>
#include <limits>
#include "dataFields3D.h"

namespace olb {

/////// Class ScalarField3D //////////////////////////////////

template<typename T>
ScalarField3D<T>::ScalarField3D(int nx_, int ny_, int nz_)
    : nx(nx_), ny(ny_), nz(nz_), rawData(0), field(0),
      xSlice(ny,nz), ySlice(nx,nz), zSlice(nx,ny)
{ }

template<typename T>
ScalarField3D<T>::~ScalarField3D() {
    deConstruct();
}

template<typename T>
ScalarField3D<T>::ScalarField3D(ScalarField3D<T> const& rhs)
    : nx(rhs.nx), ny(rhs.ny), nz(rhs.nz), rawData(0), field(0),
      xSlice(ny,nz), ySlice(nx,nz), zSlice(nx,ny)
{
    if (rhs.isConstructed()) {
        construct();
        for (int iData=0; iData<getSize(); ++iData) {
            (*this)[iData] = rhs[iData];
        }
    }
}

template<typename T>
ScalarField3D<T>& ScalarField3D<T>::operator=(ScalarField3D<T> const& rhs) {
    ScalarField3D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template<typename T>
bool ScalarField3D<T>::isConstructed() const {
    return rawData;
}

template<typename T>
void ScalarField3D<T>::construct() {
    if (!isConstructed()) {
        allocateMemory();
    }
}

template<typename T>
void ScalarField3D<T>::deConstruct() {
    if (isConstructed()) {
        releaseMemory();
    }
}

template<typename T>
void ScalarField3D<T>::reset() {
    OLB_PRECONDITION(isConstructed());
    for (int index=0; index<nx*ny; ++index) {
        (*this)[index] = T();
    }
}

template<typename T>
ScalarField2D<T>& ScalarField3D<T>::sliceX(int xVal) const {
    xSlice.construct();
    for (int iY=0; iY<ny; ++iY) {
        for (int iZ=0; iZ<nz; ++iZ) {
            xSlice.get(iY,iZ) = this->get(xVal,iY,iZ);
        }
    }
    return xSlice;
}

template<typename T>
ScalarField2D<T>& ScalarField3D<T>::sliceY(int yVal) const {
    ySlice.construct();
    for (int iX=0; iX<nx; ++iX) {
        for (int iZ=0; iZ<nz; ++iZ) {
            ySlice.get(iX,iZ) = this->get(iX,yVal,iZ);
        }
    }
    return ySlice;
}

template<typename T>
ScalarField2D<T>& ScalarField3D<T>::sliceZ(int zVal) const {
    zSlice.construct();
    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            zSlice.get(iX,iY) = this->get(iX,iY,zVal);
        }
    }
    return zSlice;
}

template<typename T>
void ScalarField3D<T>::swap(ScalarField3D<T>& rhs) {
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}

template<typename T>
void ScalarField3D<T>::allocateMemory() {
    rawData = new T [nx*ny*nz];
    field   = new T** [nx];
    for (int iX=0; iX<nx; ++iX) {
        field[iX] = new T* [ny];
        for (int iY=0; iY<ny; ++iY) {
            field[iX][iY] = rawData + nz*(iY+ny*iX);
        }
    }
}

template<typename T>
void ScalarField3D<T>::releaseMemory() {
    delete [] rawData; rawData = 0;
    for (int iX=0; iX<nx; ++iX) {
      delete [] field[iX];
    }
    delete [] field;
}


template<typename T>
T computeMin(ScalarField3D<T> const& field) {
    T minimum = std::numeric_limits<T>::max();
    for (int iData=0; iData<field.getSize(); ++iData) {
        if (field[iData] < minimum) {
            minimum = field[iData];
        }
    }
    return minimum;
}

template<typename T>
T computeMax(ScalarField3D<T> const& field) {
    T maximum = std::numeric_limits<T>::min();
    for (int iData=0; iData<field.getSize(); ++iData) {
        if (field[iData] > maximum) {
            maximum = field[iData];
        }
    }
    return maximum;
}

template<typename T>
T computeAverage(ScalarField3D<T> const& field) {
    T average = T();
    for (int iData=0; iData<field.getSize(); ++iData) {
        average += field[iData];
    }
    average /= (T)(field.getSize());
    return average;
}

template<typename T>
T computeRMS(ScalarField3D<T> const& field) {
    return sqrt(computeNormSqr(field));
}

template<typename T>
T computeNormSqr(ScalarField3D<T> const& field) {
    T normSqr = T();
    for (int iData=0; iData<field.getSize(); ++iData) {
        normSqr += field[iData]*field[iData];
    }
    normSqr /= (T)(field.getSize());
    return normSqr;
}


//////// Class TensorField3D //////////////////////////////////

template<typename T, int nDim>
TensorField3D<T,nDim>::TensorField3D(int nx_, int ny_, int nz_)
    : nx(nx_), ny(ny_), nz(nz_), rawData(0), field(0),
      xSlice(ny,nz), ySlice(nx,nz), zSlice(nx,ny)
{
    for (int iDim=0; iDim<nDim; ++iDim) {
        components[iDim] =
            new ScalarField3D<T>(nx, ny, nz);
    }
}

template<typename T, int nDim>
TensorField3D<T,nDim>::~TensorField3D() {
    deConstruct();
    for (int iDim=0; iDim<nDim; ++iDim) {
        delete components[iDim];
    }
}

template<typename T, int nDim>
TensorField3D<T,nDim>::TensorField3D(TensorField3D<T,nDim> const& rhs)
    : nx(rhs.nx), ny(rhs.ny), nz(rhs.nz), rawData(0), field(0),
      xSlice(ny,nz), ySlice(nx,nz), zSlice(nx,ny)
{
    for (int iDim=0; iDim<nDim; ++iDim) {
        components[iDim] = new ScalarField3D<T>(nx, ny, nz);
    }
    if (rhs.isConstructed()) {
        construct();
        for (int iData=0; iData<nx*ny*nz; ++iData) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                (*this)[iData][iDim] = rhs[iData][iDim];
            }
        }
    }
}

template<typename T, int nDim>
TensorField3D<T,nDim>& TensorField3D<T,nDim>::operator=(TensorField3D<T,nDim> const& rhs) {
    TensorField3D<T,nDim> tmp(rhs);
    swap(tmp);
    return *this;
}


template<typename T, int nDim>
bool TensorField3D<T,nDim>::isConstructed() const {
    return rawData;
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::construct() {
    if (!isConstructed()) {
        allocateMemory();
    }
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::deConstruct() {
    if (isConstructed()) {
        releaseMemory();
    }
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::reset() {
    OLB_PRECONDITION(isConstructed());
    for (int index=0; index<nx*ny; ++index) {
        for (int iDim=0; iDim<nDim; ++iDim) {
            (*this)[index][iDim] = T();
        }
    }
}

template<typename T, int nDim>
TensorField2D<T,nDim> const& TensorField3D<T,nDim>::sliceX(int xVal) const {
    xSlice.construct();
    for (int iY=0; iY<ny; ++iY) {
        for (int iZ=0; iZ<nz; ++iZ) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                xSlice.get(iY,iZ)[iDim] = this->get(xVal,iY,iZ)[iDim];
            }
        }
    }
    return xSlice;
}

template<typename T, int nDim>
TensorField2D<T,nDim> const& TensorField3D<T,nDim>::sliceY(int yVal) const {
    ySlice.construct();
    for (int iX=0; iX<nx; ++iX) {
        for (int iZ=0; iZ<nz; ++iZ) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                ySlice.get(iX,iZ)[iDim], this->get(iX,yVal,iZ)[iDim];
            }
        }
    }
    return ySlice;
}

template<typename T, int nDim>
TensorField2D<T,nDim> const& TensorField3D<T,nDim>::sliceZ(int zVal) const {
    zSlice.construct();
    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            for (int iDim=0; iDim<nDim; ++iDim) {
                zSlice.get(iX,iY)[iDim], this->get(iX,iY,zVal)[iDim];
            }
        }
    }
    return zSlice;
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::swap(TensorField3D<T,nDim>& rhs) {
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
    for (int iDim=0; iDim<nDim; ++iDim) {
        std::swap(components[iDim], rhs.components[iDim]);
    }
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::allocateMemory() {
    rawData = new Tensor   [nx*ny*nz];
    field   = new Tensor** [nx];
    for (int iX=0; iX<nx; ++iX) {
        field[iX] = new Tensor* [ny];
        for (int iY=0; iY<ny; ++iY) {
            field[iX][iY] = rawData + nz*(iY+ny*iX);
        }
    }
}

template<typename T, int nDim>
void TensorField3D<T,nDim>::releaseMemory() {
    delete [] rawData; rawData = 0;
    for (int iX=0; iX<nx; ++iX) {
      delete [] field[iX];
    }
    delete [] field;
}
template<typename T, int nDim>
ScalarField3D<T> const& TensorField3D<T,nDim>::extractComponent(int iDim)
    const
{
    components[iDim]->construct();
    for (int iEl=0; iEl<nx*ny*nz; ++iEl) {
        (*components[iDim])[iEl] = (*this)[iEl][iDim];
    };
    return *components[iDim];
}

}  // namespace olb

#endif  

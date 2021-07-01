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
 * Scalar, vector and tensor fields for 3D data analysis -- header file.
 */
#ifndef DATA_FIELDS_3D_H
#define DATA_FIELDS_3D_H

#include "olbDebug.h"
#include "dataFields2D.h"

namespace olb {

template<typename T>
class ScalarField3D {
public:
    ScalarField3D(int nx_, int ny_, int nz_);
    ~ScalarField3D();
    ScalarField3D(ScalarField3D<T> const& rhs);
    ScalarField3D<T>& operator=(ScalarField3D<T> const& rhs);
    bool isConstructed() const;
    void construct();
    void deConstruct();
    void reset();
    int getNx() const { return nx; }
    int getNy() const { return ny; }
    int getNz() const { return nz; }
    int getSize() const { return nx*ny*nz; }
    T& get(int iX, int iY, int iZ) {
        OLB_PRECONDITION(iX>=0 && iX<nx);
        OLB_PRECONDITION(iY>=0 && iY<ny);
        OLB_PRECONDITION(iZ>=0 && iZ<nz);
        OLB_PRECONDITION(isConstructed());
        return field[iX][iY][iZ];
    }
    T const& get(int iX, int iY, int iZ) const {
        OLB_PRECONDITION(iX>=0 && iX<nx);
        OLB_PRECONDITION(iY>=0 && iY<ny);
        OLB_PRECONDITION(iZ>=0 && iZ<nz);
        OLB_PRECONDITION(isConstructed());
        return field[iX][iY][iZ];
    }
    T& operator[] (int ind) {
        OLB_PRECONDITION(ind>=0 && ind<nx*ny*nz);
        OLB_PRECONDITION(isConstructed());
        return rawData[ind];
    }
    T const& operator[] (int ind) const {
        OLB_PRECONDITION(ind>=0 && ind<nx*ny*nz);
        OLB_PRECONDITION(isConstructed());
        return rawData[ind];
    }
    ScalarField2D<T>& sliceX(int xVal) const;
    ScalarField2D<T>& sliceY(int yVal) const;
    ScalarField2D<T>& sliceZ(int zVal) const;
    void swap(ScalarField3D<T>& rhs);
private:
    void allocateMemory();
    void releaseMemory();
private:
    int nx;
    int ny;
    int nz;
    T   *rawData;
    T   ***field;
    mutable ScalarField2D<T> xSlice;
    mutable ScalarField2D<T> ySlice;
    mutable ScalarField2D<T> zSlice;
};

template<typename T, int nDim>
class TensorField3D {
public:
    typedef T Tensor[nDim];
public:
    TensorField3D(int nx_, int ny_, int nz_);
    ~TensorField3D();
    TensorField3D(TensorField3D<T,nDim> const& rhs);
    TensorField3D<T,nDim>& operator=(TensorField3D<T,nDim> const& rhs);
    bool isConstructed() const;
    void construct();
    void deConstruct();
    void reset();
    int getNx() const { return nx; }
    int getNy() const { return ny; }
    int getNz() const { return nz; }
    Tensor& get(int iX, int iY, int iZ) {
        OLB_PRECONDITION(iX>=0 && iX<nx);
        OLB_PRECONDITION(iY>=0 && iY<ny);
        OLB_PRECONDITION(iZ>=0 && iZ<nz);
        OLB_PRECONDITION(isConstructed());
        return field[iX][iY][iZ];
    }
    Tensor const& get(int iX, int iY, int iZ) const {
        OLB_PRECONDITION(iX>=0 && iX<nx);
        OLB_PRECONDITION(iY>=0 && iY<ny);
        OLB_PRECONDITION(iZ>=0 && iZ<nz);
        OLB_PRECONDITION(isConstructed());
        return field[iX][iY][iZ];
    }
    Tensor& operator[] (int ind) {
        OLB_PRECONDITION(ind>=0 && ind<nx*ny*nz);
        OLB_PRECONDITION(isConstructed());
        return rawData[ind];
    }
    Tensor const& operator[] (int ind) const {
        OLB_PRECONDITION(ind>=0 && ind<nx*ny*nz);
        OLB_PRECONDITION(isConstructed());
        return rawData[ind];
    }
    TensorField2D<T,nDim> const& sliceX(int xVal) const;
    TensorField2D<T,nDim> const& sliceY(int yVal) const;
    TensorField2D<T,nDim> const& sliceZ(int zVal) const;
    ScalarField3D<T> const& extractComponent(int iDim) const;
    void swap(TensorField3D<T,nDim>& rhs);
private:
    void allocateMemory();
    void releaseMemory();
private:
    int nx;
    int ny;
    int nz;
    Tensor *rawData;
    Tensor ***field;
    mutable TensorField2D<T,nDim> xSlice;
    mutable TensorField2D<T,nDim> ySlice;
    mutable TensorField2D<T,nDim> zSlice;

    mutable ScalarField3D<T> *components[nDim];
};


}  // namespace olb


#endif

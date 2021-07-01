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
 * Scalar, vector and tensor fields for 2D data analysis -- header file.
 */
#ifndef DATA_FIELDS_2D_H
#define DATA_FIELDS_2D_H

#include "olbDebug.h"

namespace olb {

template<typename T>
class ScalarField2D {
public:
    ScalarField2D(int nx_, int ny_);
    ~ScalarField2D();
    ScalarField2D(ScalarField2D<T> const& rhs);
    ScalarField2D<T>& operator=(ScalarField2D<T> const& rhs);
    bool isConstructed() const;
    void construct();
    void deConstruct();
    void reset();
    int getNx() const { return nx; }
    int getNy() const { return ny; }
    int getSize() const { return nx*ny; }
    T& get(int iX, int iY) {
        OLB_PRECONDITION(iX>=0 && iX<nx);
        OLB_PRECONDITION(iY>=0 && iY<ny);
        OLB_PRECONDITION(isConstructed());
        return field[iX][iY];
    }
    T const& get(int iX, int iY) const {
        OLB_PRECONDITION(iX>=0 && iX<nx);
        OLB_PRECONDITION(iY>=0 && iY<ny);
        OLB_PRECONDITION(isConstructed());
        return field[iX][iY];
    }
    T& operator[] (int ind) {
        OLB_PRECONDITION(ind>=0 && ind<nx*ny);
        OLB_PRECONDITION(isConstructed());
        return rawData[ind];
    }
    T const& operator[] (int ind) const {
        OLB_PRECONDITION(ind>=0 && ind<nx*ny);
        OLB_PRECONDITION(isConstructed());
        return rawData[ind];
    }
    void swap(ScalarField2D<T>& rhs);
private:
    void allocateMemory();
    void releaseMemory();
private:
    int nx;
    int ny;
    T   *rawData;
    T   **field;
};

template<typename T, int nDim>
class TensorField2D {
public:
    typedef T Tensor[nDim];
public:
    TensorField2D(int nx_, int ny_);
    ~TensorField2D();
    TensorField2D(TensorField2D<T,nDim> const& rhs);
    TensorField2D<T,nDim>& operator=(TensorField2D<T,nDim> const& rhs);
    bool isConstructed() const;
    void construct();
    void deConstruct();
    void reset();
    int getNx() const { return nx; }
    int getNy() const { return ny; }
    Tensor& get(int iX, int iY) {
        OLB_PRECONDITION(iX>=0 && iX<nx);
        OLB_PRECONDITION(iY>=0 && iY<ny);
        OLB_PRECONDITION(isConstructed());
        return field[iX][iY];
    }
    Tensor const& get(int iX, int iY) const {
        OLB_PRECONDITION(iX>=0 && iX<nx);
        OLB_PRECONDITION(iY>=0 && iY<ny);
        OLB_PRECONDITION(isConstructed());
        return field[iX][iY];
    }
    Tensor& operator[] (int ind) {
        OLB_PRECONDITION(ind>=0 && ind<nx*ny);
        OLB_PRECONDITION(isConstructed());
        return rawData[ind];
    }
    Tensor const& operator[] (int ind) const {
        OLB_PRECONDITION(ind>=0 && ind<nx*ny);
        OLB_PRECONDITION(isConstructed());
        return rawData[ind];
    }
    void swap(TensorField2D<T,nDim>& rhs);
    ScalarField2D<T> const& extractComponent(int whichDim) const;
private:
    void allocateMemory();
    void releaseMemory();
private:
    int nx;
    int ny;
    Tensor *rawData;
    Tensor **field;
    mutable ScalarField2D<T> *components[nDim];
};

template<typename T, template<typename U> class ScalarField>
T computeMin(ScalarField<T> const& field);

template<typename T, template<typename U> class ScalarField>
T computeMax(ScalarField<T> const& field);

template<typename T, template<typename U> class ScalarField>
T computeAverage(ScalarField<T> const& field);

template<typename T, template<typename U> class ScalarField>
T computeRMS(ScalarField<T> const& field);

template<typename T, template<typename U> class ScalarField>
T computeNormSqr(ScalarField<T> const& field);


}


#endif

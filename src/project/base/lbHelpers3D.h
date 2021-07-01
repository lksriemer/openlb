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
 * Template specializations for some computationally intensive LB 
 * functions of the header file lbHelpers.h, for some D3Q19 grids.
 */

#ifndef LB_HELPERS_3D_H
#define LB_HELPERS_3D_H

namespace olb {

// Efficient specialization for D3Q19 lattice
template<typename T>
struct lbHelpers<T, descriptors::D3Q19Descriptor> {

    static T equilibrium( int iPop, T rho, const T u[3], const T uSqr ) {
        typedef descriptors::D3Q19Descriptor<T> L;
        T c_u = L::c[iPop][0]*u[0] + L::c[iPop][1]*u[1] + L::c[iPop][2]*u[2];
        return rho * L::t[iPop] * ( 1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr )
               - L::t[iPop];
    }

    static T incEquilibrium( int iPop, const T j[3],
                             const T jSqr, const T pressure )
    {
        typedef descriptors::D3Q19Descriptor<T> L;
        T c_j = L::c[iPop][0]*j[0] + L::c[iPop][1]*j[1] + L::c[iPop][2]*j[2];
        return L::t[iPop] * ( 3.*pressure + 3.*c_j + 4.5*c_j*c_j - 1.5*jSqr )
               - L::t[iPop];
    }

    static void computeFneq (
            Cell<T,descriptors::D3Q19Descriptor> const& cell,
            T fNeq[3], T rho, const T u[3] )
    {
        const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
        for (int iPop=0; iPop < 19; ++iPop) {
            fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
        }
    }

    static T bgkCollision (
            Cell<T,descriptors::D3Q19Descriptor>& cell,
            T rho, const T u[3], T omega)
    {
        const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
        for (int iPop=0; iPop < 19; ++iPop) {
            cell[iPop] *= (T)1-omega;
            cell[iPop] += omega *
                lbHelpers<T,descriptors::D3Q19Descriptor>::equilibrium (
                              iPop, rho, u, uSqr
                 );
        }
        return uSqr;
    }

    static T incBgkCollision (
            Cell<T,descriptors::D3Q19Descriptor>& cell,
            T pressure, const T j[2], T omega)
    {
        const T jSqr = util::normSqr<T,descriptors::D3Q19Descriptor<T>::d>(j);
        for (int iPop=0; iPop < descriptors::D3Q19Descriptor<T>::q; ++iPop) {
            cell[iPop] *= (T)1-omega;
            cell[iPop] += omega * lbHelpers<T,descriptors::D3Q19Descriptor>::incEquilibrium (
                              iPop, j, jSqr, pressure );
        }
        return jSqr;
    }

    static T constRhoBgkCollision (
            Cell<T,descriptors::D3Q19Descriptor>& cell,
            T const& rho, T u[3], T ratioRho, T omega)
    {
        const T uSqr = util::normSqr<T,descriptors::D3Q19Descriptor<T>::d>(u);
        for (int iPop=0; iPop < descriptors::D3Q19Descriptor<T>::q; ++iPop) {
            T feq = lbHelpers<T,descriptors::D3Q19Descriptor>::
                         equilibrium(iPop, rho, u, uSqr );
            cell[iPop] =
              ratioRho*(feq+descriptors::D3Q19Descriptor<T>::t[iPop])
              -descriptors::D3Q19Descriptor<T>::t[iPop] +
                  ((T)1-omega)*(cell[iPop]-feq);
        }
        return uSqr;
    }

    static void partial_rho (
        Cell<T,descriptors::D3Q19Descriptor> const& cell,
        T& surfX_M1, T& surfX_0, T& surfX_P1,
        T& surfY_M1, T& surfY_P1, T& surfZ_M1, T& surfZ_P1 )
    {
        surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
        surfX_0  = cell[0] + cell[2] + cell[3] + cell[8] +
                   cell[9] + cell[11] + cell[12] + cell[17] + cell[18];
        surfX_P1 = cell[10] + cell[13] + cell[14] + cell[15] + cell[16];

        surfY_M1 = cell[2] + cell[4] + cell[8] + cell[9] + cell[14];
        surfY_P1 = cell[5] + cell[11] + cell[13] + cell[17] + cell[18];

        surfZ_M1 = cell[3] + cell[6] + cell[8] + cell[16] + cell[18];
        surfZ_P1 = cell[7] + cell[9] + cell[12] + cell[15] + cell[17];

    }

    static void computeRhoU (
        Cell<T,descriptors::D3Q19Descriptor> const& cell, T& rho, T u[3] )
    {
        T surfX_M1, surfX_0, surfX_P1,
          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;

        u[0]  = ( surfX_P1 - surfX_M1 ) / rho;
        u[1]  = ( surfY_P1 - surfY_M1 ) / rho;
        u[2]  = ( surfZ_P1 - surfZ_M1 ) / rho;
    }

    static void computeJ (
        Cell<T,descriptors::D3Q19Descriptor> const& cell, T j[3] )
    {
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
        surfX_P1 = cell[10] + cell[13] + cell[14] + cell[15] + cell[16];

        surfY_M1 = cell[2] + cell[4] + cell[8] + cell[9] + cell[14];
        surfY_P1 = cell[5] + cell[11] + cell[13] + cell[17] + cell[18];

        surfZ_M1 = cell[3] + cell[6] + cell[8] + cell[16] + cell[18];
        surfZ_P1 = cell[7] + cell[9] + cell[12] + cell[15] + cell[17];



        j[0]  = ( surfX_P1 - surfX_M1 );
        j[1]  = ( surfY_P1 - surfY_M1 );
        j[2]  = ( surfZ_P1 - surfZ_M1 );
    }

    static void computeStress (
        Cell<T,descriptors::D3Q19Descriptor> const& cell, T rho,
        const T u[3], T pi[6] )
    {
        typedef descriptors::D3Q19Descriptor<T> L;
        // Workaround for Intel(r) compiler 9.1;
        // "using namespace util::tensorIndices3D" is not sufficient
        using util::tensorIndices3D::xx;
        using util::tensorIndices3D::yy;
        using util::tensorIndices3D::zz;
        using util::tensorIndices3D::xy;
        using util::tensorIndices3D::xz;
        using util::tensorIndices3D::yz;

        T surfX_M1, surfX_0, surfX_P1,
          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        pi[xx] = surfX_P1+surfX_M1 - 1./L::invCs2*(rho-(T)1) - rho*u[0]*u[0];
        pi[yy] = surfY_P1+surfY_M1 - 1./L::invCs2*(rho-(T)1) - rho*u[1]*u[1];
        pi[zz] = surfZ_P1+surfZ_M1 - 1./L::invCs2*(rho-(T)1) - rho*u[2]*u[2];

        pi[xy] = cell[4] - cell[5] + cell[13] - cell[14] - rho*u[0]*u[1];
        pi[xz] = cell[6] - cell[7] + cell[15] - cell[16] - rho*u[0]*u[2];
        pi[yz] = cell[8] - cell[9] + cell[17] - cell[18] - rho*u[1]*u[2];
    }

    static void computeAllMomenta (
        Cell<T,descriptors::D3Q19Descriptor> const& cell,
        T& rho, T u[3], T pi[6] )
    {
        typedef descriptors::D3Q19Descriptor<T> L;
        // Workaround for Intel(r) compiler 9.1;
        // "using namespace util::tensorIndices3D" is not sufficient
        using util::tensorIndices3D::xx;
        using util::tensorIndices3D::yy;
        using util::tensorIndices3D::zz;
        using util::tensorIndices3D::xy;
        using util::tensorIndices3D::xz;
        using util::tensorIndices3D::yz;

        T surfX_M1, surfX_0, surfX_P1,
          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;

        T rhoU0  = ( surfX_P1 - surfX_M1 ) / rho;
        T rhoU1  = ( surfY_P1 - surfY_M1 ) / rho;
        T rhoU2  = ( surfZ_P1 - surfZ_M1 ) / rho;
        u[0] = rhoU0 / rho;
        u[1] = rhoU1 / rho;
        u[2] = rhoU2 / rho;

        pi[xx] = surfX_P1+surfX_M1 - 1./L::invCs2*(rho-(T)1) - rhoU0*u[0];
        pi[yy] = surfY_P1+surfY_M1 - 1./L::invCs2*(rho-(T)1) - rhoU1*u[1];
        pi[zz] = surfZ_P1+surfZ_M1 - 1./L::invCs2*(rho-(T)1) - rhoU2*u[2];

        pi[xy] = cell[4] - cell[5] + cell[13] - cell[14] - rhoU0*u[1];
        pi[xz] = cell[6] - cell[7] + cell[15] - cell[16] - rhoU0*u[2];
        pi[yz] = cell[8] - cell[9] + cell[17] - cell[18] - rhoU1*u[2];
    }

    static void swapAndStreamCell (
          Cell<T,descriptors::D3Q19Descriptor> ***grid,
          int iX, int iY, int iZ, int nX, int nY, int nZ, int iPop, T& fTmp )
    {
        fTmp                     = grid[iX][iY][iZ][iPop];
        grid[iX][iY][iZ][iPop]   = grid[iX][iY][iZ][iPop+9];
        grid[iX][iY][iZ][iPop+9] = grid[nX][nY][nZ][iPop];
        grid[nX][nY][nZ][iPop]   = fTmp;
    }

    static void swapAndStream3D(Cell<T,descriptors::D3Q19Descriptor> ***grid,
                                int iX, int iY, int iZ)
    {
        T fTmp;
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ,   4, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ,   5, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY  , iZ-1, 6, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY  , iZ+1, 7, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX  , iY-1, iZ-1, 8, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX  , iY-1, iZ+1, 9, fTmp);
    }

    static T computeRho (
        Cell<T,descriptors::D3Q19Descriptor> const& cell)
    {
        T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4]
                        + cell[5] + cell[6] + cell[7] + cell[8]
                        + cell[9] + cell[10] + cell[11] + cell[12]
                        + cell[13] + cell[14] + cell[15] + cell[16]
                        + cell[17] + cell[18] + (T)1;
        return rho;
    }

    static void modifyVelocity (
            Cell<T,descriptors::D3Q19Descriptor>& cell, const T newU[3] )
    {
        T rho, oldU[3];
        computeRhoU(cell, rho, oldU);
        const T oldUSqr = util::normSqr<T,3>(oldU);
        const T newUSqr = util::normSqr<T,3>(newU);
        for (int iPop=0; iPop<19; ++iPop) {
            cell[iPop] = cell[iPop]
                             - equilibrium(iPop, rho, oldU, oldUSqr)
                             + equilibrium(iPop, rho, newU, newUSqr);
        }
    }

};  //struct lbHelpers<D3Q19Descriptor>


// Efficient specialization for D3Q15 lattice
template<typename T>
struct lbHelpers<T, descriptors::D3Q15Descriptor> {

    static T equilibrium( int iPop, T rho, const T u[3], const T uSqr ) {
        T c_u = descriptors::D3Q15Descriptor<T>::c[iPop][0]*u[0] +
                descriptors::D3Q15Descriptor<T>::c[iPop][1]*u[1] +
                descriptors::D3Q15Descriptor<T>::c[iPop][2]*u[2];
        return rho * descriptors::D3Q15Descriptor<T>::t[iPop] * (
                   1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr)
               - descriptors::D3Q15Descriptor<T>::t[iPop];
    }

    static T incEquilibrium( int iPop, const T j[3],
                             const T jSqr, const T pressure )
    {
        typedef descriptors::D3Q15Descriptor<T> L;
        T c_j = L::c[iPop][0]*j[0] + L::c[iPop][1]*j[1] + L::c[iPop][2]*j[2];
        return L::t[iPop] * ( 3.*pressure + 3.*c_j + 4.5*c_j*c_j - 1.5*jSqr )
               - L::t[iPop];
    }

    static void computeFneq (
            Cell<T,descriptors::D3Q15Descriptor> const& cell,
            T fNeq[3], T rho, const T u[3] )
    {
        const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
        for (int iPop=0; iPop < 15; ++iPop) {
            fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
        }
    }

    static T bgkCollision (
            Cell<T,descriptors::D3Q15Descriptor>& cell,
            T rho, const T u[3], T omega)
    {
        const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
        for (int iPop=0; iPop < 15; ++iPop) {
            cell[iPop] *= (T)1-omega;
            cell[iPop] += omega *
                lbHelpers<T,descriptors::D3Q15Descriptor>::equilibrium (
                              iPop, rho, u, uSqr
                 );
        }
        return uSqr;
    }

    static T incBgkCollision (
            Cell<T,descriptors::D3Q15Descriptor>& cell,
            T pressure, const T j[2], T omega)
    {
        const T jSqr = util::normSqr<T,descriptors::D3Q15Descriptor<T>::d>(j);
        for (int iPop=0; iPop < descriptors::D3Q19Descriptor<T>::q; ++iPop) {
            cell[iPop] *= (T)1-omega;
            cell[iPop] += omega * lbHelpers<T,descriptors::D3Q15Descriptor>::incEquilibrium (
                              iPop, j, jSqr, pressure );
        }
        return jSqr;
    }

    static T constRhoBgkCollision (
            Cell<T,descriptors::D3Q15Descriptor>& cell,
            T const& rho, T u[3], T ratioRho, T omega)
    {
        const T uSqr = util::normSqr<T,descriptors::D3Q15Descriptor<T>::d>(u);
        for (int iPop=0; iPop < descriptors::D3Q15Descriptor<T>::q; ++iPop) {
            T feq = lbHelpers<T,descriptors::D3Q15Descriptor>::
                         equilibrium(iPop, rho, u, uSqr );
            cell[iPop] =
              ratioRho*(feq+descriptors::D3Q15Descriptor<T>::t[iPop])
              -descriptors::D3Q15Descriptor<T>::t[iPop] +
                  ((T)1-omega)*(cell[iPop]-feq);
        }
        return uSqr;
    }

    static void partial_rho (
        Cell<T,descriptors::D3Q15Descriptor> const& cell,
        T& surfX_M1, T& surfX_0, T& surfX_P1,
        T& surfY_M1, T& surfY_P1, T& surfZ_M1, T& surfZ_P1 )
    {
        surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
        surfX_0  = cell[0] + cell[2] + cell[3] + cell[9] + cell[10];
        surfX_P1 = cell[8] + cell[11] + cell[12] + cell[13] + cell[14];

        surfY_M1 = cell[2] + cell[4] + cell[5] + cell[13] + cell[14];
        surfY_P1 = cell[6] + cell[7] + cell[9] + cell[11] + cell[12];

        surfZ_M1 = cell[3] + cell[4] + cell[6] + cell[12] + cell[14];
        surfZ_P1 = cell[5] + cell[7] + cell[10] + cell[11] + cell[13];
    }

    static T computeRho (
        Cell<T,descriptors::D3Q15Descriptor> const& cell)
    {
        T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4]
                        + cell[5] + cell[6] + cell[7] + cell[8]
                        + cell[9] + cell[10] + cell[11] + cell[12]
                        + cell[13] + cell[14] + (T)1;
        return rho;
    }

    static void computeRhoU (
        Cell<T,descriptors::D3Q15Descriptor> const& cell,
        T& rho, T u[3] )
    {
        T surfX_M1, surfX_0, surfX_P1,
          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;

        u[0]  = ( surfX_P1 - surfX_M1 ) / rho;
        u[1]  = ( surfY_P1 - surfY_M1 ) / rho;
        u[2]  = ( surfZ_P1 - surfZ_M1 ) / rho;
    }

    static void computeJ (
        Cell<T,descriptors::D3Q15Descriptor> const& cell, T j[3] )
    {
        T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
        surfX_P1 = cell[8] + cell[11] + cell[12] + cell[13] + cell[14];

        surfY_M1 = cell[2] + cell[4] + cell[5] + cell[13] + cell[14];
        surfY_P1 = cell[6] + cell[7] + cell[9] + cell[11] + cell[12];

        surfZ_M1 = cell[3] + cell[4] + cell[6] + cell[12] + cell[14];
        surfZ_P1 = cell[5] + cell[7] + cell[10] + cell[11] + cell[13];

        j[0]  = ( surfX_P1 - surfX_M1 );
        j[1]  = ( surfY_P1 - surfY_M1 );
        j[2]  = ( surfZ_P1 - surfZ_M1 );
    }

    static void computeStress (
        Cell<T,descriptors::D3Q15Descriptor> const& cell, T rho,
        const T u[3], T pi[6] )
    {
        typedef descriptors::D3Q15Descriptor<T> L;
        using namespace util::tensorIndices3D;

        T surfX_M1, surfX_0, surfX_P1,
          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        pi[xx] = surfX_P1+surfX_M1 - 1./L::invCs2*(rho-(T)1) - rho*u[0]*u[0];
        pi[yy] = surfY_P1+surfY_M1 - 1./L::invCs2*(rho-(T)1) - rho*u[1]*u[1];
        pi[zz] = surfZ_P1+surfZ_M1 - 1./L::invCs2*(rho-(T)1) - rho*u[2]*u[2];

        pi[xy] =   cell[4] + cell[5] - cell[6]  - cell[7]
                 + cell[11] + cell[12] - cell[13] - cell[14] - rho*u[0]*u[1];
        pi[xz] =   cell[4] - cell[5] + cell[6]  - cell[7]
                 + cell[11] - cell[12] + cell[13] - cell[14] - rho*u[0]*u[2];
        pi[yz] =   cell[4] - cell[5] - cell[6]  + cell[7]
                 + cell[11] - cell[12] - cell[13] + cell[14] - rho*u[1]*u[2];

    }

    static void computeAllMomenta (
        Cell<T,descriptors::D3Q15Descriptor> const& cell,
        T& rho, T u[3], T pi[6] )
    {
        typedef descriptors::D3Q15Descriptor<T> L;
        using namespace util::tensorIndices3D;

        T surfX_M1, surfX_0, surfX_P1,
          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

        partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                          surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

        rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;

        T rhoU0  = ( surfX_P1 - surfX_M1 ) / rho;
        T rhoU1  = ( surfY_P1 - surfY_M1 ) / rho;
        T rhoU2  = ( surfZ_P1 - surfZ_M1 ) / rho;
        u[0] = rhoU0 / rho;
        u[1] = rhoU1 / rho;
        u[2] = rhoU2 / rho;

        pi[xx] = surfX_P1+surfX_M1 - 1./L::invCs2*(rho-(T)1) - rhoU0*u[0];
        pi[yy] = surfY_P1+surfY_M1 - 1./L::invCs2*(rho-(T)1) - rhoU1*u[1];
        pi[zz] = surfZ_P1+surfZ_M1 - 1./L::invCs2*(rho-(T)1) - rhoU2*u[2];

        pi[xy] =   cell[4] + cell[5] - cell[6]  - cell[7]
                 + cell[11] + cell[12] - cell[13] - cell[14] - rhoU0*u[1];
        pi[xz] =   cell[4] - cell[5] + cell[6]  - cell[7]
                 + cell[11] - cell[12] + cell[13] - cell[14] - rhoU0*u[2];
        pi[yz] =   cell[4] - cell[5] - cell[6]  + cell[7]
                 + cell[11] - cell[12] - cell[13] + cell[14] - rhoU1*u[2];
    }

    static void swapAndStreamCell (
          Cell<T,descriptors::D3Q15Descriptor> ***grid,
          int iX, int iY, int iZ, int nX, int nY, int nZ, int iPop, T& fTmp )
    {
        fTmp                     = grid[iX][iY][iZ][iPop];
        grid[iX][iY][iZ][iPop]   = grid[iX][iY][iZ][iPop+7];
        grid[iX][iY][iZ][iPop+7] = grid[nX][nY][nZ][iPop];
        grid[nX][nY][nZ][iPop]   = fTmp;
    }

    static void swapAndStream3D(Cell<T,descriptors::D3Q15Descriptor> ***grid,
                                int iX, int iY, int iZ)
    {
        T fTmp;
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);

        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ-1, 4, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ+1, 5, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ-1, 6, fTmp);
        swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ+1, 7, fTmp);
    }

    static void modifyVelocity (
            Cell<T,descriptors::D3Q15Descriptor>& cell, const T newU[3] )
    {
        T rho, oldU[3];
        computeRhoU(cell, rho, oldU);
        const T oldUSqr = util::normSqr<T,3>(oldU);
        const T newUSqr = util::normSqr<T,3>(newU);
        for (int iPop=0; iPop<15; ++iPop) {
            cell[iPop] = cell[iPop]
                             - equilibrium(iPop, rho, oldU, oldUSqr)
                             + equilibrium(iPop, rho, newU, newUSqr);
        }
    }

};  //struct lbHelpers<D3Q15Descriptor>


}  // namespace olb

#endif

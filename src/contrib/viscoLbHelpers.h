/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Orestis Malaspinas, Jonas Latt
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef VISCO_LB_HELPERS_H
#define VISCO_LB_HELPERS_H

#include <cmath>

namespace olb {

template<typename T, template<typename U> class Lattice>
struct viscoSolventLbHelpers 
{
    /// Computation of equilibrium distribution
    static T equilibrium( int iPop, const T u[Lattice<T>::d], T uSqr, 
                          const T A[Lattice<T>::q], const T B[Lattice<T>::q], 
                          const T C[Lattice<T>::q], const T D[Lattice<T>::q],
                          const T E[Lattice<T>::q][Lattice<T>::d][Lattice<T>::d])
    {
        typedef Lattice<T> L;
        
        T c_u = T();
        T cc_uu = T();
        T cc_E;
        for (int iAlpha = 0; iAlpha < L::d; ++iAlpha)
        {
            c_u += L::c[iPop][iAlpha]*u[iAlpha];
            for (int iBeta = 0; iBeta < L::d; ++iBeta)
            {
                cc_uu += u[iAlpha]*u[iBeta]*L::c[iPop][iAlpha]*L::c[iPop][iBeta];
                cc_E += E[iPop][iAlpha][iBeta]*L::c[iPop][iAlpha]*L::c[iPop][iBeta];
            }
        }
        
        return A[iPop] + B[iPop] * c_u + C[iPop] * uSqr + D[iPop] * cc_uu + cc_E;
    }
    
    static void computeConstants(T rho, const T sigma[], 
                                 std::vector<T> &normIndex, std::vector<T> &diagIndex,
                                 const T A[Lattice<T>::q], const T B[Lattice<T>::q], 
                                 const T C[Lattice<T>::q], const T D[Lattice<T>::q],
                                 const T E[Lattice<T>::q][Lattice<T>::d][Lattice<T>::d])
    {
        using namespace util::tensorIndices2D;
        
        const T As = (sigma[xx]+sigma[yy])/(T)16;
        const T Bs = rho / (T)12;
        const T Cs = -rho/(T)16;
        const T Ds = rho/(T)8;
        const T Exxs = (sigma[xx]-sigma[yy])/(T)16;
        const T Exys = sigma[xy]/(T)8;
        
        A[0] = rho - (T)2*As;
        B[0] = T();
        C[0] = (T)12*Cs;
        D[0] = T();
        for (int iAlpha = 0; iAlpha < L::d; ++iAlpha)
        {
            for (int iBeta = 0; iBeta < L::d; ++iBeta)
            {
                E[0][iAlpha][iBeta] = T();
            }
        }
        
        for (unsigned iPop = 0; iPop < normIndex.size(); ++iPop)
        {
            A[normIndex[iPop]] = (T)2*As;
            B[normIndex[iPop]] = (T)4*Bs;
            C[normIndex[iPop]] = (T)2*Cs;
            D[normIndex[iPop]] = (T)4*Ds;
            E[normIndex[iPop]][0][0] = (T)4*Exxs;
            E[normIndex[iPop]][1][1] = -(T)4*Exxs;
            E[normIndex[iPop]][0][1] = T();
            E[normIndex[iPop]][1][0] = T();
        }
        
        for (unsigned iPop = 0; iPop < diagIndex.size(); ++iPop)
        {
            A[diagIndex[iPop]] = As;
            B[diagIndex[iPop]] = Bs;
            C[diagIndex[iPop]] = Cs;
            D[diagIndex[iPop]] = Ds;
            E[diagIndex[iPop]][0][0] = Exxs;
            E[diagIndex[iPop]][1][1] = -Exxs;
            E[diagIndex[iPop]][0][1] = Exys;
            E[diagIndex[iPop]][1][0] = Exys;
        }
    }
};

}

#endif

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
 * Helper functions for the implementation of MRT (multiple relaxation
 * time) dynamics. This file is all about efficiency.
 * The generic template code is specialized for commonly used Lattices,
 * so that a maximum performance can be taken out of each case.
 */
#ifndef MRT_HELPERS_H
#define MRT_HELPERS_H

#include "latticeDescriptors.h"
#include "cell.h"
#include "util.h"


namespace olb {

/// Helper functions for the (somewhat special) D3Q13 lattice
template<typename T>
struct d3q13Helpers {
    typedef descriptors::D3Q13Descriptor<T> Lattice;

    /// BGK collision step
    static T collision( Cell<T,descriptors::D3Q13Descriptor>& cell,
                        T rho, const T u[Lattice::d],
                        T lambda_nu, T lambda_nu_prime)
    {
        const T lambda_e = descriptors::D3Q13Descriptor<T>::lambda_e;
        const T lambda_h = descriptors::D3Q13Descriptor<T>::lambda_h;

        T uxSqr = util::sqr(u[0]);
        T uySqr = util::sqr(u[1]);
        T uzSqr = util::sqr(u[2]);
        T uSqr = uxSqr + uySqr + uzSqr;

        T s1 = cell[7] + cell[8] + cell[9] + cell[10];
        T s2 = cell[11] + cell[12];
        T s3 = cell[1] + cell[2] + cell[3] + cell[4];
        T s4 = cell[5] + cell[6];
        T sTot = s1 + s2 + s3 + s4;
        T d1 = cell[7] + cell[8] - cell[9] - cell[10];
        T d2 = cell[1] + cell[2] - cell[3] - cell[4];

        T M[13]; // The terms M[0]-M[3] are not used (they correspond
                 // to rho, rho*ux, rho*uy, rho*uz respectively),
                 // but we still use a 13-component vector to preserve
                 // the clarity of the code.
        M[4] = -(T)12*cell[0] + sTot - (T)11/(T)2;
          // The 11/2 correction term in M4 accounts for the fact that
          // cell[i] corresponds to f[i]-ti, and not f[i]
        M[5] = s1 - (T)2*s2 + s3 - (T)2*s4;
        M[6] = d1 + d2;
        M[7] = cell[7] - cell[8] + cell[1] - cell[2];
        M[8] = cell[11] - cell[12] + cell[5] - cell[6];
        M[9] = cell[9] - cell[10] + cell[3] - cell[4];
        M[10] = d1 - d2;
        M[11] = -cell[7] + cell[8] + s2 + cell[1] - cell[2] - s4;
        M[12] = cell[9] - cell[10] - cell[11] + cell[12]
                - cell[3] + cell[4] + cell[5] - cell[6];
        T Mneq[13]; // The terms Mneq[0]-Mneq[3] are not used (they are
                    // actually all zero, because of conservation laws),
                    // but we still use a 13-component vector to preserve
                    // the clarity of the code.
        Mneq[4] = M[4] + (T)11/(T)2*rho - (T)13/(T)2*rho*uSqr;
        Mneq[5] = M[5] - rho*( (T)2*uxSqr-(uySqr+uzSqr) );
        Mneq[6] = M[6] - rho*( uySqr-uzSqr );
        Mneq[7] = M[7] - rho*( u[0]*u[1] );
        Mneq[8] = M[8] - rho*( u[1]*u[2] );
        Mneq[9] = M[9] - rho*( u[0]*u[2] );
        Mneq[10] = M[10];
        Mneq[11] = M[11];
        Mneq[12] = M[12];

        Mneq[4]  *= lambda_e  / (T)156;
        Mneq[5]  *= lambda_nu / (T)24;
        Mneq[6]  *= lambda_nu / (T)8;
        Mneq[7]  *= lambda_nu_prime / (T)4;
        Mneq[8]  *= lambda_nu_prime / (T)4;
        Mneq[9]  *= lambda_nu_prime / (T)4;
        Mneq[10] *= lambda_h / (T)8;
        Mneq[11] *= lambda_h / (T)8;
        Mneq[12] *= lambda_h / (T)8;

        T F1 = Mneq[4] + Mneq[5] + Mneq[6];
        T F2 = Mneq[4] + Mneq[5] - Mneq[6];
        T F3 = Mneq[4] - (T)2*Mneq[5];

        cell[0]  -= (T)-12*Mneq[4];
        cell[1]  -= F1 + Mneq[7]                -Mneq[10]+Mneq[11];
        cell[2]  -= F1 - Mneq[7]                -Mneq[10]-Mneq[11];
        cell[3]  -= F2                  +Mneq[9]+Mneq[10]         -Mneq[12];
        cell[4]  -= F2                  -Mneq[9]+Mneq[10]         +Mneq[12];
        cell[5]  -= F3          +Mneq[8]                 -Mneq[11]+Mneq[12];
        cell[6]  -= F3          -Mneq[8]                 -Mneq[11]-Mneq[12];
        cell[7]  -= F1 + Mneq[7]                +Mneq[10]-Mneq[11];
        cell[8]  -= F1 - Mneq[7]                +Mneq[10]+Mneq[11];
        cell[9]  -= F2                  +Mneq[9]-Mneq[10]         +Mneq[12];
        cell[10] -= F2                  -Mneq[9]-Mneq[10]         -Mneq[12];
        cell[11] -= F3          +Mneq[8]                 +Mneq[11]-Mneq[12];
        cell[12] -= F3          -Mneq[8]                 +Mneq[11]+Mneq[12];

        return uSqr;
    }


    /*
    static T collision( Cell<T,descriptors::D3Q13Descriptor>& cell,
                        T rho, const T u[Lattice::d],
                        T lambda_nu, T lambda_nu_prime)
    {
        for  (int i=0; i<13; ++i) {
            cell[i] += descriptors::D3Q13Descriptor<T>::t[i];
        }
        static const int e_[13][13] = {
            { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
            { 0, 1, 1, 1, 1, 0, 0,-1,-1,-1,-1, 0, 0 },
            { 0, 1,-1, 0, 0, 1, 1,-1, 1, 0, 0,-1,-1 },
            { 0, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1 },
           {-12, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
            { 0, 1, 1, 1, 1,-2,-2, 1, 1, 1, 1,-2,-2 },
            { 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0 },
            { 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1 },
            { 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0 },
            { 0, 1, 1,-1,-1, 0, 0,-1,-1, 1, 1, 0, 0 },
            { 0,-1, 1, 0, 0, 1, 1, 1,-1, 0, 0,-1,-1 },
            { 0, 0, 0, 1,-1,-1, 1, 0, 0,-1, 1, 1,-1 } };

        int eNormSqr[13];
        for (int k=0; k<13; ++k) {
            eNormSqr[k] = 0;
            for (int i=0; i<13; ++i) {
                eNormSqr[k] += e_[k][i]*e_[k][i];
            }
        }

        T M_[13];
        for (int k=0; k<13; ++k) {
            M_[k] = 0.;
            for (int i=0; i<13; ++i) {
                M_[k] += cell[i] * (T)e_[k][i];
            }
        }

        T rho_ = M_[0];
        T jx_  = M_[1];
        T jy_  = M_[2];
        T jz_  = M_[3];
        T jxSqr_ = jx_*jx_;
        T jySqr_ = jy_*jy_;
        T jzSqr_ = jz_*jz_;
        T jSqr_  = jxSqr_+jySqr_+jzSqr_;

        static const T lambda_e = 1.5;
        static const T lambda_h = 1.8;
        static const T cs_sqr   = 1./3.;

        static const T lambda[13] = {
            0., 0., 0., 0.,
            lambda_e, lambda_nu, lambda_nu,
            lambda_nu_prime, lambda_nu_prime, lambda_nu_prime,
            lambda_h, lambda_h, lambda_h };

        T Meq_[13];

        Meq_[0] = rho_;
        Meq_[1] = jx_;
        Meq_[2] = jy_;
        Meq_[3] = jz_;
        Meq_[4] = 3./2.*(13.*cs_sqr-8.)*rho_+13./2./rho_*jSqr_;
        Meq_[5] = (2.*jxSqr_-(jySqr_+jzSqr_))/rho_;
        Meq_[6] = (jySqr_-jzSqr_)/rho_;
        Meq_[7] = jx_*jy_/rho_;
        Meq_[8] = jy_*jz_/rho_;
        Meq_[9] = jx_*jz_/rho_;
        Meq_[10] = 0.;
        Meq_[11] = 0.;
        Meq_[12] = 0.;

        T Nout[13];
        for (int i=0; i<13; ++i) {
            Nout[i] = cell[i];
        }

        for (int i=0; i<13; ++i) {
            for (int k=4; k<13; ++k) {
                Nout[i] -= 
                    lambda[k]/(T)eNormSqr[k]*(cell[k]-Meq_[k])*(T)e_[k][i];
            }
        }

        for (int i=0; i<13; ++i) {
            cell[i] = Nout[i];
        }

        for  (int i=0; i<13; ++i) {
            cell[i] -= descriptors::D3Q13Descriptor<T>::t[i];
        }

        return jSqr_/rho_/rho_;
    }
*/

    /*
    static T collision( Cell<T,descriptors::D3Q13Descriptor>& cell,
                        T rho, const T u[Lattice::d],
                        T lambda_nu, T lambda_nu_prime)
    {
        for  (int i=0; i<13; ++i) {
            cell[i] += descriptors::D3Q13Descriptor<T>::t[i];
        }
        static const int e_[13][13] = {
            { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
            { 0, 1, 1, 1, 1, 0, 0,-1,-1,-1,-1, 0, 0 },
            { 0, 1,-1, 0, 0, 1, 1,-1, 1, 0, 0,-1,-1 },
            { 0, 0, 0, 1,-1, 1,-1, 0, 0,-1, 1,-1, 1 },
           {-12, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
            { 0, 1, 1, 1, 1,-2,-2, 1, 1, 1, 1,-2,-2 },
            { 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0 },
            { 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1 },
            { 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 0, 0 },
            { 0, 1, 1,-1,-1, 0, 0,-1,-1, 1, 1, 0, 0 },
            { 0,-1, 1, 0, 0, 1, 1, 1,-1, 0, 0,-1,-1 },
            { 0, 0, 0, 1,-1,-1, 1, 0, 0,-1, 1, 1,-1 } };

        int eNormSqr[13];
        for (int k=0; k<13; ++k) {
            eNormSqr[k] = 0;
            for (int i=0; i<13; ++i) {
                eNormSqr[k] += e_[k][i]*e_[k][i];
            }
        }

        T M_[13];
        for (int k=0; k<13; ++k) {
            M_[k] = 0.;
            for (int i=0; i<13; ++i) {
                M_[k] += cell[i] * (T)e_[k][i];
            }
        }

        T rho_ = M_[0];
        T jx_  = M_[1];
        T jy_  = M_[2];
        T jz_  = M_[3];
        T jxSqr_ = jx_*jx_;
        T jySqr_ = jy_*jy_;
        T jzSqr_ = jz_*jz_;
        T jSqr_  = jxSqr_+jySqr_+jzSqr_;

        T uxSqr = util::sqr(u[0]);
        T uySqr = util::sqr(u[1]);
        T uzSqr = util::sqr(u[2]);
        T uSqr = uxSqr + uySqr + uzSqr;

        assert( abs(uxSqr-jxSqr_/rho_/rho_)<1e-10);

        if( abs(uySqr-jySqr_/rho_/rho_)>=1e-10) {
            std::cout << uySqr << std::endl;
            std::cout << jySqr_ << std::endl;
            std::cout << rho_ << std::endl;
            std::cout << jySqr_/rho_/rho_ << std::endl;
        }
        assert( abs(uySqr-jySqr_/rho_/rho_)<1e-10);
        assert( abs(uzSqr-jzSqr_/rho_/rho_)<1e-10);

        if ((abs(rho-rho_) >=1e-10)) {
            std::cout << rho << std::endl;
            std::cout << rho_ << std::endl;
        }
        assert( abs(rho-rho_)<1e-10);

        static const T lambda_e = 1.5;
        static const T lambda_h = 1.8;
        static const T cs_sqr   = 1./3.;

        const T lambda_e_ = descriptors::D3Q13Descriptor<T>::lambda_e;
        const T lambda_h_ = descriptors::D3Q13Descriptor<T>::lambda_h;

        assert( abs(lambda_e - lambda_e_) < 1e-10 );
        assert( abs(lambda_h - lambda_h_) < 1e-10);


        static const T lambda[13] = {
            0., 0., 0., 0.,
            lambda_e, lambda_nu, lambda_nu,
            lambda_nu_prime, lambda_nu_prime, lambda_nu_prime,
            lambda_h, lambda_h, lambda_h };

        T Meq_[13];

        Meq_[0] = rho_;
        Meq_[1] = jx_;
        Meq_[2] = jy_;
        Meq_[3] = jz_;
        Meq_[4] = 3./2.*(13.*cs_sqr-8.)*rho_+13./2./rho_*jSqr_;
        Meq_[5] = (2.*jxSqr_-(jySqr_+jzSqr_))/rho_;
        Meq_[6] = (jySqr_-jzSqr_)/rho_;
        Meq_[7] = jx_*jy_/rho_;
        Meq_[8] = jy_*jz_/rho_;
        Meq_[9] = jx_*jz_/rho_;
        Meq_[10] = 0.;
        Meq_[11] = 0.;
        Meq_[12] = 0.;

        T Nout[13];
        for (int i=0; i<13; ++i) {
            Nout[i] = cell[i];
        }

        for (int i=0; i<13; ++i) {
            for (int k=4; k<13; ++k) {
                Nout[i] -= 
                    lambda[k]/(T)eNormSqr[k]*(cell[k]-Meq_[k])*(T)e_[k][i];
            }
        }

        for (int i=0; i<13; ++i) {
            cell[i] = Nout[i];
        }

        for  (int i=0; i<13; ++i) {
            cell[i] -= descriptors::D3Q13Descriptor<T>::t[i];
        }

        return jSqr_/rho_/rho_;
    }
*/

    /// BGK collision step with density correction
    static T constRhoCollision( Cell<T,descriptors::D3Q13Descriptor>& cell,
                                T rho, const T u[Lattice::d],
                                T ratioRho,
                                T lambda_nu, T lambda_nu_prime)
    {
        const T uSqr = util::normSqr<T,Lattice::d>(u);

        return uSqr;
    }
}; // struct d3q13helpers

}  // namespace olb

#endif

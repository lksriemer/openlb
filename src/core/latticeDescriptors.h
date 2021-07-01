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
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef LATTICE_DESCRIPTORS_H
#define LATTICE_DESCRIPTORS_H

#include <vector>
#include "olbDebug.h"

namespace olb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {

    struct NoExternalField {
        static const int numScalars = 0;
        static const int numSpecies = 0;
    };

    struct Force2dDescriptor {
        static const int numScalars = 2;
        static const int numSpecies = 1;
        static const int forceBeginsAt = 0;
        static const int sizeOfForce   = 2;
    };

    struct Force3dDescriptor {
        static const int numScalars = 3;
        static const int numSpecies = 1;
        static const int forceBeginsAt = 0;
        static const int sizeOfForce   = 3;
    };

    template<typename T, typename ExternalField>
    class ExternalFieldArray {
    public:
        T* get(int index) {
            OLB_PRECONDITION( index < ExternalField::numScalars );
            return data+index;
        }
        T const* get(int index) const {
            OLB_PRECONDITION( index < ExternalField::numScalars );
            return data+index;
        }
    private:
        T data[ExternalField::numScalars];
    };

    template<typename T>
    class ExternalFieldArray<T,descriptors::NoExternalField> {
    public:
        T* get(unsigned index) {
            OLB_PRECONDITION( false );
            static T data = T();
            return &data;
        }
        T const* get(unsigned index) const {
            OLB_PRECONDITION( false );
            static T data = T();
            return &data;
        }
    };


    /// D2Q9 lattice
    template <typename T>
        struct D2Q9Descriptor
    {
        typedef NoExternalField ExternalField;
        enum { d = 2, q = 9 };      ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
        static const T t[q];        ///< lattice weights
        static const T invCs2;      ///< inverse square of speed of sound
    };

    /// D3Q13 lattice
    template <typename T>
        struct D3Q13Descriptor
    {
        typedef NoExternalField ExternalField;
        enum { d = 3, q = 13 };     ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
        static const T t[q];        ///< lattice weights
        static const T invCs2;      ///< inverse square of speed of sound
        static const T lambda_e;    ///< relaxation parameter for the bulk stress
        static const T lambda_h;    ///< additional relaxation parameter
    };

    /// D3Q15 lattice
    template <typename T>
        struct D3Q15Descriptor
    {
        typedef NoExternalField ExternalField;
        enum { d = 3, q = 15 };     ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
        static const T t[q];        ///< lattice weights
        static const T invCs2;      ///< inverse square of speed of sound
    };

    /// D3Q19 lattice
    template <typename T>
        struct D3Q19Descriptor
    {
        typedef NoExternalField ExternalField;
        enum { d = 3, q = 19 };     ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
        static const T t[q];        ///< lattice weights
        static const T invCs2;      ///< inverse square of speed of sound
    };

    /// D3Q27 lattice
    template <typename T>
        struct D3Q27Descriptor
    {
        typedef NoExternalField ExternalField;
        enum { d = 3, q = 27 };     ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
        static const T t[q];        ///< lattice weights
        static const T invCs2;      ///< inverse square of speed of sound
    };

    /// D2Q9 lattice with force
    template <typename T>
        struct ForcedD2Q9Descriptor
    {
        typedef Force2dDescriptor ExternalField;
        enum { d = 2, q = 9 };      ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
        static const T t[q];        ///< lattice weights
        static const T invCs2;      ///< inverse square of speed of sound
    };

    /// D3Q19 lattice with force
    template <typename T>
        struct ForcedD3Q19Descriptor
    {
        typedef Force3dDescriptor ExternalField;
        enum { d = 3, q = 19 };     ///< number of dimensions/distr. functions
        static const int c[q][d];   ///< lattice directions
        static const T t[q];        ///< lattice weights
        static const T invCs2;      ///< inverse square of speed of sound
    };




}  // namespace descriptors

}  // namespace olb

#endif

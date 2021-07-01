/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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
 * Parallel dynamics object -- generic template code
 */
#ifndef PARALLEL_DYNAMICS_HH
#define PARALLEL_DYNAMICS_HH

#include "complexGrids/mpiManager/mpiManager.h"
#include "parallelDynamics.h"

namespace olb {


#ifdef PARALLEL_MODE_MPI


////////////////////// Class ParallelDynamics /////////////////////////////

template<typename T, template<typename U> class Lattice>
ParallelDynamics<T,Lattice>::ParallelDynamics(Cell<T,Lattice>* baseCell_, int onProcessor_)
    : baseCell(baseCell_),
      onProcessor(onProcessor_)
{ }

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice>* ParallelDynamics<T,Lattice>::clone() const {
    return new ParallelDynamics(baseCell, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
                                             LatticeStatistics<T>& statistics_)
{
    if (baseCell) {
        baseCell -> collide(statistics_);
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::staticCollide (
        Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
        LatticeStatistics<T>& statistics_ )
{
    if (baseCell) {
        baseCell -> staticCollide(u, statistics_);
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::iniEquilibrium (
        Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d] )
{
    if (baseCell) {
        baseCell -> iniEquilibrium(rho, u);
    }
}

template<typename T, template<typename U> class Lattice>
T ParallelDynamics<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const {
    T rho;
    if (baseCell) {
        rho = baseCell -> computeRho();
    }
    singleton::mpi().bCast(&rho, 1, onProcessor);
    return rho;
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
    if (baseCell) {
        baseCell -> computeU(u);
    }
    singleton::mpi().bCast(u, Lattice<T>::d, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d] ) const
{
    if (baseCell) {
        baseCell -> getDynamics() -> computeJ(*baseCell, j);
    }
    singleton::mpi().bCast(j, Lattice<T>::d, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeStress (
        Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
        T pi[util::TensorVal<Lattice<T> >::n] ) const
{
    if (baseCell) {
        baseCell -> getDynamics() -> computeStress(*baseCell, rho, u, pi);
    }
    singleton::mpi().bCast(pi, util::TensorVal<Lattice<T> >::n, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeRhoU (
        Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
    if (baseCell) {
        baseCell -> getDynamics() -> computeRhoU(*baseCell, rho, u);
    }
    singleton::mpi().bCast(&rho, 1, onProcessor);
    singleton::mpi().bCast(u, Lattice<T>::d, onProcessor);
}


template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeAllMomenta (
        Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d],
        T pi[util::TensorVal<Lattice<T> >::n] ) const
{
    if (baseCell) {
        baseCell -> computeAllMomenta(rho, u, pi);
    }
    singleton::mpi().bCast(&rho, 1, onProcessor);
    singleton::mpi().bCast(u, Lattice<T>::d, onProcessor);
    singleton::mpi().bCast(pi, util::TensorVal<Lattice<T> >::n, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computePopulations(Cell<T,Lattice> const& cell, T* f) const
{
    if (baseCell) {
        baseCell -> computePopulations(f);
    }
    singleton::mpi().bCast(f, Lattice<T>::q, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeExternalField (
        Cell<T,Lattice> const& cell, int pos, int size, T* ext ) const
{
    if (baseCell) {
        baseCell -> computeExternalField(pos, size, ext);
    }
    singleton::mpi().bCast(ext, Lattice<T>::ExternalField::numScalars, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho) {
    if (baseCell) {
        baseCell -> defineRho(rho);
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]) {
    if (baseCell) {
        baseCell -> defineU(u);
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineRhoU(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d]) {
    if (baseCell) {
        baseCell -> defineRhoU(rho, u);
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineAllMomenta (
        Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
        const T pi[util::TensorVal<Lattice<T> >::n] )
{
    if (baseCell) {
        baseCell -> defineAllMomenta(rho, u, pi);
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::definePopulations(Cell<T,Lattice>& cell, const T* f)
{
    if (baseCell) {
        baseCell -> definePopulations(f);
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineExternalField (
        Cell<T,Lattice>& cell, int pos, int size, const T* ext )
{
    if (baseCell) {
        baseCell -> defineExternalField(pos, size, ext);
    }
}

template<typename T, template<typename U> class Lattice>
T ParallelDynamics<T,Lattice>::getOmega() const {
    T omega;
    if (baseCell) {
        omega = baseCell -> getDynamics() -> getOmega();
    }
    singleton::mpi().bCast(&omega, 1, onProcessor);
    return omega;
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::setOmega(T omega_) {
    if (baseCell) {
        baseCell -> getDynamics() -> setOmega(omega_);
    }
}

////////////////// Class ConstParallelDynamics /////////////////////////

template<typename T, template<typename U> class Lattice>
ConstParallelDynamics<T,Lattice>::ConstParallelDynamics (
        Cell<T,Lattice> const* baseCell_, int onProcessor_ )
    : baseCell(baseCell_),
      onProcessor(onProcessor_)
{ }

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice>* ConstParallelDynamics<T,Lattice>::clone() const {
    return new ConstParallelDynamics(baseCell, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
                                                  LatticeStatistics<T>& statistics_)
{ }

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::staticCollide (
        Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
        LatticeStatistics<T>& statistics_ )
{ }

template<typename T, template<typename U> class Lattice>
T ConstParallelDynamics<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const {
    T rho;
    if (baseCell) {
        rho = baseCell -> computeRho();
    }
    singleton::mpi().bCast(&rho, 1, onProcessor);
    return rho;
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
    if (baseCell) {
        baseCell -> computeU(u);
    }
    singleton::mpi().bCast(u, Lattice<T>::d, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d] ) const
{
    if (baseCell) {
        baseCell -> getDynamics() -> computeJ(*baseCell, j);
    }
    singleton::mpi().bCast(j, Lattice<T>::d, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeStress (
        Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
        T pi[util::TensorVal<Lattice<T> >::n] ) const
{
    if (baseCell) {
        baseCell -> getDynamics() -> computeStress(*baseCell, rho, u, pi);
    }
    singleton::mpi().bCast(pi, util::TensorVal<Lattice<T> >::n, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeRhoU (
        Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
    if (baseCell) {
        baseCell -> getDynamics() -> computeRhoU(*baseCell, rho, u);
    }
    singleton::mpi().bCast(&rho, 1, onProcessor);
    singleton::mpi().bCast(u, Lattice<T>::d, onProcessor);
}


template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeAllMomenta (
        Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d],
        T pi[util::TensorVal<Lattice<T> >::n] ) const
{
    if (baseCell) {
        baseCell -> getDynamics() -> computeAllMomenta(*baseCell, rho, u, pi);
    }
    singleton::mpi().bCast(&rho, 1, onProcessor);
    singleton::mpi().bCast(u, Lattice<T>::d, onProcessor);
    singleton::mpi().bCast(pi, util::TensorVal<Lattice<T> >::n, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computePopulations (
        Cell<T,Lattice> const& cell, T* f ) const
{
    if (baseCell) {
        baseCell -> computePopulations(f);
    }
    singleton::mpi().bCast(f, Lattice<T>::q, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeExternalField (
        Cell<T,Lattice> const& cell, int pos, int size, T* ext ) const
{
    if (baseCell) {
        baseCell -> computeExternalField(pos, size, ext);
    }
    singleton::mpi().bCast(ext, Lattice<T>::ExternalField::numScalars, onProcessor);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{ }

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d])
{ }

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::defineRhoU(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d])
{ }

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::defineAllMomenta (
        Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
        const T pi[util::TensorVal<Lattice<T> >::n] )
{ }

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::definePopulations (
        Cell<T,Lattice>& cell, const T* f)
{ }

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::defineExternalField (
        Cell<T,Lattice>& cell, int pos, int size, const T* ext )
{ }

template<typename T, template<typename U> class Lattice>
T ConstParallelDynamics<T,Lattice>::getOmega() const {
    T omega;
    if (baseCell) {
        omega = baseCell -> getDynamics() -> getOmega();
    }
    singleton::mpi().bCast(&omega, 1, onProcessor);
    return omega;
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::setOmega(T omega_)
{ }

#endif


}

#endif // defined MULTIBLOCK_DYNAMICS_H

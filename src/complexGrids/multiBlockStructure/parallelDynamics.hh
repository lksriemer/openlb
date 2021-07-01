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
ParallelDynamics<T,Lattice>::ParallelDynamics(std::vector<Cell<T,Lattice>*>& baseCell_)
    : baseCell(baseCell_)
{ }

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice>* ParallelDynamics<T,Lattice>::clone() const {
    return new ParallelDynamics(baseCell);
}

template<typename T, template<typename U> class Lattice>
T ParallelDynamics<T,Lattice>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
    T eq = T();
    if (baseCell[0]) {
        eq = baseCell[0] -> computeEquilibrium(iPop, rho, u, uSqr);
    }
    singleton::mpi().sendToMaster(&eq, 1, baseCell[0]);
    return eq;
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::iniEquilibrium(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d]) {
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> getDynamics() -> iniEquilibrium(*baseCell[iCell], rho, u);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
                                          LatticeStatistics<T>& statistics_)
{
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> collide(statistics_);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::staticCollide (
        Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
        LatticeStatistics<T>& statistics_ )
{
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> staticCollide(u, statistics_);
        }
    }
}

template<typename T, template<typename U> class Lattice>
T ParallelDynamics<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const {
    T rho = T();
    if (baseCell[0]) {
        rho = baseCell[0] -> computeRho();
    }
    singleton::mpi().sendToMaster(&rho, 1, baseCell[0]);
    return rho;
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computeU(u);
    }
    singleton::mpi().sendToMaster(u, Lattice<T>::d, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> getDynamics() -> computeJ(*baseCell[0], j);
    }
    singleton::mpi().sendToMaster(j, Lattice<T>::d, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeStress (
        Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
        T pi[util::TensorVal<Lattice<T> >::n] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> getDynamics() -> computeStress(*baseCell[0], rho, u, pi);
    }
    singleton::mpi().sendToMaster(pi, util::TensorVal<Lattice<T> >::n, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeRhoU (
        Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computeRhoU(rho, u);
    }
    singleton::mpi().sendToMaster(&rho, 1, baseCell[0]);
    singleton::mpi().sendToMaster(u, Lattice<T>::d, baseCell[0]);
}


template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeAllMomenta (
        Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d],
        T pi[util::TensorVal<Lattice<T> >::n] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computeAllMomenta(rho, u, pi);
    }
    singleton::mpi().sendToMaster(&rho, 1, baseCell[0]);
    singleton::mpi().sendToMaster(u, Lattice<T>::d, baseCell[0]);
    singleton::mpi().sendToMaster(pi, util::TensorVal<Lattice<T> >::n, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computePopulations(Cell<T,Lattice> const& cell, T* f) const
{
    if (baseCell[0]) {
        baseCell[0] -> computePopulations(f);
    }
    singleton::mpi().sendToMaster(f, Lattice<T>::q, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::computeExternalField (
        Cell<T,Lattice> const& cell, int pos, int size, T* ext ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computeExternalField(pos, size, ext);
    }
    singleton::mpi().sendToMaster(ext, Lattice<T>::ExternalField::numScalars, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho) {
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> defineRho(rho);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]) {
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> defineU(u);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineRhoU(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d]) {
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> defineU(u);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineAllMomenta (
        Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
        const T pi[util::TensorVal<Lattice<T> >::n] )
{
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> defineAllMomenta(rho, u, pi);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::definePopulations(Cell<T,Lattice>& cell, const T* f)
{
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> definePopulations(f);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::defineExternalField (
        Cell<T,Lattice>& cell, int pos, int size, const T* ext )
{
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> defineExternalField(pos, size, ext);
        }
    }
}

template<typename T, template<typename U> class Lattice>
T ParallelDynamics<T,Lattice>::getOmega() const {
    T omega;
    if (baseCell[0]) {
        omega = baseCell[0] -> getDynamics() -> getOmega();
    }
    singleton::mpi().sendToMaster(&omega, 1, baseCell[0]);
    return omega;
}

template<typename T, template<typename U> class Lattice>
void ParallelDynamics<T,Lattice>::setOmega(T omega_) {
    for (unsigned iCell=0; iCell<baseCell.size(); ++iCell) {
        if (baseCell[iCell]) {
            baseCell[iCell] -> getDynamics() -> setOmega(omega_);
        }
    }
}

////////////////// Class ConstParallelDynamics /////////////////////////

template<typename T, template<typename U> class Lattice>
ConstParallelDynamics<T,Lattice>::ConstParallelDynamics(std::vector<Cell<T,Lattice> const*>& baseCell_)
    : baseCell(baseCell_)
{ }

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice>* ConstParallelDynamics<T,Lattice>::clone() const {
    return new ConstParallelDynamics(baseCell);
}

template<typename T, template<typename U> class Lattice>
T ConstParallelDynamics<T,Lattice>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
    T eq = T();
    if (baseCell[0]) {
        eq = baseCell[0] -> computeEquilibrium(iPop, rho, u, uSqr);
    }
    singleton::mpi().sendToMaster(&eq, 1, baseCell[0]);
    return eq;
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::iniEquilibrium(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d])
{ }

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
    if (baseCell[0]) {
        rho = baseCell[0] -> computeRho();
    }
    singleton::mpi().sendToMaster(&rho, 1, baseCell[0]);
    return rho;
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computeU(u);
    }
    singleton::mpi().sendToMaster(u, Lattice<T>::d, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> getDynamics() -> computeJ(*baseCell[0], j);
    }
    singleton::mpi().sendToMaster(j, Lattice<T>::d, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeStress (
        Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
        T pi[util::TensorVal<Lattice<T> >::n] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> getDynamics() -> computeStress(*baseCell[0], rho, u, pi);
    }
    singleton::mpi().sendToMaster(pi, util::TensorVal<Lattice<T> >::n, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeRhoU (
        Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computeRhoU(rho, u);
    }
    singleton::mpi().sendToMaster(&rho, 1, baseCell[0]);
    singleton::mpi().sendToMaster(u, Lattice<T>::d, baseCell[0]);
}


template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeAllMomenta (
        Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d],
        T pi[util::TensorVal<Lattice<T> >::n] ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computeAllMomenta(rho, u, pi);
    }
    singleton::mpi().sendToMaster(&rho, 1, baseCell[0]);
    singleton::mpi().sendToMaster(u, Lattice<T>::d, baseCell[0]);
    singleton::mpi().sendToMaster(pi, util::TensorVal<Lattice<T> >::n, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computePopulations (
        Cell<T,Lattice> const& cell, T* f ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computePopulations(f);
    }
    singleton::mpi().sendToMaster(f, Lattice<T>::q, baseCell[0]);
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::computeExternalField (
        Cell<T,Lattice> const& cell, int pos, int size, T* ext ) const
{
    if (baseCell[0]) {
        baseCell[0] -> computeExternalField(pos, size, ext);
    }
    singleton::mpi().sendToMaster(ext, Lattice<T>::ExternalField::numScalars, baseCell[0]);
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
    if (baseCell[0]) {
        omega = baseCell[0] -> getDynamics() -> getOmega();
    }
    singleton::mpi().sendToMaster(&omega, 1, baseCell[0]);
    return omega;
}

template<typename T, template<typename U> class Lattice>
void ConstParallelDynamics<T,Lattice>::setOmega(T omega_)
{ }

#endif


}

#endif // defined MULTIBLOCK_DYNAMICS_H

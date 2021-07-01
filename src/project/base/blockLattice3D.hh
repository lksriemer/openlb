/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  OMP parallel code by Mathias Krause, Copyright (C) 2007
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
 * The dynamics of a 3D block lattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_3D_HH
#define BLOCK_LATTICE_3D_HH

#include <algorithm>
#include "blockLattice3D.h"
#include "dynamics.h"
#include "cell.h"
#include "lbHelpers.h"
#include "util.h"
#include "ompManager.h"
#include "loadBalancer.h"


namespace olb {

////////////////////// Class BlockLattice3D /////////////////////////

/** \param nx_ lattice width (first index)
 *  \param ny_ lattice height (second index)
 *  \param nz_ lattice depth (third index)
 */
template<typename T, template<typename U> class Lattice>
BlockLattice3D<T,Lattice>::
    BlockLattice3D(int nx_, int ny_, int nz_)
        : nx(nx_), ny(ny_), nz(nz_)
{
    allocateMemory();
    resetPostProcessors();
}

/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */
template<typename T, template<typename U> class Lattice>
BlockLattice3D<T,Lattice>::~BlockLattice3D()
{
    releaseMemory();
    clearPostProcessors();
}

/** The whole lattice is duplicated, but not the dynamics objects
 * pointed to by the cells.
 */
template<typename T, template<typename U> class Lattice>
BlockLattice3D<T,Lattice>::BlockLattice3D (
        BlockLattice3D<T,Lattice> const& rhs )
{
    nx = rhs.nx;
    ny = rhs.ny;
    nz = rhs.ny;
    allocateMemory();
    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            for (int iZ=0; iZ<nz; ++iZ) {
                grid[iX][iY][iZ] = rhs.grid[iX][iY][iZ];
            }
        }
    }
}

/** The current lattice is deallocated, then the lattice from the rhs
 * is duplicated.
 * \param rhs the lattice to be duplicated
 */
template<typename T, template<typename U> class Lattice>
BlockLattice3D<T,Lattice>& BlockLattice3D<T,Lattice>::operator= (
        BlockLattice3D<T,Lattice> const& rhs )
{
    BlockLattice3D<T,Lattice> tmp(rhs);
    swap(tmp);
    return *this;
}

/** The swap is efficient, in the sense that only pointers to the 
 * lattice are copied, and not the lattice itself.
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::swap(BlockLattice3D& rhs) {
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(rawData, rhs.rawData);
    std::swap(grid, rhs.grid);
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::initialize() {
    postProcess();
}

/** The dynamics object is not duplicated: all cells of the rectangular
 * domain point to the same dynamics.
 *
 * The dynamics object is not owned by the BlockLattice3D object, its
 * memory management must be taken care of by the user.
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::defineDynamics (
        int x0, int x1, int y0, int y1, int z0, int z1,
        Dynamics<T,Lattice>* dynamics )
{
    OLB_PRECONDITION(x0>=0 && x1<nx);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=0 && y1<ny);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=0 && z1<nz);
    OLB_PRECONDITION(z1>=z0);

    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
                grid[iX][iY][iZ].defineDynamics(dynamics);
            }
        }
    }
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::specifyStatisticsStatus (
        int x0, int x1, int y0, int y1, int z0, int z1, bool status)
{
    OLB_PRECONDITION(x0>=0 && x1<nx);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=0 && y1<ny);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=0 && z1<nz);
    OLB_PRECONDITION(z1>=z0);

    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
                grid[iX][iY][iZ].specifyStatisticsStatus(status);
            }
        }
    }
}

/** 
 * This method is automatically parallelized if your compiler understands
 * OpenMP
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::collide (
    int x0, int x1, int y0, int y1, int z0, int z1)
{
    OLB_PRECONDITION(x0>=0 && x1<nx);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=0 && y1<ny);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=0 && z1<nz);
    OLB_PRECONDITION(z1>=z0);

    int iX, iY, iZ;
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for private (iY,iZ) schedule(dynamic,1)
    #endif
    for (iX=x0; iX<=x1; ++iX) {
        for (iY=y0; iY<=y1; ++iY) {
            for (iZ=z0; iZ<=z1; ++iZ) {
                grid[iX][iY][iZ].collide();
                grid[iX][iY][iZ].revert();
            }
        }
    }
}

/** \sa collide(int,int,int,int) */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::collide() {
    collide(0, nx-1, 0, ny-1, 0, nz-1);
}

/** 
 * A useful method for initializing the flow field to a given velocity
 * profile.
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::staticCollide (
        int x0, int x1, int y0, int y1, int z0, int z1,
        TensorField3D<T,3> const& u )
{
    OLB_PRECONDITION(x0>=0 && x1<nx);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=0 && y1<ny);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=0 && z1<nz);
    OLB_PRECONDITION(z1>=z0);

    int iX, iY, iZ;
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for private (iY,iZ) schedule(dynamic,1)
    #endif
    for (iX=x0; iX<=x1; ++iX) {
        for (iY=y0; iY<=y1; ++iY) {
            for (iZ=z0; iZ<=z1; ++iZ) {
                grid[iX][iY][iZ].staticCollide(u.get(iX,iY,iZ));
                grid[iX][iY][iZ].revert();
            }
        }
    }
}

/** \sa collide(int,int,int,int) */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::staticCollide(TensorField3D<T,3> const& u) {
    staticCollide(0, nx-1, 0, ny-1, 0, nz-1, u);
}

/** The distribution function never leave the rectangular domain. On the
 * domain boundaries, the (outgoing) distribution functions that should 
 * be streamed outside are simply left untouched.
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::stream (
    int x0, int x1, int y0, int y1, int z0, int z1, bool periodic )
{
    OLB_PRECONDITION(x0>=0 && x1<nx);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=0 && y1<ny);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=0 && z1<nz);
    OLB_PRECONDITION(z1>=z0);

    bulkStream(x0+1,x1-1, y0+1,y1-1, z0+1,z1-1);

    boundaryStream(x0,x1,y0,y1,z0,z1, x0,x0, y0,y1, z0,z1);
    boundaryStream(x0,x1,y0,y1,z0,z1, x1,x1, y0,y1, z0,z1);
    boundaryStream(x0,x1,y0,y1,z0,z1, x0+1,x1-1, y0,y0, z0,z1);
    boundaryStream(x0,x1,y0,y1,z0,z1, x0+1,x1-1, y1,y1, z0,z1);
    boundaryStream(x0,x1,y0,y1,z0,z1, x0+1,x1-1, y0+1,y1-1, z0,z0);
    boundaryStream(x0,x1,y0,y1,z0,z1, x0+1,x1-1, y0+1,y1-1, z1,z1);

    if (periodic) {
        makePeriodic();
    }
}

/** Post-processing steps are called at the end of this method.
 * \sa stream(int,int,int,int,int,int) */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::stream(bool periodic) {
    stream(0, nx-1, 0, ny-1, 0, nz-1);

    if (periodic) {
        makePeriodic();
    }

    postProcess();
}

/** This operation is more efficient than a successive application of
 * collide(int,int,int,int,int,int) and stream(int,int,int,int,int,int),
 * because memory is traversed only once instead of twice.
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::collideAndStream (
    int x0, int x1, int y0, int y1, int z0, int z1, bool periodic )
{
    OLB_PRECONDITION(x0>=0 && x1<nx);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=0 && y1<ny);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=0 && z1<nz);
    OLB_PRECONDITION(z1>=z0);

    collide(x0,x0, y0,y1, z0,z1);
    collide(x1,x1, y0,y1, z0,z1);
    collide(x0+1,x1-1, y0,y0, z0,z1);
    collide(x0+1,x1-1, y1,y1, z0,z1);
    collide(x0+1,x1-1, y0+1,y1-1, z0,z0);
    collide(x0+1,x1-1, y0+1,y1-1, z1,z1);

    bulkCollideAndStream(x0+1,x1-1, y0+1,y1-1, z0+1,z1-1);

    boundaryStream(x0,x1,y0,y1,z0,z1, x0,x0, y0,y1, z0,z1);
    boundaryStream(x0,x1,y0,y1,z0,z1, x1,x1, y0,y1, z0,z1);
    boundaryStream(x0,x1,y0,y1,z0,z1, x0+1,x1-1, y0,y0, z0,z1);
    boundaryStream(x0,x1,y0,y1,z0,z1, x0+1,x1-1, y1,y1, z0,z1);
    boundaryStream(x0,x1,y0,y1,z0,z1, x0+1,x1-1, y0+1,y1-1, z0,z0);
    boundaryStream(x0,x1,y0,y1,z0,z1, x0+1,x1-1, y0+1,y1-1, z1,z1);

    if (periodic) {
        makePeriodic();
    }
}

/** Post-processing steps are called at the end of this method.
 * \sa collideAndStream(int,int,int,int,int,int) */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::collideAndStream(bool periodic) {
    collideAndStream(0, nx-1, 0, ny-1, 0, nz-1);

    if (periodic) {
        makePeriodic();
    }

    postProcess();
}

template<typename T, template<typename U> class Lattice>
T BlockLattice3D<T,Lattice>::computeAverageDensity (
        int x0, int x1, int y0, int y1, int z0, int z1) const
{
    T sumRho = T();
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
                T rho, u[Lattice<T>::d];
                get(iX,iY,iZ).computeRhoU(rho, u);
                sumRho += rho;
            }
        }
    }
    return sumRho / (T)(x1-x0+1) / (T)(y1-y0+1) / (T)(z1-z0+1);
}

template<typename T, template<typename U> class Lattice>
T BlockLattice3D<T,Lattice>::computeAverageDensity() const {
    return computeAverageDensity(0, nx-1, 0, ny-1, 0, nz-1);
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::stripeOffDensityOffset (
        int x0, int x1, int y0, int y1, int z0, int z1, T offset )
{
    for (int iX=x0; iX<=x1; ++iX) {
        for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
                for (int iPop=0; iPop<Lattice<T>::q; ++iPop) {
                    get(iX,iY,iZ)[iPop] -= Lattice<T>::t[iPop] * offset;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::stripeOffDensityOffset(T offset) {
    stripeOffDensityOffset(0, nx-1, 0, ny-1, 0, nz-1, offset);
}


template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::addPostProcessor (
    PostProcessorGenerator3D<T,Lattice> const& ppGen )
{
    postProcessors.push_back(ppGen.generate(*this));
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::resetPostProcessors() {
    clearPostProcessors();
    StatPPGenerator3D<T,Lattice> statPPGenerator;
    addPostProcessor(statPPGenerator);
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::clearPostProcessors() {
    typename std::vector<PostProcessor3D<T,Lattice>*>::iterator ppIt
        = postProcessors.begin();
    for (; ppIt != postProcessors.end(); ++ppIt) {
        delete *ppIt;
    }
    postProcessors.clear();
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::postProcess() {
    for (unsigned iPr=0; iPr<postProcessors.size(); ++iPr) {
        postProcessors[iPr] -> process(*this);
    }
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::postProcess (
        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
    for (unsigned iPr=0; iPr<postProcessors.size(); ++iPr) {
        postProcessors[iPr] -> process(*this, x0_, x1_, y0_, y1_, z0_, z1_);
    }
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::subscribeReductions(Reductor<T>& reductor) {
    for (unsigned iPr=0; iPr<postProcessors.size(); ++iPr) {
        postProcessors[iPr] -> subscribeReductions(*this, &reductor);
    }
}

template<typename T, template<typename U> class Lattice>
LatticeStatistics<T>& BlockLattice3D<T,Lattice>::getStatistics() {
    return statistics;
}

template<typename T, template<typename U> class Lattice>
LatticeStatistics<T> const&
    BlockLattice3D<T,Lattice>::getStatistics() const
{
    return statistics;
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::allocateMemory() {
    rawData = new Cell<T,Lattice> [nx*ny*nz];
    grid    = new Cell<T,Lattice>** [nx];
    for (int iX=0; iX<nx; ++iX) {
        grid[iX] = new Cell<T,Lattice>* [ny];
        for (int iY=0; iY<ny; ++iY) {
            grid[iX][iY] = rawData + nz*(iY+ny*iX);
        }
    }
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::releaseMemory() {
    delete [] rawData;
    for (int iX=0; iX<nx; ++iX) {
      delete [] grid[iX];
    }
    delete [] grid;
}

/** This method is slower than bulkStream(int,int,int,int), because it must
 * be verified which distribution functions are to be kept from leaving
 * the domain.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::boundaryStream (
    int lim_x0, int lim_x1, int lim_y0, int lim_y1, int lim_z0, int lim_z1,
    int x0, int x1, int y0, int y1, int z0, int z1 )
{
    OLB_PRECONDITION(lim_x0>=0 && lim_x1<nx);
    OLB_PRECONDITION(lim_x1>=lim_x0);
    OLB_PRECONDITION(lim_y0>=0 && lim_y1<ny);
    OLB_PRECONDITION(lim_y1>=lim_y0);
    OLB_PRECONDITION(lim_z0>=0 && lim_z1<nz);
    OLB_PRECONDITION(lim_z1>=lim_z0);

    OLB_PRECONDITION(x0>=lim_x0 && x1<=lim_x1);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=lim_y0 && y1<=lim_y1);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=lim_z0 && z1<=lim_z1);
    OLB_PRECONDITION(z1>=z0);

    int iX, iY, iZ, iPop;

    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for private(iY,iZ,iPop)
    #endif
    for (iX=x0; iX<=x1; ++iX) {
        for (iY=y0; iY<=y1; ++iY) {
            for (iZ=z0; iZ<=z1; ++iZ) {
                for (iPop=1; iPop<=Lattice<T>::q/2; ++iPop) {
                    int nextX = iX + Lattice<T>::c[iPop][0];
                    int nextY = iY + Lattice<T>::c[iPop][1];
                    int nextZ = iZ + Lattice<T>::c[iPop][2];
                    if ( nextX>=lim_x0 && nextX<=lim_x1 &&
                         nextY>=lim_y0 && nextY<=lim_y1 &&
                         nextZ>=lim_z0 && nextZ<=lim_z1 )
                    {
                        std::swap(grid[iX][iY][iZ][iPop+Lattice<T>::q/2],
                                  grid[nextX][nextY][nextZ][iPop]);
                    }
                }
            }
        }
    }
}

/** This method is faster than boundaryStream(int,int,int,int,int,int), but it
 * is erroneous when applied to boundary cells.
 * \sa stream(int,int,int,int,int,int)
 * \sa stream()
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::bulkStream (
    int x0, int x1, int y0, int y1, int z0, int z1 )
{
    OLB_PRECONDITION(x0>=0 && x1<nx);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=0 && y1<ny);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=0 && z1<nz);
    OLB_PRECONDITION(z1>=z0);

    int iX, iY, iZ, iPop;
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for private(iY,iZ,iPop)
    #endif
    for (iX=x0; iX<=x1; ++iX) {
        for (iY=y0; iY<=y1; ++iY) {
            for (iZ=z0; iZ<=z1; ++iZ) {
                for (iPop=1; iPop<=Lattice<T>::q/2; ++iPop) {
                    int nextX = iX + Lattice<T>::c[iPop][0];
                    int nextY = iY + Lattice<T>::c[iPop][1];
                    int nextZ = iZ + Lattice<T>::c[iPop][2];
                    std::swap(grid[iX][iY][iZ][iPop+Lattice<T>::q/2],
                              grid[nextX][nextY][nextZ][iPop]);
                }
            }
        }
    }
}

/** This method is fast, but it is erroneous when applied to boundary
 * cells.
 * \sa collideAndStream(int,int,int,int,int,int)
 * \sa collideAndStream()
 */
template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::bulkCollideAndStream (
    int x0, int x1, int y0, int y1, int z0, int z1 )
{
    OLB_PRECONDITION(x0>=0 && x1<nx);
    OLB_PRECONDITION(x1>=x0);
    OLB_PRECONDITION(y0>=0 && y1<ny);
    OLB_PRECONDITION(y1>=y0);
    OLB_PRECONDITION(z0>=0 && z1<nz);
    OLB_PRECONDITION(z1>=z0);

    #ifdef PARALLEL_MODE_OMP
        if (omp.get_size() <= x1-x0+1) {
            #pragma omp parallel
                {
                loadBalancer loadbalance(omp.get_rank(), omp.get_size(), x1-x0+1, x0);
                int iX, iY, iZ, iPop;

                iX=loadbalance.get_firstGlobNum();
                for (int iY=y0; iY<=y1; ++iY) {
                    for (int iZ=z0; iZ<=z1; ++iZ) {
                        grid[iX][iY][iZ].collide();
                        grid[iX][iY][iZ].revert();
                        }
                }

                for (iX=loadbalance.get_firstGlobNum()+1; iX<=loadbalance.get_lastGlobNum(); ++iX) {
                    for (iY=y0; iY<=y1; ++iY) {
                        for (iZ=z0; iZ<=z1; ++iZ) {
                            grid[iX][iY][iZ].collide();
                            /** The method beneath doesnt work with Intel compiler 9.1044 and 9.1046 for Itanium prozessors
                             *    lbHelpers<T,Lattice>::swapAndStream3D(grid, iX, iY, iZ);
                             *  Therefore we use:
                             */
                                int half = Lattice<T>::q/2;
                                for (int iPop=1; iPop<=half; ++iPop) {
                                    int nextX = iX + Lattice<T>::c[iPop][0];
                                    int nextY = iY + Lattice<T>::c[iPop][1];
                                    int nextZ = iZ + Lattice<T>::c[iPop][2];
                                    T fTmp                          = grid[iX][iY][iZ][iPop];
                                    grid[iX][iY][iZ][iPop]          = grid[iX][iY][iZ][iPop+half];
                                    grid[iX][iY][iZ][iPop+half]     = grid[nextX][nextY][nextZ][iPop];
                                    grid[nextX][nextY][nextZ][iPop] = fTmp;
                                }
                        }
                    }
                }

                #pragma omp barrier

                iX=loadbalance.get_firstGlobNum();
                for (iY=y0; iY<=y1; ++iY) {
                    for (iZ=z0; iZ<=z1; ++iZ) {
                        for (iPop=1; iPop<=Lattice<T>::q/2; ++iPop) {
                            int nextX = iX + Lattice<T>::c[iPop][0];
                            int nextY = iY + Lattice<T>::c[iPop][1];
                            int nextZ = iZ + Lattice<T>::c[iPop][2];
                            std::swap(grid[iX][iY][iZ][iPop+Lattice<T>::q/2],
                                grid[nextX][nextY][nextZ][iPop]);
                        }
                    }
                }
                }
            }
        else {
            for (int iX=x0; iX<=x1; ++iX) {
                for (int iY=y0; iY<=y1; ++iY) {
                    for (int iZ=z0; iZ<=z1; ++iZ) {
                        grid[iX][iY][iZ].collide();
                        lbHelpers<T,Lattice>::swapAndStream3D(grid, iX, iY, iZ);
                    }
                }
            }
        }
    #else
        for (int iX=x0; iX<=x1; ++iX) {
            for (int iY=y0; iY<=y1; ++iY) {
                for (int iZ=z0; iZ<=z1; ++iZ) {
                    grid[iX][iY][iZ].collide();
                    lbHelpers<T,Lattice>::swapAndStream3D(grid, iX, iY, iZ);
                }
            }
        }
    #endif
}

template<typename T, template<typename U> class Lattice>
void BlockLattice3D<T,Lattice>::makePeriodic()
{
    for (int iX = 0; iX < getNx(); ++iX)
    {
        for (int iY = 0; iY < getNy(); ++iY)
        {
            for (int iPop=1; iPop<=Lattice<T>::q/2; ++iPop)
            {
                int nextX = iX + Lattice<T>::c[iPop][0];
                int nextY = iY + Lattice<T>::c[iPop][1];

                int nextZ = Lattice<T>::c[iPop][2];
                if ( nextX<0 || nextX>=getNx() ||
                     nextY<0 || nextY>=getNy() ||
                     nextZ<0 || nextZ>=getNz() )
                {
                    nextX = (nextX+getNx())%getNx();
                    nextY = (nextY+getNy())%getNy();
                    nextZ = (nextZ+getNz())%getNz();
                    std::swap (
                        grid[iX]   [iY]   [0]    [iPop+Lattice<T>::q/2],
                        grid[nextX][nextY][nextZ][iPop] );
                }

                nextX = iX + Lattice<T>::c[iPop][0];
                nextY = iY + Lattice<T>::c[iPop][1];
                nextZ = getNz()-1 + Lattice<T>::c[iPop][2];
                if ( nextX<0 || nextX>=getNx() ||
                     nextY<0 || nextY>=getNy() ||
                     nextZ<0 || nextZ>=getNz() )
                {
                    nextX = (nextX+getNx())%getNx();
                    nextY = (nextY+getNy())%getNy();
                    nextZ = (nextZ+getNz())%getNz();
                    std::swap (
                        grid[iX]   [iY]   [getNz()-1][iPop+Lattice<T>::q/2],
                        grid[nextX][nextY][nextZ]    [iPop] );
                }
            }
        }
    }

    for (int iX = 0; iX < getNx(); ++iX)
    {
        for (int iZ = 0; iZ < getNz(); ++iZ)
        {
            for (int iPop=1; iPop<=Lattice<T>::q/2; ++iPop)
            {
                int nextX = iX + Lattice<T>::c[iPop][0];
                int nextZ = iZ + Lattice<T>::c[iPop][2];

                int nextY = Lattice<T>::c[iPop][1];
                if ( nextX<0 || nextX>=getNx() ||
                     nextY<0 || nextY>=getNy() ||
                     nextZ<0 || nextZ>=getNz() )
                {
                    nextX = (nextX+getNx())%getNx();
                    nextY = (nextY+getNy())%getNy();
                    nextZ = (nextZ+getNz())%getNz();
                    std::swap (
                        grid[iX]   [0]    [iZ]   [iPop+Lattice<T>::q/2],
                        grid[nextX][nextY][nextZ][iPop] );
                }

                nextX = iX + Lattice<T>::c[iPop][0];
                nextZ = iZ + Lattice<T>::c[iPop][2];

                nextY = getNy() - 1 + Lattice<T>::c[iPop][1];
                if ( nextX<0 || nextX>=getNx() ||
                     nextY<0 || nextY>=getNy() ||
                     nextZ<0 || nextZ>=getNz() )
                {
                    nextX = (nextX+getNx())%getNx();
                    nextY = (nextY+getNy())%getNy();
                    nextZ = (nextZ+getNz())%getNz();
                    std::swap (
                        grid[iX]   [getNy()-1][iZ]   [iPop+Lattice<T>::q/2],
                        grid[nextX][nextY]    [nextZ][iPop] );
                }
            }
        }
    }

    for (int iY = 0; iY < getNy(); ++iY)
    {
        for (int iZ = 0; iZ < getNz(); ++iZ)
        {
            for (int iPop=1; iPop<=Lattice<T>::q/2; ++iPop)
            {
                int nextY = iY + Lattice<T>::c[iPop][1];
                int nextZ = iZ + Lattice<T>::c[iPop][2];

                int nextX = Lattice<T>::c[iPop][0];
                if ( nextX<0 || nextX>=getNx() ||
                     nextY<0 || nextY>=getNy() ||
                     nextZ<0 || nextZ>=getNz() )
                {
                    nextX = (nextX+getNx())%getNx();
                    nextY = (nextY+getNy())%getNy();
                    nextZ = (nextZ+getNz())%getNz();
                    std::swap (
                        grid[0]    [iY]   [iZ]   [iPop+Lattice<T>::q/2],
                        grid[nextX][nextY][nextZ][iPop] );
                }

                nextY = iY + Lattice<T>::c[iPop][1];
                nextZ = iZ + Lattice<T>::c[iPop][2];

                nextX = getNx() - 1 + Lattice<T>::c[iPop][0];
                if ( nextX<0 || nextX>=getNx() ||
                     nextY<0 || nextY>=getNy() ||
                     nextZ<0 || nextZ>=getNz() )
                {
                    nextX = (nextX+getNx())%getNx();
                    nextY = (nextY+getNy())%getNy();
                    nextZ = (nextZ+getNz())%getNz();
                    std::swap (
                        grid[getNx()-1][iY]   [iZ]   [iPop+Lattice<T>::q/2],
                        grid[nextX]    [nextY][nextZ][iPop] );
                }
            }
        }
    }
}

}  // namespace olb

#endif

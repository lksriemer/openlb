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
 * Interface for post-processing steps -- header file.
 */
#ifndef N_BLOCK_LATTICE_POST_PROCESSING_H
#define N_BLOCK_LATTICE_POST_PROCESSING_H

#include <vector>
#include "viscoPostProcessing2D.h"

namespace olb {

/////////////////// 2D Postprocessing ///////////////////////////////


template<typename T, template<typename U> class Lattice>
ViscoelasticCoupler2D<T,Lattice>::ViscoelasticCoupler2D(std::vector<BlockStructure2D<T,Lattice>*> partners_) :
        partners(partners_)
{
}

template<typename T, template<typename U> class Lattice>
void ViscoelasticCoupler2D<T,Lattice>::process(BlockLattice2D<T,Lattice>& blockLattice) {
    processSubDomain(blockLattice, 0, blockLattice.getNx()-1, 0, blockLattice.getNy()-1);
}

template<typename T, template<typename U> class Lattice>
void ViscoelasticCoupler2D<T,Lattice>::processSubDomain (
        BlockLattice2D<T,Lattice>& currentBlock,
        int x0_, int x1_, int y0_, int y1_)
{
    OLB_PRECONDITION(x0_ <= x1_);
    OLB_PRECONDITION(x1_ < partner->getNx());
    OLB_PRECONDITION(x1_ < currentBlock.getNx());
    OLB_PRECONDITION(y0_ <= y1_);
    OLB_PRECONDITION(y1_ < partner->getNy());
    OLB_PRECONDITION(y1_ < currentBlock.getNy());

    for (int iX=x0_; iX<=x1_; ++iX)
    {
        for (int iY=y0_; iY<=y1_; ++iY)
        {
            T u[Lattice<T>::d];
            partner -> get(iX,iY).computeU(u);
            T* externalVelocity = currentBlock.get(iX,iY).getExternal (
                                      Lattice<T>::ExternalField::velocityBeginsAt );
            for (int iD=0; iD<Lattice<T>::d; ++iD)
            {
                *(externalVelocity+iD) = u[iD];
            }
            
            T ux
        }
    }
}

template<typename T, template<typename U> class Lattice>
int ViscoelasticCoupler2D<T,Lattice>::extent() const {
    return 1; 
}

template<typename T, template<typename U> class Lattice>
int ViscoelasticCoupler2D<T,Lattice>::extent(int direction) const {
    return 1; 
}

/// Solvent coupling with other lattices

template<typename T, template<typename U> class Lattice>
SolventCoupler2D<T,Lattice>::SolventCoupler2D(std::vector<BlockStructure2D<T,Lattice>*> partners_) :
        partners(partners_)
{
}

template<typename T, template<typename U> class Lattice>
void SolventCoupler2D<T,Lattice>::process(BlockLattice2D<T,Lattice>& blockLattice) {
    processSubDomain(blockLattice, 0, blockLattice.getNx()-1, 0, blockLattice.getNy()-1);
}

template<typename T, template<typename U> class Lattice>
void SolventCoupler2D<T,Lattice>::processSubDomain (
        BlockLattice2D<T,Lattice>& currentBlock,
        int x0_, int x1_, int y0_, int y1_)
{
//     OLB_PRECONDITION(x0_ <= x1_);
//     OLB_PRECONDITION(x1_ < partner->getNx());
//     OLB_PRECONDITION(x1_ < currentBlock.getNx());
//     OLB_PRECONDITION(y0_ <= y1_);
//     OLB_PRECONDITION(y1_ < partner->getNy());
//     OLB_PRECONDITION(y1_ < currentBlock.getNy());

    for (int iX=x0_; iX<=x1_; ++iX) 
    {
        for (int iY=y0_; iY<=y1_; ++iY) 
        {
            T Q[3];
            Q[xx] = partner[xx] -> get(iX,iY).computeRho()+(T)1;
            Q[xy] = partner[xy] -> get(iX,iY).computeRho()+(T)1;
            Q[yy] = partner[yy] -> get(iX,iY).computeRho()+(T)1;
            
            T* externalConformation = currentBlock.get(iX,iY).getExternal (
                    Lattice<T>::ExternalField::conformationBeginsAt );

            for (int iD=0; iD<3; ++iD) {
                *(externalConformation+iD) = Q[iD];
            }
        }
    }
}

template<typename T, template<typename U> class Lattice>
int SolventCoupler2D<T,Lattice>::extent() const {
    return 0; 
}

template<typename T, template<typename U> class Lattice>
int SolventCoupler2D<T,Lattice>::extent(int direction) const {
    return 0; 
}


/// LatticeCouplingGenerator for viscoelastic coupling

template<typename T, template<typename U> class Lattice>
ViscoelasticCouplerGenerator2D<T,Lattice>::ViscoelasticCouplerGenerator2D(int x0_, int x1_, int y0_, int y1_)
    : LatticeCouplingGenerator2D<T,Lattice>(x0_, x1_, y0_, y1_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>* ViscoelasticCouplerGenerator2D<T,Lattice>::generate (
                                std::vector<BlockStructure2D<T,Lattice>*> partners) const
{
    return new ViscoelasticCoupler2D<T,Lattice>(partners);
}

template<typename T, template<typename U> class Lattice>
LatticeCouplingGenerator2D<T,Lattice>* ViscoelasticCouplerGenerator2D<T,Lattice>::clone() const {
    return new ViscoelasticCouplerGenerator2D<T,Lattice>(*this);
}




}  // namespace olb

#endif

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
#include "core/ompManager.h"
#include "core/postProcessing.h"

namespace olb {

/////////////////// Forward Declarations /////////////////////////////

template<typename T, template<typename U> class Lattice>
class BlockLattice2D;

// template<typename T, template<typename U> class Lattice>
// class BlockLattice3D;


/////////////////// 2D Postprocessing ///////////////////////////////

//=================================================================//
//=========================Visco Post-proc=========================//
//=================================================================//

/// Viscoelastic coupler post processor
template<typename T, template<typename U> class Lattice>
class ViscoelasticCoupler2D : public NblockLatticePostProcessor2D<T,Lattice>
{
public:
    ViscoelasticCoupler2D(std::vector<BlockStructure2D<T,Lattice>*> partners_);
    virtual void process(BlockLattice2D<T,Lattice>& blockLattice);
    virtual void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_);
    virtual int extent() const;
    virtual int extent(int direction) const;
private:
    std::vector<BlockStructure2D<T,Lattice>* > partners;
};

/// Viscoelastic Lattice Coupling Generator
template<typename T, template<typename U> class Lattice>
class ViscoelasticCouplerGenerator2D : public LatticeCouplingGenerator2D<T,Lattice>
{
public:
    ViscoelasticCouplerGenerator2D(int x0_, int x1_, int y0_, int y1_);
    virtual PostProcessor2D<T,Lattice>* generate(std::vector<BlockStructure2D<T,Lattice>*> partners) const;
    virtual LatticeCouplingGenerator2D<T,Lattice>* clone() const;
};


//=================================================================//
//=========================Solvent Post-proc=========================//
//=================================================================//
/// Solvent coupler post processor
template<typename T, template<typename U> class Lattice>
class SolventCoupler2D : public NblockLatticePostProcessor2D<T,Lattice>
{
public:
    SolventCoupler2D(std::vector<BlockStructure2D<T,Lattice>*> partners_);
    virtual void process(BlockLattice2D<T,Lattice>& blockLattice);
    virtual void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_);
    virtual int extent() const;
    virtual int extent(int direction) const;
private:
    std::vector<BlockStructure2D<T,Lattice>* > partners;
};

/// Solvent Lattice Coupling Generator
template<typename T, template<typename U> class Lattice>
class SolventCouplerGenerator2D : public LatticeCouplingGenerator2D<T,Lattice>
{
public:
    SolventCouplerGenerator2D(int x0_, int x1_, int y0_, int y1_);
    virtual PostProcessor2D<T,Lattice>* generate(std::vector<BlockStructure2D<T,Lattice>*> partners) const;
    virtual LatticeCouplingGenerator2D<T,Lattice>* clone() const;
};


}  // namespace olb

#endif

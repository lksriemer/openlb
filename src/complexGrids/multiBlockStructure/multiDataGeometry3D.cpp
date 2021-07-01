/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt and Bernd Stahl
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
 * Geometry specifications for 3D multiblocks -- implementation.
 */

#include "multiDataGeometry3D.h"
#include "core/olbDebug.h"

namespace olb {

    

////////////////////// Class BlockParameters3D /////////////////////////////

BlockParameters3D::BlockParameters3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                     int envelopeWidth_, int procId_,
                                     bool leftX, bool rightX, bool leftY, bool rightY, bool leftZ, bool rightZ)
    : envelopeWidth(envelopeWidth_),
      procId(procId_),
      bulk(x0_, x1_, y0_, y1_, z0_, z1_),
      envelope(x0_-envelopeWidth, x1_+envelopeWidth, y0_-envelopeWidth, y1_+envelopeWidth, z0_-envelopeWidth, z1_+envelopeWidth),
      nonPeriodicEnvelope(envelope)
{
    if (leftX) {
        nonPeriodicEnvelope.x0 += envelopeWidth;
    }
    if (rightX) {
        nonPeriodicEnvelope.x1 -= envelopeWidth;
    }

    if (leftY) {
        nonPeriodicEnvelope.y0 += envelopeWidth;
    }
    if (rightY) {
        nonPeriodicEnvelope.y1 -= envelopeWidth;
    }
    if (leftZ) {
        nonPeriodicEnvelope.z0 += envelopeWidth;
    }
    if (rightZ) {
        nonPeriodicEnvelope.z1 -= envelopeWidth;
    }
}


////////////////////// Class MultiDataDistribution3D /////////////////////

MultiDataDistribution3D::MultiDataDistribution3D(int nx_, int ny_, int nz_)
    : nx(nx_), ny(ny_), nz(nz_)
{ }

MultiDataDistribution3D& MultiDataDistribution3D::operator=(MultiDataDistribution3D const& rhs) {
    nx = rhs.nx;
    ny = rhs.ny;
    nz = rhs.nz;
    blocks = rhs.blocks;
    normalOverlaps = rhs.normalOverlaps;
    periodicOverlaps = rhs.periodicOverlaps;
    return (*this);
}


BlockParameters3D const& MultiDataDistribution3D::getBlockParameters(int whichBlock) const {
    OLB_PRECONDITION( whichBlock < getNumBlocks() );
    return blocks[whichBlock];
}
Overlap3D const& MultiDataDistribution3D::getNormalOverlap(int whichOverlap) const {
    OLB_PRECONDITION( whichOverlap < getNumNormalOverlaps() );
    return normalOverlaps[whichOverlap];
}
Overlap3D const& MultiDataDistribution3D::getPeriodicOverlap(int whichOverlap) const {
    OLB_PRECONDITION( whichOverlap < getNumPeriodicOverlaps() );
    return periodicOverlaps[whichOverlap];
}

void MultiDataDistribution3D::addBlock(int x0, int x1, int y0, int y1, int z0, int z1, int envelopeWidth, int procId) {
    OLB_PRECONDITION( x0>=0 && y0>=0 && z0>=0);
    OLB_PRECONDITION( x1<nx && y1<ny && z1<nz);
    OLB_PRECONDITION( x0 <= x1 && y0 <= y1 && z0 <= z1);

    BlockParameters3D newBlock(x0, x1, y0, y1, z0, z1, envelopeWidth, procId,
                               x0==0, x1==nx-1, y0==0, y1==ny-1, z0==0, z1==nz-1);
    computeNormalOverlaps(newBlock);
    blocks.push_back(newBlock);
    computePeriodicOverlaps();
}

int MultiDataDistribution3D::locate(int x, int y, int z, int guess) const {
    OLB_PRECONDITION( x>=0 && x < nx );
    OLB_PRECONDITION( y>=0 && y < ny );
    OLB_PRECONDITION( z>=0 && z < nz );
    OLB_PRECONDITION( guess < getNumBlocks() );

    // Check first whether (x,y,z) is contained in the guessed block
    BlockCoordinates3D const& guessCoord = blocks[guess].getBulk();
    if (util::contained(x, y, z, guessCoord.x0, guessCoord.x1, guessCoord.y0, guessCoord.y1, guessCoord.z0, guessCoord.z1)) {
        return guess;
    }
    // If not, try all the blocks
    for (int iBlock=0; iBlock<(int)blocks.size(); ++iBlock) {
        if (iBlock != guess) {
            BlockCoordinates3D const& coord = blocks[iBlock].getBulk();
            if (util::contained(x, y, z, coord.x0, coord.x1, coord.y0, coord.y1, coord.z0, coord.z1)) {
                return iBlock;
            }
        }
    }
    return -1;
}
//// MODIF COMPLETE UNTIL HERE
void MultiDataDistribution3D::computeNormalOverlaps(BlockParameters3D const& newBlock) {
    BlockCoordinates3D intersection;
    int iNew = getNumBlocks();
    for (int iBlock=0; iBlock<getNumBlocks(); ++iBlock) {
        if (util::intersect(blocks[iBlock].getBulk(), newBlock.getNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap3D(iBlock, iNew, intersection));
        }
        if (util::intersect(newBlock.getBulk(), blocks[iBlock].getNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap3D(iNew, iBlock, intersection));
        }
    }
}

void MultiDataDistribution3D::computePeriodicOverlaps() {
    int iNew = getNumBlocks()-1;
    BlockParameters3D const& newBlock = blocks[iNew];
    BlockCoordinates3D intersection;
    for (int dx=-1; dx<=+1; dx+=1) {
        for (int dy=-1; dy<=+1; dy+=1) {
            for (int dz=-1; dz<=+1; dz+=1) {
                if (dx!=0 || dy!=0 || dz!=0) {
                    int shiftX = dx*getNx();
                    int shiftY = dy*getNy();
                    int shiftZ = dz*getNz();
                    BlockCoordinates3D newBulk(newBlock.getBulk().shift(shiftX,shiftY,shiftZ));
                    BlockCoordinates3D newEnvelope(newBlock.getEnvelope().shift(shiftX,shiftY,shiftZ));
                    for (int iBlock=0; iBlock<getNumBlocks(); ++iBlock) {
                        if (util::intersect(blocks[iBlock].getBulk(), newEnvelope, intersection)) {
                            periodicOverlaps.push_back( Overlap3D(iBlock, iNew, intersection, shiftX, shiftY, shiftZ) );
                        }
                        if (!(iBlock==iNew) &&
                            util::intersect(newBulk, blocks[iBlock].getEnvelope(), intersection))
                        {
                            intersection = intersection.shift(-shiftX,-shiftY, -shiftZ);
                            periodicOverlaps.push_back( Overlap3D(iNew, iBlock, intersection, -shiftX, -shiftY, -shiftZ) );
                        }
                    }
                }
            }
        }
    }
}

int MultiDataDistribution3D::getNumAllocatedBulkCells() const {
    int numCells = 0;
    for (unsigned iBlock=0; iBlock<blocks.size(); ++iBlock) {
        numCells += blocks[iBlock].getBulkLx() * blocks[iBlock].getBulkLy() * blocks[iBlock].getBulkLz();
    }
    return numCells;
}

bool MultiDataDistribution3D::getNextChunkX(int iX, int iY, int iZ, int& nextLattice, int& nextChunkSize) const {
    nextLattice = locate(iX,iY,iZ);
    if (nextLattice == -1) {
        int exploreX = iX+1;
        while(exploreX<getNx() && locate(exploreX,iY,iZ)==-1) {
            ++exploreX;
        }
        nextChunkSize = exploreX-iX;
        return false;
    }
    else {
        nextChunkSize = blocks[nextLattice].getBulk().x1-iX+1;
        return true;
    }
}

bool MultiDataDistribution3D::getNextChunkY(int iX, int iY, int iZ, int& nextLattice, int& nextChunkSize) const {
    nextLattice = locate(iX,iY,iZ);
    if (nextLattice == -1) {
        int exploreY = iY+1;
        while(exploreY<getNy() && locate(iX,exploreY,iZ)==-1) {
            ++exploreY;
        }
        nextChunkSize = exploreY-iY;
        return false;
    }
    else {
        nextChunkSize = blocks[nextLattice].getBulk().y1-iY+1;
        return true;
    }
}

bool MultiDataDistribution3D::getNextChunkZ(int iX, int iY, int iZ, int& nextLattice, int& nextChunkSize) const {
    nextLattice = locate(iX,iY,iZ);
    if (nextLattice == -1) {
        int exploreZ = iZ+1;
        while(exploreZ<getNz() && locate(iX,iY,exploreZ)==-1) {
            ++exploreZ;
        }
        nextChunkSize = exploreZ-iZ;
        return false;
    }
    else {
        nextChunkSize = blocks[nextLattice].getBulk().z1-iZ+1;
        return true;
    }
}

}  // namespace olb

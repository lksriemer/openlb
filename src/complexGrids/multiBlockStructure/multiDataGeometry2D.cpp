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
 * Geometry specifications for 2D multiblocks -- implementation.
 */

#include "multiDataGeometry2D.h"
#include "core/olbDebug.h"


namespace olb {

////////////////////// Class BlockParameters2D /////////////////////////////

BlockParameters2D::BlockParameters2D(int x0_, int x1_, int y0_, int y1_,
                                     int envelopeWidth_, int procId_,
                                     bool leftX, bool rightX, bool leftY, bool rightY)
    : envelopeWidth(envelopeWidth_),
      procId(procId_),
      bulk(x0_, x1_, y0_, y1_),
      envelope(x0_-envelopeWidth, x1_+envelopeWidth, y0_-envelopeWidth, y1_+envelopeWidth),
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
}


////////////////////// Class MultiDataDistribution2D /////////////////////

MultiDataDistribution2D::MultiDataDistribution2D(int nx_, int ny_)
    : nx(nx_), ny(ny_)
{ }

MultiDataDistribution2D& MultiDataDistribution2D::operator=(MultiDataDistribution2D const& rhs) {
    nx = rhs.nx;
    ny = rhs.ny;
    blocks = rhs.blocks;
    normalOverlaps = rhs.normalOverlaps;
    periodicOverlaps = rhs.periodicOverlaps;
    return *this;
}


BlockParameters2D const& MultiDataDistribution2D::getBlockParameters(int whichBlock) const {
    OLB_PRECONDITION( whichBlock < getNumBlocks() );
    return blocks[whichBlock];
}
Overlap2D const& MultiDataDistribution2D::getNormalOverlap(int whichOverlap) const {
    OLB_PRECONDITION( whichOverlap < getNumNormalOverlaps() );
    return normalOverlaps[whichOverlap];
}
Overlap2D const& MultiDataDistribution2D::getPeriodicOverlap(int whichOverlap) const {
    OLB_PRECONDITION( whichOverlap < getNumPeriodicOverlaps() );
    return periodicOverlaps[whichOverlap];
}

void MultiDataDistribution2D::addBlock(int x0, int x1, int y0, int y1, int envelopeWidth, int procId) {
    OLB_PRECONDITION( x0>=0 && y0>=0 );
    OLB_PRECONDITION( x1<nx && y1<ny );
    OLB_PRECONDITION( x0 <= x1 && y0 <= y1 );

    BlockParameters2D newBlock(x0, x1, y0, y1, envelopeWidth, procId,
                               x0==0, x1==nx-1, y0==0, y1==ny-1);
    computeNormalOverlaps(newBlock);
    blocks.push_back(newBlock);
    computePeriodicOverlaps();
}

int MultiDataDistribution2D::locate(int x, int y, int guess) const {
    OLB_PRECONDITION( x>=0 && x < nx );
    OLB_PRECONDITION( y>=0 && y < ny );
    OLB_PRECONDITION( guess < getNumBlocks() );

    // Check first whether (x,y) is contained in the guessed block
    BlockCoordinates2D const& guessCoord = blocks[guess].getBulk();
    if (util::contained(x, y, guessCoord.x0, guessCoord.x1, guessCoord.y0, guessCoord.y1)) {
        return guess;
    }
    // If not, try all the blocks
    for (int iBlock=0; iBlock<(int)blocks.size(); ++iBlock) {
        if (iBlock != guess) {
            BlockCoordinates2D const& coord = blocks[iBlock].getBulk();
            if (util::contained(x, y, coord.x0, coord.x1, coord.y0, coord.y1)) {
                return iBlock;
            }
        }
    }
    return -1;
}

void MultiDataDistribution2D::computeNormalOverlaps(BlockParameters2D const& newBlock) {
    BlockCoordinates2D intersection;
    int iNew = getNumBlocks();
    for (int iBlock=0; iBlock<getNumBlocks(); ++iBlock) {
        if (util::intersect(blocks[iBlock].getBulk(), newBlock.getNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap2D(iBlock, iNew, intersection));
        }
        if (util::intersect(newBlock.getBulk(), blocks[iBlock].getNonPeriodicEnvelope(), intersection)) {
            normalOverlaps.push_back(Overlap2D(iNew, iBlock, intersection));
        }
    }
}

void MultiDataDistribution2D::computePeriodicOverlaps() {
    int iNew = getNumBlocks()-1;
    BlockParameters2D const& newBlock = blocks[iNew];
    BlockCoordinates2D intersection;
    for (int dx=-1; dx<=+1; dx+=1) {
        for (int dy=-1; dy<=+1; dy+=1) {
            if (dx!=0 || dy!=0) {
                int shiftX = dx*getNx();
                int shiftY = dy*getNy();
                BlockCoordinates2D newBulk(newBlock.getBulk().shift(shiftX,shiftY));
                BlockCoordinates2D newEnvelope(newBlock.getEnvelope().shift(shiftX,shiftY));
                for (int iBlock=0; iBlock<getNumBlocks(); ++iBlock) {
                    if (util::intersect(blocks[iBlock].getBulk(), newEnvelope, intersection)) {
                        periodicOverlaps.push_back( Overlap2D(iBlock, iNew, intersection, shiftX, shiftY) );
                    }
                    if (!(iBlock==iNew) &&
                        util::intersect(newBulk, blocks[iBlock].getEnvelope(), intersection))
                    {
                        intersection = intersection.shift(-shiftX,-shiftY);
                        periodicOverlaps.push_back( Overlap2D(iNew, iBlock, intersection, -shiftX, -shiftY) );
                    }
                }
            }
        }
    }
}

int MultiDataDistribution2D::getNumAllocatedBulkCells() const {
    int numCells = 0;
    for (unsigned iBlock=0; iBlock<blocks.size(); ++iBlock) {
        numCells += blocks[iBlock].getBulkLx() * blocks[iBlock].getBulkLy();
    }
    return numCells;
}

bool MultiDataDistribution2D::getNextChunkX(int iX, int iY, int& nextLattice, int& nextChunkSize) const {
    nextLattice = locate(iX,iY);
    if (nextLattice == -1) {
        int exploreX = iX+1;
        while(exploreX<getNx() && locate(exploreX,iY)==-1) {
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

bool MultiDataDistribution2D::getNextChunkY(int iX, int iY, int& nextLattice, int& nextChunkSize) const {
    nextLattice = locate(iX,iY);
    if (nextLattice == -1) {
        int exploreY = iY+1;
        while(exploreY<getNy() && locate(iX,exploreY)==-1) {
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

}  // namespace olb

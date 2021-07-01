/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias Krause
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


#ifdef PARALLEL_MODE_OMP

    #include "loadBalancer.h"
    #include "olbDebug.h"

    loadBalancer::loadBalancer(int rank, int size, int globChunkSize, int offset) { 

        OLB_PRECONDITION(rank>=0 && size>1 && offset>=0)
        OLB_PRECONDITION(size<=globChunkSize &&  rank<size);

        locChunkSize = (globChunkSize+size-rank-1)/size;
        if (rank+1 <= globChunkSize-(globChunkSize/size)*size) {
            firstGlobNum = globChunkSize/size * rank + rank + offset;
            lastGlobNum  = firstGlobNum + locChunkSize - 1;
        }
        else {
            firstGlobNum = globChunkSize/size * rank + globChunkSize - (globChunkSize/size)*size + offset;
            lastGlobNum  = firstGlobNum + locChunkSize - 1;
        }
    }

    int loadBalancer::get_locChunkSize() const {
        return locChunkSize;
    }

    int loadBalancer::get_firstGlobNum() const {
        return firstGlobNum;
    }

    int loadBalancer::get_lastGlobNum() const {
        return lastGlobNum;
    }

#endif

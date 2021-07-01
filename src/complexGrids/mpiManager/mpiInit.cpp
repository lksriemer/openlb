/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 The OpenLB project
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
 * Wrapper functions that simplify the use of MPI, template instatiations
 */

#include "mpiManager.h"
#include "io/parallelIO.h"

namespace olb {

namespace singleton {
    MpiManager& mpi() {
        static MpiManager instance;
        return instance;
}

}

void olbInit(int *argc, char ***argv, bool verbous) {
    singleton::mpi().init(argc, argv, verbous);
#ifdef PARALLEL_MODE_MPI
    ParBuf *newCoutBuf = new ParBuf(std::cout.rdbuf());
    ParBuf *newCerrBuf = new ParBuf(std::cerr.rdbuf());
    ParBuf *newClogBuf = new ParBuf(std::clog.rdbuf());
    ParBuf *newCinBuf  = new ParBuf(std::cin.rdbuf());

    std::cout.rdbuf(newCoutBuf);
    std::cerr.rdbuf(newCerrBuf);
    std::clog.rdbuf(newClogBuf);
    std::cin. rdbuf(newCinBuf);
#endif
}

}


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
 * Groups all the include files for basic 2D dynamics.
 */
#include "singleton.h"
#include "units.h"
#include "latticeDescriptors.h"
#include "dynamics.h"
#include "cell.h"
#include "boundaries.h"
#include "boundaries2D.h"
#include "postProcessing.h"
#include "dataFields2D.h"
#include "blockStructure2D.h"
#include "fdBoundaries2D.h"
#include "blockLattice2D.h"
#include "blockLatticeView2D.h"
#include "blockStatistics2D.h"
#include "dynamicsHelpers2D.h"
#include "simulationSetup2D.h"
#include "imageCreator.h"
#include "olbAddOn2D.h"

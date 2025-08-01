/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
 *                2024 Dennis Teutscher
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
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
 * Groups all the generic 3D template files in the directory utilities.
 */

#include "benchmarkUtil.hh"
#include "timer.hh"
#include "functorPtr.hh"
#include "functorDsl3D.hh"
#include "hyperplane3D.hh"
#include "hyperplaneLattice3D.hh"
#include "line3D.hh"
#include "lineLattice3D.hh"
#include "random.hh"
#ifdef FEATURE_PROJ
#include "osmParser.hh"
#endif

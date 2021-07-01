/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias Krause
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

#ifndef BLOCK_GIF_WRITER_3D_H
#define BLOCK_GIF_WRITER_3D_H

//#include <sstream>
#include <iomanip>
#include <vector>

#include "colormaps.h"
#include "io/ostreamManager.h"
#include "functors/blockLatticeBaseF3D.h"


namespace olb {

// to start with, there are only XY,XZ,YZ plane implemented
// thus origin is not necessary
//
// normal with (1,0,0) is a YZ plane
// vector in plane is given by s*vecU +t*vecV

template<typename T, template <typename U> class DESCRIPTOR>
class BlockGifWriter3D {
public:
  BlockGifWriter3D(std::string const& map);
  void write( BlockLatticeF3D<T,DESCRIPTOR>& f, std::vector<int> normal );

private:
//  // get u which spans the plane
//  std::vector<int> vecU( std::vector<int> normal );
//  // get v which spans the plane
//  std::vector<int> vecV( std::vector<int> normal );
  // get origin O, set _imageOriginX and _imageX
  void setFrame( BlockLatticeF3D<T,DESCRIPTOR>& f, std::vector<int> normal );
  //   
  void writeGif( BlockLatticeF3D<T,DESCRIPTOR>& f, std::vector<int> normal );

  int _extNx;
  int _extNy;
  int _extNz;

  mutable OstreamManager clout;
  int _colorRange;
  int _numColors;
  graphics::ColorMap<T> _colorMap;
};

}  // namespace olb

#endif

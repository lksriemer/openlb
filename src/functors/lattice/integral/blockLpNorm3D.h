/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef BLOCK_LP_NORM_3D_H
#define BLOCK_LP_NORM_3D_H

namespace olb {

template<typename T> class BlockF3D;
template<typename T> class BlockIndicatorF3D;

/// Block level functor that returns the Lp norm over omega of the euklid norm of the input block functor.
template <typename T, typename W, int P>
class BlockLpNorm3D final : public BlockF3D<W> {
protected:
  BlockF3D<W>&          _f;
  BlockIndicatorF3D<T>& _indicatorF;
public:
  /**
   * \param f          data functor
   * \param indicatorF indicator functor describing the subset to be integrated
   **/
  BlockLpNorm3D(BlockF3D<W>& f, BlockIndicatorF3D<T>& indicatorF);
  bool operator() (W output[], const int input[]) override;
};

}

#endif

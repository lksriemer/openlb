/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Albert Mink, Lukas Baron, Mathias J. Krause
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
#ifndef BLOCK_LATTICE_CALC_F_2D_HH
#define BLOCK_LATTICE_CALC_F_2D_HH

#include<vector>
#include "functors/blockLatticeCalcF2D.h"
#include "functors/genericF.h"


namespace olb {

////////////////////////////// BlockLatticeCalc2D //////////////////////////////

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeCalc2D<T,DESCRIPTOR>::BlockLatticeCalc2D
  (BlockLatticeF2D<T,DESCRIPTOR>& f, BlockLatticeF2D<T,DESCRIPTOR>& g)
  : BlockLatticeF2D<T,DESCRIPTOR>( f.getBlockLattice2D(), f.getTargetDim() ),
    _f(f), _g(g) { }

template <typename T, template <typename U> class DESCRIPTOR>
void BlockLatticeCalc2D<T,DESCRIPTOR>::myErase(GenericF<T,int>* ptr)
{
  this->removeChild(ptr);
  delete ptr;
  // additional delete recursively this from father, father from grandfather for all calc classes
  if (this->_pointerVec.size() == 0)  this->_f.myErase(this);
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePlus2D<T,DESCRIPTOR>::BlockLatticePlus2D
  (BlockLatticeF2D<T,DESCRIPTOR>& f, BlockLatticeF2D<T,DESCRIPTOR>& g)
  : BlockLatticeCalc2D<T,DESCRIPTOR>(f,g)
{
  this->_name = "(" + f.getName() + "+" + g.getName() + ")";
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePlus2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] + this->_g(input)[i] );
  }
  // 'NULL' initiates the recursion for delete
  this->myErase(NULL);
  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeMinus2D<T,DESCRIPTOR>::BlockLatticeMinus2D
  (BlockLatticeF2D<T,DESCRIPTOR>& f, BlockLatticeF2D<T,DESCRIPTOR>& g)
  : BlockLatticeCalc2D<T,DESCRIPTOR>(f,g)
{
  this->_name = "(" + f.getName() + "-" + g.getName() + ")";
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeMinus2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] - this->_g(input)[i] );
  }
  // 'NULL' initiates the recursion for delete
  this->myErase(NULL);
  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeMultiplication2D<T,DESCRIPTOR>::BlockLatticeMultiplication2D
  (BlockLatticeF2D<T,DESCRIPTOR>& f, BlockLatticeF2D<T,DESCRIPTOR>& g)
  : BlockLatticeCalc2D<T,DESCRIPTOR>(f,g)
{
  this->_name = f.getName() + "*" + g.getName();
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeMultiplication2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] * this->_g(input)[i] );
  }
  // 'NULL' initiates the recursion for delete
  this->myErase(NULL);
  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeDivision2D<T,DESCRIPTOR>::BlockLatticeDivision2D
  (BlockLatticeF2D<T,DESCRIPTOR>& f, BlockLatticeF2D<T,DESCRIPTOR>& g)
  : BlockLatticeCalc2D<T,DESCRIPTOR>(f,g)
{
  this->_name = f.getName() + "/" + g.getName();
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeDivision2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] / this->_g(input)[i] );
  }
  // 'NULL' initiates the recursion for delete
  this->myErase(NULL);
  return output;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>& BlockLatticeF2D<T,DESCRIPTOR>::operator+(BlockLatticeF2D<T,DESCRIPTOR>& rhs)
{
  BlockLatticeF2D<T,DESCRIPTOR>* tmp = new BlockLatticePlus2D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>& BlockLatticeF2D<T,DESCRIPTOR>::operator-(BlockLatticeF2D<T,DESCRIPTOR>& rhs)
{
  BlockLatticeF2D<T,DESCRIPTOR>* tmp = new BlockLatticeMinus2D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>& BlockLatticeF2D<T,DESCRIPTOR>::operator*(BlockLatticeF2D<T,DESCRIPTOR>& rhs)
{
  BlockLatticeF2D<T,DESCRIPTOR>* tmp = new BlockLatticeMultiplication2D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>& BlockLatticeF2D<T,DESCRIPTOR>::operator/(BlockLatticeF2D<T,DESCRIPTOR>& rhs)
{
  BlockLatticeF2D<T,DESCRIPTOR>* tmp = new BlockLatticeDivision2D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}



} // end namespace olb

#endif

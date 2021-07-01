/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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
#ifndef BLOCK_LATTICE_CALC_F_3D_HH
#define BLOCK_LATTICE_CALC_F_3D_HH

#include<vector>
#include "functors/genericF.h"
#include "functors/blockLatticeCalcF3D.h"


namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeCalc3D<T,DESCRIPTOR>::BlockLatticeCalc3D
  (BlockLatticeF3D<T,DESCRIPTOR>& f, BlockLatticeF3D<T,DESCRIPTOR>& g)
  : BlockLatticeF3D<T,DESCRIPTOR>( f.getBlockLattice3D(), f.getTargetDim() ),
    _f(f), _g(g) { }

template <typename T, template <typename U> class DESCRIPTOR>
void BlockLatticeCalc3D<T,DESCRIPTOR>::myErase(GenericF<T,int>* ptr)
{
  this->removeChild(ptr);
  delete ptr;
  // additional delete recursively this from father, father from grandfather for all calc classes
  if (this->_pointerVec.size() == 0)  this->_f.myErase(this);
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePlus3D<T,DESCRIPTOR>::BlockLatticePlus3D
  (BlockLatticeF3D<T,DESCRIPTOR>& f, BlockLatticeF3D<T,DESCRIPTOR>& g)
  : BlockLatticeCalc3D<T,DESCRIPTOR>(f,g)
{
  this->_name = "(" + f.getName() + "+" + g.getName() + ")";
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticePlus3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
BlockLatticeMinus3D<T,DESCRIPTOR>::BlockLatticeMinus3D
  (BlockLatticeF3D<T,DESCRIPTOR>& f, BlockLatticeF3D<T,DESCRIPTOR>& g)
  : BlockLatticeCalc3D<T,DESCRIPTOR>(f,g)
{
  this->_name = "(" + f.getName() + "-" + g.getName() + ")";
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeMinus3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
BlockLatticeMultiplication3D<T,DESCRIPTOR>::BlockLatticeMultiplication3D
  (BlockLatticeF3D<T,DESCRIPTOR>& f, BlockLatticeF3D<T,DESCRIPTOR>& g)
  : BlockLatticeCalc3D<T,DESCRIPTOR>(f,g)
{
  this->_name = f.getName() + "*" + g.getName();
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeMultiplication3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
BlockLatticeDivision3D<T,DESCRIPTOR>::BlockLatticeDivision3D
  (BlockLatticeF3D<T,DESCRIPTOR>& f, BlockLatticeF3D<T,DESCRIPTOR>& g)
  : BlockLatticeCalc3D<T,DESCRIPTOR>(f,g)
{
  this->_name = f.getName() + "/" + g.getName();
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeDivision3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
BlockLatticeF3D<T,DESCRIPTOR>& BlockLatticeF3D<T,DESCRIPTOR>::operator+(BlockLatticeF3D<T,DESCRIPTOR>& rhs)
{
  BlockLatticeF3D<T,DESCRIPTOR>* tmp = new BlockLatticePlus3D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>& BlockLatticeF3D<T,DESCRIPTOR>::operator-(BlockLatticeF3D<T,DESCRIPTOR>& rhs)
{
  BlockLatticeF3D<T,DESCRIPTOR>* tmp = new BlockLatticeMinus3D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>& BlockLatticeF3D<T,DESCRIPTOR>::operator*(BlockLatticeF3D<T,DESCRIPTOR>& rhs)
{
  BlockLatticeF3D<T,DESCRIPTOR>* tmp = new BlockLatticeMultiplication3D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeF3D<T,DESCRIPTOR>& BlockLatticeF3D<T,DESCRIPTOR>::operator/(BlockLatticeF3D<T,DESCRIPTOR>& rhs)
{
  BlockLatticeF3D<T,DESCRIPTOR>* tmp = new BlockLatticeDivision3D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}



} // end namespace olb

#endif

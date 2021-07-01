/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#ifndef SUPER_LATTICE_CALC_F_2D_HH
#define SUPER_LATTICE_CALC_F_2D_HH

#include<vector>
#include "functors/superLatticeCalcF2D.h"
#include "functors/genericF.h"




namespace olb {


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeCalc2D<T,DESCRIPTOR>::SuperLatticeCalc2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f, SuperLatticeF2D<T,DESCRIPTOR>& g)
  : SuperLatticeF2D<T,DESCRIPTOR>( f.getSuperLattice2D(), f.getTargetDim() ),
    _f(f), _g(g) { }

template <typename T, template <typename U> class DESCRIPTOR>
void SuperLatticeCalc2D<T,DESCRIPTOR>::myErase(GenericF<T,int>* ptr) {
  this->removeChild(ptr);
  delete ptr;
  // additional delete recursively this from father, father from grandfather for all calc classes
  if (this->_pointerVec.size() == 0)  this->_f.myErase(this);
}


// addition
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePlus2D<T,DESCRIPTOR>::SuperLatticePlus2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f, SuperLatticeF2D<T,DESCRIPTOR>& g)
  : SuperLatticeCalc2D<T,DESCRIPTOR>(f,g)
{ this->_name = "(" + f.getName() + "+" + g.getName() + ")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePlus2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] + this->_g(input)[i] );
  }
  // 'NULL' initiates the recursion for delete
  this->myErase(NULL);
  return output;
}


// subtraction
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeMinus2D<T,DESCRIPTOR>::SuperLatticeMinus2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f, SuperLatticeF2D<T,DESCRIPTOR>& g)
  : SuperLatticeCalc2D<T,DESCRIPTOR>(f,g)
{ this->_name = "(" + f.getName() + "-" + g.getName() + ")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeMinus2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] - this->_g(input)[i] );
  }
  // 'NULL' initiates the recursion for delete
  this->myErase(NULL);
  return output;
}


// multiplication
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeMultiplication2D<T,DESCRIPTOR>::SuperLatticeMultiplication2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f, SuperLatticeF2D<T,DESCRIPTOR>& g)
  : SuperLatticeCalc2D<T,DESCRIPTOR>(f,g)
{ this->_name = f.getName() + "*" + g.getName(); }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeMultiplication2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
{
  std::vector<T> output;
  for(int i = 0; i < this->_f.getTargetDim(); i++) {
    output.push_back( this->_f(input)[i] * this->_g(input)[i] );
  }
  // 'NULL' initiates the recursion for delete
  this->myErase(NULL);
  return output;
}


// division
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeDivision2D<T,DESCRIPTOR>::SuperLatticeDivision2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f, SuperLatticeF2D<T,DESCRIPTOR>& g)
  : SuperLatticeCalc2D<T,DESCRIPTOR>(f,g)
{ this->_name = f.getName() + "/" + g.getName(); }


template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeDivision2D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
SuperLatticeF2D<T,DESCRIPTOR>& SuperLatticeF2D<T,DESCRIPTOR>::operator+(SuperLatticeF2D<T,DESCRIPTOR>& rhs)
{
  SuperLatticeF2D<T,DESCRIPTOR>* tmp = new SuperLatticePlus2D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF2D<T,DESCRIPTOR>& SuperLatticeF2D<T,DESCRIPTOR>::operator-(SuperLatticeF2D<T,DESCRIPTOR>& rhs)
{
  SuperLatticeF2D<T,DESCRIPTOR>* tmp = new SuperLatticeMinus2D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF2D<T,DESCRIPTOR>& SuperLatticeF2D<T,DESCRIPTOR>::operator*(SuperLatticeF2D<T,DESCRIPTOR>& rhs)
{
  SuperLatticeF2D<T,DESCRIPTOR>* tmp = new SuperLatticeMultiplication2D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF2D<T,DESCRIPTOR>& SuperLatticeF2D<T,DESCRIPTOR>::operator/(SuperLatticeF2D<T,DESCRIPTOR>& rhs)
{
  SuperLatticeF2D<T,DESCRIPTOR>* tmp = new SuperLatticeDivision2D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}



} // end namespace olb

#endif

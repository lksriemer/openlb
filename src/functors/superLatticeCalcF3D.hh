/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef SUPER_LATTICE_CALC_F_3D_HH
#define SUPER_LATTICE_CALC_F_3D_HH

#include<vector>
#include "functors/superLatticeCalcF3D.h"
#include "functors/genericF.h"


namespace olb {



template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeCalc3D<T,DESCRIPTOR>::SuperLatticeCalc3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, SuperLatticeF3D<T,DESCRIPTOR>& g)
  : SuperLatticeF3D<T,DESCRIPTOR>( f.getSuperLattice3D(), f.getTargetDim() ),
    _f(f), _g(g) { }

template <typename T, template <typename U> class DESCRIPTOR>
void SuperLatticeCalc3D<T,DESCRIPTOR>::myErase(GenericF<T,int>* ptr) {
  this->removeChild(ptr);
  delete ptr;
  // additional delete recursively this from father, father from grandfather for all calc classes
  if (this->_pointerVec.size() == 0)  this->_f.myErase(this);
}


// addition
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePlus3D<T,DESCRIPTOR>::SuperLatticePlus3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, SuperLatticeF3D<T,DESCRIPTOR>& g)
  : SuperLatticeCalc3D<T,DESCRIPTOR>(f,g)
{ this->_name = "(" + f.getName() + "+" + g.getName() + ")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticePlus3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
SuperLatticeMinus3D<T,DESCRIPTOR>::SuperLatticeMinus3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, SuperLatticeF3D<T,DESCRIPTOR>& g)
  : SuperLatticeCalc3D<T,DESCRIPTOR>(f,g)
{ this->_name = "(" + f.getName() + "-" + g.getName() + ")"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeMinus3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
SuperLatticeMultiplication3D<T,DESCRIPTOR>::SuperLatticeMultiplication3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, SuperLatticeF3D<T,DESCRIPTOR>& g)
  : SuperLatticeCalc3D<T,DESCRIPTOR>(f,g)
{ this->_name = f.getName() + "*" + g.getName(); }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeMultiplication3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
SuperLatticeDivision3D<T,DESCRIPTOR>::SuperLatticeDivision3D(
  SuperLatticeF3D<T,DESCRIPTOR>& f, SuperLatticeF3D<T,DESCRIPTOR>& g)
  : SuperLatticeCalc3D<T,DESCRIPTOR>(f,g)
{ this->_name = f.getName() + "/" + g.getName(); }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeDivision3D<T,DESCRIPTOR>::operator()(std::vector<int> input)
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
SuperLatticeF3D<T,DESCRIPTOR>& SuperLatticeF3D<T,DESCRIPTOR>::operator+(SuperLatticeF3D<T,DESCRIPTOR>& rhs)
{
  SuperLatticeF3D<T,DESCRIPTOR>* tmp = new SuperLatticePlus3D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF3D<T,DESCRIPTOR>& SuperLatticeF3D<T,DESCRIPTOR>::operator-(SuperLatticeF3D<T,DESCRIPTOR>& rhs)
{
  SuperLatticeF3D<T,DESCRIPTOR>* tmp = new SuperLatticeMinus3D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF3D<T,DESCRIPTOR>& SuperLatticeF3D<T,DESCRIPTOR>::operator*(SuperLatticeF3D<T,DESCRIPTOR>& rhs)
{
  SuperLatticeF3D<T,DESCRIPTOR>* tmp = new SuperLatticeMultiplication3D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeF3D<T,DESCRIPTOR>& SuperLatticeF3D<T,DESCRIPTOR>::operator/(SuperLatticeF3D<T,DESCRIPTOR>& rhs)
{
  SuperLatticeF3D<T,DESCRIPTOR>* tmp = new SuperLatticeDivision3D<T,DESCRIPTOR>(*this,rhs);
  this->addChild(tmp);
  return *tmp;
}



} // end namespace olb

#endif

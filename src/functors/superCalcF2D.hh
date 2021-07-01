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

#ifndef SUPER_CALC_F_2D_HH
#define SUPER_CALC_F_2D_HH


#include "functors/superCalcF2D.h"



namespace olb {


template <typename T>
SuperCalc2D<T>::SuperCalc2D(SuperF2D<T>& f, SuperF2D<T>& g)
  : SuperF2D<T>( f.getSuperStructure(), f.getTargetDim() ), _f(f), _g(g)
{
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

// addition
template <typename T>
SuperPlus2D<T>::SuperPlus2D(SuperF2D<T>& f, SuperF2D<T>& g) : SuperCalc2D<T>(f,g)
{
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";
}

template <typename T>
bool SuperPlus2D<T>::operator()(T output[], const int input[])
{
  this->_f(output,input);
  T tmp[this->_g.getTargetDim()];
  this->_g(tmp,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]+=tmp[i];
  }
  return true;
}


// subtraction
template <typename T>
SuperMinus2D<T>::SuperMinus2D(SuperF2D<T>& f, SuperF2D<T>& g) : SuperCalc2D<T>(f,g)
{
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
}

template <typename T>
bool SuperMinus2D<T>::operator()(T output[], const int input[])
{
  this->_f(output,input);
  T tmp[this->_g.getTargetDim()];
  this->_g(tmp,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]-=tmp[i];
  }
  return true;
}


// multiplication
template <typename T>
SuperMultiplication2D<T>::SuperMultiplication2D(SuperF2D<T>& f, SuperF2D<T>& g)
  : SuperCalc2D<T>(f,g)
{
  this->getName() = f.getName() + "*" + g.getName();
}

template <typename T>
bool SuperMultiplication2D<T>::operator()(T output[], const int input[])
{
  this->_f(output,input);
  T tmp[this->_g.getTargetDim()];
  this->_g(tmp,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]*=tmp[i];
  }
  return true;
}


// division
template <typename T>
SuperDivision2D<T>::SuperDivision2D(SuperF2D<T>& f, SuperF2D<T>& g)
  : SuperCalc2D<T>(f,g)
{
  this->getName() = f.getName() + "/" + g.getName();
}


template <typename T>
bool SuperDivision2D<T>::operator()(T output[], const int input[])
{
  this->_f(output,input);
  T tmp[this->_g.getTargetDim()];
  this->_g(tmp,input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i]/=tmp[i];
  }
  return true;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T>
SuperF2D<T>& SuperF2D<T>::operator+(SuperF2D<T>& rhs)
{
  auto tmp = std::make_shared< SuperPlus2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
SuperF2D<T>& SuperF2D<T>::operator-(SuperF2D<T>& rhs)
{
  auto tmp = std::make_shared< SuperMinus2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
SuperF2D<T>& SuperF2D<T>::operator*(SuperF2D<T>& rhs)
{
  auto tmp = std::make_shared< SuperMultiplication2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
SuperF2D<T>& SuperF2D<T>::operator/(SuperF2D<T>& rhs)
{
  auto tmp = std::make_shared< SuperDivision2D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}



} // end namespace olb

#endif

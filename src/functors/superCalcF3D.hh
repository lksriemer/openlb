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

#ifndef SUPER_CALC_F_3D_HH
#define SUPER_CALC_F_3D_HH


#include "functors/superCalcF3D.h"


namespace olb {


template <typename T>
SuperCalc3D<T>::SuperCalc3D(SuperF3D<T>& f, SuperF3D<T>& g)
  : SuperF3D<T>( f.getSuperStructure(), f.getTargetDim() ), _f(f), _g(g)
{
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

// addition
template <typename T>
SuperPlus3D<T>::SuperPlus3D(SuperF3D<T>& f, SuperF3D<T>& g) : SuperCalc3D<T>(f,g)
{
  this->getName() = "(" + f.getName() + "+" + g.getName() + ")";
}

template <typename T>
bool SuperPlus3D<T>::operator()(T output[], const int input[])
{
  this->_f(output,input);
  T tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] += tmp[i];
  }
  return true;
}

// subtraction
template <typename T>
SuperMinus3D<T>::SuperMinus3D(SuperF3D<T>& f, SuperF3D<T>& g) : SuperCalc3D<T>(f,g)
{
  this->getName() = "(" + f.getName() + "-" + g.getName() + ")";
}

template <typename T>
bool SuperMinus3D<T>::operator()(T output[], const int input[])
{
  this->_f(output,input);
  T tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] -= tmp[i];
  }
  return true;
}


// multiplication
template <typename T>
SuperMultiplication3D<T>::SuperMultiplication3D(SuperF3D<T>& f, SuperF3D<T>& g)
  : SuperCalc3D<T>(f,g)
{
  this->getName() = f.getName() + "*" + g.getName();
}

template <typename T>
bool SuperMultiplication3D<T>::operator()(T output[], const int input[])
{
  this->_f(output,input);
  T tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] *= tmp[i];
  }
  return true;
}


// division
template <typename T>
SuperDivision3D<T>::SuperDivision3D(SuperF3D<T>& f, SuperF3D<T>& g)
  : SuperCalc3D<T>(f,g)
{
  this->getName() = f.getName() + "/" + g.getName();
}

template <typename T>
bool SuperDivision3D<T>::operator()(T output[], const int input[])
{
  this->_f(output,input);
  T tmp[this->_g.getTargetDim()];
  this->_g(tmp, input);
  for (int i = 0; i < this->_f.getTargetDim(); i++) {
    output[i] /= tmp[i];
  }
  return true;
}


/////////////////////////////////operator()/// ////////////////////////////////
template <typename T>
SuperF3D<T>& SuperF3D<T>::operator+(SuperF3D<T>& rhs)
{
  auto tmp = std::make_shared< SuperPlus3D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
SuperF3D<T>& SuperF3D<T>::operator-(SuperF3D<T>& rhs)
{
  auto tmp = std::make_shared< SuperMinus3D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
SuperF3D<T>& SuperF3D<T>::operator*(SuperF3D<T>& rhs)
{
  auto tmp = std::make_shared< SuperMultiplication3D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T>
SuperF3D<T>& SuperF3D<T>::operator/(SuperF3D<T>& rhs)
{
  auto tmp = std::make_shared< SuperDivision3D<T> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}



} // end namespace olb

#endif

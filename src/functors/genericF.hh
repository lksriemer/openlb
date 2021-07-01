 /*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2013 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef GENERIC_F_HH
#define GENERIC_F_HH


/** \file
 * The description of a generic interface for all functor classes
 * -- generic implementation.
 */

#include"genericF.h"
#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm

namespace olb {

template <typename T, typename S>
int GenericF<T,S>::getSourceDim() const {
  return _m;
}

template <typename T, typename S>
int GenericF<T,S>::getTargetDim() const {
  return _n;
}

template <typename T, typename S>
std::string& GenericF<T,S>::getName()
{ return _name; }

template <typename T, typename S>
std::string const& GenericF<T,S>::getName() const
{ return _name; }

template <typename T, typename S>
std::vector<T> GenericF<T,S>::operator() () {
  std::vector<S> v;
  return operator()(v);
}

template <typename T, typename S>
std::vector<T> GenericF<T,S>::operator() (S input0) {
  std::vector<S> v;
  v.push_back(input0);
  return operator()(v);
}

template <typename T, typename S>
std::vector<T> GenericF<T,S>::operator() (S input0, S input1) {
  std::vector<S> v;
  v.push_back(input0);
  v.push_back(input1);
  return operator()(v);
}

template <typename T, typename S>
std::vector<T> GenericF<T,S>::operator() (S input0, S input1, S input2) {
  std::vector<S> v;
  v.push_back(input0);
  v.push_back(input1);
  v.push_back(input2);
  return operator()(v);
}

template <typename T, typename S>
std::vector<T> GenericF<T,S>::operator() (S input0, S input1, S input2, S input3) {
  std::vector<S> v;
  v.push_back(input0);
  v.push_back(input1);
  v.push_back(input2);
  v.push_back(input3);
  return operator()(v);
}

template <typename T, typename S>
void GenericF<T,S>::myErase(GenericF<T,S>* ptr) {
  // remove ptr from child list
  this->removeChild(ptr);
  // delete object
  delete ptr;
}

template <typename T, typename S>
void GenericF<T,S>::addChild(GenericF<T,S>* ptr) {
  _pointerVec.push_back(ptr);
}

template <typename T, typename S>
void GenericF<T,S>::removeChild(GenericF<T,S>* ptr) {
  for ( unsigned i = 0; i < _pointerVec.size(); i++ ) {
    if ( _pointerVec[i] == ptr ) {
      _pointerVec.erase(_pointerVec.begin() + i);
      // jump out of for-loop
      // break;
    }
  }
}


} // end namespace olb

#endif

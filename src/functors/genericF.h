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

/** \file
 * The description of a generic interface for all functor classes
 * -- header file.
 */


#ifndef GENERIC_F_H
#define GENERIC_F_H

#include<vector>
#include<string>

namespace olb {

/**
 *  GenericF is a base class, that can represent continuous as well as discrete
 *  functions.
 *  Please take care about the source and target dimensions in the constructor.
 *                F: S^m -> T^n (S=source, T=target)
 *
 *  \param _m     source dimension
 *  \param _n     target dimension
 *  \param _name  is functor name e.g. velocity, pressure
 */
template <typename T, typename S>
class GenericF {
protected:
  GenericF() {};
  GenericF(int targetDim, int sourceDim) : _n(targetDim), _m(sourceDim) { };

  int _n;
  int _m;
  std::vector< GenericF<T,S>* > _pointerVec;
  std::string _name;

public:
  virtual ~GenericF() {};
  int getSourceDim() const;
  int getTargetDim() const;

  /// read and write access to name
  std::string& getName();
  /// read only access to name
  std::string const& getName() const;

  virtual std::vector<T> operator() (std::vector<S> input) = 0;

  virtual std::vector<T> operator() ();
  virtual std::vector<T> operator() (S input0);
  virtual std::vector<T> operator() (S input0, S input1);
  virtual std::vector<T> operator() (S input0, S input1, S input2);
  virtual std::vector<T> operator() (S input0, S input1, S input2, S input3);

  // memory management
  /// adds ptr to the child list
  virtual void addChild(GenericF<T,S>* ptr);
  /// remove ptr from the child list
  virtual void removeChild(GenericF<T,S>* ptr);
  /// deletes ptr from child list and frees the memory for the object *ptr
  /// CALC CLASSES additionally delete recursively father from grandfather
  /// and free the memory (if possible)
  virtual void myErase(GenericF<T,S>* ptr);
};

} // end namespace olb

#endif

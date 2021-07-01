/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerländer
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

#ifndef FUNCTOR_PTR_HH
#define FUNCTOR_PTR_HH

#include "functorPtr.h"

namespace olb {

template<typename F>
FunctorPtr<F>::FunctorPtr(F& f):
  _ownF(),
  _f(&f)
{ }

template<typename F>
FunctorPtr<F>::FunctorPtr(F* f):
  _ownF(),
  _f(f)
{ }

template<typename F>
FunctorPtr<F>::FunctorPtr(std::unique_ptr<F>&& f):
  _ownF(f.release()),
  _f(_ownF.get())
{ }

template<typename F>
template<typename... Args>
bool FunctorPtr<F>::operator()(Args... args)
{
  return _f->operator()(std::forward<Args>(args)...);
}

template<typename F>
typename std::add_lvalue_reference<F>::type FunctorPtr<F>::operator*() const
{
  return *_f;
}

template<typename F>
typename std::add_pointer<F>::type FunctorPtr<F>::operator->() const noexcept
{
  return _f;
}

template<typename F>
FunctorPtr<F>::operator bool() const
{
  return _f != nullptr;
}

}

#endif

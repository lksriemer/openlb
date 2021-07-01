/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias Krause, Peter Weisbrod
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


#ifndef LOAD_BALANCER_HH
#define LOAD_BALANCER_HH

#include <vector>
#include <map>
#include "communication/loadBalancer.h"
#include "core/olbDebug.h"

namespace olb {

template<typename T>
int LoadBalancer<T>::loc(const int& glob) {
  return _loc[glob];
}

template<typename T>
int LoadBalancer<T>::loc(int glob) const {
  std::map<int,int>::const_iterator iter = _loc.find(glob);
  return iter->second;
}

template<typename T>
int LoadBalancer<T>::glob(int loc) const {
  return _glob[loc];
}

template<typename T>
int LoadBalancer<T>::rank(const int& glob) {
  return _rank[glob];
}

template<typename T>
int LoadBalancer<T>::rank(int glob) const {
  std::map<int,int>::const_iterator iter = _rank.find(glob);
  return iter->second;
}

template<typename T>
int LoadBalancer<T>::size() const {
  return _size;
}

template<typename T>
void LoadBalancer<T>::print() const {
  OstreamManager clout(std::cout,"LoadBalancer");
  for(unsigned i = 0; i < this->_glob.size(); i++) {
    clout << "glob[" << i << "]=" << this->_glob[i] << std::endl;
  }
  for(std::map<int,int>::const_iterator it = this->_loc.begin(); it != this->_loc.end();
      it++) {
    clout << "loc[" << (*it).first << "]=" << (*it).second << std::endl;
  }
  for(std::map<int,int>::const_iterator it = this->_rank.begin(); it != this->_rank.end();
      it++) {
    clout << "rank[" << (*it).first << "]=" << (*it).second << std::endl;
  }
}


}  // namespace olb
#endif

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Thomas Henn
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

#ifndef CONTACTDETECTION_H
#define CONTACTDETECTION_H

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

/*
 * Prototype for future contact detection algorithms
 **/

template<typename T, template<typename U> class PARTICLETYPE>
class ContactDetection {
public:
  ContactDetection(ParticleSystem3D<T, PARTICLETYPE>& pSys) : _pSys(pSys) {};

  virtual void sort() {};
  virtual int getMatches(int pInt, std::vector<std::pair<size_t, T> >& matches) {
    return 0;
  };

  virtual ~ContactDetection() {};
  //TODO =0 is needed at some time
  virtual ContactDetection<T, PARTICLETYPE>* generate(ParticleSystem3D<T, PARTICLETYPE>& pSys) {
    return this;
  };


protected:
  ParticleSystem3D<T, PARTICLETYPE>& _pSys;
};

} // namespace olb

#endif /*CONTACTDETECTION_H*/

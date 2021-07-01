/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause
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

#ifndef PARTICLE_3D_H
#define PARTICLE_3D_H

#include <string>
#include <vector>

namespace olb {

template<typename T>
class Particle3D {
public:
  Particle3D();
  Particle3D(std::vector<T> pos, T mas = 1., T rad = 1., int id=0);
  Particle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1., int id=0);
  Particle3D(const Particle3D<T>& p);

  inline void setPos(std::vector<T> pos) {
    _pos = pos;
  }
  inline std::vector<T>& getPos() {
    return _pos;
  }
  inline const std::vector<T>& getPos() const {
    return _pos;
  }
  inline void setVel(std::vector<T> vel) {
    _vel = vel;
  }

  inline int getID() {
    return _id;
  }

  inline std::vector<T>& getVel() {
    return _vel;
  }
  inline const std::vector<T>& getVel() const {
    return _vel;
  }
  //  void rotAndResetForce();
  inline void addForce(std::vector<T>& frc);
  inline void setForce(std::vector<T>& frc);
  inline std::vector<T>& getForce();
  inline const std::vector<T>& getForce() const;
  inline void resetForce();

  void serialize(T serial[]);
  void unserialize(T*);
  void print();

  inline const T& getMass() {
    return _mas;
  }

  inline const T& getInvMass() {
    return _invMas;
  }

  inline const T& getMass() const {
    return _mas;
  }
  inline void setMass(T m) {
    _mas = m;
    _invMas = 1. / _mas;
  }
  inline const T& getRad() {
    return _rad;
  }
  inline const T& getRad() const {
    return _rad;
  }
  inline void setRad(T r) {
    _rad = r;
  }
  inline const int& getCuboid() {
    return _cuboid;
  }
  inline void setCuboid(int c) {
    _cuboid = c;
  }
  inline const bool& getActive() {
    return _active;
  }
  inline const bool& getActive() const {
    return _active;
  }
  inline void setActive(bool act) {
    _active = act;
    if (!act) {
      _vel[0] = 0;
      _vel[1] = 0;
      _vel[2] = 0;
    }
  }

  static const int serialPartSize = 14;
  std::vector<std::pair<size_t, T> > _verletList;

protected:
  std::vector<T> _pos;
  std::vector<T> _vel;
  std::vector<T> _force;
  T _invMas;
  T _mas;
  T _rad;
  int _cuboid;
  int _id;
  bool _active;
};

/*
 * Electric Particles
 */
template<typename T>
class ElParticle3D : public Particle3D<T> {
public:
  ElParticle3D();
  ElParticle3D(std::vector<T> pos, T mas = 1., T rad = 1., T charge = 1.);
  ElParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1.,
               T charge = 1.);
  ElParticle3D(const ElParticle3D<T>& p);
  virtual ~ElParticle3D() {
  }
  ;
  virtual void serialize(T serial[]);
  virtual void unserialize(T*);
  static const int serialPartSize = 10;
  T _charge;
};

/*
 * Particles for Agglomeration
 */

template<typename T>
class AggParticle3D : public Particle3D<T> {
public:
  AggParticle3D();
  AggParticle3D(std::vector<T> pos, T mas = 1., T rad = 1.);
  AggParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1.);
  AggParticle3D(const Particle3D<T>& p);

  inline void setMass(T mas) {
    this->_mas = mas;
  }
  inline void setRad(T rad) {
    this->_rad = rad;
  }
  inline const bool& getAggl() {
    return _aggl;
  }
  inline const bool& getAggl() const {
    return _aggl;
  }
  inline void setAggl(bool aggl) {
    _aggl = aggl;
    if (aggl) {
      this->setActive(false);
    }
  }

  static const int serialPartSize = 14;
  void serialize(T serial[]);
  void unserialize(T*);

private:
  bool _aggl;
};

/*
 * Rotating Particles
 */
template<typename T>
class RotatingParticle3D : public Particle3D<T> {
public:
  RotatingParticle3D();
  RotatingParticle3D(std::vector<T> pos, T mas = 1., T rad = 1.);
  RotatingParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1.);
  RotatingParticle3D(const RotatingParticle3D<T>& p);

  static const int serialPartSize = 19;
  void serialize(T serial[]);
  void unserialize(T*);

  inline std::vector<T>& getAVel() {
    return _aVel;
  }
  inline const std::vector<T>& getAVel() const {
    return _aVel;
  }
  inline std::vector<T>& getTorque() {
    return _torque;
  }
  inline const std::vector<T>& getTorque() const {
    return _torque;
  }
private:
  std::vector<T> _aVel;
  std::vector<T> _torque;
};

}

#endif /* PARTICLE_3D_H */

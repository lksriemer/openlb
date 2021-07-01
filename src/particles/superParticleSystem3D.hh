/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn
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

#ifndef SUPERPARTICLESYSTEM_3D_HH
#define SUPERPARTICLESYSTEM_3D_HH

#define shadows

#include <utility>
#include <string>
#include "io/fileName.h"
#include "superParticleSystem3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  CuboidGeometry3D<T>& cuboidGeometry, LoadBalancer<T>& loadBalancer,
  SuperGeometry3D<T>& superGeometry, LBconverter<T>& conv)
  : SuperStructure3D<T>(cuboidGeometry, loadBalancer,
                        superGeometry.getOverlap()),
  clout(std::cout, "SuperParticleSystem3d"),
  _superGeometry(superGeometry),
  _conv(conv)
{
  init();
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  SuperGeometry3D<T>& superGeometry, LBconverter<T>& conv)
  : SuperStructure3D<T>(superGeometry.getCuboidGeometry(),
                        superGeometry.getLoadBalancer(),
                        superGeometry.getOverlap()),
  clout(std::cout, "SuperParticleSystem3d"),
  _superGeometry(superGeometry),
  _conv(conv)
{
  init();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::init()
{
  int rank = 0;
  int size = 0;
#if PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
  size = singleton::mpi().getSize();
#endif
  for (int i = 0; i < this->_cuboidGeometry.getNc(); ++i) {
    if (this->_loadBalancer.rank(i) == rank) {
      auto dummy = new ParticleSystem3D<T, PARTICLETYPE>(_superGeometry, _conv);
      this->_cuboidGeometry.get(i).getOrigin().toStdVector()[0];
      std::vector<T> physPos = this->_cuboidGeometry.get(i).getOrigin().toStdVector();
      std::vector<T> physExtend(3, 0);
      T physR = this->_cuboidGeometry.get(i).getDeltaR();
      for (int j = 0; j < 3; j++) {
        physPos[j] -= .5 * physR;
        physExtend[j] = (this->_cuboidGeometry.get(i).getExtend()[j] + 1)
                        * physR;
      }
      dummy->setPosExt(physPos, physExtend);
      _pSystems.push_back(dummy);
    }
  }

#if PARALLEL_MODE_MPI
  for (int i=0; i<this->_cuboidGeometry.getNc(); ++i) {
    if (this->_loadBalancer.rank(i) == rank) {
      std::vector<int> dummy;
      this->getCuboidGeometry().getNeighbourhood(i, dummy, 3);
      _rankNeighbours.insert(_rankNeighbours.end(), dummy.begin(), dummy.end());
      _cuboidNeighbours.push_back(dummy);
    }
  }
  for (auto& N : _rankNeighbours) {
    N = this->_loadBalancer.rank(N);
  }
#endif

  /* Ein jeder ist sein eigener Nachbar*/
  if (rank == 0) {
    _rankNeighbours.push_back(size-1);
  }
  if (rank == size-1) {
    _rankNeighbours.push_back(0);
  }
  _rankNeighbours.push_back(rank);
  _rankNeighbours.sort();
  _rankNeighbours.unique();
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  SuperParticleSystem3D<T, PARTICLETYPE>& spSys)
  : SuperStructure3D<T>(spSys._cuboidGeometry, spSys._loadBalancer,
                        spSys._overlap),
  clout(std::cout, "SuperParticleSystem3d"),
  _pSystems(spSys._pSystems),
  _superGeometry(spSys._superGeometry),
  _conv(spSys._conv),
  _rankNeighbours(spSys._rankNeighbours),
  _overlap(spSys._overlap)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  SuperParticleSystem3D<T, PARTICLETYPE> const& spSys)
  : SuperStructure3D<T>(spSys._cuboidGeometry, spSys._loadBalancer,
                        spSys._overlap),
  clout(std::cout, "SuperParticleSystem3d"),
  _pSystems(spSys._pSystems),
  _superGeometry(spSys._superGeometry),
  _conv(spSys._conv),
  _rankNeighbours(spSys._rankNeighbours),
  _overlap(spSys._overlap)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  SuperParticleSystem3D<T, PARTICLETYPE> && spSys)
  : SuperStructure3D<T>(spSys._cuboidGeometry, spSys._loadBalancer,
                        spSys._overlap),
  clout(std::cout, "SuperParticleSystem3d"),
  _superGeometry(spSys._superGeometry),
  _conv(spSys._conv),
  _rankNeighbours(spSys._rankNeighbours),
  _overlap(spSys._overlap)
{
  _pSystems = std::move(spSys._pSystems);
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::print()
{
  int no = globalNumOfParticles();
  int active = globalNumOfActiveParticles();
  clout << "activeParticles= " << active << " (" << no << ") " << std::endl;
  //  cout << "[SuperParticleSystem3D] " << _pSystems.size()
  //      << " pSystems on rank " << singleton::mpi().getRank() << "\n";
  //  for (auto pS : _pSystems) {
  //    cout << pS->_particles.size() << " ";
  //  }
  //  std::cout << std::endl;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::print(std::list<int> mat)
{
  std::list<int>::iterator _matIter;
  int no = globalNumOfParticles();
  clout << "globalNumOfParticles=" << no;
  int active = globalNumOfActiveParticles();
  clout << "; activeParticles=" << active;
  for (_matIter = mat.begin(); _matIter != mat.end(); _matIter++) {
    clout << "; material" << *_matIter << "=" << countMaterial(*_matIter);
  }
  clout << std::endl;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::captureEscapeRate(
  std::list<int> mat)
{
  std::list<int>::iterator _matIter;
//  int no = globalNumOfParticles();
  T sum = T();
  // capture rate
  for (_matIter = mat.begin(); _matIter != mat.end(); _matIter++) {
    sum += (T) countMaterial(*_matIter);
  }
  clout << "captureRate=" << 1. - sum / globalNumOfParticles()
        << "; escapeRate=" << sum / globalNumOfParticles() << std::endl;
}

template<typename T, template<typename U> class PARTICLETYPE>
std::vector<ParticleSystem3D<T, PARTICLETYPE>*>& SuperParticleSystem3D<T,
    PARTICLETYPE>::getPSystems()    //+*
{
  return _pSystems;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::simulate(T dT)
{

  //  for (auto pS : _pSystems) {
  //      pS->velocityVerlet1(dT);
  //      if (_overlap > 0) {
  //        pS->makeTree();
  //        cout << "makeTree" << endl;
  //      }
  //  }
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //      pS->velocityVerlet2(dT);
  //  }

  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->computeForce();
  //    pS->predictorCorrector1(dT);
  //  }
  //
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->computeForce();
  //    pS->predictorCorrector2(dT);
  //  }

  //  for (auto pS : _pSystems) {
  //    pS->rungeKutta4_1(dT);
  //  }
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->rungeKutta4_2(dT);
  //  }
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->rungeKutta4_3(dT);
  //  }
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->rungeKutta4_4(dT);
  //  }
  //  updateParticleDistribution();

  for (auto pS : _pSystems) {
    pS->_contactDetection->sort();

    pS->simulate(dT);
    pS->computeBoundary();

  }
  updateParticleDistribution();

}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setOverlap(T overlap)
{
  _overlap = overlap;
  //  for (auto pS : _pSystems) {
  //    pS->setOverlap(_overlap);
  //  }
}

template<typename T, template<typename U> class PARTICLETYPE>
T SuperParticleSystem3D<T, PARTICLETYPE>::getOverlap()
{
  return _overlap;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::numOfPSystems()
{
  return _pSystems.size();
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::globalNumOfParticles()
{
#if PARALLEL_MODE_MPI
  int buffer = rankNumOfParticles();
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  return buffer;
#else
  return rankNumOfParticles();
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::globalNumOfShadowParticles()
{
#if PARALLEL_MODE_MPI
  int buffer = rankNumOfShadowParticles();
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  return buffer;
#else
  return rankNumOfShadowParticles();
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::rankNumOfParticles()
{
  int num = 0;
  for (auto pS : _pSystems) {
    num += pS->size();
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::rankNumOfShadowParticles()
{
  int num = 0;
  for (auto pS : _pSystems) {
    num += pS->_shadowParticles.size();
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::rankNumOfActiveParticles()
{
  int num = 0;
  for (auto pS : _pSystems) {
    num += pS->numOfActiveParticles();
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::countLocMaterial(int mat)
{
  int num = 0;
  for (auto pS : _pSystems) {
    num += pS->countMaterial(mat);
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::countMaterial(int mat)
{
#if PARALLEL_MODE_MPI
  int buffer = countLocMaterial(mat);
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  return buffer;
#else
  return countLocMaterial(mat);
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::globalNumOfActiveParticles()
{
#if PARALLEL_MODE_MPI
  int buffer = rankNumOfActiveParticles();
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  return buffer;
#else
  return rankNumOfActiveParticles();
#endif
}

// TODO class olb::ParticleSystem3D<double, olb::Particle3D>’ has no member named ‘radiusDistribution
/*
 template<typename T, template<typename U> class PARTICLETYPE>
 std::map<T, int> SuperParticleSystem3D<T, PARTICLETYPE>::rankRadiusDistribution()
 {
 std::map<T, int> distrAll;
 typename std::map<T, int>::iterator ita;
 for (auto pS : _pSystems) {
 std::map<T, int> distrPs;
 distrPs = pS->radiusDistribution();
 for (typename std::map<T, int>::iterator itp = distrPs.begin();
 itp != distrPs.end(); ++itp) {
 T radKeyP = itp->first;
 int countP = itp->second;
 if (distrAll.count(radKeyP) > 0) {
 ita = distrAll.find(radKeyP);
 ita->second = ita->second + countP;
 } else {
 distrAll.insert(std::pair<T, int>(radKeyP, countP));
 }
 }
 }
 return distrAll;
 }
 */

template<typename T, template<typename U> class PARTICLETYPE>
std::vector<int> SuperParticleSystem3D<T, PARTICLETYPE>::numOfForces()
{
  std::vector<int> dummy;
  for (auto pS : _pSystems) {
    dummy.push_back(pS->numOfForces());
  }
  return dummy;
}

template<typename T, template<typename U> class PARTICLETYPE>
template<template<typename V> class DESCRIPTOR>
void SuperParticleSystem3D<T, PARTICLETYPE>::setVelToFluidVel(
  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setVelToAnalyticalVel(
  AnalyticalConst3D<T,T>& aVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToAnalyticalVel(aVel);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
bool SuperParticleSystem3D<T, PARTICLETYPE>::findCuboid(PARTICLETYPE<T> &p)
{
  return findCuboid(p, _overlap);
}

template<typename T, template<typename U> class PARTICLETYPE>
bool SuperParticleSystem3D<T, PARTICLETYPE>::findCuboid(PARTICLETYPE<T> &p,
    T overlap)
{
  int C = this->_cuboidGeometry.get_iC(p.getPos()[0], p.getPos()[1],
                                       p.getPos()[2], overlap);
  if (C != this->_cuboidGeometry.getNc()) {
    p.setCuboid(C);
    return 1;
  }
  clout << "Lost Particle! Pos: " << p.getPos()[0] << " " << p.getPos()[1]
        << " " << p.getPos()[2] << " Vel: " << p.getVel()[0] << " "
        << p.getVel()[1] << " " << p.getVel()[2] << " " << std::endl;
  p.setActive(false);
  return 0;
}

template<typename T, template<typename U> class PARTICLETYPE>
bool SuperParticleSystem3D<T, PARTICLETYPE>::checkCuboid(PARTICLETYPE<T>& p,
    T overlap)
{
  return checkCuboid(p, overlap, p.getCuboid());
}

template<typename T, template<typename U> class PARTICLETYPE>
bool SuperParticleSystem3D<T, PARTICLETYPE>::checkCuboid(PARTICLETYPE<T>& p,
    T overlap, int iC)
{
  return this->_cuboidGeometry.get(iC).physCheckPoint(p.getPos()[0],
         p.getPos()[1],
         p.getPos()[2], overlap);
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::updateParticleDistribution()
{
  /* Find particles on wrong cuboid, store in relocate and delete*/
  //maps particles to their new rank
  //  std::multimap<int, PARTICLETYPE<T> > relocate;
  _relocate.clear();
  //#define shadows
#ifdef shadows
  _relocateShadow.clear();
#endif
  for (unsigned int pS = 0; pS < _pSystems.size(); ++pS) {
    _pSystems[pS]->_shadowParticles.clear();
    auto par = _pSystems[pS]->_particles.begin();
    while (par != _pSystems[pS]->_particles.end()) {
      //      if (par->getActive()) {
      //Check if particle is still in his cuboid
      if (checkCuboid(*par, 0)) {
#ifdef shadows
        //If yes --> check if its inside boundary layer
        if (!(checkCuboid(*par, -_overlap))) {
          std::set<int> sendTo;
          for (int iC = 0; iC < this->_cuboidGeometry.getNc(); iC++) {
            if (par->getCuboid() != iC && checkCuboid(*par, _overlap, iC)) {
              int rank = this->_loadBalancer.rank(iC);
              if (!sendTo.count(rank)) {
                _relocateShadow.insert(std::make_pair(rank, (*par)));
                sendTo.insert(rank);
              }
            }
          }
        }
#endif
        //If not --> find new cuboid
      } else {
        findCuboid(*par, 0);
        _relocate.insert(
          std::make_pair(this->_loadBalancer.rank(par->getCuboid()), (*par)));
        par = _pSystems[pS]->_particles.erase(par);
        par--;
      }
      //      }
      par++;
    }
  }

  /* Communicate number of Particles per cuboid*/
#ifdef PARALLEL_MODE_MPI
  singleton::MpiNonBlockingHelper mpiNbHelper;
  //mpiNbHelper.allocate(_rankNeighbours.size());
  int k=0;

  /* Serialize particles */
  _send_buffer.clear();
  T buffer[PARTICLETYPE<T>::serialPartSize];
  for (auto rN : _relocate) {
    rN.second.serialize(buffer);
    _send_buffer[rN.first].insert(_send_buffer[rN.first].end(), buffer, buffer+PARTICLETYPE<T>::serialPartSize);
  }

  /*Send Particles */
  k=0;
  int noSends = 0;
  for (auto rN : _rankNeighbours) {
    if (_send_buffer[rN].size() > 0) {
      ++noSends;
    }
  }
//  int noSends = _send_buffer.size();
  if (noSends > 0) {
    mpiNbHelper.allocate(noSends);
//    cout << mpiNbHelper.get_size() << std::endl;
    for (auto rN : _rankNeighbours) {
      if (_send_buffer[rN].size() > 0) {
        singleton::mpi().iSend<double>(&_send_buffer[rN][0], _relocate.count(rN)*PARTICLETYPE<T>::serialPartSize, rN,&mpiNbHelper.get_mpiRequest()[k++], 1);
      }
    }
  }

  /*Receive and add particles*/
  singleton::mpi().barrier();
  k=0;
  int flag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
  if (flag) {
    for (auto rN : _rankNeighbours) {
      MPI_Status status;
      int flag = 0;
      MPI_Iprobe(rN, 1, MPI_COMM_WORLD, &flag, &status);
      if (flag) {
        int number_amount = 0;
        MPI_Get_count(&status, MPI_DOUBLE, &number_amount);
        T recv_buffer[number_amount];
        singleton::mpi().receive(recv_buffer, number_amount, rN, 1);

        for (int iPar=0; iPar<number_amount; iPar+=PARTICLETYPE<T>::serialPartSize) {
          PARTICLETYPE<T> p;
          p.unserialize(&recv_buffer[iPar]);
          if (singleton::mpi().getRank() == this->_loadBalancer.rank(p.getCuboid())) {
            _pSystems[this->_loadBalancer.loc(p.getCuboid())]->addParticle(p);
          }
        }
      }
    }
  }
  if (noSends > 0) {
    singleton::mpi().waitAll(mpiNbHelper);
  }

#ifdef shadows
  /**************************************************************************************************************/
  //Same Again for shadowParticles
  mpiNbHelper.allocate(_rankNeighbours.size());
  k=0;
  /* Serialize particles */
  _send_buffer.clear();
  //  T buffer[PARTICLETYPE<T>::serialPartSize];
  for (auto rN : _relocateShadow) {
    //    std::vector<T> buffer = rN.second.serialize();
    rN.second.serialize(buffer);
    _send_buffer[rN.first].insert(_send_buffer[rN.first].end(), buffer, buffer+PARTICLETYPE<T>::serialPartSize);
  }

  /*Send Particles */
  k=0;
  noSends = _send_buffer.size();
  if (noSends > 0) {
    mpiNbHelper.allocate(noSends);
    for (auto rN : _rankNeighbours) {
      if (_send_buffer[rN].size() > 0) {
        singleton::mpi().iSend<double>(&_send_buffer[rN][0], _relocateShadow.count(rN)*PARTICLETYPE<T>::serialPartSize, rN,&mpiNbHelper.get_mpiRequest()[k++], 4);
      }
    }
  }

  singleton::mpi().barrier();
  k=0;
  flag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
  if (flag) {
    for (auto rN : _rankNeighbours) {
      MPI_Status status;
      int flag = 0;
      MPI_Iprobe(rN, 4, MPI_COMM_WORLD, &flag, &status);
      //      cout << "Message from " << rN << " found on " << singleton::mpi().getRank() << std::endl;
      if (flag) {
        int number_amount = 0;
        MPI_Get_count(&status, MPI_DOUBLE, &number_amount);
        //        cout << "Contains " << number_amount << " infos" << std::endl;
        T recv_buffer[number_amount];
        singleton::mpi().receive(recv_buffer, number_amount, rN, 4);
        for (int iPar=0; iPar<number_amount; iPar+=PARTICLETYPE<T>::serialPartSize) {
          //          std::cout << "Particle unserialized" << std::endl;
          PARTICLETYPE<T> p;
          p.unserialize(&recv_buffer[iPar]);
          addShadowParticle(p);
        }
      }
    }
  }
  if (noSends > 0) {
    singleton::mpi().waitAll(mpiNbHelper);
  }
#endif
  mpiNbHelper.free();
#else

  for (auto& par : _relocate) {
    addParticle(par.second);
  }
  for (auto& par : _relocateShadow) {
    addShadowParticle(par.second);
  }
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticle(PARTICLETYPE<T>& p)
{
  if (findCuboid(p)) {
#if PARALLEL_MODE_MPI
    if (singleton::mpi().getRank() == this->_loadBalancer.rank(p.getCuboid())) {
      _pSystems[this->_loadBalancer.loc(p.getCuboid())]->addParticle(p);
    } else {
      clout << "Particle not found on Cuboid: " << p.getCuboid() << std::endl;
      clout << "Ppos: " << p.getPos()[0] << " " << p.getPos()[1] << " " << p.getPos()[2]<< std::endl;
    }
#else
    _pSystems[this->_loadBalancer.loc(p.getCuboid())]->addParticle(p);
#endif
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticle(
  IndicatorF3D<T>& ind, T mas, T rad, int no, std::vector<T> vel)
{
  srand(clock());
  //  srand(rand());
  std::vector<T> pos(3, 0.);
  bool indic[1] = { false };

  no += globalNumOfParticles();
  while (globalNumOfParticles() < no) {
    pos[0] = ind.getMin()[0]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
    pos[1] = ind.getMin()[1]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
    pos[2] = ind.getMin()[2]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#if PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    int x0, y0, z0, C;
    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      C = locLat[0];
      if (this->_loadBalancer.rank(C) == singleton::mpi().getRank()) {
        x0 = locLat[1];
        y0 = locLat[2];
        z0 = locLat[3];
        if (_superGeometry.get(C, x0, y0, z0) == 1
            && _superGeometry.get(C, x0, y0 + 1, z0) == 1
            && _superGeometry.get(C, x0, y0, z0 + 1) == 1
            && _superGeometry.get(C, x0, y0 + 1, z0 + 1) == 1
            && _superGeometry.get(C, x0 + 1, y0, z0) == 1
            && _superGeometry.get(C, x0 + 1, y0 + 1, z0) == 1
            && _superGeometry.get(C, x0 + 1, y0, z0 + 1) == 1
            && _superGeometry.get(C, x0 + 1, y0 + 1, z0 + 1) == 1
            && ind(indic, &pos[0])) {
          PARTICLETYPE<T> p(pos, vel, mas, rad);
          addParticle(p);
        }
      }
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticle(int material, T mas,
    T rad, int no)
{
  srand(time(0));
  std::vector<T> pos(3, 0.), min(3, 0.), max(3, 0.);
  min = _superGeometry.getStatistics().getMinPhysR(material);
  max = _superGeometry.getStatistics().getMaxPhysR(material);
  for (int i = 0; i < 3; i++) {
    min[i] -= _conv.getLatticeL();
    max[i] += 2 * _conv.getLatticeL();
  }
  no += globalNumOfParticles();
  while (globalNumOfParticles() < no) {
    pos[0] = min[0] + (T) (rand() % 100000) / 100000. * (max[0] - min[0]);
    pos[1] = min[1] + (T) (rand() % 100000) / 100000. * (max[1] - min[1]);
    pos[2] = min[2] + (T) (rand() % 100000) / 100000. * (max[2] - min[2]);

#if PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      if (this->_loadBalancer.rank(locLat[0]) == singleton::mpi().getRank()) {
        if ((_superGeometry.get(locLat[0], locLat[1], locLat[2], locLat[3])
             == material
             || _superGeometry.get(locLat[0], locLat[1], locLat[2], locLat[3])
             == 1)
            && (_superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                   locLat[3]) == material
                || _superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                      locLat[3]) == 1)
            && (_superGeometry.get(locLat[0], locLat[1], locLat[2],
                                   locLat[3] + 1) == material
                || _superGeometry.get(locLat[0], locLat[1], locLat[2],
                                      locLat[3] + 1) == 1)
            && (_superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                   locLat[3] + 1) == material
                || _superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                      locLat[3] + 1) == 1)
            && (_superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                   locLat[3]) == material
                || _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                      locLat[3]) == 1)
            && (_superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                   locLat[3]) == material
                || _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                      locLat[3]) == 1)
            && (_superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                   locLat[3] + 1) == material
                || _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                      locLat[3] + 1) == 1)
            && (_superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                   locLat[3] + 1) == material
                || _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                      locLat[3] + 1) == 1)) {
          PARTICLETYPE<T> p(pos, mas, rad);
          addParticle(p);
        }
      }
    }
  }
  singleton::mpi().barrier();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addShadowParticle(
  PARTICLETYPE<T>& p)
{
  for (unsigned int pS = 0; pS < _pSystems.size(); ++pS) {
    int globIC = this->_loadBalancer.glob(pS);
    if (globIC != p.getCuboid() && checkCuboid(p, _overlap, globIC)
        && !checkCuboid(p, 0, globIC)) {
      _pSystems[pS]->addShadowParticle(p);
      //      cout << "ShadowParticle added" << std::endl;
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
std::vector<ParticleSystem3D<T, PARTICLETYPE>*> SuperParticleSystem3D<T,
    PARTICLETYPE>::getParticleSystems()
{
  return _pSystems;
}

template<typename T, template<typename U> class PARTICLETYPE>
ParticleSystem3D<T, PARTICLETYPE>& SuperParticleSystem3D<T, PARTICLETYPE>::operator[](
  int i)
{
  return *(_pSystems[i]);
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addForce(
  std::shared_ptr<Force3D<T, PARTICLETYPE> > f)
{
  for (auto pS : _pSystems) {
    pS->addForce(f);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addBoundary(
  std::shared_ptr<Boundary3D<T, PARTICLETYPE> > b)
{
  for (auto pS : _pSystems) {
    pS->addBoundary(b);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::saveToFile(std::string name)
{
  std::string fullName = createFileName(name) + ".particles";

  int rank = 0;

#ifdef PARALLEL_MODE_MPI
  int size = 1;
  size = singleton::mpi().getSize();
  rank = singleton::mpi().getRank();
#endif

  if (rank == 0) {
    std::ofstream fout(fullName.c_str(), std::ios::trunc);
    if (!fout) {
      clout << "Error: could not open " << fullName << std::endl;
    }
    fout.close();
  }
#if PARALLEL_MODE_MPI
  if (rank > 0) {
    int prev = rank - 1;
    int buffer = 0;
    MPI_Status status;
    MPI_Recv(&buffer, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &status);
  }
#endif
  for (auto pS : _pSystems) {
    pS->saveToFile(fullName);
  }
#if PARALLEL_MODE_MPI
  if (rank < size - 1) {
    int next = rank + 1;
    int buffer = 0;
    MPI_Send(&buffer, 1, MPI_INT, next, 0, MPI_COMM_WORLD);
  }
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticlesFromFile(
  std::string name, T mass, T radius)
{
  std::string fullName = createFileName(name) + ".particles";
  std::ifstream fin(fullName.c_str());

  std::string line;
  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    T buffer[PARTICLETYPE<T>::serialPartSize];
    for (int i = 0; i < PARTICLETYPE<T>::serialPartSize; i++) {
      iss >> buffer[i];
    }
    PARTICLETYPE<T> p;
    p.unserialize(buffer);
    if (radius != 0) {
      p.setRad(radius);
    }
    if (mass != 0) {
      p.setMass(mass);
    }
    addParticle(p);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setContactDetection(
  ContactDetection<T, PARTICLETYPE>& contactDetection)
{
  clout << "Setting ContactDetectionAlgorithm" << std::endl;

  for (auto pS : _pSystems) {
    pS->setContactDetection(contactDetection);
  }
}

//TODO move to contrib
/*
 template<typename T, template<typename U> class PARTICLETYPE>
 void SuperParticleSystem3D<T, PARTICLETYPE>::setContactDetectionAlgorithmNanoflann(
 T sRad2)
 {
 clout << "Setting ContactDetectionAlgorithm = Nanoflann" << std::endl;

 for (auto pS : _pSystems) {
 pS->setContactDetectionAlgorithmNanoflann(sRad2);
 }
 }

 template<typename T, template<typename U> class PARTICLETYPE>
 void SuperParticleSystem3D<T, PARTICLETYPE>::setContactDetectionAlgorithmMunjizaNBS(T spacing)
 {
 int i = 0;
 if (spacing == 0) {
 spacing = 2*this->_conv.getLatticeL();
 }
 clout << "Setting ContactDetectionAlgorithm = MunjizaNBS" << std::endl;
 clout << "Spacing: "  << spacing << std::endl;
 for (auto pS : _pSystems) {
 pS->setContactDetectionAlgorithmMunjizaNBS(
 _overlap, spacing);
 i++;
 }
 }

 template<typename T, template<typename U> class PARTICLETYPE>
 void SuperParticleSystem3D<T, PARTICLETYPE>::setContactDetectionAlgorithmDCell(T spacing)
 {
 int i = 0;
 if (spacing == 0) {
 spacing == this->_conv.getLatticeL();
 }

 clout << "Setting ContactDetectionAlgorithm = DCell" << std::endl;

 for (auto pS : _pSystems) {
 pS->setContactDetectionAlgorithmDCell(
 _overlap, spacing);
 i++;
 }
 }

 template<typename T, template<typename U> class PARTICLETYPE>
 void SuperParticleSystem3D<T, PARTICLETYPE>::setContactDetectionAlgorithmPLattice()
 {
 int i = 0;
 clout << "Setting ContactDetectionAlgorithm = PLattice" << std::endl;

 for (auto pS : _pSystems) {
 pS->setContactDetectionAlgorithmPLattice(
 _overlap, this->_cuboidGeometry.get(i).getDeltaR());
 i++;
 }
 }
 */

//template<typename T, template<typename U> class PARTICLETYPE>
//template<template<typename V> class DESCRIPTOR>
//void SuperParticleSystem3D<T, PARTICLETYPE>::particleOnFluid(
//    SuperLattice3D<T, DESCRIPTOR>& sLattice, T eps,
//    SuperGeometry3D<T>& sGeometry) {
//  for (unsigned int i = 0; i < _pSystems.size(); ++i) {
//    _pSystems[i]->particleOnFluid(
//        //        sLattice.getExtendedBlockLattice(i),
//        sLattice.getBlockLattice(i),
//        sLattice.getCuboidGeometry().get(this->_loadBalancer.glob(i)),
//        sLattice.getOverlap(), eps, sGeometry.getExtendedBlockGeometry(i));
//  }
//}
//
//template<typename T, template<typename U> class PARTICLETYPE>
//template<template<typename V> class DESCRIPTOR>
//void SuperParticleSystem3D<T, PARTICLETYPE>::resetFluid(
//    SuperLattice3D<T, DESCRIPTOR>& sLattice) {
//  for (unsigned int i = 0; i < _pSystems.size(); ++i) {
//    _pSystems[i]->resetFluid(
//        //        sLattice.getExtendedBlockLattice(i),
//        sLattice.getBlockLattice(i),
//        sLattice.getCuboidGeometry().get(this->_loadBalancer.glob(i)),
//        sLattice.getOverlap());
//  }
//}

//template<typename T>
//class SuperRotParticleSystem3D : public SuperParticleSystem3D<T, RotatingParticle3D> {
// public:
//  SuperRotParticleSystem3D(SuperGeometry3D<T>& sg, LBconverter<T>& conv): SuperParticleSystem3D<T, RotatingParticle3D>(sg, conv)
//  {};
//  void simulate(T dT) {
//   for (auto pS : this->_pSystems) {
//     pS->_contactDetection->sort();
//     pS->simulate(dT);
//     pS->computeTorque();
//     pS->computeBoundary();
//   }
//   this->updateParticleDistribution();
// }
//};

}//namespace olb

#endif /* SUPERPARTICLESYSTEM_3D_HH */

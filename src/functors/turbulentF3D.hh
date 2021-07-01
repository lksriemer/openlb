/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Patrick Nathen, Mathias J. Krause
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

#ifndef TURBULENT_F_3D_HH
#define TURBULENT_F_3D_HH

#include<vector>
#include<cmath>
#include<string>

#include "functors/turbulentF3D.h"
#include "functors/blockLatticeBaseF3D.h"
#include "functors/blockLatticeLocalF3D.h"
#include "functors/superLatticeBaseF3D.h"
#include "functors/indicatorBaseF.h"
#include "core/superLattice3D.h"
#include "utilities/vectorHelpers.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity


namespace olb {


///////////////////////////// SuperLatticeYplus3D //////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeYplus3D<T,DESCRIPTOR>::SuperLatticeYplus3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
  const LBconverter<T>& converter, SuperGeometry3D<T>& superGeometry,
  IndicatorF3D<bool,T>& indicator, const int material )
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1),
    _superGeometry(superGeometry), _indicator(indicator), _material(material)
{
  this->_name = "yPlus";
}

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeYplus3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];
  int lociz = input[3];

  std::vector<T> output(1,T());

  if( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() )
  {
    std::vector<T> normalTemp(3,T());
    std::vector<T> normal(3,T());       // normalized
    T counter = T();
    T distance = T();
    if(_superGeometry.get(input) == 1)
    {
      for(int iPop = 1; iPop < DESCRIPTOR<T>::q; iPop++)
      {
        if(_superGeometry.get(input[0],
                              input[1] + DESCRIPTOR<T>::c[iPop][0],
                              input[2] + DESCRIPTOR<T>::c[iPop][1],
                              input[3] + DESCRIPTOR<T>::c[iPop][2]) == _material)
        {
          counter++;
          normalTemp[0] += DESCRIPTOR<T>::c[iPop][0];
          normalTemp[1] += DESCRIPTOR<T>::c[iPop][1];
          normalTemp[2] += DESCRIPTOR<T>::c[iPop][2];
        }
      }
      if(counter != 0)
      {
        // get physical Coordinates at intersection
        std::vector<T> physR = _superGeometry.getCuboidGeometry().getPhysR(input);
        T voxelSize = _superGeometry.getCuboidGeometry().get(globIC).getDeltaR();

        normal = util::normalize(normalTemp);

        std::vector<T> direction(normal);
        direction[0] = voxelSize*normal[0]*2.;
        direction[1] = voxelSize*normal[1]*2.;
        direction[2] = voxelSize*normal[2]*2.;

        // calculate distance to STL file
        if( _indicator.distance(distance, physR, direction) )
        {
          // call stress at this point
          T rho;
          T u[3];
          T pi[6];
          this->_sLattice.getBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix, lociy, lociz).computeRhoU(rho, u);
          this->_sLattice.getBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).get(locix, lociy, lociz).computeStress(pi);

          // Compute phys stress tau = mu*du/dx
          T omega = this->_converter.getOmega();
          T dt = this->_converter.physTime();
          T physFactor = -omega*DESCRIPTOR<T>::invCs2/rho/2./dt*this->_converter.physRho(rho)*this->_converter.getCharNu();

          //  Totel Stress projected from cell in normal direction on obstacle
          T Rx = pi[0]*physFactor*normal[0] + pi[1]*physFactor*normal[1] + pi[2]*physFactor*normal[2];
          T Ry = pi[1]*physFactor*normal[0] + pi[3]*physFactor*normal[1] + pi[4]*physFactor*normal[2];
          T Rz = pi[2]*physFactor*normal[0] + pi[4]*physFactor*normal[1] + pi[5]*physFactor*normal[2];

          // Stress appearing as pressure in corresponding direction is calculated and substracted
          T R_res_pressure = normal[0]*pi[0]*physFactor*normal[0] + normal[0]*pi[1]*physFactor*normal[1] + normal[0]*pi[2]*physFactor*normal[2]
                            +normal[1]*pi[1]*physFactor*normal[0] + normal[1]*pi[3]*physFactor*normal[1] + normal[1]*pi[4]*physFactor*normal[2]
                            +normal[2]*pi[2]*physFactor*normal[0] + normal[2]*pi[4]*physFactor*normal[1] + normal[2]*pi[5]*physFactor*normal[2];

          Rx -= R_res_pressure *normal[0];
          Ry -= R_res_pressure *normal[1];
          Rz -= R_res_pressure *normal[2];

          T tau_wall = sqrt(Rx*Rx+Ry*Ry+Rz*Rz);
          T u_tau = sqrt(tau_wall/this->_converter.physRho(rho));
          //y_plus
          output[0] = distance*u_tau / this->_converter.getCharNu(); 
        } // if 4
      }
    }
  }
  return output;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticeADM3D<T,DESCRIPTOR>::BlockLatticeADM3D
  (BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, T sigma, int order )
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,4), _sigma(sigma), _order(order)
{ this->_name = "ADMfilter"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> BlockLatticeADM3D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  std::vector<T> output(4,T());

  //////Declaration of all variables needed for filtering //////
 
  int globX = input[0];
  int globY = input[1];
  int globZ = input[2];

  std::vector<int> u_000(3,int());
  u_000[0]=globX, u_000[1]=globY, u_000[2]=globZ;

  ///x filter ini
  std::vector<int> u_p00(3,int());
  u_p00[0]=globX+1, u_p00[1]=globY, u_p00[2]=globZ;

  std::vector<int> u_2p00(3,int());
  u_2p00[0]=globX+2, u_2p00[1]=globY, u_2p00[2]=globZ;

  std::vector<int> u_3p00(3,int());
  u_3p00[0]=globX+3, u_3p00[1]=globY, u_3p00[2]=globZ;

  std::vector<int> u_n00(3,int());
  u_n00[0]=globX-1, u_n00[1]=globY, u_n00[2]=globZ;

  std::vector<int> u_2n00(3,int());
  u_2n00[0]=globX-2, u_2n00[1]=globY, u_2n00[2]=globZ;

  std::vector<int> u_3n00(3,int());
  u_3n00[0]=globX-3, u_3n00[1]=globY, u_3n00[2]=globZ;

  ///y filter ini
  std::vector<int> u_0p0(3,int());
  u_0p0[0]=globX, u_0p0[1]=globY+1, u_0p0[2]=globZ;

  std::vector<int> u_02p0(3,int());
  u_02p0[0]=globX, u_02p0[1]=globY+2, u_02p0[2]=globZ;

  std::vector<int> u_03p0(3,int());
  u_03p0[0]=globX, u_03p0[1]=globY+3, u_03p0[2]=globZ;

  std::vector<int> u_0n0(3,int());
  u_0n0[0]=globX, u_0n0[1]=globY-1, u_0n0[2]=globZ;

  std::vector<int> u_02n0(3,int());
  u_02n0[0]=globX, u_02n0[1]=globY-2, u_02n0[2]=globZ;

  std::vector<int> u_03n0(3,int());
  u_03n0[0]=globX, u_03n0[1]=globY-3, u_03n0[2]=globZ;
   
  ///z filter ini
  std::vector<int> u_00p(3,int());
  u_00p[0]=globX, u_00p[1]=globY, u_00p[2]=globZ+1;

  std::vector<int> u_002p(3,int());
  u_002p[0]=globX, u_002p[1]=globY, u_002p[2]=globZ+2;

  std::vector<int> u_003p(3,int());
  u_003p[0]=globX, u_003p[1]=globY, u_003p[2]=globZ+3;

  std::vector<int> u_00n(3,int());
  u_00n[0]=globX, u_00n[1]=globY, u_00n[2]=globZ-1;

  std::vector<int> u_002n(3,int());
  u_002n[0]=globX, u_002n[1]=globY, u_002n[2]=globZ-2;

  std::vector<int> u_003n(3,int());
  u_003n[0]=globX, u_003n[1]=globY, u_003n[2]=globZ-3;
   
  BlockLatticeDensity3D<T, DESCRIPTOR> density(this-> _blockLattice);
  BlockLatticeVelocity3D<T, DESCRIPTOR> velocity(this-> _blockLattice);
  //////// end of declaration 

  switch(_order) {
    case(2): //////// second order//////// 
    {
      T d_0 = 6./16.; 
      T d_1 = -4./16.;
      T d_2 = 1./16.;
      output[0] = density(u_000)[0] - _sigma*(d_2*(density(u_2n00)[0]+density(u_02n0)[0]+density(u_002n)[0]) +
                                    d_1*(density(u_n00)[0]+density(u_0n0)[0]+density(u_00n)[0])+
                                    d_0*(density(u_000)[0]+density(u_000)[0]+density(u_000)[0])+
                                    d_1*(density(u_p00)[0]+density(u_0p0)[0]+density(u_00p)[0])+
                                    d_2*(density(u_2p00)[0]+density(u_02p0)[0]+density(u_002p)[0]) );
      for (int i = 0; i < 3; ++i ) {        
        output[i+1] = velocity(u_000)[i] - _sigma*(d_2*(velocity(u_2n00)[i] + velocity(u_02n0)[i] + velocity(u_002n)[i]) +
                                                  d_1*(velocity(u_n00)[i] + velocity(u_0n0)[i] + velocity(u_00n)[i])+
                                                  d_0*(velocity(u_000)[i] + velocity(u_000)[i] + velocity(u_000)[i])+
                                                  d_1*(velocity(u_p00)[i] + velocity(u_0p0)[i] + velocity(u_00p)[i])+
                                                  d_2*(velocity(u_2p00)[i] + velocity(u_02p0)[i] + velocity(u_002p)[i]) );
      }
    }
    case(3): ////// third order ///////////
    { 
      T d_0 = 5./16.;  
      T d_1 = -15./64.;
      T d_2 = 3./32.;
      T d_3 = -1./64.;
      output[0] = density(u_000)[0]- _sigma*(d_3*(density(u_3n00)[0]+density(u_03n0)[0]+density(u_003n)[0])+
                                    d_2*(density(u_2n00)[0]+density(u_02n0)[0]+density(u_002n)[0]) +
                                    d_1*(density(u_n00)[0]+density(u_0n0)[0]+density(u_00n)[0])+
                                    d_0*(density(u_000)[0]+density(u_000)[0]+density(u_000)[0])+
                                    d_1*(density(u_p00)[0]+density(u_0p0)[0]+density(u_00p)[0])+
                                    d_2*(density(u_2p00)[0]+density(u_02p0)[0]+density(u_002p)[0])+
                                    d_3*(density(u_3p00)[0]+density(u_03p0)[0]+density(u_003p)[0]) );
      for (int i = 0; i < 3; ++i ) {        
        output[i+1] = velocity(u_000)[i] - _sigma*(d_3*(velocity(u_3n00)[i] + velocity(u_03n0)[i]+velocity(u_003n)[i])+
                                                  d_2*(velocity(u_2n00)[i] + velocity(u_02n0)[i]+velocity(u_002n)[i])+
                                                  d_1*(velocity(u_n00)[i] + velocity(u_0n0)[i]+velocity(u_00n)[i])+
                                                  d_0*(velocity(u_000)[i] + velocity(u_000)[i]+velocity(u_000)[i])+
                                                  d_1*(velocity(u_p00)[i] + velocity(u_0p0)[i]+velocity(u_00p)[i])+
                                                  d_2*(velocity(u_2p00)[i] + velocity(u_02p0)[i]+velocity(u_002p)[i])+
                                                  d_3*(velocity(u_3p00)[i] + velocity(u_03p0)[i]+velocity(u_003p)[i]) );
      }
    }
  }
  return output;
}

template <typename T, template <typename U> class DESCRIPTOR>
void BlockLatticeADM3D<T,DESCRIPTOR>::execute(std::vector<int> input) {

  this-> _blockLattice.get(input[0],input[1],input[2]).defineExternalField( DESCRIPTOR<T>::ExternalField::rhoIsAt, 1, &this->operator()(input)[0] );
  this-> _blockLattice.get(input[0],input[1],input[2]).defineExternalField( DESCRIPTOR<T>::ExternalField::velocityBeginsAt, 3, &this->operator()(input)[0] );
}

template <typename T, template <typename U> class DESCRIPTOR>
void BlockLatticeADM3D<T,DESCRIPTOR>::execute() {

  int nX = this-> _blockLattice.getNx();
  int nY = this-> _blockLattice.getNy();
  int nZ = this-> _blockLattice.getNz();
  std::vector<int> i(3,int());
  for (i[0]=0; i[0]<nX; i[0]++) {
    for (i[1]=0; i[1]<nY; i[1]++) {
      for (i[2]=0; i[2]<nZ; i[2]++) {
        execute(i);
      }
    }
  }
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeADM3D<T,DESCRIPTOR>::SuperLatticeADM3D
  (SuperLattice3D<T,DESCRIPTOR>& sLattice, T sigma, int order) : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,4),
  _sigma(sigma), _order(order)
{ this->_name = "ADMfilter"; }

template <typename T, template <typename U> class DESCRIPTOR>
std::vector<T> SuperLatticeADM3D<T,DESCRIPTOR>::operator() (std::vector<int> input)
{
  int globIC = input[0];
  if ( this->_sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
    std::vector<int> inputLocal(3,T());
    T overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    BlockLatticeADM3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)), _sigma, _order );

    return blockLatticeF(inputLocal);
  } else {
    return std::vector<T>();
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
void SuperLatticeADM3D<T,DESCRIPTOR>::execute(SuperGeometry3D<T>& superGeometry, const int material) {

  this->_sLattice.communicate();
  CuboidGeometry3D<T>& cGeometry =  this->_sLattice.getCuboidGeometry();
  LoadBalancer<T>& load = this->_sLattice.getLoadBalancer();
  int overlap = this->_sLattice.getOverlap();

  for (int iC=0; iC<load.size(); iC++) {
    BlockLatticeADM3D<T,DESCRIPTOR> blockLatticeF( this->_sLattice.getExtendedBlockLattice(iC), _sigma, _order );
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    int nZ = cGeometry.get(load.glob(iC)).getNz();

    std::vector<int> i(3,overlap);
    for (; i[0]<nX+overlap; i[0]++) {
      for (; i[1]<nY+overlap; i[1]++) {
        for (; i[2]<nZ+overlap; i[2]++) {
          if (superGeometry.getExtendedBlockGeometry(iC).get(i) == material)
            blockLatticeF.execute(i);
        }
      }
    }
  }
}


} // end namespace olb
#endif



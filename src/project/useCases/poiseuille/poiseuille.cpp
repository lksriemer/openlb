/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006 Jonas Latt
 *  Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
 *  E-mail: jonas.latt@gmail.com
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

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef enum {velocity, pressure} BoundaryType;

typedef double T;

T poiseuilleVelocity(int iY, UnitConverter<T> const& converter) {
    T y = (T)iY / converter.getN();
    return 4.*converter.getU() * (y-y*y);
}

T poiseuillePressure(int iX, UnitConverter<T> const& converter) {
    T Lx = converter.getNx()-1;
    T Ly = converter.getNy()-1;
    return 8.*converter.getNu()*converter.getU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}


void iniGeometry( BlockLattice2D<T, D2Q9Descriptor>& lattice,
                  UnitConverter<T> const& converter,
                  Dynamics<T, D2Q9Descriptor>& bulkDynamics,
                  BoundaryCondition2D<T,D2Q9Descriptor>& boundaryCondition,
                  BoundaryType boundaryType )
{
    int nx = converter.getNx();
    int ny = converter.getNy();
    T   omega = converter.getOmega();
    lattice.defineDynamics(0,nx-1, 0,ny-1, &bulkDynamics);

    boundaryCondition.addVelocityBoundary1P(1,nx-2,ny-1,ny-1, omega);
    boundaryCondition.addVelocityBoundary1N(1,nx-2,   0,   0, omega);

    if (boundaryType==velocity) {
        boundaryCondition.addVelocityBoundary0N(0,0, 1, ny-2, omega);
        boundaryCondition.addVelocityBoundary0P(nx-1,nx-1, 1, ny-2, omega);
    }
    else {
        boundaryCondition.addPressureBoundary0N(0,0, 1, ny-2, omega);
        boundaryCondition.addPressureBoundary0P(nx-1,nx-1, 1, ny-2, omega);
    }

    boundaryCondition.addExternalVelocityCornerNN(0,0, omega);
    boundaryCondition.addExternalVelocityCornerNP(0,ny-1, omega);
    boundaryCondition.addExternalVelocityCornerPN(nx-1,0, omega);
    boundaryCondition.addExternalVelocityCornerPP(nx-1,ny-1, omega);

    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            T u[2] = {poiseuilleVelocity(iY, converter),0.};
            T rho = (T)1 + poiseuillePressure(iX, converter) *
                           D2Q9Descriptor<T>::invCs2;
                lattice.get(iX,iY).defineRhoU(rho, u);
                lattice.get(iX,iY).iniEquilibrium(rho, u);
        }
    }

    lattice.initialize();
}

void plotStatistics(int iT,
                    BlockLattice2D<T, D2Q9Descriptor>&    lattice,
                    UnitConverter<T>& converter,
                    BlockStatistics2D<T, D2Q9Descriptor>& statistics,
                    T                                     middleU)
{
    T dx = converter.getDeltaX();
    T dt = converter.getDeltaT();
    cout << "Iteration " << setw(5) << iT;
    cout << "; t=" << setprecision(3) << setw(6)
         << iT*converter.getDeltaT();
    cout << "; E=" << setprecision(10) << setw(15)
         << lattice.getStatistics().getAverageEnergy()*dx*dx/dt/dt;
    cout << "; rho=" << setprecision(8) << setw(11)
         << lattice.getStatistics().getAverageRho();
    cout << "; uMax= " << setprecision(8) << setw(11)
         << middleU*dx/dt;
    cout << endl;

    ImageCreator<T> imageCreator("jet.map");
    imageCreator.writeScaledGif(createFileName("p", iT, 6),
                                statistics.getPressure());
    imageCreator.writeScaledGif(createFileName("u", iT, 6),
                                statistics.getVelocityNorm());
    statistics.reset();
}

void writeProfile ( string fName,
                    BlockLattice2D<T, D2Q9Descriptor>& lattice,
                    UnitConverter<T>& converter )
{
    ofstream ofile((singleton::directories().getLogOutDir()+fName).c_str());
    for (int iY=0; iY<converter.getNy(); ++iY) {
        T dx = converter.getDeltaX();
        T dt = converter.getDeltaT();
        T analytical = poiseuilleVelocity(iY, converter);
        T numerical[2];
        lattice.get(converter.getNx()/2, iY).computeU(numerical);
        ofile << iY*dx << " " << analytical*dx/dt
                       << " " << numerical[0]*dx/dt << "\n";
    }
}

int main() {
    singleton::directories().setOlbDir("../../../../");
    singleton::directories().setOutputDir("./tmp/");

    const BoundaryType boundaryType = pressure;
    const T uMax = 0.02;
    const T Re   = 10.;
    const int N  = 20;

    const T lx   = 2.;
    const T ly   = 1.;

    const int maxIter  = 20000;
    const int saveIter = 500;

    UnitConverter<T> converter(uMax, Re, N, lx, ly);
    writeLogFile(converter, "2D Poiseuille flow");
    BlockLattice2D<T, D2Q9Descriptor> lattice(converter.getNx(),
                                              converter.getNy() );

    BoundaryCondition2D<T,D2Q9Descriptor>* boundaryCondition =
        createLocalBoundaryCondition2D(lattice);

    BGKdynamics<T, D2Q9Descriptor> bulkDynamics (
                      converter.getOmega(),
                      instances::getBulkMomenta<T,D2Q9Descriptor>(),
                      lattice.getStatistics()
    );

    BlockStatistics2D<T,D2Q9Descriptor> statistics(lattice);

    iniGeometry(lattice, converter,
                bulkDynamics, *boundaryCondition, boundaryType);

    for (int iT=0; iT<maxIter; ++iT) {
        T middleU[2];
        lattice.get( converter.getNx()/2,
                     converter.getNy()/2 ).computeU(middleU);

        lattice.collideAndStream();

        if (iT%saveIter==0) {
            plotStatistics(iT, lattice, converter, statistics, middleU[0]);
            writeProfile("centerVel.dat", lattice, converter);
        }
    }

    delete boundaryCondition;
}


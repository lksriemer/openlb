/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
#include <fstream>
#include "olb2D.h"
//#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

T poiseuilleVelocity(int iY, int N, T u) {
    T y = (T)iY / (T)N;
    return 4.*u * (y-y*y);
}


void iniGeometry( BlockLattice2D<T,D2Q9Descriptor>& lattice,
                  UnitConverter<T> const& converter,
                  Dynamics<T, D2Q9Descriptor>& bulkDynamics,
                  BoundaryCondition2D<T,D2Q9Descriptor>& boundaryCondition )
{
    const T lx1   = 5.0;
    const T ly1   = 0.75;
    const T omega = converter.getOmega();

    const int nx = converter.getNx();
    const int ny = converter.getNy();
    const int nx1 = converter.nCell(lx1);
    const int ny1 = converter.nCell(ly1);

    lattice.defineDynamics(0,nx-1,0,ny-1,    &bulkDynamics);
    lattice.defineDynamics(0,nx1-1,0,ny1-1,
                           &instances::getNoDynamics<T,D2Q9Descriptor>());

    boundaryCondition.addVelocityBoundary0N(0,0,        ny1+1,ny-2, omega);
    boundaryCondition.addVelocityBoundary0N(nx1,nx1,    1,ny1-1, omega);
    boundaryCondition.addVelocityBoundary0P(nx-1,nx-1,  1, ny-2, omega);
    boundaryCondition.addVelocityBoundary1P(1,nx-2,     ny-1,ny-1, omega);
    boundaryCondition.addVelocityBoundary1N(1,nx1-1,    ny1,ny1, omega);
    boundaryCondition.addVelocityBoundary1N(nx1+1,nx-2, 0,0, omega);

    boundaryCondition.addExternalVelocityCornerNN(0,ny1, omega);
    boundaryCondition.addExternalVelocityCornerNN(nx1,0, omega);

    boundaryCondition.addExternalVelocityCornerNP(0,ny-1, omega);
    boundaryCondition.addExternalVelocityCornerPN(nx-1,0, omega);
    boundaryCondition.addExternalVelocityCornerPP(nx-1,ny-1, omega);

    boundaryCondition.addInternalVelocityCornerNN(nx1,ny1, omega);

    for (int iX=0; iX<=nx1; ++iX) {
        for (int iY=ny1; iY<ny; ++iY) {
            T vel[] = {
                poiseuilleVelocity(iY-ny1, ny-ny1-1, converter.getU()),
                T()
            };
            lattice.get(iX,iY).defineRhoU((T)1, vel);
            lattice.get(iX,iY).iniEquilibrium((T)1, vel);
        }
    }

    for (int iX=nx1+1; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            T vel[] = {
                poiseuilleVelocity( iY, ny-1,
                                    converter.getU()*((T)1-(T)ny1/(T)ny) ),
                T()
            };
            lattice.get(iX,iY).defineRhoU((T)1, vel);
            lattice.get(iX,iY).iniEquilibrium((T)1, vel);
        }
    }

    lattice.initialize();
}

void writeGifs(BlockLattice2D<T,D2Q9Descriptor>& lattice,
               UnitConverter<T> const& converter, int iter)
{
    const int imSize = 600;

    BlockLatticeView2D<T,D2Q9Descriptor> lView( lattice,
            converter.nCell(2.5), converter.nCell(20.),
            0, converter.getNy()-1 );
    BlockStatistics2D<T,D2Q9Descriptor> subStatistics(lView);
    ImageCreator<T> imageCreator("jet.map");
    imageCreator.writeScaledGif(createFileName("u", iter, 6),
                                subStatistics.getVelocityNorm(), imSize,imSize);
    imageCreator.writeScaledGif(createFileName("omega", iter, 6),
                                subStatistics.getVorticity(), imSize,imSize);
    //imageCreator.writeVtkData(createFileName("vtk", iter, 6),
    //                          subStatistics.getVorticity(),
    //                          subStatistics.getVelocity(),
    //                          converter.getDeltaX(), converter.getDeltaT() );

}


int main() {
    singleton::directories().setOlbDir("../../../../");
    singleton::directories().setOutputDir("./tmp/");

    UnitConverter<T> converter(
            (T) 2e-2,  // uMax
            (T) 500.,  // Re
            60,        // N
            40.,       // lx
            1.5        // ly
    );
    writeLogFile(converter, "2D Backward facing step");

    const T maxT        = (T)100000;
    const int  iterSave = 100;

    BlockLattice2D<T, D2Q9Descriptor> lattice(converter.getNx(),
                                              converter.getNy());

    ConstRhoBGKdynamics<T, D2Q9Descriptor> bulkDynamics (
                      converter.getOmega(),
                      instances::getBulkMomenta<T,D2Q9Descriptor>(),
                      lattice.getStatistics()
    );

    BoundaryCondition2D<T,D2Q9Descriptor>* boundaryCondition =
         createInterpBoundaryCondition2D(lattice);

    iniGeometry(lattice, converter, bulkDynamics, *boundaryCondition);


    for (int iT=0; iT*converter.getDeltaT()<maxT; ++iT) {
        if (iT%iterSave==0) {
            cout << iT << endl;
            cout << "step " << iT
                 << "; t=" << iT*converter.getDeltaT()
                 << "; av energy="
                 << lattice.getStatistics().getAverageEnergy()
                 << "; av rho="
                 << lattice.getStatistics().getAverageRho() << endl;
            writeGifs(lattice, converter, iT);
        }

        lattice.collideAndStream();
    }

    delete boundaryCondition;
}

/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Mathias Krause, Jonas Latt, Vincent Heuveline
 *  Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
 *  E-mail: Jonas.Latt@cui.unige.ch
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


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

void iniGeometry( BlockLattice2D<T,D2Q9Descriptor>& lattice,
                  UnitConverter<T> const& converter,
                  Dynamics<T, D2Q9Descriptor>& bulkDynamics,
                  BoundaryCondition2D<T,D2Q9Descriptor>& bc )
{
    const int nx = converter.getNx();
    const int ny = converter.getNy();

    const T omega = converter.getOmega();

    lattice.defineDynamics(0,nx-1, 0,ny-1, &bulkDynamics);

    bc.addVelocityBoundary0N(   0,   0,   1,ny-2, omega);
    bc.addVelocityBoundary0P(nx-1,nx-1,   1,ny-2, omega);
    bc.addVelocityBoundary1N(   1,nx-2,   0,   0, omega);
    bc.addVelocityBoundary1P(   1,nx-2,ny-1,ny-1, omega);

    bc.addExternalVelocityCornerNN(   0,   0, omega);
    bc.addExternalVelocityCornerNP(   0,ny-1, omega);
    bc.addExternalVelocityCornerPN(nx-1,   0, omega);
    bc.addExternalVelocityCornerPP(nx-1,ny-1, omega);

    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            T vel[] = { T(), T()};
            lattice.get(iX,iY).defineRhoU((T)1, vel);
            lattice.get(iX,iY).iniEquilibrium((T)1, vel);
        }
    }

    for (int iX=1; iX<nx-1; ++iX) {
        T u = converter.getU();
        T vel[] = { u, T() };
        lattice.get(iX,ny-1).defineRhoU((T)1, vel);
        lattice.get(iX,ny-1).iniEquilibrium((T)1, vel);
    }

    lattice.initialize();
}

void writeGifs(BlockLattice2D<T,D2Q9Descriptor>& lattice,
               BlockStatistics2D<T,D2Q9Descriptor>& statistics,
               UnitConverter<T> const& converter, int iter)
{
    const int imSize = 600;

    ImageCreator<T> imageCreator("jet.map");
    imageCreator.writeScaledGif(createFileName("uz", iter, 6),
                                statistics.getVelocityNorm(),
                                imSize, imSize );
    statistics.reset();
}

void writeVTK(BlockStatistics2D<T,D2Q9Descriptor>& statistics,
              UnitConverter<T> const& converter, int iter)
{
    VTKOut2D<T>:: writeFlowField (
            createFileName("vtk", iter, 6),
            "Pressure", statistics.getPressure(),
            "Velocity", statistics.getVelocity(),
            converter.getDeltaX(), converter.getDeltaT() );
   statistics.reset();
}


int main() {

    #ifdef PARALLEL_MODE_OMP
        #pragma omp parallel
            omp.init();
    #endif

    singleton::directories().setOlbDir("../../../../");
    singleton::directories().setOutputDir("./tmp/");

    UnitConverter<T> converter(
            (T) 1e-2,  // uMax
            (T) 1000.,  // Re
            128,        // N
            1.,        // lx
            1.         // ly 
    );
    const T logT     = (T)0.01;
    const T imSave   = (T)10.;
    const T vtkSave  = (T)10.;
    const T maxT     = (T)1.;

    writeLogFile(converter, "2D cavity");

    BlockLattice2D<T, D2Q9Descriptor> lattice(converter.getNx(),
                                               converter.getNy());

    ConstRhoBGKdynamics<T, D2Q9Descriptor> bulkDynamics (
                      converter.getOmega(),
                      instances::getBulkMomenta<T,D2Q9Descriptor>(),
                      lattice.getStatistics()
    );

    BoundaryCondition2D<T,D2Q9Descriptor>* boundaryCondition
        = createLocalBoundaryCondition2D(lattice);

    iniGeometry(lattice, converter, bulkDynamics, *boundaryCondition);

    BlockStatistics2D<T,D2Q9Descriptor> statistics(lattice);


    int iT=0;
    for (iT=0; iT*converter.getDeltaT()<maxT; ++iT) {
        if (iT%converter.nStep(logT)==0) {
            cout << "step " << iT
                 << "; t=" << iT*converter.getDeltaT()
                 << "; av energy="
                 << lattice.getStatistics().getAverageEnergy()
                 << "; av rho="
                 << lattice.getStatistics().getAverageRho() << endl;
        }

        if (iT%converter.nStep(imSave)==0 && iT>0) {
            cout << "Saving Gif ..." << endl;
            writeGifs(lattice, statistics, converter, iT);
        }

        if (iT%converter.nStep(vtkSave)==0 && iT>0) {
            cout << "Saving VTK file ..." << endl;
            writeVTK(statistics, converter, iT);
        }

        lattice.collideAndStream();

        //lattice.collide();
        //lattice.stream();

    }
    cout << iT << endl;
    cout << lattice.getStatistics().getAverageEnergy() << endl;

    delete boundaryCondition;

   return 0;
}

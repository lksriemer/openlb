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
#include "olb3D.h"
//#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

void iniGeometry( BlockLattice3D<T,D3Q19Descriptor>& lattice,
                  UnitConverter<T> const& converter,
                  Dynamics<T, D3Q19Descriptor>& bulkDynamics,
                  BoundaryCondition3D<T,D3Q19Descriptor>& bc )
{
    const T lx1   = 5.0;
    const T ly1   = 0.75;
    const T omega = converter.getOmega();

    const int nx = converter.getNx();
    const int ny = converter.getNy();
    const int nz = converter.getNz();
    const int nx1 = converter.nCell(lx1);
    const int ny1 = converter.nCell(ly1);

    lattice.defineDynamics(0,nx-1, 0,ny-1, 0,nz-1, &bulkDynamics);
    lattice.defineDynamics(0,nx1-1,0,ny1-1,0,nz-1,
                           &instances::getNoDynamics<T,D3Q19Descriptor>());

    bc.addVelocityBoundary0N(    0,    0,ny1+1, ny-2,   1,nz-2, omega);
    bc.addVelocityBoundary0N(  nx1,  nx1,    1,ny1-1,   1,nz-2, omega);
    bc.addVelocityBoundary0P(nx-1 , nx-1,    1, ny-2,   1,nz-2, omega);
    bc.addVelocityBoundary1N(    1,nx1-1,  ny1,  ny1,   1,nz-2, omega);
    bc.addVelocityBoundary1N(nx1+1, nx-2,    0,    0,   1,nz-2, omega);
    bc.addVelocityBoundary1P(    1, nx-2, ny-1, ny-1,   1,nz-2, omega);
    bc.addVelocityBoundary2N(    1, nx-2,ny1+1, ny-2,   0,   0, omega);
    bc.addVelocityBoundary2N(nx1+1, nx-2,    1,  ny1,   0,   0, omega);
    bc.addVelocityBoundary2P(    1, nx-2,ny1+1, ny-2,nz-1,nz-1, omega);
    bc.addVelocityBoundary2P(nx1+1, nx-2,    1,  ny1,nz-1,nz-1, omega);
    
    bc.addExternalVelocityEdge0NN(    1,nx1-1,  ny1,  ny1,    0,    0, omega);
    bc.addExternalVelocityEdge0NN(nx1+1, nx-2,    0,    0,    0,    0, omega);
    bc.addExternalVelocityEdge0NP(    1,nx1-1,  ny1,  ny1, nz-1, nz-1, omega);
    bc.addExternalVelocityEdge0NP(nx1+1, nx-2,    0,    0, nz-1, nz-1, omega);
    bc.addExternalVelocityEdge0PN(    1, nx-2, ny-1, ny-1,    0,    0, omega);
    bc.addExternalVelocityEdge0PP(    1, nx-2, ny-1, ny-1, nz-1, nz-1, omega);
 
    bc.addExternalVelocityEdge1NN(    0,    0,ny1+1, ny-2,    0,    0, omega);
    bc.addExternalVelocityEdge1NN(  nx1,  nx1,    1,ny1-1,    0,    0, omega);
    bc.addExternalVelocityEdge1NP( nx-1, nx-1,    1, ny-2,    0,    0, omega);
    bc.addExternalVelocityEdge1PN(    0,    0,ny1+1, ny-2, nz-1, nz-1, omega);
    bc.addExternalVelocityEdge1PN(  nx1,  nx1,    1,ny1-1, nz-1, nz-1, omega);
    bc.addExternalVelocityEdge1PP( nx-1, nx-1,    1, ny-2, nz-1, nz-1, omega);

    bc.addExternalVelocityEdge2NN(   0,   0, ny1, ny1,   1,nz-2, omega);
    bc.addExternalVelocityEdge2NN( nx1, nx1,   0,   0,   1,nz-2, omega);
    bc.addExternalVelocityEdge2NP(   0,   0,ny-1,ny-1,   1,nz-2, omega);
    bc.addExternalVelocityEdge2PN(nx-1,nx-1,   0,   0,   1,nz-2, omega);
    bc.addExternalVelocityEdge2PP(nx-1,nx-1,ny-1,ny-1,   1,nz-2, omega);

    bc.addExternalVelocityCornerNNN(   0, ny1,   0, omega);
    bc.addExternalVelocityCornerNNN( nx1,   0,   0, omega);
    bc.addExternalVelocityCornerNNP(   0, ny1,nz-1, omega);
    bc.addExternalVelocityCornerNNP( nx1,   0,nz-1, omega);
    bc.addExternalVelocityCornerNPN(   0,ny-1,   0, omega);
    bc.addExternalVelocityCornerNPN( nx1, ny1,   0, omega);
    bc.addExternalVelocityCornerNPP(   0,ny-1,nz-1, omega);
    bc.addExternalVelocityCornerNPP( nx1, ny1,nz-1, omega);
    bc.addExternalVelocityCornerPNN(nx-1,   0,   0, omega);
    bc.addExternalVelocityCornerPNP(nx-1,   0,nz-1, omega);
    bc.addExternalVelocityCornerPPN(nx-1,ny-1,   0, omega);
    bc.addExternalVelocityCornerPPP(nx-1,ny-1,nz-1, omega);

    bc.addInternalVelocityEdge2NN(  nx1,  nx1,  ny1,  ny1,  1, nz-2, omega);

    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            for (int iZ=0; iZ<nz; ++iZ) {
                T vel[] = { T(), T(), T() };
                lattice.get(iX,iY,iZ).defineRhoU((T)1, vel);
                lattice.get(iX,iY,iZ).iniEquilibrium((T)1, vel);
            }
        }
    }

    for (int iX=0; iX<=nx1; ++iX) {
        for (int iY=ny1+1; iY<ny-1; ++iY) {
            for (int iZ=1; iZ<nz-1; ++iZ) {
                T vel[] = { converter.getU(), T(), T() };
                lattice.get(iX,iY,iZ).defineRhoU((T)1, vel);
                lattice.get(iX,iY,iZ).iniEquilibrium((T)1, vel);
            }
        }
    }

    for (int iX=nx1+1; iX<nx; ++iX) {
        for (int iY=1; iY<ny-1; ++iY) {
            for (int iZ=1; iZ<nz-1; ++iZ) {
                T u = converter.getU()*ly1/converter.getLy();
                T vel[] = { u, T(), T() };
                lattice.get(iX,iY,iZ).defineRhoU((T)1, vel);
                lattice.get(iX,iY,iZ).iniEquilibrium((T)1, vel);
            }
        }
    }

    lattice.initialize();

}

void writeGifs(BlockLattice3D<T,D3Q19Descriptor>& lattice,
               UnitConverter<T> const& converter, int iter)
{
    const int imSize = 600;
    const int nx = converter.getNx();
    const int ny = converter.getNy();
    const int nz = converter.getNz();

    BlockLatticeView3D<T,D3Q19Descriptor> lView( lattice,
            converter.nCell(2.5), converter.nCell(18.),
            0, ny-1, 0, nx-1 );
    BlockStatistics3D<T,D3Q19Descriptor> subStatistics(lView);
    ImageCreator<T> imageCreator("jet.map");

    //imageCreator.writeScaledGif(createFileName("pz", iter, 6),
    //                            subStatistics.getPressure().sliceZ(nz/2),
    //                            imSize, imSize);
    imageCreator.writeScaledGif(createFileName("uz", iter, 6),
                                subStatistics.getVelocityNorm().sliceZ(nz/2),
                                imSize, imSize );
}

void writeVTK(BlockLattice3D<T,D3Q19Descriptor>& lattice,
              UnitConverter<T> const& converter, int iter)
{
    BlockLatticeView3D<T,D3Q19Descriptor> lView( lattice,
            converter.nCell(2.5), converter.nCell(18.),
            0, converter.getNy()-1,
            0, converter.getNz()-1 );
    BlockStatistics3D<T,D3Q19Descriptor> subStatistics(lView);

    VTKOut3D<T>:: writeFlowField (
            createFileName("vtk", iter, 6),
            "Vorticity", subStatistics.getVorticityNorm(),
            "Velocity", subStatistics.getVelocity(),
            converter.getDeltaX(), converter.getDeltaT() );
}

int main() {
    singleton::directories().setOlbDir("../../../../");
    singleton::directories().setOutputDir("./tmp/");

    UnitConverter<T> converter(
            (T) 4e-2,  // uMax
            (T) 100.,  // Re
            20,        // N
            30.,       // lx
            1.5,       // ly
            1.5        // lz
    );

    const T logT     = (T)1/(T)100;
    const T imSave   = (T)1/(T)10;
    const T vtkSave  = (T)10.;
    const T maxT     = (T)100.1;

    writeLogFile(converter, "3D backward facing step");

    BlockLattice3D<T, D3Q19Descriptor> lattice(converter.getNx(),
                                               converter.getNy(),
                                               converter.getNz());

    ConstRhoBGKdynamics<T, D3Q19Descriptor> bulkDynamics (
                      converter.getOmega(),
                      instances::getBulkMomenta<T,D3Q19Descriptor>(),
                      lattice.getStatistics()
    );

    BoundaryCondition3D<T,D3Q19Descriptor>* boundaryCondition
        = createLocalBoundaryCondition3D(lattice);

    iniGeometry(lattice, converter, bulkDynamics, *boundaryCondition);

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
            writeGifs(lattice, converter, iT);
        }

        if (iT%converter.nStep(vtkSave)==0 && iT>0) {
            cout << "Saving VTK file ..." << endl;
            writeVTK(lattice, converter, iT);
        }

        lattice.collideAndStream();
    }

    delete boundaryCondition;
}

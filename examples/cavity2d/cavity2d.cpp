/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006 - 2011 Mathias J. Krause, Jonas Fietz, Jonas Latt
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


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
  #include "olb2D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <io/xmlReader.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor


void iniGeometry( BlockStructure2D<T,DESCRIPTOR>& lattice,
                  LBunits<T> const& converter,
                  Dynamics<T, DESCRIPTOR>& bulkDynamics,
                  OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& bc )
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
        T u = converter.getLatticeU();
        T vel[] = { u, T() };
        lattice.get(iX,ny-1).defineRhoU((T)1, vel);
        lattice.get(iX,ny-1).iniEquilibrium((T)1, vel);
    }

    lattice.initialize();
}

void writeGifs(std::string filename, BlockStructure2D<T,DESCRIPTOR>& lattice,
               LBunits<T> const& converter, int iter)
{
    const int imSize = 600;

    DataAnalysisBase2D<T,DESCRIPTOR> const& analysis = lattice.getDataAnalysis();

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName(filename, iter, 6),
                               analysis.getVelocityNorm(),
                               imSize, imSize );
}

void writeVTK(std::string filename, BlockStructure2D<T,DESCRIPTOR>& lattice,
              LBunits<T> const& converter, int iter)
{
    DataAnalysisBase2D<T,DESCRIPTOR> const& analysis = lattice.getDataAnalysis();

    T dx = converter.getDeltaX();
    T dt = converter.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName(filename, iter, 6), dx);
    vtkOut.writeData<double>(analysis.getVorticity(), "vorticity", (T)1/dt);
    vtkOut.writeData<2,double>(analysis.getVelocity(), "velocity", dx/dt);
}


int main(int argc, char* argv[]) {
    olbInit(&argc, &argv);

    string fName("demo.xml");
    XMLreader config(fName);

    std::string olbdir, outputdir;
    config["Application"]["OlbDir"].read(olbdir);
    config["Application"]["OutputDir"].read(outputdir);
    singleton::directories().setOlbDir(olbdir);
    singleton::directories().setOutputDir(outputdir);

    LBunits<T> converter(
    	config["Application"]["MaxU"].get<T>(),
    	config["Application"]["Re"].get<T>(), 
    	config["Mesh"]["Refinement"].get<int>(),
        config["Mesh"]["lx"].get<T>(),
        config["Mesh"]["ly"].get<T>()
    );

    T logT;
    T imSave;
    T vtkSave;
    T maxT;

    config["Log"]["SaveTime"].read(logT);
    config["VisualizationImages"]["SaveTime"].read(imSave);
    config["VisualizationVTK"]["SaveTime"].read(vtkSave);
    config["Application"]["MaxTime"].read(maxT);


    writeLogFile(converter, config["Log"]["Filename"].get<std::string>()); 


#ifndef PARALLEL_MODE_MPI  // sequential program execution
    BlockLattice2D<T, DESCRIPTOR> lattice(converter.getNx(), converter.getNy() );
#else                      // parallel program execution
    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            createRegularDataDistribution( converter.getNx(), converter.getNy() ) );
#endif

    ConstRhoBGKdynamics<T, DESCRIPTOR> bulkDynamics (
            converter.getOmega(),
            instances::getBulkMomenta<T,DESCRIPTOR>()
            );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createInterpBoundaryCondition2D(lattice);

    iniGeometry(lattice, converter, bulkDynamics, *boundaryCondition);

    int iT=0;
    for (iT=0; iT*converter.getDeltaT()<maxT; ++iT) {
        if (iT%converter.nStep(logT)==0) {
            cout << "step " << iT
                << "; lattice time=" << lattice.getStatistics().getTime()
                << "; t=" << iT*converter.getDeltaT()
                << "; av energy="
                << lattice.getStatistics().getAverageEnergy()
                << "; av rho="
                << lattice.getStatistics().getAverageRho() << endl;
        }

        if (iT%converter.nStep(imSave)==0 && iT>0) {
            cout << "Saving Gif ..." << endl;
            writeGifs(config["VisualizationImages"]["Filename"].get<std::string>(), lattice, converter, iT);
        }
        if (iT%converter.nStep(vtkSave)==0 && iT>0) {
            cout << "Saving VTK file ..." << endl;
            writeVTK(config["VisualizationVTK"]["Filename"].get<std::string>(), lattice, converter, iT);
        }

        lattice.collideAndStream();
    }
    cout << iT << endl;
    cout << lattice.getStatistics().getAverageEnergy() << endl;

    delete boundaryCondition;
}

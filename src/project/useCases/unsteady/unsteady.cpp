/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
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

/* unsteady.c:
 * This example examines an unsteady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, where as the outlet implements an outflow condition:
 * grad_x p = 0. At Reynolds numbers around 100, an unstable periodic pattern is
 * created, the Karman vortex street.
 */


#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

// Parameters for the simulation setup
const int nx     = 336;
const int ny     = 64;
const T uMax     = 0.02;
const T Re       = 100.;
const int obst_x = nx/5-1;
const int obst_y = ny/2-1;
const int obst_r = ny/10+1;
const T nu       = uMax * 2.*obst_r / Re;
const T omega    = 1. / (3.*nu+1./2.);
const int maxT   = 50000;
const int tSave  = 100;

// Initialize a nx-by-ny lattice
BlockLattice2D<T, D2Q9Descriptor> lattice(nx,ny);

// Initialize an object that describes the BGK
// dynamics in the bulk
ConstRhoBGKdynamics<T, D2Q9Descriptor> bulkDynamics (
                  omega,
                  instances::getBulkMomenta<T,D2Q9Descriptor>(),
                  lattice.getStatistics()
);

// Computation of a Poiseuille velocity profile
T poiseuilleVelocity(int iY) {
    T y = (T)iY;
    T L = (T)(ny-1);
    return 4.*uMax / (L*L) * (L*y-y*y);
}

// Computation of a Poiseuille pressure profile
T poiseuillePressure(int iX) {
    T x = (T)iX;
    T Lx = (T)(nx-1);
    T Ly = (T)(ny-1);
    return 8.*nu*uMax / (Ly*Ly) * (Lx/2.-x);
}

// Set up the geometry of the simulation
void iniGeometry(BoundaryCondition2D<T,D2Q9Descriptor>& boundaryCondition) {
    // Bulk dynamics
    lattice.defineDynamics(0,nx-1,0,ny-1,    &bulkDynamics);

    // top boundary
    boundaryCondition.addVelocityBoundary1P(1,nx-2,ny-1,ny-1, omega); 
    // bottom boundary
    boundaryCondition.addVelocityBoundary1N(1,nx-2,   0,   0, omega);
    // left boundary
    boundaryCondition.addVelocityBoundary0N(0,0, 1, ny-2, omega);
    // right boundary
    // Note that the right boundary implements a Neumann condition
    // on the pressure. Try to set up a VelocityBoundary
    // (only this line needs to be changed) to get
    // a Neumann condition of the velocity instead
    boundaryCondition.addPressureBoundary0P(nx-1,nx-1, 1, ny-2, omega);

    // Corner nodes
    boundaryCondition.addExternalVelocityCornerNN(0,0, omega);
    boundaryCondition.addExternalVelocityCornerNP(0,ny-1, omega);
    boundaryCondition.addExternalVelocityCornerPN(nx-1,0, omega);
    boundaryCondition.addExternalVelocityCornerPP(nx-1,ny-1, omega);

    // Definition of the obstacle (bounce-back nodes)
    for (int iX=0; iX<nx; ++iX) {
        for (int iY=0; iY<ny; ++iY) {
            T u[2] = {poiseuilleVelocity(iY),0.};
            T rho = (T)1 + poiseuillePressure(iX) *
                               D2Q9Descriptor<T>::invCs2;
            if ( (iX-obst_x)*(iX-obst_x) +
                 (iY-obst_y)*(iY-obst_y) <= obst_r*obst_r )
            {
                lattice.defineDynamics( iX,iX,iY,iY,
                        &instances::getBounceBack<T,D2Q9Descriptor>() );
            }
            else {
                lattice.get(iX,iY).defineRhoU(rho, u);
                lattice.get(iX,iY).iniEquilibrium(rho, u);
            }
        }
    }

    // Make the lattice ready for simulation
    lattice.initialize();
}

int main() {
    singleton::directories().setOlbDir("../../../../");
    singleton::directories().setOutputDir("./tmp/");

    // Initialize an object that produces the boundary condition.
    // createLocalBoundaryCondition2D: local, regularized boundary condition
    // createInterpBoundaryCondition2D: non-local boundary, based on an
    //                                  interpolation of the stress tensor.

    BoundaryCondition2D<T,D2Q9Descriptor>* boundaryCondition =
        createInterpBoundaryCondition2D(lattice);


    iniGeometry(*boundaryCondition);

    // Computation of simulation results (e.g. the velocity field)
    BlockStatistics2D<T,D2Q9Descriptor> statistics(lattice);
    // Creation of images representing intermediate simulation results
    ImageCreator<T> imageCreator("jet.map");

    // Main loop over time
    for (int iT=0; iT<maxT; ++iT) {
        if (iT%tSave==0 && iT !=0) {
            cout << iT << endl;
            cout << "av energy="
                 << lattice.getStatistics().getAverageEnergy()
                 << "; av rho="
                 << lattice.getStatistics().getAverageRho() << endl;

            // Creation of gif images. This works only on systems on
            // which ImageMagick is installed. If you have the
            // program gifmerge, you can create an animation with the help
            // of a command of the kind "gifmerge -5 u*.gif > anim_u"
            imageCreator.writeScaledGif(createFileName("p", iT, 6),
                                        statistics.getPressure());
            imageCreator.writeScaledGif(createFileName("u", iT, 6),
                                        statistics.getVelocityNorm());
            statistics.reset();
        }

        lattice.collideAndStream();

        // Neumann (zero-gradient) condition for the right boundary:
        // x-component of the velocity and pressure are copied from
        // the previous cell.
        if (iT>0) {
            for (int iY=1; iY<ny-1; ++iY) {
                T rho2, u2[2];
                lattice.get(nx-2,iY).computeRhoU(rho2, u2);
                u2[1] = T(); // Outflow is orthogonal to boundary
                lattice.get(nx-1,iY).defineRhoU(rho2, u2);
            }
        }
    }

    delete boundaryCondition;
}

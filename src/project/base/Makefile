# This file is part of the OpenLB library
#
# Copyright (C) 2006, 2007 Jonas Latt
# Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
# E-mail: jonas.latt@gmail.com
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public 
# License along with this program; if not, write to the Free 
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.

## This file is used to compile all .cpp source files to .o object
## files for the precompiled version of the library.
## By default, invoking the command "make" in the current directory
## lists all cpp files.

OLB_BASE_PATH = .

include dependencies.mk

default:
	echo $(DEP__BOUNDARIES_2D_CPP)
	echo $(DEP__BOUNDARIES_3D_CPP)

$(OLB_BASE_PATH)/loadBalancer.o: \
     $(OLB_BASE_PATH)/loadBalancer.cpp \
     $(DEP__LOADBALANCER_CPP)

$(OLB_BASE_PATH)/ompManager.o: \
     $(OLB_BASE_PATH)/ompManager.cpp \
     $(DEP__OMPMANAGER_CPP)

$(OLB_BASE_PATH)/latticeDescriptors.o: \
     $(OLB_BASE_PATH)/latticeDescriptors.cpp \
     $(DEP__LATTICE_DESCRIPTORS_CPP)

$(OLB_BASE_PATH)/dynamics.o: \
     $(OLB_BASE_PATH)/dynamics.cpp \
     $(DEP__DYNAMICS_CPP)

$(OLB_BASE_PATH)/cell.o: \
     $(OLB_BASE_PATH)/cell.cpp \
     $(DEP__CELL_CPP)

$(OLB_BASE_PATH)/imageCreator.o: \
     $(OLB_BASE_PATH)/imageCreator.cpp \
     $(DEP__IMAGE_CREATOR_CPP)

$(OLB_BASE_PATH)/dataFields2D.o: \
     $(OLB_BASE_PATH)/dataFields2D.cpp \
     $(DEP__DATA_FIELDS_2D_CPP)

$(OLB_BASE_PATH)/blockLattice2D.o: \
     $(OLB_BASE_PATH)/blockLattice2D.cpp \
     $(DEP__BLOCK_LATTICE_2D_CPP)

$(OLB_BASE_PATH)/blockLatticeView2D.o: \
     $(OLB_BASE_PATH)/blockLatticeView2D.cpp \
     $(DEP__BLOCK_STRUCTURE_2D_CPP)

$(OLB_BASE_PATH)/blockStatistics2D.o: \
     $(OLB_BASE_PATH)/blockStatistics2D.cpp \
     $(DEP__BLOCK_STATISTICS_2D_CPP)

$(OLB_BASE_PATH)/boundaries2D.o: \
     $(OLB_BASE_PATH)/boundaries2D.cpp \
     $(DEP__BOUNDARIES_2D_CPP)

$(OLB_BASE_PATH)/fdBoundaries2D.o: \
     $(OLB_BASE_PATH)/fdBoundaries2D.cpp \
     $(DEP__FD_BOUNDARIES_2D_CPP)

$(OLB_BASE_PATH)/dynamicsHelpers2D.o: \
     $(OLB_BASE_PATH)/dynamicsHelpers2D.cpp \
     $(DEP__DYNAMICS_HELPERS_2D_CPP)

$(OLB_BASE_PATH)/simulationSetup2D.o: \
     $(OLB_BASE_PATH)/simulationSetup2D.cpp \
     $(DEP__SIMULATION_SETUP_2D_CPP)

$(OLB_BASE_PATH)/postProcessing2D.o: \
     $(OLB_BASE_PATH)/postProcessing2D.cpp \
     $(DEP__POST_PROCESSING_2D_CPP)

$(OLB_BASE_PATH)/dataFields3D.o: \
     $(OLB_BASE_PATH)/dataFields3D.cpp \
     $(DEP__DATA_FIELDS_3D_CPP)

$(OLB_BASE_PATH)/blockLattice3D.o: \
     $(OLB_BASE_PATH)/blockLattice3D.cpp \
     $(DEP__BLOCK_LATTICE_3D_CPP)

$(OLB_BASE_PATH)/blockLatticeView3D.o: \
     $(OLB_BASE_PATH)/blockLatticeView3D.cpp \
     $(DEP__BLOCK_STRUCTURE_3D_CPP)

$(OLB_BASE_PATH)/blockStatistics3D.o: \
     $(OLB_BASE_PATH)/blockStatistics3D.cpp \
     $(DEP__BLOCK_STATISTICS_3D_CPP)

$(OLB_BASE_PATH)/boundaries3D.o: \
     $(OLB_BASE_PATH)/boundaries3D.cpp \
     $(DEP__BOUNDARIES_3D_CPP)

$(OLB_BASE_PATH)/fdBoundaries3D.o: \
     $(OLB_BASE_PATH)/fdBoundaries3D.cpp \
     $(DEP__FD_BOUNDARIES_3D_CPP)

$(OLB_BASE_PATH)/dynamicsHelpers3D.o: \
     $(OLB_BASE_PATH)/dynamicsHelpers3D.cpp \
     $(DEP__DYNAMICS_HELPERS_3D_CPP)

$(OLB_BASE_PATH)/simulationSetup3D.o: \
     $(OLB_BASE_PATH)/simulationSetup3D.cpp \
     $(DEP__SIMULATION_SETUP_3D_CPP)

$(OLB_BASE_PATH)/postProcessing3D.o: \
     $(OLB_BASE_PATH)/postProcessing3D.cpp \
     $(DEP__POST_PROCESSING_3D_CPP)


clean:
	/bin/rm -f *.o *~

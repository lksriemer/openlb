# This file is part of the OpenLB library
#
# Copyright (C) 2006 Jonas Latt, Vincent Heuveline
# Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
# E-mail: Jonas.Latt@cui.unige.ch
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

## This file contains the compilation directives for all .cpp
## source files. They are required only when the precompiled
## version of the library is created.

BASIC_OBJECTS = \
     $(OLB_BASE_PATH)/latticeDescriptors.o \
     $(OLB_BASE_PATH)/dynamics.o \
     $(OLB_BASE_PATH)/cell.o \
     $(OLB_BASE_PATH)/imageCreator.o \
     $(OLB_BASE_PATH)/dataFields2D.o \
     $(OLB_BASE_PATH)/ompManager.o \
     $(OLB_BASE_PATH)/loadBalancer.o

OBJECTS_2D = \
     $(OLB_BASE_PATH)/blockLattice2D.o \
     $(OLB_BASE_PATH)/blockLatticeView2D.o \
     $(OLB_BASE_PATH)/blockStatistics2D.o \
     $(OLB_BASE_PATH)/boundaries2D.o \
     $(OLB_BASE_PATH)/fdBoundaries2D.o \
     $(OLB_BASE_PATH)/dynamicsHelpers2D.o \
     $(OLB_BASE_PATH)/simulationSetup2D.o \
     $(OLB_BASE_PATH)/postProcessing2D.o

OBJECTS_3D = \
     $(OLB_BASE_PATH)/dataFields3D.o \
     $(OLB_BASE_PATH)/blockLattice3D.o \
     $(OLB_BASE_PATH)/blockLatticeView3D.o \
     $(OLB_BASE_PATH)/blockStatistics3D.o \
     $(OLB_BASE_PATH)/boundaries3D.o \
     $(OLB_BASE_PATH)/fdBoundaries3D.o \
     $(OLB_BASE_PATH)/dynamicsHelpers3D.o \
     $(OLB_BASE_PATH)/simulationSetup3D.o \
     $(OLB_BASE_PATH)/postProcessing3D.o

ALL_OBJECTS = $(BASIC_OBJECTS) $(OBJECTS_2D) $(OBJECTS_3D)

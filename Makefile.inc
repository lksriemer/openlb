#  This file is part of the OpenLB library
#
#  Copyright (C) 2007 Mathias Krause
#  E-mail contact: info@openlb.net
#  The most recent release of OpenLB can be downloaded at
#  <http://www.openlb.net/>
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public 
#  License along with this program; if not, write to the Free 
#  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.

###########################################################################
###########################################################################
## DEFINITIONS TO BE CHANGED

CXX             := g++
#CXX             := icpc -D__aligned__=ignored
#CXX             := mpiCC
#CXX             := mpic++

OPTIM           := -O3 -Wall
DEBUG           := -g -DOLB_DEBUG

CXXFLAGS        := $(OPTIM)
#CXXFLAGS        := $(DEBUG)

#CXXFLAGS        += -fdiagnostics-color=auto
#CXXFLAGS        += -std=c++0x                  

ARPRG           := ar
#ARPRG           := xiar                  # mandatory for intel compiler

LDFLAGS         :=

PARALLEL_MODE   := OFF
#PARALLEL_MODE   := MPI
#PARALLEL_MODE   := OMP
#PARALLEL_MODE   := HYBRID

MPIFLAGS        :=
OMPFLAGS        := -fopenmp

BUILDTYPE       := precompiled
#BUILDTYPE       := generic


###########################################################################
## conditional settings

ifeq ($(BUILDTYPE), precompiled)
   CXXFLAGS := -DOLB_PRECOMPILED $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), MPI)
   CXXFLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS) $(CXXFLAGS)
endif

ifeq ($(PARALLEL_MODE), OMP)
   CXXFLAGS := -DPARALLEL_MODE_OMP $(OMPFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

ifeq ($(PARALLEL_MODE), HYBRID)
   CXXFLAGS := -DPARALLEL_MODE_OMP -DPARALLEL_MODE_MPI $(OMPFLAGS) $(MPIFLAGS) $(CXXFLAGS)
   LDFLAGS  := $(OMPFLAGS) $(LDFLAGS)
endif

###########################################################################
## defines shell

SHELL           := /bin/sh

###########################################################################
## dependencies, object, library directory and library name

DEPENDDIR       := build/$(BUILDTYPE)/dep
OBJDIR          := build/$(BUILDTYPE)/obj
LIBDIR          := build/$(BUILDTYPE)/lib
LIB             := olb

###########################################################################
## search directories

SUBDIRS         := src/boundary \
                   src/communication \
                   src/dynamics \
                   src/core \
                   src/geometry \
                   src/external/tinyxml \
                   src/functors \
                   src/io \
                   src/utilities

EXAMPLEDIRS     := examples/aorta3d \
                   examples/bstep2d \
                   examples/bstep3d \
                   examples/cavity2d/sequential \
                   examples/cavity2d/parallel \
                   examples/cavity3d/sequential \
                   examples/cavity3d/parallel \
                   examples/cylinder2d \
                   examples/cylinder3d \
                   examples/d3q13cavity \
                   examples/multiComponent2d \
                   examples/multiComponent3d \
                   examples/nozzle3d \
                   examples/phaseSeparation2d \
                   examples/phaseSeparation3d \
                   examples/poiseuille2d/bgk \
                   examples/poiseuille2d/forced \
                   examples/poiseuille2d/mrt \
                   examples/thermal2d \
                   examples/thermal3d \
                   examples/venturi3d

INCLUDEDIRS     := src \
                   src/ \
                   src/external
			
BUILDTYPEDIRS   := build/precompiled \
                   build/generic

SUBDIRSH        := $(foreach d,$(SUBDIRS),$(ROOT)/$(d))
IDIR            := $(foreach d,$(INCLUDEDIRS),-I$(ROOT)/$(d))


###########################################################################
###########################################################################

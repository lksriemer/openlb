# This file is part of the OpenLB library
#
# Copyright (C) 2006 Jonas Latt
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

## Makefile for bstep.cpp, generic version.
## In this version, all OpenLB source code is recompiled when
## the source bstep.cpp is changed. For this to work, the source file
## must #include the file olb2D.h, as well as olb2D.hh.
## If you want to restrict yourself to the D2Q9 lattice and speed up
## the compilation process, use the other Makefile, the one for
## precompiled code.

default: bstep

CXX      = g++
OPTIM    = -O3
DEBUG    = -g -DOLB_DEBUG
OLB_BASE_PATH = ../../base
OLB_ADDON_PATH = ../../addon

COMPFLAGS = $(OPTIM)
#COMPFLAGS = $(DEBUG)
CXXFLAGS  = -I$(OLB_BASE_PATH) -I$(OLB_ADDON_PATH) $(COMPFLAGS)
LOPT      = #$(DEBUG)

include $(OLB_BASE_PATH)/dependencies.mk
include $(OLB_BASE_PATH)/compilation.mk

include $(OLB_ADDON_PATH)/dependencies.mk
include $(OLB_ADDON_PATH)/compilation.mk

bstep: bstep.cpp $(DEP__OLB_2D_H) $(DEP__OLB_2D_HH)
	$(CXX) $(CXXFLAGS) $(LOPT) -o bstep bstep.cpp


clean:
	rm -f *~ *.o *.exe bstep; $(MAKE) clean_addon; $(MAKE) clean_base

clean_addon:
	cd $(OLB_ADDON_PATH); $(MAKE) clean;

clean_base:
	cd $(OLB_BASE_PATH); $(MAKE) clean;

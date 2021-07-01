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

## This file specifies all compilation dependencies between
## the source files of the OpenLB project


################  .h files #########################################

DEP__LOADBALANCER_H = $(OLB_BASE_PATH)/loadBalancer.h

DEP__OMPMANAGER_H = $(OLB_BASE_PATH)/ompManager.h

DEP__OLB_DEBUG_H = $(OLB_BASE_PATH)/olbDebug.h

DEP__SINGLETON_H = $(OLB_BASE_PATH)/singleton.h

DEP__LATTICE_DESCRIPTORS_H = $(OLB_BASE_PATH)/latticeDescriptors.h

DEP__UNITS_H = \
     $(OLB_BASE_PATH)/units.h \
     $(DEP__SINGLETON_H)

DEP__UTIL_H = $(OLB_BASE_PATH)/util.h

DEP__DYNAMICS_H = \
     $(OLB_BASE_PATH)/dynamics.h \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__UTIL_H) \
     $(DEP__POST_PROCESSING_H)

DEP__CELL_H = \
     $(OLB_BASE_PATH)/cell.h \
     $(DEP__OLB_DEBUG_H) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__DYNAMICS_H)

DEP__POST_PROCESSING_H = $(OLB_BASE_PATH)/postProcessing.h

DEP__DATA_FIELDS_2D_H = \
     $(OLB_BASE_PATH)/dataFields2D.h \
     $(DEP__OLB_DEBUG_H)

DEP__DATA_FIELDS_3D_H = \
     $(OLB_BASE_PATH)/dataFields3D.h \
     $(DEP__OLB_DEBUG_H) \
     $(DEP__DATA_FIELDS_2D_H)

DEP__BLOCK_STRUCTURE_2D_H = \
     $(OLB_BASE_PATH)/blockStructure2D.h \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__DATA_FIELDS_2D_H)

DEP__BLOCK_STRUCTURE_3D_H = \
     $(OLB_BASE_PATH)/blockStructure3D.h \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__DATA_FIELDS_3D_H)

DEP__BLOCK_LATTICE_2D_H = \
     $(OLB_BASE_PATH)/blockLattice2D.h \
     $(DEP__OLB_DEBUG_H) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__DATA_FIELDS_2D_H) \
     $(DEP__BLOCK_STRUCTURE_2D_H)

DEP__BLOCK_LATTICE_3D_H = \
     $(OLB_BASE_PATH)/blockLattice3D.h \
     $(DEP__OLB_DEBUG_H) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__DATA_FIELDS_3D_H) \
     $(DEP__BLOCK_STRUCTURE_3D_H)

DEP__BLOCK_LATTICE_VIEW_2D_H = \
     $(OLB_BASE_PATH)/blockLatticeView2D.h \
     $(DEP__BLOCK_LATTICE_2D_H)

DEP__BLOCK_LATTICE_VIEW_3D_H = \
     $(OLB_BASE_PATH)/blockLatticeView3D.h \
     $(DEP__BLOCK_LATTICE_3D_H)

DEP__BLOCK_STATISTICS_2D_H = \
     $(OLB_BASE_PATH)/blockStatistics2D.h \
     $(DEP__BLOCK_LATTICE_2D_H) \
     $(DEP__BLOCK_LATTICE_VIEW_2D_H) \
     $(DEP__DATA_FIELDS_2D_H)

DEP__BLOCK_STATISTICS_3D_H = \
     $(OLB_BASE_PATH)/blockStatistics3D.h \
     $(DEP__BLOCK_LATTICE_3D_H) \
     $(DEP__BLOCK_LATTICE_VIEW_3D_H) \
     $(DEP__DATA_FIELDS_3D_H)

DEP__BOUNDARIES_H = \
     $(OLB_BASE_PATH)/boundaries.h \
     $(DEP__DYNAMICS_H) \
     $(DEP__UTIL_H) \
     $(DEP__CELL_H)

DEP__BOUNDARIES_2D_H = \
     $(OLB_BASE_PATH)/boundaries2D.h \
     $(DEP__BOUNDARIES_H)

DEP__BOUNDARIES_3D_H = \
     $(OLB_BASE_PATH)/boundaries3D.h \
     $(DEP__BOUNDARIES_H)

DEP__FD_BOUNDARIES_2D_H = \
     $(OLB_BASE_PATH)/fdBoundaries2D.h \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__BOUNDARIES_H)

DEP__FD_BOUNDARIES_3D_H = \
     $(OLB_BASE_PATH)/fdBoundaries3D.h \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__BOUNDARIES_H)

DEP__DYNAMICS_HELPERS_2D_H = \
     $(OLB_BASE_PATH)/dynamicsHelpers2D.h \
     $(DEP__BOUNDARIES_2D_H) \
     $(DEP__FD_BOUNDARIES_2D_H)

DEP__DYNAMICS_HELPERS_3D_H = \
     $(OLB_BASE_PATH)/dynamicsHelpers3D.h \
     $(DEP__BOUNDARIES_3D_H) \
     $(DEP__FD_BOUNDARIES_3D_H)

DEP__SIMULATION_SETUP_2D_H = \
    $(OLB_BASE_PATH)/simulationSetup2D.h \
    $(DEP__BLOCK_LATTICE_2D_H) \
    $(DEP__BLOCK_LATTICE_VIEW_2D_H) 

DEP__SIMULATION_SETUP_3D_H = \
    $(OLB_BASE_PATH)/simulationSetup3D.h \
    $(DEP__BLOCK_LATTICE_3D_H) \
    $(DEP__BLOCK_LATTICE_VIEW_3D_H) 

DEP__IMAGE_CREATOR_H = \
     $(OLB_BASE_PATH)/imageCreator.h \
     $(DEP__DATA_FIELDS_2D_H) \
     $(DEP__DATA_FIELDS_3D_H) 

DEP__FINITE_DIFFERENCE_H = $(OLB_BASE_PATH)/finiteDifference.h

DEP__LB_HELPERS_H = \
    $(OLB_BASE_PATH)/lbHelpers.h \
    $(OLB_BASE_PATH)/lbHelpers2D.h \
    $(OLB_BASE_PATH)/lbHelpers3D.h \
    $(DEP__LATTICE_DESCRIPTORS_H) \
    $(DEP__CELL_H) \
    $(DEP__UTIL_H)

DEP__MRT_HELPERS_H = \
    $(OLB_BASE_PATH)/mrtHelpers.h \
    $(DEP__LATTICE_DESCRIPTORS_H) \
    $(DEP__CELL_H) \
    $(DEP__UTIL_H)

DEP__FIRST_ORDER_LB_HELPERS_H = \
    $(OLB_BASE_PATH)/firstOrderLbHelpers.h \
    $(OLB_BASE_PATH)/firstOrderLbHelpers2D.h \
    $(OLB_BASE_PATH)/firstOrderLbHelpers3D.h \
    $(DEP__LB_HELPERS_H) \
    $(DEP__LATTICE_DESCRIPTORS_H) \
    $(DEP__CELL_H) \
    $(DEP__UTIL_H)

DEP__OLB_2D_H = \
     $(OLB_BASE_PATH)/olb2D.h \
     $(DEP__SINGLETON_H) \
     $(DEP__UNITS_H) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__DYNAMICS_H) \
     $(DEP__CELL_H) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__DATA_FIELDS_2D_H) \
     $(DEP__BLOCK_STRUCTURE_2D_H) \
     $(DEP__BLOCK_LATTICE_2D_H) \
     $(DEP__BLOCK_LATTICE_VIEW_2D_H) \
     $(DEP__BLOCK_STATISTICS_2D_H) \
     $(DEP__BOUNDARIES_H) \
     $(DEP__FD_BOUNDARIES_2D_H) \
     $(DEP__DYNAMICS_HELPERS_2D_H) \
     $(DEP__SIMULATION_SETUP_2D_H) \
     $(DEP__IMAGE_CREATOR_H) \
     $(DEP__OLB_ADDON_2D_H)

DEP__OLB_3D_H = \
     $(OLB_BASE_PATH)/olb3D.h \
     $(DEP__SINGLETON_H) \
     $(DEP__UNITS_H) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__DYNAMICS_H) \
     $(DEP__CELL_H) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__DATA_FIELDS_3D_H) \
     $(DEP__BLOCK_STRUCTURE_3D_H) \
     $(DEP__BLOCK_LATTICE_3D_H) \
     $(DEP__BLOCK_LATTICE_VIEW_3D_H) \
     $(DEP__BLOCK_STATISTICS_3D_H) \
     $(DEP__BOUNDARIES_H) \
     $(DEP__FD_BOUNDARIES_3D_H) \
     $(DEP__DYNAMICS_HELPERS_3D_H) \
     $(DEP__SIMULATION_SETUP_3D_H) \
     $(DEP__IMAGE_CREATOR_H)
     $(DEP__OLB_ADDON_3D_H)

     
################  .hh files #########################################

DEP__LATTICE_DESCRIPTORS_HH = $(OLB_BASE_PATH)/latticeDescriptors.hh

DEP__DYNAMICS_HH = \
     $(OLB_BASE_PATH)/dynamics.hh \
     $(DEP__DYNAMICS_H) \
     $(DEP__CELL_H) \
     $(DEP__LB_HELPERS_H) \
     $(DEP__FIRST_ORDER_LB_HELPERS_H) \
     $(DEP__MRT_HELPERS_H)

DEP__CELL_HH = \
     $(OLB_BASE_PATH)/cell.hh \
     $(DEP__CELL_H)

DEP__POST_PROCESSING_HH = \
     $(OLB_BASE_PATH)/postProcessing.hh \
     $(DEP__BLOCK_LATTICE_2D_H) \
     $(DEP__BLOCK_LATTICE_3D_H)

DEP__DATA_FIELDS_2D_HH = \
     $(OLB_BASE_PATH)/dataFields2D.hh \
     $(DEP__DATA_FIELDS_2D_H)

DEP__DATA_FIELDS_3D_HH = \
     $(OLB_BASE_PATH)/dataFields3D.hh \
     $(DEP__DATA_FIELDS_3D_H)

DEP__BLOCK_LATTICE_2D_HH = \
     $(OLB_BASE_PATH)/blockLattice2D.hh \
     $(DEP__BLOCK_LATTICE_2D_H) \
     $(DEP__DYNAMICS_H) \
     $(DEP__CELL_H) \
     $(DEP__LB_HELPERS_H) \
     $(DEP__UTIL_H)

DEP__BLOCK_LATTICE_3D_HH = \
     $(OLB_BASE_PATH)/blockLattice3D.hh \
     $(DEP__BLOCK_LATTICE_3D_H) \
     $(DEP__DYNAMICS_H) \
     $(DEP__CELL_H) \
     $(DEP__LB_HELPERS_H) \
     $(DEP__UTIL_H)

DEP__BLOCK_LATTICE_VIEW_2D_HH = \
     $(OLB_BASE_PATH)/blockLatticeView2D.hh \
     $(DEP__BLOCK_LATTICE_VIEW_2D_H) \
     $(DEP__CELL_H)


DEP__BLOCK_LATTICE_VIEW_3D_HH = \
     $(OLB_BASE_PATH)/blockLatticeView3D.hh \
     $(DEP__BLOCK_LATTICE_VIEW_3D_H) \
     $(DEP__CELL_H)

DEP__BLOCK_STATISTICS_2D_HH = \
     $(OLB_BASE_PATH)/blockStatistics2D.h \
     $(DEP__BLOCK_STATISTICS_2D_H) \
     $(DEP__FINITE_DIFFERENCE_H) \
     $(DEP__UTIL_H)

DEP__BLOCK_STATISTICS_3D_HH = \
     $(OLB_BASE_PATH)/blockStatistics3D.hh \
     $(DEP__BLOCK_STATISTICS_3D_H) \
     $(DEP__FINITE_DIFFERENCE_H) \
     $(DEP__UTIL_H)

DEP__BOUNDARIES_HH = \
     $(OLB_BASE_PATH)/boundaries.hh \
     $(DEP__BOUNDARIES_H) \
     $(DEP__LB_HELPERS_H) \
     $(DEP__FIRST_ORDER_LB_HELPERS_H)

DEP__BOUNDARIES_2D_HH = \
     $(OLB_BASE_PATH)/boundaries2D.hh \
     $(DEP__BOUNDARIES_2D_H) \
     $(DEP__LB_HELPERS_H) \
     $(DEP__FIRST_ORDER_LB_HELPERS_H)

DEP__BOUNDARIES_3D_HH = \
     $(OLB_BASE_PATH)/boundaries3D.hh \
     $(DEP__BOUNDARIES_3D_H) \
     $(DEP__LB_HELPERS_H) \
     $(DEP__FIRST_ORDER_LB_HELPERS_H)

DEP__FD_BOUNDARIES_2D_HH = \
     $(OLB_BASE_PATH)/fdBoundaries2D.hh \
     $(DEP__FD_BOUNDARIES_2D_H) \
     $(DEP__FINITE_DIFFERENCE_H) \
     $(DEP__BLOCK_LATTICE_2D_H) \
     $(DEP__UTIL_H)

DEP__FD_BOUNDARIES_3D_HH = \
     $(OLB_BASE_PATH)/fdBoundaries3D.hh \
     $(DEP__FD_BOUNDARIES_3D_H) \
     $(DEP__FINITE_DIFFERENCE_H) \
     $(DEP__BLOCK_LATTICE_3D_H) \
     $(DEP__UTIL_H)

DEP__DYNAMICS_HELPERS_2D_HH = \
     $(OLB_BASE_PATH)/dynamicsHelpers2D.hh \
     $(DEP__DYNAMICS_HELPERS_2D_H) \
     $(DEP__BLOCK_LATTICE_2D_H)

DEP__DYNAMICS_HELPERS_3D_HH = \
     $(OLB_BASE_PATH)/dynamicsHelpers3D.hh \
     $(DEP__DYNAMICS_HELPERS_3D_H) \
     $(DEP__BLOCK_LATTICE_3D_H)

DEP__SIMULATION_SETUP_2D_HH = \
     $(OLB_BASE_PATH)/simulationSetup2D.hh \
     $(DEP__SIMULATION_SETUP_2D_H) \
     $(DEP__BLOCK_LATTICE_2D_H) \
     $(DEP__BLOCK_STATISTICS_2D_H) \
     $(DEP__LB_HELPERS_H) \
     $(DEP__FIRST_ORDER_LB_HELPERS_H) \
     $(DEP__UTIL_H)

DEP__SIMULATION_SETUP_3D_HH = \
     $(OLB_BASE_PATH)/simulationSetup3D.hh \
     $(DEP__SIMULATION_SETUP_3D_H) \
     $(DEP__BLOCK_LATTICE_3D_H) \
     $(DEP__BLOCK_STATISTICS_3D_H) \
     $(DEP__LB_HELPERS_H) \
     $(DEP__FIRST_ORDER_LB_HELPERS_H) \
     $(DEP__UTIL_H)

DEP__IMAGE_CREATOR_HH = \
     $(OLB_BASE_PATH)/imageCreator.hh \
     $(DEP__IMAGE_CREATOR_H) \
     $(DEP__SINGLETON_H)

DEP__OLB_2D_HH = \
     $(DEP__LATTICE_DESCRIPTORS_HH) \
     $(DEP__CELL_HH) \
     $(DEP__DYNAMICS_HH) \
     $(DEP__BOUNDARIES_HH) \
     $(DEP__BOUNDARIES_2D_HH) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__DATA_FIELDS_2D_HH) \
     $(DEP__FD_BOUNDARIES_2D_HH) \
     $(DEP__BLOCK_LATTICE_2D_HH) \
     $(DEP__BLOCK_LATTICE_VIEW_2D_HH) \
     $(DEP__BLOCK_STATISTICS_2D_HH) \
     $(DEP__DYNAMICS_HELPERS_2D_HH) \
     $(DEP__SIMULATION_SETUP_2D_HH) \
     $(DEP__IMAGE_CREATOR_HH) \
     $(DEP__OLB_ADDON_2D_HH)

DEP__OLB_3D_HH = \
     $(DEP__LATTICE_DESCRIPTORS_HH) \
     $(DEP__CELL_HH) \
     $(DEP__DYNAMICS_HH) \
     $(DEP__BOUNDARIES_HH) \
     $(DEP__BOUNDARIES_3D_HH) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__DATA_FIELDS_3D_HH) \
     $(DEP__FD_BOUNDARIES_3D_HH) \
     $(DEP__BLOCK_LATTICE_3D_HH) \
     $(DEP__BLOCK_LATTICE_VIEW_3D_HH) \
     $(DEP__BLOCK_STATISTICS_3D_HH) \
     $(DEP__DYNAMICS_HELPERS_3D_HH) \
     $(DEP__SIMULATION_SETUP_3D_HH) \
     $(DEP__IMAGE_CREATOR_HH) \
     $(DEP__OLB_ADDON_3D_HH)


################  .cpp files #########################################

DEP__LOADBALANCER_CPP = \
     $(OLB_BASE_PATH)/loadBalancer.cpp \
     $(DEP__LOADBALANCER_H)

DEP__OMPMANAGER_CPP = \
     $(OLB_BASE_PATH)/ompManager.cpp \
     $(DEP__OMPMANAGER_H)

DEP__LATTICE_DESCRIPTORS_CPP = \
     $(OLB_BASE_PATH)/latticeDescriptors.cpp \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__DYNAMICS_CPP = \
     $(OLB_BASE_PATH)/dynamics.cpp \
     $(DEP__DYNAMICS_H) \
     $(DEP__DYNAMICS_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__CELL_CPP = \
     $(OLB_BASE_PATH)/cell.cpp \
     $(DEP__CELL_H) \
     $(DEP__CELL_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__POST_PROCESSING_2D_CPP = \
     $(OLB_BASE_PATH)/postProcessing2D.cpp \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__POST_PROCESSING_3D_CPP = \
     $(OLB_BASE_PATH)/postProcessing3D.cpp \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__DATA_FIELDS_2D_CPP = \
     $(OLB_BASE_PATH)/dataFields2D.cpp \
     $(DEP__DATA_FIELDS_2D_H) \
     $(DEP__DATA_FIELDS_2D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__DATA_FIELDS_3D_CPP = \
     $(OLB_BASE_PATH)/dataFields3D.cpp \
     $(DEP__DATA_FIELDS_3D_H) \
     $(DEP__DATA_FIELDS_3D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__BLOCK_LATTICE_2D_CPP = \
     $(OLB_BASE_PATH)/blockLattice2D.cpp \
     $(DEP__BLOCK_LATTICE_2D_H) \
     $(DEP__BLOCK_LATTICE_2D_HH) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__BLOCK_LATTICE_3D_CPP = \
     $(OLB_BASE_PATH)/blockLattice3D.cpp \
     $(DEP__BLOCK_LATTICE_3D_H) \
     $(DEP__BLOCK_LATTICE_3D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__BLOCK_LATTICE_VIEW_2D_CPP = \
     $(OLB_BASE_PATH)/blockLatticeView2D.cpp \
     $(DEP__BLOCK_LATTICE_VIEW_2D_H) \
     $(DEP__BLOCK_LATTICE_VIEW_2D_HH) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__BLOCK_LATTICE_VIEW_3D_CPP = \
     $(OLB_BASE_PATH)/blockLatticeView3D.cpp \
     $(DEP__BLOCK_LATTICE_VIEW_3D_H) \
     $(DEP__BLOCK_LATTICE_VIEW_3D_HH) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__BLOCK_STATISTICS_2D_CPP = \
     $(OLB_BASE_PATH)/blockStatistics2D.cpp \
     $(DEP__BLOCK_STATISTICS_2D_H) \
     $(DEP__BLOCK_STATISTICS_2D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__BLOCK_STATISTICS_3D_CPP = \
     $(OLB_BASE_PATH)/blockStatistics3D.cpp \
     $(DEP__BLOCK_STATISTICS_3D_H) \
     $(DEP__BLOCK_STATISTICS_3D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__BOUNDARIES_2D_CPP = \
     $(OLB_BASE_PATH)/boundaries2D.cpp \
     $(DEP__BOUNDARIES_H) \
     $(DEP__BOUNDARIES_HH) \
     $(DEP__BOUNDARIES_2D_H) \
     $(DEP__BOUNDARIES_2D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__BOUNDARIES_3D_CPP = \
     $(OLB_BASE_PATH)/boundaries3D.cpp \
     $(DEP__BOUNDARIES_H) \
     $(DEP__BOUNDARIES_HH) \
     $(DEP__BOUNDARIES_3D_H) \
     $(DEP__BOUNDARIES_3D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__FD_BOUNDARIES_2D_CPP = \
     $(OLB_BASE_PATH)/fdBoundaries2D.cpp \
     $(DEP__BOUNDARIES_2D_H) \
     $(DEP__BOUNDARIES_2D_HH) \
     $(DEP__FD_BOUNDARIES_2D_H) \
     $(DEP__FD_BOUNDARIES_2D_HH) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__FD_BOUNDARIES_3D_CPP = \
     $(OLB_BASE_PATH)/fdBoundaries3D.cpp \
     $(DEP__BOUNDARIES_3D_H) \
     $(DEP__BOUNDARIES_3D_HH) \
     $(DEP__FD_BOUNDARIES_3D_H) \
     $(DEP__FD_BOUNDARIES_3D_HH) \
     $(DEP__POST_PROCESSING_H) \
     $(DEP__POST_PROCESSING_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__DYNAMICS_HELPERS_2D_CPP = \
     $(OLB_BASE_PATH)/dynamicsHelpers2D.cpp \
     $(DEP__DYNAMICS_HELPERS_2D_H) \
     $(DEP__DYNAMICS_HELPERS_2D_HH) \
     $(DEP__BOUNDARIES_H) \
     $(DEP__BOUNDARIES_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__DYNAMICS_HELPERS_3D_CPP = \
     $(OLB_BASE_PATH)/dynamicsHelpers3D.cpp \
     $(DEP__DYNAMICS_HELPERS_3D_H) \
     $(DEP__DYNAMICS_HELPERS_3D_HH) \
     $(DEP__BOUNDARIES_H) \
     $(DEP__BOUNDARIES_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__SIMULATION_SETUP_2D_CPP = \
     $(OLB_BASE_PATH)/simulationSetup2D.cpp \
     $(DEP__SIMULATION_SETUP_2D_H) \
     $(DEP__SIMULATION_SETUP_2D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__SIMULATION_SETUP_3D_CPP = \
     $(OLB_BASE_PATH)/simulationSetup3D.cpp \
     $(DEP__SIMULATION_SETUP_3D_H) \
     $(DEP__SIMULATION_SETUP_3D_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)

DEP__IMAGE_CREATOR_CPP = \
     $(OLB_BASE_PATH)/imageCreator.cpp \
     $(DEP__IMAGE_CREATOR_H) \
     $(DEP__IMAGE_CREATOR_HH) \
     $(DEP__LATTICE_DESCRIPTORS_H) \
     $(DEP__LATTICE_DESCRIPTORS_HH)


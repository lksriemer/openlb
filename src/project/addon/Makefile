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

## This file is used to compile all .cpp source files to .o object
## files for the precompiled version of the library.
## By default, invoking the command "make" in the current directory
## lists all cpp files.

OLB_BASE_PATH = .

include dependencies.mk

default:
	echo $(DEP__ADDON_BOUNDARIES_2D_CPP)
	echo $(DEP__ADDON_BOUNDARIES_3D_CPP)

clean:
	/bin/rm -f *.o *~

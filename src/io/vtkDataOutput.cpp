/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Bernd Stahl, Jonas Latt
 *  Address: Battelle Batiment A, Route de Drize 7, 1227 Carouge, Switzerland
 *  E-mail: bernd.stahl@cui.unige.ch
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

#include "complexGrids/mpiManager/mpiManager.h"
#include "vtkDataOutput.h"
#include "vtkDataOutput.hh"
#include "serializerIO.h"
#include "serializerIO.hh"
#include "base64.h"
#include "base64.hh"

namespace olb {

template class VtkDataWriter3D<double>;

template void writeVTKData3D<double> (
        std::string const& fName,
        std::string const& scalarFieldName,
        ScalarFieldBase3D<double> const& scalarField,
        std::string const& vectorFieldName,
        TensorFieldBase3D<double,3> const& vectorField,
        double deltaX, double deltaT );

template<>
std::string VtkTypeNames<bool>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<char>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned char>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<short int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned short int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<long int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned long int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<float>::getBaseName() {
    return "Float";
}

template<>
std::string VtkTypeNames<double>::getBaseName() {
    return "Float";
}

template<>
std::string VtkTypeNames<long double>::getBaseName() {
    return "Float";
}

}

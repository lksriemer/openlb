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

#ifndef VTK_DATA_OUTPUT_HH
#define VTK_DATA_OUTPUT_HH

#include "complexGrids/mpiManager/mpiManager.h"
#include "vtkDataOutput.h"
#include "serializerIO.h"
#include "core/singleton.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace olb {


////////// struct VtkTypeNames ///////////////////////////////////

template<typename T>
class VtkTypeNames {
public:
    static std::string getName();
private:
    static std::string getBaseName();
};

template<typename T>
std::string VtkTypeNames<T>::getName() {
    std::stringstream sstream;
    sstream << getBaseName();
    sstream << 8 * sizeof(T);

    std::string tName;
    sstream >> tName;
    return tName;
}


////////// class VtkDataWriter3D ////////////////////////////////////////

template<typename T>
void VtkDataWriter3D::writeDataField(DataSerializer<T> const& serializer,
                                     std::string const& name, T scalingFactor, int nDim)
{
    if (singleton::mpi().isMainProcessor()) {
        (*ostr) << "<DataArray type=\"" << VtkTypeNames<T>::getName()
                << "\" Name=\"" << name
                << "\" format=\"binary\" encoding=\"base64";
        if (nDim>1) {
            (*ostr) << "\" NumberOfComponents=\"" << nDim;
        }
        (*ostr) << "\">\n";
    }

    // undocumented requirement of the vtk xml file format:
    // in front of every binary blob, base64 or raw-binary, appended or not, 
    // there is an UInt32 length indicator, giving the size of the binary blob in bytes;
    // when using base64 encoding, that length header must be encoded separately;
    // there must be no newline between the encoded length indicator and the encoded data block.
    //
    // those properties are properly handled by the serializer2ostr function

    ScalingSerializer<T> scaledSerializer(serializer, scalingFactor);
    serializer2ostr(scaledSerializer, ostr);

    if (singleton::mpi().isMainProcessor()) {
        (*ostr) << "\n</DataArray>\n";
    }
}



////////// Free Functions //////////////////////////////////////////////

template<typename T> void writeVTKData3D (
        std::string const& fName,
        std::string const& scalarFieldName,
        ScalarFieldBase3D<T> const& scalarField,
        std::string const& vectorFieldName,
        TensorFieldBase3D<T,3> const& vectorField,
        T deltaX, T deltaT )
{
    std::string fullName = singleton::directories().getVtkOutDir() + fName+".vti";
    VtkDataWriter3D vtiOut(fullName);
    int nx = scalarField.getNx();
    int ny = scalarField.getNy();
    int nz = scalarField.getNz();
    vtiOut.writeHeader(0,nx-1,0,ny-1,0,nz-1,0,0,0,deltaX);

    vtiOut.startPiece(0,nx-1,0,ny-1,0,nz-1);
    vtiOut.writeDataField(scalarField.getSerializer(IndexOrdering::backward), scalarFieldName, 1./deltaT, 1);

    vtiOut.writeDataField(vectorField.getSerializer(IndexOrdering::backward), vectorFieldName, deltaX/deltaT, 3);
    vtiOut.endPiece();

    vtiOut.writeFooter();
}


}  // namespace olb

#endif



/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef IMAGE_CREATOR_HH
#define IMAGE_CREATOR_HH

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include "imageCreator.h"
#include "singleton.h"
#include "dataFields3D.hh"

namespace olb {

namespace graphics {


////////// class ImageCreator ////////////////////////////////////////

template<typename T>
ImageCreator<T>::ImageCreator(std::string const& map) {
    setMap(map);
}

template<typename T>
void ImageCreator<T>::setMap(std::string const& map) {
    mapColors.resize(256);
    std::string fullName = singleton::directories().getOlbDir() +
                      "src/project/graphics/colormaps/" +
                      map;
    std::ifstream mapFile(fullName.c_str());
    if (!mapFile) {
        std::cout << "could not open " << fullName << std::endl;
    }
    for (int iMap=0; iMap<256; ++iMap) {
        mapFile >> mapColors[iMap].r >> mapColors[iMap].g
                >> mapColors[iMap].b;
    }
}

template<typename T>
void ImageCreator<T>::writePpm(std::string const& fName,
                               ScalarField2D<T> const& field,
                               T minVal, T maxVal) const
{
    std::string fullName = singleton::directories().getImageOutDir() +
                      fName+".ppm";
    std::ofstream fout(fullName.c_str());
    fout << "P3\n";
    fout << field.getNx() << " " << field.getNy() << "\n";
    fout << "255\n";

    for (int iY=field.getNy()-1; iY>=0; --iY) {
        for (int iX=0; iX<field.getNx(); ++iX) {
            int outputValue = (int)
                ( (field.get(iX,iY)-minVal)/(maxVal-minVal)*(T)255 );
            if (outputValue<0) outputValue = 0;
            if (outputValue>255) outputValue = 255;
            fout << mapColors[outputValue].r << " "
                 << mapColors[outputValue].g << " "
                 << mapColors[outputValue].b << "\n";
        }
    }
}

template<typename T>
void ImageCreator<T>::writeGif(std::string const& fName,
                               ScalarField2D<T> const& field,
                               T minVal, T maxVal) const
{
    writePpm(fName, field, minVal, maxVal);

    std::string convCommand =
        std::string("convert ") +
        singleton::directories().getImageOutDir() + fName + ".ppm " +
        singleton::directories().getImageOutDir() + fName + ".gif ";
        
    std::string rmCommand =
        std::string("/bin/rm ") +
        singleton::directories().getImageOutDir() + fName + ".ppm";

    system(convCommand.c_str());
    system(rmCommand.c_str());
}

template<typename T>
void ImageCreator<T>::writeGif(std::string const& fName,
                               ScalarField2D<T> const& field,
                               T minVal, T maxVal,
                               T sizeX, T sizeY) const
{
    writePpm(fName, field, minVal, maxVal);
    std::stringstream imStream;
    imStream << "convert -resize "
             << sizeX << "x" << sizeY << " "
             << singleton::directories().getImageOutDir()
             << fName << ".ppm "
             << singleton::directories().getImageOutDir()
             << fName << ".gif";
    system(imStream.str().c_str());

    std::string rmCommand =
        std::string("/bin/rm ") +
        singleton::directories().getImageOutDir() + fName + ".ppm";
    system(rmCommand.c_str());
}

template<typename T>
void ImageCreator<T>::writeScaledGif(std::string const& fName,
                                     ScalarField2D<T> const& field) const
{
    T minVal = computeMin(field);
    T maxVal = computeMax(field);
    writeGif(fName, field, minVal, maxVal);
}

template<typename T>
void ImageCreator<T>::writeScaledGif(std::string const& fName,
                                     ScalarField2D<T> const& field,
                                     T sizeX, T sizeY) const
{
    T minVal = computeMin(field);
    T maxVal = computeMax(field);
    writeGif(fName, field, minVal, maxVal, sizeX, sizeY);
}

template<typename T>
void ImageCreator<T>::writeScaledPpm(std::string const& fName,
                                     ScalarField2D<T> const& field) const
{
    T minVal = computeMin(field);
    T maxVal = computeMax(field);
    writePpm(fName, field, minVal, maxVal);
}

template<typename T>
void ImageCreator<T>::writeText(std::string const& fName,
                                ScalarField2D<T> const& field) const
{
    std::string fullName = singleton::directories().getImageOutDir() +
                      fName+".txt";
    std::ofstream fout(fullName.c_str());
    if (!fout) {
        std::cout << "could not open " << fullName << std::endl;
    }

    for (int iY=field.getNy()-1; iY>=0; --iY) {
        for (int iX=0; iX<field.getNx(); ++iX) {
            fout << field.get(iX,iY) << " ";
        }
        fout << "\n";
    }
}


////////// class VTKOut2D ////////////////////////////////////////

template<typename T>
void VTKOut2D<T>::writeFlowField (
        std::string const& fName,
        std::string const& scalarFieldName,
        ScalarField2D<T> const& scalarField,
        std::string const& vectorFieldName,
        TensorField2D<T,2> const& vectorField,
        T deltaX, T deltaT )
{
    OLB_PRECONDITION( scalarField.getNx() == vectorField.getNx() );
    OLB_PRECONDITION( scalarField.getNy() == vectorField.getNy() );

    int nx = scalarField.getNx();
    int ny = scalarField.getNy();

    std::string fullName = singleton::directories().getVtkOutDir() +
                      fName+".vti";
    std::ofstream fout(fullName.c_str());
    if (!fout) {
        std::cout << "could not open " << fullName << std::endl;
    }

    writePreamble(fout, nx, ny, deltaX);

    fout << "<PointData Scalars=\""
         << scalarFieldName << "\" "
         <<            "Vectors=\""
         << vectorFieldName << "\">\n";

    fout << "<DataArray type=\"Float32\" Name=\""
         << scalarFieldName << "\">\n";
    for (int iY=0; iY<ny; ++iY) {
        for (int iX=0; iX<nx; ++iX) {
            fout << scalarField.get(iX,iY)/deltaT << " ";
        }
        fout << "\n";
    }
    fout << "</DataArray>\n";

    fout << "<DataArray type=\"Float32\" Name=\""
         << vectorFieldName << "\" "
         << "NumberOfComponents=\"2\">\n";
    for (int iY=0; iY<ny; ++iY) {
        for (int iX=0; iX<nx; ++iX) {
            fout << vectorField.get(iX,iY)[0]*deltaX/deltaT << " ";
            fout << vectorField.get(iX,iY)[1]*deltaX/deltaT << " 0 ";
        }
        fout << "\n";
    }
    fout << "</DataArray>\n";

    fout << "</PointData>\n";

    writePostScript(fout);
}

template<typename T>
void VTKOut2D<T>::writePreamble(std::ofstream& fout, int nx, int ny,
                                T deltaX)
{
    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
         << "byte_order=\"LittleEndian\">\n";
    fout << "<ImageData WholeExtent=\"0 "
         << nx-1 << " 0 "
         << ny-1 << " 0 0\" "
         << "Origin=\"0 0 0\" Spacing=\""
         << deltaX << " " << deltaX << " " << deltaX << "\">\n";
    fout << "<Piece Extent=\"0 "
         << nx-1 << " 0 "
         << ny-1 << " 0 0\">\n";
}

template<typename T>
void VTKOut2D<T>::writePostScript(std::ofstream& fout) {
    fout << "</Piece>\n";
    fout << "</ImageData>\n";
    fout << "</VTKFile>\n";
}


////////// class VTKOut3D ////////////////////////////////////////

template<typename T>
    void VTKOut3D<T>::writeFlowField (
            std::string const& fName,
            std::string const& scalarFieldName,
            ScalarField3D<T> const& scalarField,
            std::string const& vectorFieldName,
            TensorField3D<T,3> const& vectorField,
            T deltaX, T deltaT )
{
    OLB_PRECONDITION( scalarField.getNx() == vectorField.getNx() );
    OLB_PRECONDITION( scalarField.getNy() == vectorField.getNy() );
    OLB_PRECONDITION( scalarField.getNz() == vectorField.getNz() );

    int nx = scalarField.getNx();
    int ny = scalarField.getNy();
    int nz = scalarField.getNz();

    std::string fullName = singleton::directories().getVtkOutDir() +
                      fName+".vti";
    std::ofstream fout(fullName.c_str());
    if (!fout) {
        std::cout << "could not open " << fullName << std::endl;
    }

    writePreamble(fout, nx, ny, nz, deltaX);

    fout << "<PointData Scalars=\""
         << scalarFieldName << "\" "
         <<            "Vectors=\""
         << vectorFieldName << "\">\n";

    fout << "<DataArray type=\"Float32\" Name=\""
         << scalarFieldName << "\">\n";
    for (int iZ=0; iZ<nz; ++iZ) {
        for (int iY=0; iY<ny; ++iY) {
            for (int iX=0; iX<nx; ++iX) {
                fout << scalarField.get(iX,iY,iZ)/deltaT << " ";
            }
        }
        fout << "\n";
    }    
    fout << "</DataArray>\n";

    fout << "<DataArray type=\"Float32\" Name=\""
         << vectorFieldName << "\" "
         << "NumberOfComponents=\"3\">\n";
    for (int iZ=0; iZ<nz; ++iZ) {
        for (int iY=0; iY<ny; ++iY) {
            for (int iX=0; iX<nx; ++iX) {
                fout << vectorField.get(iX,iY,iZ)[0]*deltaX/deltaT << " ";
                fout << vectorField.get(iX,iY,iZ)[1]*deltaX/deltaT << " ";
                fout << vectorField.get(iX,iY,iZ)[2]*deltaX/deltaT << " ";
            }
        }
        fout << "\n";
    }    

    fout << "</DataArray>\n";

    fout << "</PointData>\n";

    writePostScript(fout);

}

template<typename T>
void VTKOut3D<T>::writePreamble(std::ofstream& fout, int nx, int ny, int nz,
                                T deltaX)
{
    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
         << "byte_order=\"LittleEndian\">\n";
    fout << "<ImageData WholeExtent=\"0 "
         << nx-1 << " 0 "
         << ny-1 << " 0 "
         << nz-1 << " \" "
         << "Origin=\"0 0 0\" Spacing=\""
         << deltaX << " " << deltaX << " " << deltaX << "\">\n";
    fout << "<Piece Extent=\"0 "
         << nx-1 << " 0 "
         << ny-1 << " 0 "
         << nz-1 << " \">\n";
}

template<typename T>
void VTKOut3D<T>::writePostScript(std::ofstream& fout) {
    fout << "</Piece>\n";
    fout << "</ImageData>\n";
    fout << "</VTKFile>\n";
}

}  // namespace graphics

}  // namespace olb

#endif

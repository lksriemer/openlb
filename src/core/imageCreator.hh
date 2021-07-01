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

#include "complexGrids/mpiManager/mpiManager.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "imageCreator.h"
#include "singleton.h"
#include "dataFields2D.h"

namespace olb {

namespace graphics {

////////// class ImageCreator ////////////////////////////////////////

template<typename T>
ImageCreator<T>::ImageCreator(std::string const& map) {
    setMap(map);
}

template<typename T>
void ImageCreator<T>::setMap(std::string const& map) {
    if (singleton::mpi().isMainProcessor()) {
        mapColors.resize(256);
        std::string fullName = singleton::directories().getOlbDir() + "src/io/colormaps/" + map;
        std::ifstream mapFile(fullName.c_str());
        if (!mapFile) {
            std::cout << "could not open " << fullName << std::endl;
        }
        for (int iMap=0; iMap<256; ++iMap) {
            mapFile >> mapColors[iMap].r >> mapColors[iMap].g
                    >> mapColors[iMap].b;
        }
    }
}

template<typename T>
void ImageCreator<T>::writePpm(std::string const& fName,
                               ScalarFieldBase2D<T> const& field,
                               T minVal, T maxVal) const
{
    ScalarField2D<T> localField(field.getNx(), field.getNy());
    localField.construct();
    copyDataBlock(field, localField);
    if (singleton::mpi().isMainProcessor()) {
        if (std::fabs(minVal-maxVal)<1.e-12) {
            minVal = computeMin(localField);
            maxVal = computeMax(localField);
        }
        std::string fullName = singleton::directories().getImageOutDir() + fName+".ppm";
        std::ofstream fout(fullName.c_str());
        fout << "P3\n";
        fout << localField.getNx() << " " << localField.getNy() << "\n";
        fout << "255\n";

        for (int iY=localField.getNy()-1; iY>=0; --iY) {
            for (int iX=0; iX<localField.getNx(); ++iX) {
                int outputValue = (int)
                    ( (localField.get(iX,iY)-minVal)/(maxVal-minVal)*(T)255 );
                if (outputValue<0) outputValue = 0;
                if (outputValue>255) outputValue = 255;
                fout << mapColors[outputValue].r << " "
                     << mapColors[outputValue].g << " "
                     << mapColors[outputValue].b << "\n";
            }
        }
    }
}

template<typename T>
void ImageCreator<T>::writeGif(std::string const& fName,
                               ScalarFieldBase2D<T> const& field,
                               T minVal, T maxVal) const
{
    writePpm(fName, field, minVal, maxVal);

    if (singleton::mpi().isMainProcessor()) {
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
}

template<typename T>
void ImageCreator<T>::writeGif(std::string const& fName,
                               ScalarFieldBase2D<T> const& field,
                               T minVal, T maxVal,
                               T sizeX, T sizeY) const
{
    writePpm(fName, field, minVal, maxVal);
    if (singleton::mpi().isMainProcessor()) {
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
}

template<typename T>
void ImageCreator<T>::writeScaledGif(std::string const& fName,
                                     ScalarFieldBase2D<T> const& field) const
{
    writeGif(fName, field, T(), T());
}

template<typename T>
void ImageCreator<T>::writeScaledGif(std::string const& fName,
                                     ScalarFieldBase2D<T> const& field,
                                     T sizeX, T sizeY) const
{
    writeGif(fName, field, T(), T(), sizeX, sizeY);
}

template<typename T>
void ImageCreator<T>::writeScaledPpm(std::string const& fName,
                                     ScalarFieldBase2D<T> const& field) const
{
    writePpm(fName, field, T(), T());
}

template<typename T>
void ImageCreator<T>::writeText(std::string const& fName,
                                ScalarFieldBase2D<T> const& field) const
{
    ScalarField2D<T> localField(field.getNx(), field.getNy());
    localField.construct();
    copyDataBlock(field, localField);
    if (singleton::mpi().isMainProcessor()) {
        std::string fullName = singleton::directories().getImageOutDir() + fName+".txt";
        std::ofstream fout(fullName.c_str());
        if (!fout) {
            std::cout << "could not open " << fullName << std::endl;
        }

        for (int iY=localField.getNy()-1; iY>=0; --iY) {
            for (int iX=0; iX<localField.getNx(); ++iX) {
                fout << localField.get(iX,iY) << " ";
            }
            fout << "\n";
        }
    }
}

}  // namespace graphics

}  // namespace olb

#endif

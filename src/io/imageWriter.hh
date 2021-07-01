/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
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

#ifndef IMAGE_WRITER_HH
#define IMAGE_WRITER_HH

#include "communication/mpiManager.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "imageWriter.h"
#include "colormaps.h"
#include "core/singleton.h"
#include "core/dataFields2D.h"


namespace olb {

namespace graphics {

////////// class ImageWriter ////////////////////////////////////////

template<typename T>
ImageWriter<T>::ImageWriter(std::string const& map)
  : clout(std::cout,"ImageWriter"), _colorRange(1024), _numColors(1024),
    _colorMap( mapGenerators::generateMap<T>(map) )
{ }

template<typename T>
ImageWriter<T>::ImageWriter(std::string const& map, int colorRange, int numColors)
  : clout(std::cout,"ImageWriter"), _colorRange(colorRange), _numColors(numColors),
    _colorMap( mapGenerators::generateMap<T>(map) )
{ }

template<typename T>
void ImageWriter<T>::setMap(std::string const& map, int colorRange, int numColors) {
  _colorRange = colorRange;
  _numColors  = numColors;
  _colorMap   = mapGenerators::generateMap<T>(map);
}

/// the whole work of .gif is here
/// fist create a .ppm file the via imagemagick .gif
template<typename T>
void ImageWriter<T>::writePpm(std::string const& fName,
  ScalarFieldBase2D<T> const& field, T minVal, T maxVal) const
{
  ScalarField2D<T> localField(field.getNx(), field.getNy());
  localField.construct();
  copyDataBlock(field, localField);
  if (singleton::mpi().isMainProcessor()) {
    // wenn uebergebene min-max Werte sinnlos sind.
    if (std::fabs(minVal-maxVal)<1.e-12) {
      minVal = computeMin(localField);
      maxVal = computeMax(localField);
    }
    std::string fullName = singleton::directories().getImageOutDir() + fName+".ppm";
    std::ofstream fout(fullName.c_str());
    // write header
    fout << "P3\n";
    // dimension of image
    fout << localField.getNx() << " " << localField.getNy() << "\n";
    // dynamic range
    fout << (_colorRange-1) << "\n";

    for (int iY=localField.getNy()-1; iY>=0; --iY)
    { // laeuft von oben links los, und dann x0, x1, ...
      for (int iX=0; iX<localField.getNx(); ++iX)
      {
        T outputValue = ( (localField.get(iX,iY)-minVal)/(maxVal-minVal)/(T)_numColors*(T)(_numColors-1) );
        if (outputValue <   T()) outputValue = T();
        if (outputValue >= (T)1) outputValue = (T)(_numColors-1)/(T)_numColors;
        rgb<T> color = _colorMap.get(outputValue);
        fout << (int) (color.r*(_colorRange-1)) << " "
             << (int) (color.g*(_colorRange-1)) << " "
             << (int) (color.b*(_colorRange-1)) << "\n";
      }
    }
  }
}

template<typename T>
void ImageWriter<T>::writeGif(std::string const& fName,
  ScalarFieldBase2D<T> const& field, T minVal, T maxVal) const
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

    int retConv = system(convCommand.c_str());
    if (retConv!=0) clout << "Error: ImageMagick convert command failed" << std::endl;
//    int retRm = system(rmCommand.c_str());
//    if (retRm!=0) clout << "Error: removing temporary ppm file failed" << std::endl;
  }
}

template<typename T>
void ImageWriter<T>::writeGif(std::string const& fName,
  ScalarFieldBase2D<T> const& field, T minVal, T maxVal, T sizeX, T sizeY) const
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
    int retResize = system(imStream.str().c_str());
    if (retResize!=0) clout << "Error: ImageMagick resize command failed" << std::endl;

    std::string rmCommand =
      std::string("/bin/rm ") +
      singleton::directories().getImageOutDir() + fName + ".ppm";
//    int retRm = system(rmCommand.c_str());
//    if (retRm!=0) clout << "Error: removing temporary ppm file failed" << std::endl;
  }
}

template<typename T>
void ImageWriter<T>::writeScaledGif(std::string const& fName,
                                    ScalarFieldBase2D<T> const& field) const
{
  writeGif(fName, field, T(), T());
}

template<typename T>
void ImageWriter<T>::writeScaledGif(std::string const& fName,
                                    ScalarFieldBase2D<T> const& field,
                                    T sizeX, T sizeY) const
{
  writeGif(fName, field, T(), T(), sizeX, sizeY);
}

template<typename T>
void ImageWriter<T>::writeScaledPpm(std::string const& fName,
  ScalarFieldBase2D<T> const& field) const
{
  writePpm(fName, field, T(), T());
}

template<typename T>
void ImageWriter<T>::writeText(std::string const& fName,
  ScalarFieldBase2D<T> const& field) const
{
  ScalarField2D<T> localField(field.getNx(), field.getNy());
  localField.construct();
  copyDataBlock(field, localField);
  if (singleton::mpi().isMainProcessor()) {
    std::string fullName = singleton::directories().getImageOutDir() + fName+".txt";
    std::ofstream fout(fullName.c_str());
    if (!fout) {
      clout << "Error: could not open " << fullName << std::endl;
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

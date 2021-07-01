/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

#ifndef SERIALIZER_IO_HH
#define SERIALIZER_IO_HH

#include "complexGrids/mpiManager/mpiManager.h"
#include "serializerIO.h"
#include "base64.h"
#include "core/olbDebug.h"
#include <istream>
#include <ostream>
#include <fstream>

namespace olb {

template<typename T>
void serializer2ostr(DataSerializer<T> const& serializer, std::ostream* ostr) {
    int fullSize = 0;
    if (singleton::mpi().isMainProcessor()) {
        Base64Encoder<unsigned int> sizeEncoder(*ostr, 1);
        fullSize = serializer.getSize();
        unsigned int binarySize = (unsigned int) (fullSize * sizeof(T));
        sizeEncoder.encode(&binarySize, 1);
    }
    Base64Encoder<T>* dataEncoder = 0;
    if (singleton::mpi().isMainProcessor()) {
        dataEncoder = new Base64Encoder<T>(*ostr, fullSize);
    }
    while (!serializer.isEmpty()) {
        int bufferSize;
        const T* dataBuffer = serializer.getNextDataBuffer(bufferSize);
        if (singleton::mpi().isMainProcessor()) {
            dataEncoder->encode(dataBuffer, bufferSize);
        }
    }
    delete dataEncoder;
}

template<typename T>
void saveData(Serializable<T> const& object, std::string fName) {
    std::ofstream* ostr = 0;
    if (singleton::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fName.c_str());
        OLB_PRECONDITION( *ostr );
    }
    serializer2ostr(object.getSerializer(IndexOrdering::memorySaving), ostr);
    delete ostr;
}

template<typename T>
void istr2unSerializer(DataUnSerializer<T>& unSerializer, std::istream* istr) {
    int fullSize = 0;
    if (singleton::mpi().isMainProcessor()) {
        Base64Decoder<unsigned int> sizeDecoder(*istr, 1);
        unsigned int binarySize;
        sizeDecoder.decode(&binarySize, 1);
        fullSize = (int)(binarySize / sizeof(T));

        OLB_PRECONDITION((int)fullSize == unSerializer.getSize());
    }
    Base64Decoder<T>* dataDecoder = 0;
    if (singleton::mpi().isMainProcessor()) {
        dataDecoder = new Base64Decoder<T>(*istr, fullSize);
    }
    while (!unSerializer.isFull()) {
        int bufferSize = 0;
        T* dataBuffer = unSerializer.getNextDataBuffer(bufferSize);
        if (singleton::mpi().isMainProcessor()) {
            dataDecoder->decode(dataBuffer, bufferSize);
        }
        unSerializer.commitData();
    }
    delete dataDecoder;
}

template<typename T>
void loadData(Serializable<T>& object, std::string fName) {
    std::ifstream* istr = 0;
    if (singleton::mpi().isMainProcessor()) {
        istr = new std::ifstream(fName.c_str());
        OLB_PRECONDITION( *istr );
    }
    istr2unSerializer(object.getUnSerializer(IndexOrdering::memorySaving), istr);
    delete istr;
}

} // namespace olb

#endif

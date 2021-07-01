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


#ifndef SERIALIZER_HH
#define SERIALIZER_HH

#include "complexGrids/mpiManager/mpiManager.h"
#include "serializer.h"
#include "olbDebug.h"
#include <algorithm>

namespace olb {

////////// class ScalingSerializer ////////////////////////////

template<typename T>
ScalingSerializer<T>::ScalingSerializer(DataSerializer<T> const& baseSerializer_, T scalingFactor_)
    : baseSerializer(baseSerializer_),
      scalingFactor(scalingFactor_)
{ }

template<typename T>
int ScalingSerializer<T>::getSize() const {
    return baseSerializer.getSize();
}

template<typename T>
const T* ScalingSerializer<T>::getNextDataBuffer(int& bufferSize) const {
    const T* unscaledBuffer = baseSerializer.getNextDataBuffer(bufferSize);
    scaledBuffer.resize(bufferSize);
    for (int iBuffer=0; iBuffer<bufferSize; ++iBuffer) {
        scaledBuffer[iBuffer] = unscaledBuffer[iBuffer] * scalingFactor;
    }
    return &scaledBuffer[0];
}

template<typename T>
bool ScalingSerializer<T>::isEmpty() const {
    return baseSerializer.isEmpty();
}

template<typename T>
void copySerializedData(DataSerializer<T> const& serializer, DataUnSerializer<T>& unSerializer) {
    OLB_PRECONDITION( serializer.getSize() == unSerializer.getSize() );
    int writePos = 0, readPos = 0;
    int serializerBufferSize =0, unSerializerBufferSize =0;
    const T* serializerBuffer =0;
    T* unSerializerBuffer =0;
    while (!serializer.isEmpty()) {
        if (readPos==serializerBufferSize) {
            serializerBuffer = serializer.getNextDataBuffer(serializerBufferSize);
            readPos = 0;
        }
        if (writePos==unSerializerBufferSize) {
            unSerializerBuffer = unSerializer.getNextDataBuffer(unSerializerBufferSize);
            writePos = 0;
        }

        int remainToRead = serializerBufferSize - readPos;
        int remainToWrite = unSerializerBufferSize - writePos;
        int nextChunk = std::min(remainToRead, remainToWrite);
        if (singleton::mpi().isMainProcessor()) {
            for (int iChunk=0; iChunk<nextChunk; ++iChunk, ++readPos, ++writePos) {
                unSerializerBuffer[writePos] = serializerBuffer[readPos];
            }
        }
        if (writePos==unSerializerBufferSize) {
            unSerializer.commitData();
        }
    }
}

template<typename T>
void copyDataBlock(Serializable<T> const& from, Serializable<T>& to, IndexOrdering::OrderingT ordering) {
    copySerializedData(from.getSerializer(ordering), to.getUnSerializer(ordering));
}

}  // namespace olb

#endif



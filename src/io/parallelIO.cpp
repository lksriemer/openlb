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

#include "complexGrids/mpiManager/mpiManager.h"
#include "complexGrids/mpiManager/mpiManager.hh"
#include "parallelIO.h"

namespace olb {

///////////////////////////////////////////////////////////////////
// Class ParBuf
///////////////////////////////////////////////////////////////////

ParBuf::ParBuf(std::streambuf* _originalBuf)
  : originalBuf(_originalBuf), mode(normal)
{ }

std::streambuf::int_type
ParBuf::overflow (std::streambuf::int_type c) {
  int_type returnVal = c;
  if (c != EOF) {
#ifdef PARALLEL_MODE_MPI
    if (singleton::mpi().isMainProcessor()) {
#endif
      returnVal = originalBuf->sputc((char)c);
#ifdef PARALLEL_MODE_MPI
    }
    if (mode==normal) {
      singleton::mpi().bCast(&returnVal, 1);
    }
#endif
  }
  return returnVal;
}

ParBuf::Modes
ParBuf::getMode() const {
  return mode;
}

void
ParBuf::setMode(ParBuf::Modes _mode) {
  mode = _mode;
}

std::streamsize
ParBuf::xsputn(const char* s, std::streamsize num) {
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    return originalBuf->sputn(s,num);
#ifdef PARALLEL_MODE_MPI
  }
  else {
    return num;
  }
#endif
}

std::streambuf::int_type
ParBuf::uflow() {
  int_type value;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    value = originalBuf->sbumpc();
#ifdef PARALLEL_MODE_MPI
  }
  if (mode==normal) {
      singleton::mpi().bCast(&value, 1);
  }
#endif
  return value;
}

std::streambuf::int_type
ParBuf::underflow() {
  int_type value;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    value = originalBuf->sgetc();
#ifdef PARALLEL_MODE_MPI
  }
  if (mode==normal) {
    singleton::mpi().bCast(&value, 1);
  }
#endif
  return value;
}

std::streamsize
ParBuf::xsgetn (char* s, std::streamsize num) {
    std::streamsize sizeRead;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    sizeRead = originalBuf->sgetn(s, num);
#ifdef PARALLEL_MODE_MPI
  }
  if (mode==normal) {
      singleton::mpi().bCast(&sizeRead, 1);
      singleton::mpi().bCast(s, sizeRead);
  }
#endif
  return sizeRead;
}

///////////////////////////////////////////////////////////////////
// Class vofstream
///////////////////////////////////////////////////////////////////

vofstream::vofstream() : std::ostream(NULL), fbuf(), mybuf(&fbuf) {
  this->init(&mybuf);
}

vofstream::vofstream(const char * filename, openmode mode)
  : std::ostream(NULL), fbuf(), mybuf(&fbuf)
{ 
  init(&mybuf);
  open(filename, mode);
}

vofstream::~vofstream()
{ }

std::streambuf*
vofstream::rdbuf() const {
  return const_cast<ParBuf*>(&mybuf);
}

bool
vofstream::is_open() {
  return fbuf.is_open();
}

void
vofstream::open(const char* filename, openmode mode) {
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.open(filename, mode | ios_base::out);
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    this->setstate(ios_base::failbit);
  }
}

void
vofstream::close() {
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.close();
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    setstate(ios_base::failbit);
  }
}



///////////////////////////////////////////////////////////////////
// Class vifstream
///////////////////////////////////////////////////////////////////

vifstream::vifstream() : std::istream(NULL), fbuf(), mybuf(&fbuf) {
  init(&mybuf);
}

vifstream::vifstream(const char * filename, openmode mode)
  : std::istream(NULL), fbuf(), mybuf(&fbuf)
{ 
  init(&mybuf);
  open(filename, mode);
}

vifstream::~vifstream()
{ }

std::streambuf*
vifstream::rdbuf() const {
  return const_cast<ParBuf*>(&mybuf);
}

bool
vifstream::is_open() {
  return fbuf.is_open();
}


void
vifstream::open(const char* filename, openmode mode) {
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.open(filename, mode | ios_base::in);
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    this->setstate(ios_base::failbit);
  }
}

void
vifstream::close() {
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.close();
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    setstate(ios_base::failbit);
  }
}


///////////////////////////////////////////////////////////////////
// Class vfstream
///////////////////////////////////////////////////////////////////

vfstream::vfstream() : std::iostream(NULL), fbuf(), mybuf(&fbuf) {
  this->init(&mybuf);
}

vfstream::vfstream(const char * filename, openmode mode)
  : std::iostream(NULL), fbuf(), mybuf(&fbuf)
{ 
  init(&mybuf);
  open(filename, mode);
}

vfstream::~vfstream()
{ }

std::streambuf*
vfstream::rdbuf() const {
  return const_cast<ParBuf*>(&mybuf);
}

bool
vfstream::is_open() {
  return fbuf.is_open();
}

void
vfstream::open(const char* filename, openmode mode) {
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.open(filename, mode);
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    this->setstate(ios_base::failbit);
  }
}

void
vfstream::close() {
  int ok;
#ifdef PARALLEL_MODE_MPI
  if (singleton::mpi().isMainProcessor()) {
#endif
    ok = (bool) fbuf.close();
#ifdef PARALLEL_MODE_MPI
  }
  singleton::mpi().bCast(&ok, 1);
#endif
  if (!ok) {
    setstate(ios_base::failbit);
  }
}

}

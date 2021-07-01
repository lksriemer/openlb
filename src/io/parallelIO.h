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

#ifndef PARALLEL_IO_H
#define PARALLEL_IO_H

#include <streambuf>
#include <istream>
#include <ostream>
#include <fstream>
#include <iostream>

namespace olb {

class ParBuf : public std::streambuf {
public:
  typedef enum{normal, serial} Modes;
  // There are two modes for proceeding of parallel i/o:
  // - normal: the default mode. all nodes write data, and all
  //           nodes receive data; this mode requests broadcasting.
  // - serial: only node 0 writes and receives data.
public:
  ParBuf(std::streambuf* _originalBuf);
  Modes   getMode() const;
  void    setMode(Modes _mode);
protected:
  virtual int_type overflow (int_type c);
  virtual std::streamsize xsputn(const char* s, std::streamsize num);

  virtual int_type uflow();
  virtual int_type underflow();
  virtual std::streamsize xsgetn (char* s, std::streamsize num);
private:
  std::streambuf* originalBuf;
  Modes      mode;
};

class vofstream : public std::ostream {
public:
  vofstream();
  explicit vofstream(const char * filename,
                     openmode mode = out | trunc );
  ~vofstream();

  std::streambuf* rdbuf() const;
  bool is_open();
  void open(const char* filename, openmode mode = out | trunc);
  void close();
private:
  std::filebuf fbuf;
  ParBuf  mybuf;
};

class vifstream : public std::istream {
public:
  vifstream();
  explicit vifstream(const char * filename,
                     openmode mode = in );
  ~vifstream();

  std::streambuf* rdbuf() const;
  bool is_open();
  void open(const char* filename, openmode mode = in);
  void close();
private:
  std::filebuf fbuf;
  ParBuf  mybuf;
};

class vfstream : public std::iostream {
public:
  vfstream();
  explicit vfstream(const char * filename,
                    openmode mode = in | out );
  ~vfstream();

  std::streambuf* rdbuf() const;
  bool is_open();
  void open(const char* filename, openmode mode = in | out);
  void close();
private:
  std::filebuf fbuf;
  ParBuf  mybuf;
};

} // namespace olb

#endif

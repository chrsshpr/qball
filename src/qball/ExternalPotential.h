////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory.
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008
// LLNL-CODE-635376. All rights reserved.
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// ExternalPotential.h
//
////////////////////////////////////////////////////////////////////////////////
#ifndef EXTERNALPOTENTIAL_H
#define EXTERNALPOTENTIAL_H

#include <vector>
#include <string>
#include "Sample.h"
#include "ChargeDensity.h"

class Sample;

class ExternalPotential
{
  private:

  Sample& s_;
  int n_[3];               // real space grid size in 3 dimensions
                           // read from cube file in cube file mode,
                           // otherwise must be given in constructor
  double ecut_;
  double magnitude_;       // the magnitude of external potential, defined as
                           // the average of its largest 0.1% (absolute) values
  double amplitude_;       // overall scaling factor of external potential
  vector<double> vext_r_;  // vext in real space
  std::string filename_;   // file name for external potential
  std::string fmt_;        // file format: "cube" or "xml"

  public:

  ExternalPotential(Sample& s, std::string name, std::string fmt="xml"):
    s_(s), filename_(name), ecut_(0.0), amplitude_(1.0), magnitude_(0.0),
    fmt_(fmt)
  {
    assert( fmt_ == "cube" || fmt_ == "xml" );
  }
  ~ExternalPotential() {}

  int n(int i) const { return n_[i]; }
  double ecut(void) const { return ecut_; }
  double magnitude(void) const { return magnitude_; }
  double amplitude(void) const { return amplitude_; }
  std::string filename(void) const { return filename_; }
  double v(size_t i) const { return amplitude_ * vext_r_[i]; }
  void update(const ChargeDensity& cd, const Context& my_col_ctxt);
  void set_amplitude(double a) { amplitude_ = a; }
  void reverse(void) {amplitude_ *= -1; }
  double compute_eext(const ChargeDensity& cd);
};
#endif

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
// TDNaturalOrbital.h 
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef TDNATURALORBITAL_H
#define TDNATURALORBITAL_H

#include "Sample.h"
#include "FourierTransform.h"
#include <vector>
#include <complex>
#include <math/matrix.h>
class TDNaturalOrbital
{
  private:
  const Sample & s_;
  const ComplexMatrix& ref_;
  const Context& ctxt_;
  const Basis& basis;
  int n,nb,m,mb;
  int np0,np1,np2;
  valarray<double> nto_; 
  ComplexMatrix update_;
  ComplexMatrix nto_coeff;
  ComplexMatrix hole_coeff;
  ComplexMatrix elec_coeff; 
  FourierTransform* ft;

  public:
  void update(const ComplexMatrix& mat);
  void update_elec();
  void update_hole();
  void update_NTO();
  void print_hole_orbital(int m,string filename);
  void print_elec_orbital(int m,string filename);
  void print_nto_orbital(int m,string filename);
  void print_orbital(double* wftmp,string filename);

  double nto(int n) {return nto_[n];};
  TDNaturalOrbital(const Sample& s);
  ~TDNaturalOrbital(void);
};
#endif

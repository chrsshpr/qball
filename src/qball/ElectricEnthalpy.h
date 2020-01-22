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
// ElectricEnthalpy.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ELECTRICENTHALPY_H
#define ELECTRICENTHALPY_H

#include <vector>
#include <complex>
#include <cassert>
#include <string>
#include <math/blas.h>
#include <math/matrix.h>
#include <math/d3vector.h>
#include "Wavefunction.h"
#include "SlaterDet.h"
#include "Context.h"
#include "Sample.h"
#include "Timer.h"
#include "Basis.h"

class Sample;
class MLWFTransform;
class TDMLWFTransform;

class ElectricEnthalpy
{
  private:

  const Sample& s_;
  const Wavefunction& wf_;
  Wavefunction* dwf_;
  SlaterDet& sd_;
  const Context& ctxt_;
  const Basis& basis_;

  bool onpe0_;
  bool finite_field_;

  enum { off, berry, mlwf, tdmlwf, mlwf_ref, mlwf_ref_q } pol_type_;
  bool compute_quadrupole_;

  // electric field
  D3vector e_field_;

  Wavefunction* rwf_[3];

  // MLWFtransform is used to compute S matrix
  TDMLWFTransform* tdmlwft_;
  MLWFTransform* mlwft_;

  // s matrices
  ComplexMatrix* smat_[3];

  // total, ionic and electronic part of macroscopic polarization
  D3vector dipole_total_, dipole_ion_, dipole_el_;

  // electric enthalpy
  double enthalpy_;

  std::vector <D3vector> mlwfc_;
  std::vector <double> mlwfs_;
  std::vector <D3vector> correction_;
  //std::vector <D3tensor> quad_;

  void compute_correction(void);
  double vsst(double x) const;

  public:

  mutable TimerMap tmap;

  D3vector e_field(void) const { return e_field_; }
  D3vector dipole_total(void) const { return dipole_total_; }
  D3vector dipole_ion(void) const { return dipole_ion_; }
  D3vector dipole_el(void) const { return dipole_el_; }

  double enthalpy(Wavefunction& dwf, bool compute_hpsi);

  void update(void);
  void print(std::ostream& os) const;

  ElectricEnthalpy( const Sample& s);
  ~ElectricEnthalpy(void);
};
std::ostream& operator << ( std::ostream& os, const ElectricEnthalpy& e );
#endif

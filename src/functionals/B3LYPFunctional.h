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
////////////////////////////////////////////////////////////////////////////////
//
// B3LYPFunctional.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef B3LYPFUNCTIONAL_H
#define B3LYPFUNCTIONAL_H

#include "XCFunctional.h"
#include <vector>
using namespace std; 

class B3LYPFunctional : public XCFunctional
{
  B3LYPFunctional();

  std::vector<double> _exc, _exc_up, _exc_dn;
  std::vector<double> _vxc1, _vxc1_up, _vxc1_dn,
                 _vxc2, _vxc2_upup, _vxc2_updn, _vxc2_dnup, _vxc2_dndn;
  std::vector<double> _grad_rho[3], _grad_rho_up[3], _grad_rho_dn[3];

  void excb3lyp(double rho, double grad,
    double *exc, double *vxc1, double *vxc2);

  void excb3lyp_sp(double rho_up, double rho_dn,
    double grad_up2, double grad_dn2, double grad_up_grad_dn,
    double *exc_up, double *exc_dn, double *vxc1_up, double *vxc1_dn,
    double *vxc2_upup, double *vxc2_dndn, double *vxc2_updn, double *vxc2_dnup);


  // BLYP ex and ec functions 
  void exb88(double rho, double grad,
    double *ex, double *vx1, double *vx2);
  void eclyp(double rho, double grad,
    double *ec, double *vc1, double *vc2);

  void exb88_sp(double rho_up, double rho_dn,
    double grad_up2, double grad_dn2, double grad_up_grad_dn,
    double *ex_up, double *ex_dn, double *vx1_up, double *vx1_dn,
    double *vx2_upup, double *vx2_dndn, double *vx2_updn, double *vx2_dnup);
   void eclyp_sp(double rho_up, double rho_dn,
    double grad_up2, double grad_dn2, double grad_up_grad_dn,
    double *ec_up, double *ec_dn, double *vc1_up, double *vc1_dn,
    double *vc2_upup, double *vc2_dndn, double *vc2_updn, double *vc2_dnup);
   
   // VWN ex and ex functions 
   void alpha_c(const double rh, double &a, double &da);
   void x_unpolarized(const double rh, double &ex, double &vx);
   void c_unpolarized(const double rh, double &ec, double &vc);
   void x_polarized(const double rh, double &ex, double &vx);
   void c_polarized(const double rh, double &ec, double &vc);

   void exvwn(const double rh, double &ex, double &vx);
   void ecvwn(const double rh, double &ec, double &vc);

   void exvwn_sp(double roe_up, double roe_dn,
    double &ex, double &vx_up, double &vx_dn);
   void ecvwn_sp(double roe_up, double roe_dn,
    double &ec, double &vc_up, double &vc_dn);

  public:

  B3LYPFunctional(const std::vector<std::vector<double> > &rhoe);

  bool isGGA() { return true; };
  std::string name() { return "B3LYP"; };
  void setxc(void);
};
#endif

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
// BLYPFunctional.C
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <cmath>
#include <cassert>
#include "BHLYPFunctional.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
BHLYPFunctional::BHLYPFunctional(const vector<vector<double> > &rhoe) {
  _nspin = rhoe.size();
  if ( _nspin > 1 ) assert(rhoe[0].size() == rhoe[1].size());
  _np = rhoe[0].size();
 
  if ( _nspin == 1 )
  {
    _exc.resize(_np);
    _vxc1.resize(_np);
    _vxc2.resize(_np);
    _grad_rho[0].resize(_np);
    _grad_rho[1].resize(_np);
    _grad_rho[2].resize(_np);
    rho = &rhoe[0][0];
    grad_rho[0] = &_grad_rho[0][0];
    grad_rho[1] = &_grad_rho[1][0];
    grad_rho[2] = &_grad_rho[2][0];
    exc = &_exc[0];
    vxc1 = &_vxc1[0];
    vxc2 = &_vxc2[0];
  }
  else
  {
    _exc_up.resize(_np);
    _exc_dn.resize(_np);
    _vxc1_up.resize(_np);
    _vxc1_dn.resize(_np);
    _vxc2_upup.resize(_np);
    _vxc2_updn.resize(_np);
    _vxc2_dnup.resize(_np);
    _vxc2_dndn.resize(_np);
    _grad_rho_up[0].resize(_np);
    _grad_rho_up[1].resize(_np);
    _grad_rho_up[2].resize(_np);
    _grad_rho_dn[0].resize(_np);
    _grad_rho_dn[1].resize(_np);
    _grad_rho_dn[2].resize(_np);

    rho_up = &rhoe[0][0];
    rho_dn = &rhoe[1][0];
    grad_rho_up[0] = &_grad_rho_up[0][0];
    grad_rho_up[1] = &_grad_rho_up[1][0];
    grad_rho_up[2] = &_grad_rho_up[2][0];
    grad_rho_dn[0] = &_grad_rho_dn[0][0];
    grad_rho_dn[1] = &_grad_rho_dn[1][0];
    grad_rho_dn[2] = &_grad_rho_dn[2][0];
    exc_up = &_exc_up[0];
    exc_dn = &_exc_dn[0];
    vxc1_up = &_vxc1_up[0];
    vxc1_dn = &_vxc1_dn[0];
    vxc2_upup = &_vxc2_upup[0];
    vxc2_updn = &_vxc2_updn[0];
    vxc2_dnup = &_vxc2_dnup[0];
    vxc2_dndn = &_vxc2_dndn[0];
  }
}

////////////////////////////////////////////////////////////////////////////////
void BHLYPFunctional::setxc(void) 
{
  if ( _np == 0 ) return;
  if ( _nspin == 1 )
  {
    assert( rho != 0 );
    assert( grad_rho[0] != 0 && grad_rho[1] != 0 && grad_rho[2] != 0 );
    assert( exc != 0 );
    assert( vxc1 != 0 );
    assert( vxc2 != 0 );

    for ( int i = 0; i < _np; i++ )
    {
      double grad = sqrt(grad_rho[0][i]*grad_rho[0][i] +
                         grad_rho[1][i]*grad_rho[1][i] +
                         grad_rho[2][i]*grad_rho[2][i] );

      excbhlyp(rho[i],grad,&exc[i],&vxc1[i],&vxc2[i]);

    }
  }
  else
  {
    assert( rho_up != 0 );
    assert( rho_dn != 0 );
    assert( grad_rho_up[0] != 0 && grad_rho_up[1] != 0 && grad_rho_up[2] != 0 );
    assert( grad_rho_dn[0] != 0 && grad_rho_dn[1] != 0 && grad_rho_dn[2] != 0 );
    assert( exc_up != 0 );
    assert( exc_dn != 0 );
    assert( vxc1_up != 0 );
    assert( vxc1_dn != 0 );
    assert( vxc2_upup != 0 );
    assert( vxc2_updn != 0 );
    assert( vxc2_dnup != 0 );
    assert( vxc2_dndn != 0 );

    for ( int i = 0; i < _np; i++ )
    {
      double grx_up = grad_rho_up[0][i];
      double gry_up = grad_rho_up[1][i];
      double grz_up = grad_rho_up[2][i];
      double grx_dn = grad_rho_dn[0][i];
      double gry_dn = grad_rho_dn[1][i];
      double grz_dn = grad_rho_dn[2][i];
      double grx = grx_up + grx_dn;
      double gry = gry_up + gry_dn;
      double grz = grz_up + grz_dn;
      double grad_up = sqrt(grx_up*grx_up + gry_up*gry_up + grz_up*grz_up);
      double grad_dn = sqrt(grx_dn*grx_dn + gry_dn*gry_dn + grz_dn*grz_dn);
      double grad    = sqrt(grx*grx + gry*gry + grz*grz);
      excbhlyp(rho_up[i],rho_dn[i],grad_up,grad_dn,grad,&exc_up[i],&exc_dn[i],
                &vxc1_up[i],&vxc1_dn[i],&vxc2_upup[i],&vxc2_dndn[i],
                &vxc2_updn[i], &vxc2_dnup[i]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void BHLYPFunctional::excbhlyp(double rho, double grad,
  double *exc, double *vxc1, double *vxc2)
{
  // BHandHLYP unpolarized
  // Coefficients of the BHandHLYP functional
  // A. Becke, JCP 98, 5648 (1993)
  // See also X.Xu and W. Goddard, J.Phys.Chem. A108, 2305 (2004)
  // EcLSDA is the Vosko-Wilk-Nusair correlation energy
  // dExBecke88 is the difference ExB88 - ExLDA
  // B3LYP: 0.2 ExHF + 0.8 ExSlater + 0.19 EcLSDA + 0.81 EcLYP + 0.72 dExBecke88
  // BHandHLYP: 0.5 ExHF + 0.5 ExB88 + 1.0 EcLYP
  const double xb88_coeff = 0.5; // Becke88 exchange
  const double clyp_coeff = 1.0 ;

  double ex_b88,vx1_b88,vx2_b88;
  double ec_lyp,vc1_lyp,vc2_lyp;

  exb88(rho,grad,&ex_b88,&vx1_b88,&vx2_b88);
  eclyp(rho,grad,&ec_lyp,&vc1_lyp,&vc2_lyp);

  *exc = xb88_coeff  * ex_b88 +
         clyp_coeff  * ec_lyp;

  *vxc1 = xb88_coeff * vx1_b88 +
          clyp_coeff * vc1_lyp;

  *vxc2 = xb88_coeff * vx2_b88 +
          clyp_coeff * vc2_lyp;
}

////////////////////////////////////////////////////////////////////////////////
void BHLYPFunctional::excbhlyp(double rho_up, double rho_dn,
  double grad_up2, double grad_dn2, double grad_up_grad_dn,
  double *exc_up, double *exc_dn, double *vxc1_up, double *vxc1_dn,
  double *vxc2_upup, double *vxc2_dndn, double *vxc2_updn, double *vxc2_dnup)
{
  // BHandHLYP spin-polarized
  // Coefficients of the BHandHLYP functional
  // A. Becke, JCP 98, 5648 (1993)
  // See also X.Xu and W. Goddard, J.Phys.Chem. A108, 2305 (2004)
  // EcLSDA is the Vosko-Wilk-Nusair correlation energy
  // dExBecke88 is the difference ExB88 - ExLDA
  // B3LYP: 0.2 ExHF + 0.8 ExSlater + 0.19 EcLSDA + 0.81 EcLYP + 0.72 dExBecke88
  // BHandHLYP: 0.5 ExHF + 0.5 ExB88 + 1.0 EcLYP
  const double xb88_coeff = 0.5; // Becke88 exchange
  const double clyp_coeff = 1.0;

  double ex_b88_up,ex_b88_dn,vx1_b88_up,vx1_b88_dn,
         vx2_b88_upup,vx2_b88_dndn,vx2_b88_updn,vx2_b88_dnup;
  double ec_lyp_up,ec_lyp_dn,vc1_lyp_up,vc1_lyp_dn,
         vc2_lyp_upup,vc2_lyp_dndn,vc2_lyp_updn,vc2_lyp_dnup;

  exb88_sp(rho_up,rho_dn,grad_up2,grad_dn2,grad_up_grad_dn,
    &ex_b88_up,&ex_b88_dn,&vx1_b88_up,&vx1_b88_dn,
    &vx2_b88_upup,&vx2_b88_dndn,&vx2_b88_updn,&vx2_b88_dnup);
  eclyp_sp(rho_up,rho_dn,grad_up2,grad_dn2,grad_up_grad_dn,
    &ec_lyp_up,&ec_lyp_dn,&vc1_lyp_up,&vc1_lyp_dn,
    &vc2_lyp_upup,&vc2_lyp_dndn,&vc2_lyp_updn,&vc2_lyp_dnup);

  *exc_up = xb88_coeff  * ex_b88_up +
            clyp_coeff  * ec_lyp_up;

  *exc_dn = xb88_coeff  * ex_b88_dn +
            clyp_coeff  * ec_lyp_dn;

  *vxc1_up = xb88_coeff * vx1_b88_up +
             clyp_coeff * vc1_lyp_up;

  *vxc1_dn = xb88_coeff * vx1_b88_dn +
             clyp_coeff * vc1_lyp_dn;

  *vxc2_upup = xb88_coeff * vx2_b88_upup + clyp_coeff * vc2_lyp_upup;
  *vxc2_dndn = xb88_coeff * vx2_b88_dndn + clyp_coeff * vc2_lyp_dndn;
  *vxc2_updn = xb88_coeff * vx2_b88_updn + clyp_coeff * vc2_lyp_updn;
  *vxc2_dnup = xb88_coeff * vx2_b88_dnup + clyp_coeff * vc2_lyp_dnup;
}

////////////////////////////////////////////////////////////////////////////////
void BHLYPFunctional::exb88(double rho, double grad,
  double *ex, double *vx1, double *vx2)
{
  // Becke exchange constants
  const double fourthirds = 4.0 / 3.0;
  const double beta=0.0042;
  const double axa = -0.9305257363490999;   // -1.5*pow(3.0/(4*pi),third)

  *ex = 0.0;
  *vx1 = 0.0;
  *vx2 = 0.0;

  if ( rho < 1.e-18 ) return;

  // Becke's exchange
  // A.D.Becke, Phys.Rev. B38, 3098 (1988)

  const double rha = 0.5 * rho;
  const double grada = 0.5 * grad;

  const double rha13 = cbrt(rha);
  const double rha43 = rha * rha13;
  const double xa = grada / rha43;
  const double xa2 = xa*xa;
  const double asinhxa = asinh(xa);
  const double frac = 1.0 / ( 1.0 + 6.0 * beta * xa * asinhxa );
  const double ga = axa - beta * xa2 * frac;
  // in next line, ex is the energy density, hence rh13
  *ex = rha13 * ga;

  // potential
  const double gpa = ( 6.0*beta*beta*xa2 * ( xa/sqrt(xa2+1.0) - asinhxa )
                     - 2.0*beta*xa ) * frac*frac;
  *vx1 = rha13 * fourthirds * ( ga - xa * gpa );
  *vx2 = - 0.5 * gpa / grada;
}

////////////////////////////////////////////////////////////////////////////////
void BHLYPFunctional::eclyp(double rho, double grad,
  double *ec, double *vc1, double *vc2)
{
  // LYP constants
  const double a = 0.04918;
  const double b = 0.132;
  const double ab36 = a * b / 36.0;
  const double c = 0.2533;
  const double c_third = c / 3.0;
  const double d = 0.349;
  const double d_third = d / 3.0;
  const double cf = 2.87123400018819; // (3/10)*pow(3*pi*pi,2/3)
  const double cfb = cf * b;

  *ec = 0.0;
  *vc1 = 0.0;
  *vc2 = 0.0;

  if ( rho < 1.e-18 ) return;

  // LYP correlation
  // Phys. Rev. B 37, 785 (1988).
  // next lines specialized to the unpolarized case
  const double rh13 = cbrt(rho);
  const double rhm13 = 1.0 / rh13;
  const double rhm43 = rhm13 / rho;
  const double e = exp ( - c * rhm13 );
  const double num = 1.0 + cfb * e;
  const double den = 1.0 + d * rhm13;
  const double deninv = 1.0 / den;
  const double cfrac = num * deninv;

  const double delta = rhm13 * ( c + d * deninv );
  const double ddelta = - (1.0/3.0) * ( c * rhm43
                      + d * rhm13 * rhm13 / ((d+rh13)*(d+rh13)) );
  const double rhm53 = rhm43 * rhm13;
  const double t1 = e * deninv;
  const double t2 = rhm53;
  const double t3 = 6.0 + 14.0 * delta;

  const double g = ab36 * t1 * t2 * t3;

  // ec is the energy density, hence divide the energy by rho
  *ec = - a * cfrac + 0.25 * g * grad * grad / rho;

  // potential
  const double de = c_third * rhm43 * e;
  const double dnum = cfb * de;
  const double dden = - d_third * rhm43;
  const double dfrac = ( dnum * den - dden * num ) * deninv * deninv;

  const double dt1 = de * deninv - e * dden * deninv * deninv;
  const double dt2 = - (5.0/3.0) * rhm53/rho;
  const double dt3 = 14.0 * ddelta;

  const double dg = ab36 * ( dt1 * t2 * t3 + t1 * dt2 * t3 + t1 * t2 * dt3 );

  *vc1 = - a * ( cfrac + rho * dfrac ) + 0.25 * dg * grad * grad;
  *vc2 = -0.5 * g;
}

////////////////////////////////////////////////////////////////////////////////
void BHLYPFunctional::exb88_sp(double rho_up, double rho_dn,
  double grad_up2, double grad_dn2, double grad_up_grad_dn,
  double *ex_up, double *ex_dn, double *vx1_up, double *vx1_dn,
  double *vx2_upup, double *vx2_dndn, double *vx2_updn, double *vx2_dnup)
{
  *ex_up = 0.0;
  *ex_dn = 0.0;
  *vx1_up = 0.0;
  *vx1_dn = 0.0;
  *vx2_upup = 0.0;
  *vx2_updn = 0.0;
  *vx2_dnup = 0.0;
  *vx2_dndn = 0.0;

  if ( rho_up < 1.e-18 && rho_dn < 1.e-18 ) return;

  // Becke exchange constants
  const double fourthirds = 4.0 / 3.0;
  const double beta = 0.0042;
  const double cx = -0.9305257363490999;    // -1.5*pow(3.0/(4*pi),third)

  // Becke's exchange
  // A.D.Becke, Phys.Rev. B38, 3098 (1988)

  double ex1a = 0.0;
  double ex1b = 0.0;
  double vx1a = 0.0;
  double vx1b = 0.0;
  double vx2a = 0.0;
  double vx2b = 0.0;

  if ( rho_up > 1.e-18 )
  {
    const double& rha = rho_up;
    const double rha13 = cbrt(rha);
    const double rha43 = rha * rha13;
    const double grada = sqrt(grad_up2);
    const double xa = grada / rha43;
    const double xa2 = xa*xa;
    const double asinhxa = asinh(xa);
    const double fraca = 1.0 / ( 1.0 + 6.0 * beta * xa * asinhxa );
    const double ga = cx - beta * xa2 * fraca;
    // next line, ex is the energy density, hence rh13
    ex1a = rha13 * ga;
    const double gpa = ( 6.0*beta*beta*xa2 * ( xa/sqrt(xa2+1.0) - asinhxa )
                       - 2.0*beta*xa ) * fraca * fraca;
    vx1a = rha13 * fourthirds * ( ga - xa * gpa );
    vx2a = - gpa / grada;
  }

  if ( rho_dn > 1.e-18 )
  {
    const double& rhb = rho_dn;
    const double rhb13 = cbrt(rhb);
    const double rhb43 = rhb * rhb13;
    const double gradb = sqrt(grad_dn2);
    const double xb = gradb / rhb43;
    const double xb2 = xb*xb;
    const double asinhxb = asinh(xb);
    const double fracb = 1.0 / ( 1.0 + 6.0 * beta * xb * asinhxb );
    const double gb = cx - beta * xb2 * fracb;
    // next line, ex is the energy density, hence rh13
    ex1b = rhb13 * gb;
    const double gpb = ( 6.0*beta*beta*xb2 * ( xb/sqrt(xb2+1.0) - asinhxb )
                       - 2.0*beta*xb ) * fracb * fracb;
    vx1b = rhb13 * fourthirds * ( gb - xb * gpb );
    vx2b = - gpb / gradb;
  }

  *ex_up = ex1a;
  *ex_dn = ex1b;
  *vx1_up = vx1a;
  *vx1_dn = vx1b;
  *vx2_upup = vx2a;
  *vx2_updn = 0.0;
  *vx2_dnup = 0.0;
  *vx2_dndn = vx2b;
}

////////////////////////////////////////////////////////////////////////////////
void BHLYPFunctional::eclyp_sp(double rho_up, double rho_dn,
  double grad_up2, double grad_dn2, double grad_up_grad_dn,
  double *ec_up, double *ec_dn, double *vc1_up, double *vc1_dn,
  double *vc2_upup, double *vc2_dndn, double *vc2_updn, double *vc2_dnup)
{
  *ec_up = 0.0;
  *ec_dn = 0.0;
  *vc1_up = 0.0;
  *vc1_dn = 0.0;
  *vc2_upup = 0.0;
  *vc2_updn = 0.0;
  *vc2_dnup = 0.0;
  *vc2_dndn = 0.0;

  if ( rho_up < 1.e-18 && rho_dn < 1.e-18 ) return;
  if ( rho_up + rho_dn < 1.e-18 ) return;
  if ( rho_up < 0.0 ) rho_up = 0.0;
  if ( rho_dn < 0.0 ) rho_dn = 0.0;

  // LYP constants
  const double a = 0.04918;
  const double b = 0.132;
  const double c = 0.2533;
  const double d = 0.349;
  const double cf = 2.87123400018819; // (3/10)*pow(3*pi*pi,2/3)
  const double ninth = 1.0 / 9.0;
  const double two113 = pow(2.0,11.0/3.0);
  const double fourthirds = 4.0 / 3.0;

  const double& rha = rho_up;
  const double& rhb = rho_dn;

  const double rha13 = cbrt(rha);
  const double rhb13 = cbrt(rhb);
  const double rha43 = rha * rha13;
  const double rhb43 = rhb * rhb13;
  const double rha83 = rha43 * rha43;
  const double rhb83 = rhb43 * rhb43;
  const double rha113 = rha83 * rha;
  const double rhb113 = rhb83 * rhb;

  // LYP correlation
  // Phys. Rev. B 37, 785 (1988).

  const double rho = rha + rhb;
  const double rhoinv = 1.0/rho;
  const double rhab = rha * rhb;
  const double rhabrhm = rhab * rhoinv;
  const double rh13 = cbrt(rho);
  const double rhm13 = 1.0 / rh13;
  const double rhm43 = rhm13 / rho;
  const double rhm113 = rhm43 * rhm43 * rhm43 * rh13;
  const double rhm143 = rhm113 * rhoinv;
  const double e = exp ( - c * rhm13 );
  const double den = 1.0 + d * rhm13;
  const double deninv = 1.0 / den;

  const double w = e * deninv * rhm113;
  const double dw = (1.0/3.0)*(c*d+(c-10.0*d)*rh13-11.0*rh13*rh13)*e*rhm143/
                    ((d+rh13)*(d+rh13));
  const double abw = a * b * w;
  const double f1 = -4.0 * a * deninv * rhabrhm;
  const double f2 = - two113 * cf * abw * rha * rhb  * ( rha83 + rhb83 );

  const double delta = rhm13 * ( c + d * deninv );
  const double ddelta = - (1.0/3.0) * ( c * rhm43
                      + d * rhm13 * rhm13 / ((d+rh13)*(d+rh13)) );
  const double taa = 1.0 - 3.0 * delta + ( 11.0 - delta ) * rha * rhoinv;
  const double tbb = 1.0 - 3.0 * delta + ( 11.0 - delta ) * rhb * rhoinv;

  const double dtaa_drha = -3.0*ddelta - ddelta*rha*rhoinv +
                           (11.0-delta)*rhoinv*(1.0-rha*rhoinv);
  const double dtaa_drhb = -3.0*ddelta - ddelta*rha*rhoinv -
                           (11.0-delta)*rhoinv*rha*rhoinv;
  const double dtbb_drha = -3.0*ddelta - ddelta*rhb*rhoinv -
                           (11.0-delta)*rhoinv*rhb*rhoinv;
  const double dtbb_drhb = -3.0*ddelta - ddelta*rhb*rhoinv +
                           (11.0-delta)*rhoinv*(1.0-rhb*rhoinv);
  const double gaa = - abw * ( rhab * ninth * taa - rhb * rhb );
  const double gbb = - abw * ( rhab * ninth * tbb - rha * rha );
  const double gab = - abw * ( rhab * ninth * ( 47.0 - 7.0 * delta )
                               - fourthirds * rho * rho );

  // next line, ec is the energy density, hence divide the energy by rho
  const double ec1a = ( f1 + f2 + gaa * grad_up2
                        + gab * grad_up_grad_dn
                        + gbb * grad_dn2 ) / rho;
  const double ec1b = ( f1 + f2 + gaa * grad_up2
                        + gab * grad_up_grad_dn
                        + gbb * grad_dn2 ) / rho;

  const double A = -two113*cf*a*b;
  const double df1_drha = -fourthirds*a*d*deninv*deninv*rhm43*rhabrhm
                        - 4.0*a*deninv*rhoinv*(rhb-rhabrhm);
  const double df2_drha = A*dw*(rha113*rhb+rha*rhb113)
                        + A*w*((11.0/3.0)*rha83*rhb+rhb113);

  const double df1_drhb = -fourthirds*a*d*deninv*deninv*rhm43*rhabrhm
                        - 4.0*a*deninv*rhoinv*(rha-rhabrhm);
  const double df2_drhb = A*dw*(rha113*rhb+rha*rhb113)
                        + A*w*(rha113+rha*(11.0/3.0)*rhb83);

  const double dgaa_drha = -a*b*dw*(rhab*ninth*taa-rhb*rhb)
                         - abw*(rhb*ninth*taa+rhab*ninth*dtaa_drha);
  const double dgaa_drhb = -a*b*dw*(rhab*ninth*taa-rhb*rhb)
                         - abw*(rha*ninth*taa+rhab*ninth*dtaa_drhb-2.0*rhb);

  const double dgbb_drha = -a*b*dw*(rhab*ninth*tbb-rha*rha)
                         - abw*(rhb*ninth*tbb+rhab*ninth*dtbb_drha-2.0*rha);
  const double dgbb_drhb = -a*b*dw*(rhab*ninth*tbb-rha*rha)
                         - abw*(rha*ninth*tbb+rhab*ninth*dtbb_drhb);

  const double dgab_drha = -a*b*dw*(rhab*ninth*(47.0-7.0*delta)
                                    -fourthirds*rho*rho)
                           -abw*(rhb*ninth*(47.0-7.0*delta)
                                 -rhab*ninth*7.0*ddelta-2.0*fourthirds*rho);
  const double dgab_drhb = -a*b*dw*(rhab*ninth*(47.0-7.0*delta)
                                    -fourthirds*rho*rho)
                           -abw*(rha*ninth*(47.0-7.0*delta)
                                 -rhab*ninth*7.0*ddelta-2.0*fourthirds*rho);

  const double vc1a = df1_drha + df2_drha
                    + dgaa_drha * grad_up2
                    + dgab_drha * grad_up_grad_dn
                    + dgbb_drha * grad_dn2;
  const double vc1b = df1_drhb + df2_drhb
                    + dgaa_drhb * grad_up2
                    + dgab_drhb * grad_up_grad_dn
                    + dgbb_drhb * grad_dn2;

  *ec_up = ec1a;
  *ec_dn = ec1b;
  *vc1_up = vc1a;
  *vc1_dn = vc1b;
  *vc2_upup = - 2.0 * gaa;
  *vc2_updn = - gab;
  *vc2_dnup = - gab;
  *vc2_dndn = - 2.0 * gbb;
}

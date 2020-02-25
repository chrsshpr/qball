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
// ElectricEnthalpy.C
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
#include <cmath>
#include <algorithm> // fill

#include "Timer.h"
#include "Context.h"
#include <math/matrix.h>
#include "Sample.h"
#include <math/d3vector.h>
#include "ElectricEnthalpy.h"
#include "MLWFTransform.h"
#include "TDMLWFTransform.h"
#include "Wavefunction.h"
#include "Basis.h"
#include "SlaterDet.h"
#include "ChargeDensity.h"
#include "FourierTransform.h"
#include "UnitCell.h"
using namespace std;
///////////////////////////////////////////////////////////////////////////////
double ElectricEnthalpy::vsst(double x) const
{
  // smooth sawtooth periodic potential function
  // x in [-1/2, 1/2]
  // The slope of vsst is 1 at x=0
  //
  // The function vsst approximates the identity function x->x
  // in the interval [-1/2+xcut, 1/2-xcut]
  const double xcut = 0.05;
  const double xcut2 = xcut*xcut;
  // The function vsst is well represented by a
  // discrete Fourier transform of length np = 2*ng
  // Some aliasing error will arise if np < 2*ng
  // The error is determined by the product xcut*ng
  const int ng = 32;
  double v = 0.0;
  for ( int ig = 1; ig < ng; ig++ )
  {
    const double g = 2 * M_PI * ig;
    const double arg = -0.25 * g * g * xcut2;
    // next line: factor sgn to shift origin by 0.5
    const int sgn = 1 - 2*(ig%2);
    const double c = -2.0 * sgn * exp ( arg ) / g;
    v += c * sin(x*g);
  }
  return v;
}
///////////////////////////////////////////////////////////////////////////////
ElectricEnthalpy::ElectricEnthalpy(const Sample& s): s_(s), wf_(s.wf), 
  sd_(*(s.wf.sd(0,0))), ctxt_(s.wf.sd(0,0)->context()),
  basis_(s.wf.sd(0,0)->basis())
{
  onpe0_ = ctxt_.onpe0();
  e_field_ = s.ctrl.e_field;
  finite_field_ = norm(e_field_) != 0.0;


  if ( s.ctrl.polarization == "OFF" )
    pol_type_ = off;
  else if ( s.ctrl.polarization == "BERRY" )
    pol_type_ = berry;
  else if ( s.ctrl.polarization == "MLWF" )
    {
     pol_type_ = mlwf;
    }
  else if ( s.ctrl.polarization == "TDMLWF" )
    {
     pol_type_ = tdmlwf;
    }
  else if ( s.ctrl.polarization == "MLWF_REF" )
    pol_type_ = mlwf_ref;
  else if ( s.ctrl.polarization == "TDMLWF_REF")
    pol_type_ = tdmlwf_ref;
  else if ( s.ctrl.polarization == "MLWF_REF_Q" )
  {
    pol_type_ = mlwf_ref;
  }
  else
  {
    cerr << "ElectricEnthalpy: invalid polarization option" << endl;
    ctxt_.abort(1);
  }

  // do not allocate further objects if polarization is off
  if ( pol_type_ == off ) return;

  assert(wf_.nkp()==1);
  assert(wf_.nspin()==1);

  dwf_ = new Wavefunction(s.wf); 
  tdmlwft_ = new TDMLWFTransform(sd_);
  mlwft_ = new MLWFTransform(sd_);

  smat_[0] = smat_[1] = smat_[2] = 0;
  rwf_[0] = rwf_[1] = rwf_[2] = 0;
  int nst = sd_.nst();

  if ( pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q || pol_type_ == tdmlwf_ref )
  {
    // allocate real space wf arrays for MLWF refinement  
    for ( int idir = 0; idir < 3; idir++ )
      rwf_[idir] = new Wavefunction(wf_);
    correction_.resize(nst);
  }
  else if ( pol_type_ == berry )
  {
    // allocate complex Berry phase matrix
    int n = sd_.c().n();
    int nb = sd_.c().nb();
    for ( int idir = 0; idir < 3; idir++ )
      smat_[idir] = new ComplexMatrix(ctxt_,n,n,nb,nb);
  }

  if ( (pol_type_ != off) && onpe0_ )
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout.setf(ios::right,ios::adjustfield);
    cout.precision(8);
    cout << "<e_field> " << e_field_ << " </e_field>" << endl;
  }

  mlwfc_.resize(nst);
  mlwfs_.resize(nst);
}

///////////////////////////////////////////////////////////////////////////////
ElectricEnthalpy::~ElectricEnthalpy(void)
{
  if ( pol_type_ == off ) return;

  delete dwf_;
  delete mlwft_;
  for ( int idir = 0; idir < 3; idir++ )
  {
    delete rwf_[idir];
    delete smat_[idir];
  }

  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    s_.ctxt_.dmin(1,1,&tmin,1);
    s_.ctxt_.dmax(1,1,&tmax,1);
    if ( pol_type_ != off && s_.ctxt_.myproc()==0 )
    {
      string s = "name=\"" + (*i).first + "\"";
      cout << "<timing " << left << setw(22) << s
           << " min=\"" << setprecision(3) << tmin << "\""
           << " max=\"" << setprecision(3) << tmax << "\"/>"
           << endl;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::update(void)
{
  if ( pol_type_ == off ) return;

  const UnitCell& cell = sd_.basis().cell();
  // compute cos and sin matrices
  tmap["mlwf_update"].start();
  if (pol_type_ == tdmlwf || pol_type_ == tdmlwf_ref) tdmlwft_->update(); //CS
  if (pol_type_ ==mlwf)  mlwft_->update();
  tmap["mlwf_update"].stop();
  vector<SlaterDet*> sdcos(3), sdsin(3);
  if (pol_type_ == tdmlwf || pol_type_ == tdmlwf_ref) //CS
  {
  sdcos[0] = tdmlwft_->sdcosx();
  sdcos[1] = tdmlwft_->sdcosy();
  sdcos[2] = tdmlwft_->sdcosz();
  sdsin[0] = tdmlwft_->sdsinx();
  sdsin[1] = tdmlwft_->sdsiny();
  sdsin[2] = tdmlwft_->sdsinz();
  }
  else
  {
  sdcos[0] = mlwft_->sdcosx();
  sdcos[1] = mlwft_->sdcosy();
  sdcos[2] = mlwft_->sdcosz();
  sdsin[0] = mlwft_->sdsinx();
  sdsin[1] = mlwft_->sdsiny();
  sdsin[2] = mlwft_->sdsinz();
  }
  dipole_ion_ = s_.atoms.dipole();
  dipole_el_ = D3vector(0,0,0);

  if ( pol_type_ == mlwf || pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q || pol_type_ == tdmlwf || pol_type_ == tdmlwf_ref)
  {
    tmap["mlwf_trans"].start();
    if (pol_type_ == tdmlwf || pol_type_ == tdmlwf_ref) 
    { 
      tdmlwft_->compute_transform();
      tdmlwft_->apply_transform(sd_); //works w/ ETRS only. Comment out for FORKTD
    }
    else
    {
      mlwft_->compute_transform();
      mlwft_->apply_transform(sd_);
    }
    tmap["mlwf_trans"].stop();
    for ( int i = 0; i < sd_.nst(); i++ )
    {
      if (pol_type_ == tdmlwf || pol_type_ == tdmlwf_ref)
      {
        mlwfc_[i] = tdmlwft_->center(i);
        mlwfs_[i] = tdmlwft_->spread(i);
      }
      else
      {
        mlwfc_[i] = mlwft_->center(i);
        mlwfs_[i] = mlwft_->spread(i);
      }
    }

    if ( pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q || pol_type_ == tdmlwf_ref)
    {
      tmap["correction"].start();
      compute_correction();
      tmap["correction"].stop();
    }

    for ( int i = 0; i < sd_.nst(); i++ )
    {
      dipole_el_ -= 2.0 * mlwfc_[i];
      if ( pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q || pol_type_ == tdmlwf_ref) 
        dipole_el_ -= 2.0 * correction_[i]; //Use refined centers for dipole 
    }

    // compute gradient
    if ( finite_field_ )
    {
      dwf_->clear();
      for ( int idir = 0; idir < 3; idir++ )
      {
        if ( e_field_[idir] != 0.0 )
        {
          // MLWF part, modified for complex wavefunctions - DCY
	    if ( pol_type_ == mlwf )
	    {
              const double nst = sd_.nst();
              std::vector<double> adiag_inv_real(nst,0),adiag_inv_imag(nst,0);
              for ( int ist = 0; ist < nst; ist ++ )
              {
                const std::complex<double>
                adiag( mlwft_->adiag(idir*2,ist),mlwft_->adiag(idir*2+1,ist) );
                adiag_inv_real[ist] = real( std::complex<double>(1,0) / adiag );
                adiag_inv_imag[ist] = imag( std::complex<double>(1,0) / adiag );
              }

              DoubleMatrix ccos(sdcos[idir]->c());
              DoubleMatrix csin(sdsin[idir]->c());
              DoubleMatrix cp(dwf_->sd(0,0)->c());

              int nloc = cp.nloc();
              int mloc = cp.mloc();
	      
              int ione = 1;

              const double fac = length(cell.a(idir))
                               * e_field_[idir] / ( 2.0 * M_PI );

              for (int in = 0; in < nloc; in++)
              {
                int ist = cp.jglobal(in);
                double fac1 = adiag_inv_real[ist] * fac;
                double fac2 = adiag_inv_imag[ist] * fac;

                double *ptr1 = &cp[in*mloc],
                       *ptrcos = &ccos[in*mloc],
                       *ptrsin = &csin[in*mloc];

                daxpy(&mloc, &fac2, ptrcos, &ione, ptr1, &ione);
                daxpy(&mloc, &fac1, ptrsin, &ione, ptr1, &ione);

              }

	    } 
	    else if ( pol_type_ == tdmlwf)
	    {
              const double nst = sd_.nst();
              std::vector< complex<double> > adiag_inv_real(nst,0),adiag_inv_imag(nst,0);
              for ( int ist = 0; ist < nst; ist ++ )
              {
	         const std::complex<double>
                 adiag( real(tdmlwft_->adiag(idir*2,ist)),real(tdmlwft_->adiag(idir*2+1,ist)) );
	         const std::complex<double>
		 adiag2( imag(tdmlwft_->adiag(idir*2,ist)), imag(tdmlwft_->adiag(idir*2+1,ist)) );

                 const std::complex<double> adiag_real(real(tdmlwft_->adiag(idir*2,ist)));
	         const std::complex<double> adiag_imag(real(tdmlwft_->adiag(idir*2+1,ist)));
                 adiag_inv_real[ist] = real( std::complex<double>(1,0) / adiag );
                 adiag_inv_imag[ist] = imag( std::complex<double>(1,0) / adiag );
              }

              ComplexMatrix& ccos(sdcos[idir]->c());
              ComplexMatrix& csin(sdsin[idir]->c());
              ComplexMatrix& cp(dwf_->sd(0,0)->c());
	    
              int nloc = cp.nloc();
              int mloc = cp.mloc();
              int ione = 1;

              const complex<double> fac = length(cell.a(idir))
                               * e_field_[idir] / ( 2.0 * M_PI );

              for (int in = 0; in < nloc; in++)
              {
                int ist = cp.jglobal(in);
                std::complex<double> fac1 = adiag_inv_real[ist] * fac;
                std::complex<double> fac2 = adiag_inv_imag[ist] * fac;

                std::complex<double> *ptr1 = &cp[in*mloc],
                       *ptrcos = &ccos[in*mloc],
                       *ptrsin = &csin[in*mloc];

	        for (int ii=0; ii<mloc; ii++)
		  {
		   ptr1[ii] += (fac2 * ptrcos[ii]);
		   ptr1[ii] += (fac1 * ptrsin[ii]);
		  }
               }
          } 
          else if ( pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q )
          {
            // MLWF_REF part: real-space correction
            DoubleMatrix cc(rwf_[idir]->sd(0,0)->c());
            DoubleMatrix cp(dwf_->sd(0,0)->c());

            int size = cc.size();
            double alpha = e_field_[idir];
            int ione = 1;
            //cp.axpy(alpha,cc);
            daxpy (&size, &alpha, cc.valptr(), &ione, cp.valptr(), &ione);
          } 
	  else if (pol_type_ == tdmlwf_ref) //CS
	  { 
	  // TDMLWF refinement for complex wavefxn
	    ComplexMatrix& cc(rwf_[idir]->sd(0,0)->c());
            ComplexMatrix& cp(dwf_->sd(0,0)->c());

            int size = cc.size();
            complex<double> alpha = complex<double>(e_field_[idir],0);
            int ione = 1;
            //cp.axpy(alpha,cc);
            zaxpy (&size, &alpha, cc.valptr(), &ione, cp.valptr(), &ione);
          } // if pol_type_
        } // if e_field_[idir]
      } // for idir
    } // if finite_field_
  }
  dipole_total_ = dipole_ion_ + dipole_el_;
  cell.fold_in_ws(dipole_ion_);
  cell.fold_in_ws(dipole_el_);
  cell.fold_in_ws(dipole_total_);
}

///////////////////////////////////////////////////////////////////////////////
double ElectricEnthalpy::enthalpy(Wavefunction& dwf, bool compute_hpsi)
{
  // return zero if polarization is off, or field is zero
  if ( pol_type_ == off || !finite_field_ )
    return 0.0;

  enthalpy_ = - e_field_ * dipole_total_;
  //cout<<"dipole:\t"<<dipole_el_<<endl;
  if ( compute_hpsi )
  {
    // assert gamma-only and no spin
    assert(dwf.nkp()==1 && dwf.nspin()==1);
    dwf.sd(0,0)->c() += dwf_->sd(0,0)->c();
  }
  return enthalpy_;
}

///////////////////////////////////////////////////////////////////////////////
// Correction scheme by M. Stengel and N. Spaldin,
// Phys. Rev. B 73, 075121 (2006)
// Calculate correction in real space and derivatives of the correction
///////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::compute_correction(void)  
{
  int np0v = basis_.np(0);
  int np1v = basis_.np(1);
  int np2v = basis_.np(2);
  const ComplexMatrix& c = sd_.c();
  //DoubleMatrix cp(c);

  FourierTransform ft(basis_,np0v,np1v,np2v);

  int np012v = ft.np012();
  int np012loc = ft.np012loc();
  int nst = sd_.nst();
  int nloc = c.nloc();
  int mloc = c.mloc();

  // store (x-x0)|psi> in rwfs
  rwf_[0]->clear();
  rwf_[1]->clear();
  rwf_[2]->clear();

  ComplexMatrix& cx = rwf_[0]->sd(0,0)->c();
  ComplexMatrix& cy = rwf_[1]->sd(0,0)->c();
  ComplexMatrix& cz = rwf_[2]->sd(0,0)->c();

  //DoubleMatrix cpx(cx);
  //DoubleMatrix cpy(cy);
  //DoubleMatrix cpz(cz);

  // calculate refinements
  // ref is scaled by np012v
  vector<double> ref(nst*3);
  //if ( compute_quadrupole_ ) ref.resize(nst*9);

  // cell size;
  const UnitCell& cell = sd_.basis().cell();
  const double ax = cell.amat(0);
  const double ay = cell.amat(4);
  const double az = cell.amat(8);

  // half cell size;
  const double ax2 = ax / 2.0;
  const double ay2 = ay / 2.0;
  const double az2 = az / 2.0;

  // grid size;
  const double dx = ax / np0v;
  const double dy = ay / np1v;
  const double dz = az / np2v;

  for ( int i = 0; i < nst; i++ )
    correction_[i] = D3vector(0,0,0);

  int niter = 1;
  // check if override from the debug variable
  // use: set debug MLWF_REF_NITER <value>
  if ( s_.ctrl.debug.find("MLWF_REF_NITER") != string::npos )
  {
    istringstream is(s_.ctrl.debug);
    string s;
    is >> s >> niter;
    assert(s=="MLWF_REF_NITER");
    if ( onpe0_ )
      cout << " ElectricEnthalpy: override niter value: niter= "
           << niter << endl;
    assert(niter > 0);
  }
  for ( int iter = 0; iter < niter; iter++ )
  {
    fill(ref.begin(),ref.end(),0.0);

    // wavefunctions in real space
    vector<complex<double>>  wftmp(np012loc);
    vector<complex<double>>  xwftmp(np012loc);
    vector<complex<double>>  ywftmp(np012loc);
    vector<complex<double>>  zwftmp(np012loc);

    for ( int in = 0; in < nloc; in++ )
    {
      int n = c.jglobal(in);
      double* pref;
      pref = &ref[3*n];

      // real space wavefunction in wftmp
      tmap["ft"].start();
      ft.backward(c.cvalptr(mloc*in),&wftmp[0]);
      tmap["ft"].stop();

      tmap["real"].start();
      double x0 = mlwfc_[n][0] + correction_[n][0];
      double y0 = mlwfc_[n][1] + correction_[n][1];
      double z0 = mlwfc_[n][2] + correction_[n][2];

      // compute shifted sawtooth potentials vx, vy, vz
      vector<double> vx(ft.np0()), vy(ft.np1()), vz(ft.np2());
      for ( int i = 0; i < vx.size(); i++ )
      {
        double x = i * dx - x0;
        if ( x >  ax2 ) x -= ax;
        if ( x < -ax2 ) x += ax;
#ifdef NO_VSST
        vx[i] = x;
#else
        vx[i] = ax * vsst(x/ax);
#endif
      }
      for ( int j = 0; j < vy.size(); j++ )
      {
        double y = j * dy - y0;
        if ( y >  ay2 ) y -= ay;
        if ( y < -ay2 ) y += ay;
#ifdef NO_VSST
        vy[j] = y;
#else
        vy[j] = ay * vsst(y/ay);
#endif
      }
      for ( int k = 0; k < vz.size(); k++ )
      {
        double z = k * dz - z0;
        if ( z >  az2 ) z -= az;
        if ( z < -az2 ) z += az;
#ifdef NO_VSST
        vz[k] = z;
#else
        vz[k] = az * vsst(z/az);
#endif
      }

      for ( int i = 0; i < np012loc; i++ )
      {
        int ix = ft.i(i);
        int iy = ft.j(i);
        int iz = ft.k(i);
        const complex<double> wft = wftmp[i];
        const complex<double> xwft = wft * vx[ix];
        const complex<double> ywft = wft * vy[iy];
        const complex<double> zwft = wft * vz[iz];

        pref[0] += real(conj(wft) * wft)*vx[ix];
        pref[1] += real(conj(wft) * wft)*vy[iy];
        pref[2] += real(conj(wft) * wft)*vz[iz];

        xwftmp[i] = xwft;
        ywftmp[i] = ywft;
        zwftmp[i] = zwft;

      } // for i
      tmap["real"].stop();

      // ft to get xwf in reciprocal space at the last iteration
      if ( iter == niter - 1 )
      {
        tmap["ft"].start();
        if ( e_field_[0] != 0.0 )
          ft.forward(&xwftmp[0],cx.valptr(mloc*in));
        if ( e_field_[1] != 0.0 )
          ft.forward(&ywftmp[0],cy.valptr(mloc*in));
        if ( e_field_[2] != 0.0 )
          ft.forward(&zwftmp[0],cz.valptr(mloc*in));
        tmap["ft"].stop();
      } // if
    } //for in

    ctxt_.barrier();
    tmap["dsum"].start();
      ctxt_.dsum(3*nst,1,&ref[0],3*nst);
    tmap["dsum"].stop();

    tmap["real"].start();
    {  
      for ( int ist = 0; ist < nst; ist++ )
      {
        D3vector& pcor = correction_[ist];
        pcor[0] += ref[ist*3]/np012v;
        pcor[1] += ref[ist*3+1]/np012v;
        pcor[2] += ref[ist*3+2]/np012v;
      }
    }
    tmap["real"].stop();
  } // for iter
}

////////////////////////////////////////////////////////////////////////////////
void ElectricEnthalpy::print(ostream& os) const
{
  if ( pol_type_ == off ) return;
  os << fixed << right << setprecision(8);
  // print MLWF centers if pol_type_ == MLWF or MLWF_REF or MLWF_REF_Q or TDMLWF or TDMLWF_ref//CCS
  if ( pol_type_ == mlwf || pol_type_ == mlwf_ref || pol_type_ == mlwf_ref_q || pol_type_ == tdmlwf || pol_type_ == tdmlwf_ref)
  {
    int nst = sd_.nst();
    os << " <mlwf_set size=\"" << nst << "\">" << endl;
    for ( int i = 0; i < nst; i++ )
    {
      os << " <mlwf center=\"" << setprecision(8)
         << setw(12) << mlwfc_[i].x << " "
         << setw(12) << mlwfc_[i].y << " "
         << setw(12) << mlwfc_[i].z << " \""
         << "       spread=\" " << mlwfs_[i] << " \"/>" << endl;
      if ( pol_type_ == mlwf_ref || pol_type_ == tdmlwf_ref)
      {
        os << " <mlwf_ref center=\"" << setprecision(8)
           << setw(12) << mlwfc_[i].x + correction_[i].x << " "
           << setw(12) << mlwfc_[i].y + correction_[i].y << " "
           << setw(12) << mlwfc_[i].z + correction_[i].z << " \"";
        os << "/>" << endl;
      }
    }
    os << " </mlwf_set>" << endl;
  }

  // print dipole
  os << setprecision(10) << fixed << right;
  os << " <dipole>\n";
  os << "   <dipole_ion>   "
     << setw(14) << dipole_ion_.x << " "
     << setw(14) << dipole_ion_.y << " "
     << setw(14) << dipole_ion_.z << " </dipole_ion>\n";
  os << "   <dipole_el>    "
     << setw(14) << dipole_el_.x << " "
     << setw(14) << dipole_el_.y << " "
     << setw(14) << dipole_el_.z << " </dipole_el>\n";
  os << "   <dipole_total> "
     << setw(14) << dipole_total_.x << " "
     << setw(14) << dipole_total_.y << " "
     << setw(14) << dipole_total_.z << " </dipole_total>\n";
  os << " </dipole>\n";
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<< ( ostream& os, const ElectricEnthalpy& e )
{
  e.print(os);
  return os;
}

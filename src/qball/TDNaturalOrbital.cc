////////////////////////////////////////////////////////////////////////////////
//// Copyright (c) 2013, Lawrence Livermore National Security, LLC.
//// qb@ll:  Qbox at Lawrence Livermore
////
//// This file is part of qb@ll.
////
//// Produced at the Lawrence Livermore National Laboratory.
//// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
//// Based on the Qbox code by Francois Gygi Copyright (c) 2008
//// LLNL-CODE-635376. All rights reserved.
////
//// qb@ll is free software: you can redistribute it and/or modify
//// it under the terms of the GNU General Public License as published by
//// the Free Software Foundation, either version 3 of the License, or
//// (at your option) any later version.
////
//// This program is distributed in the hope that it will be useful,
//// but WITHOUT ANY WARRANTY; without even the implied warranty of
//// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//// GNU General Public License for more details, in the file COPYING in the
//// root directory of this distribution or <http://www.gnu.org/licenses/>.
////
//////////////////////////////////////////////////////////////////////////////////
////
//// TDNaturalOrbital.cc
////
//////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
#include <math.h>
#include "Sample.h"
#include "TDNaturalOrbital.h"
#include "SlaterDet.h"
#include "Wavefunction.h"
#include "Context.h"
#include <math/matrix.h>
#include <math/blas.h>
#include <fstream>
#include "FourierTransform.h"
using namespace std;

TDNaturalOrbital::TDNaturalOrbital(const Sample& s):
s_(s),ref_((*s.proj_wf).sd(0,0)->c()),basis( (*s.proj_wf).sd(0,0)->basis()),ctxt_(ref_.context()),update_((*s.proj_wf).sd(0,0)->c()),
nto_coeff((*s.proj_wf).sd(0,0)->c()),elec_coeff((*s.proj_wf).sd(0,0)->c()),hole_coeff((*s.proj_wf).sd(0,0)->c())
{ 
  m = ref_.m();
  n = ref_.n();
  mb= ref_.mb();
  nb= ref_.nb();
  nto_.resize(n);
  //valarray<double> nto_(n);
  //ComplexMatrix  nto_coeff(ctxt_,m,n,mb,nb);
  //ComplexMatrix  elec_coeff(ctxt_,m,n,mb,nb);
  //ComplexMatrix  hole_coeff(ctxt_,m,n,mb,nb);
  
  int np0 = basis.np(0);
  int np1 = basis.np(1);
  int np2 = basis.np(2); 
  ft = new FourierTransform(basis,np0,np1,np2);
}
TDNaturalOrbital::~TDNaturalOrbital()
{
  //delete nto_coeff;
  //delete elec_coeff;
  //delete hole_coeff;
  delete ft;
}
void TDNaturalOrbital::update(const ComplexMatrix& mat)
{
  update_=mat;
}
void TDNaturalOrbital::update_NTO()
{
  ComplexMatrix z(ctxt_,n,n,nb,nb);
  ComplexMatrix ortho(z),tmp(z);
  
  tmp.gemm('c','n',1.0,ref_,update_,0.0);
  ortho.gemm('c','n',1.0,tmp,tmp,0.0);
  valarray<double> w(n);
  ortho.heev('l',w,z);
  if (  ctxt_.oncoutpe() )
  {
     for (int ii =0; ii<n; ii++)
             nto_[ii]=w[ii];
  }
  //(*nto_coeff).print(cout);
  nto_coeff.gemm('n','n',1.0,update_,z,0.0); 
}
void TDNaturalOrbital::update_elec()
{

  ComplexMatrix z(ctxt_,n,n,nb,nb);
  z.gemm('c','n',1.0,ref_,nto_coeff,0.0);
  ComplexMatrix u(z);  
  int nloc = z.nloc();
  int mloc = z.mloc();
  

  for (int in = 0; in < nloc; in++)
  {
      int ist = z.jglobal(in);
      double occ = pow(nto(ist),-0.5);
      std::complex<double> *ptrz = &z[in*mloc],*ptru = &u[in*mloc];
      for (int ii =0; ii<mloc; ii++)
      {
          ptru[ii] = occ* ptrz[ii];
      }
  }
  elec_coeff.gemm('n','n',1.0,ref_,u,0.0);
}
void TDNaturalOrbital::update_hole()
{
   for (int in = 0; in < nb; in++)
   {
      int ist = nto_coeff.jglobal(in);
      double occ = pow(nto(ist),0.5);
      double hole = pow(1-nto(ist),0.5);
      std::complex<double> *ptrnto = nto_coeff.valptr(in*mb),*ptrelec = elec_coeff.valptr(in*mb),*ptrhole = hole_coeff.valptr(in*mb);
      for (int ii =0; ii<mb; ii++)
      {
            ptrhole[ii]= (ptrnto[ii]-occ*ptrelec[ii])/hole; 
      }
   } 
}
void TDNaturalOrbital::print_hole_orbital(int m,string filename)
{
   vector<complex<double>> wftmp(ft->np012());
   vector<double> wftmpr(ft->np012());
   int nloc = hole_coeff.y(m);
   ft->backward(hole_coeff.valptr(mb*nloc),&wftmp[0]);
   double *a = (double*) &wftmp[0];
   for ( int i = 0; i < ft->np012loc(); i++ ) 
          wftmpr[i] = sqrt(a[2*i]*a[2*i] + a[2*i+1]*a[2*i+1]);
   print_orbital(&wftmpr[0],filename);

}

void TDNaturalOrbital::print_elec_orbital(int m,string filename)
{
   vector<complex<double>> wftmp(ft->np012());
   vector<double> wftmpr(ft->np012());
   int nloc = elec_coeff.y(m);
   ft->backward(elec_coeff.cvalptr(mb*nloc),&wftmp[0]);
   double *a = (double*) &wftmp[0];
   for ( int i = 0; i < ft->np012loc(); i++ )
          wftmpr[i] = sqrt(a[2*i]*a[2*i] + a[2*i+1]*a[2*i+1]);
   print_orbital(&wftmpr[0],filename);

}

void TDNaturalOrbital::print_nto_orbital(int m,string filename)
{
   vector<complex<double>> wftmp(ft->np012());
   vector<double> wftmpr(ft->np012());
   int nloc = nto_coeff.y(m);
   ft->backward(nto_coeff.cvalptr(mb*nloc),&wftmp[0]);
   double *a = (double*) &wftmp[0];
   for ( int i = 0; i < ft->np012loc(); i++ )
          wftmpr[i] = sqrt(a[2*i]*a[2*i] + a[2*i+1]*a[2*i+1]);
   print_orbital(&wftmpr[0],filename);

}

void TDNaturalOrbital::print_orbital(double * wftmp,string filename)
{ 
  
  for ( int i = 0; i < ctxt_.nprow(); i++ )
  {
          bool iamsending = ref_.pc(n) == ctxt_.mycol() && i == ctxt_.myrow();
          int size=-1;
          if (  ctxt_.oncoutpe() )
          {
               if ( iamsending ){}
               else ctxt_.irecv(1,1,&size,1,i,ref_.pc(n));
           }
           else
           {
               if ( iamsending )
               {
                  size = ft->np012loc();
                  ctxt_.isend(1,1,&size,1,0,0);
               }
            }
            if ( ctxt_.oncoutpe() )
            {
               if ( iamsending ) {}
               else
               {
                   int istart = ft->np0() * ft->np1() * ft->np2_first(i);
                   ctxt_.drecv(size,1,&wftmp[istart],1,i,ref_.pc(n));
               }
            }
            else
            {
                 if ( iamsending )
                 {
                     ctxt_.dsend(size,1,&wftmp[0],1,0,0);
                 }
            }
            if ( ctxt_.oncoutpe() )
            {
                vector<double> tmpr; 
                tmpr.resize(ft->np012());
                for ( int i = 0; i < ft->np012(); i++ )
                {
                    tmpr[i] = wftmp[i];
                }
                ofstream os;
                ostringstream oss;
                oss.width(7);  
                oss.fill('0');  
                oss << endl;
                os.open(filename);
                os <<endl;
                os << endl;
                int natoms = s_.atoms.size();
                D3vector a0 = s_.atoms.cell().a(0);
                D3vector a1 = s_.atoms.cell().a(1);
                D3vector a2 = s_.atoms.cell().a(2);
                os << natoms << " " << -0.5*(a0+a1+a2) << endl;
                os << np0 << " " << a0/np0 << endl;
                os << np1 << " " << a1/np1 << endl;
                os << np2 << " " << a2/np2 << endl;
                const int nsp = s_.atoms.nsp();
                for ( int is = 0; is < nsp; is++ )
                {
                    Species* sp = s_.atoms.species_list[is];
                    const int z = sp->atomic_number();
                    const int na = s_.atoms.na(is);
                    for ( int ia = 0; ia < na; ia++ )
                    {
                         Atom *ap = s_.atoms.atom_list[is][ia];
                         os << setprecision(5);
                         os << z << " " << ((double) z) << " " << ap->position() << endl;

                     }
                        
                }
                os.setf(ios::scientific,ios::floatfield);
                os << setprecision(5);
                for ( int i = 0; i < np0; i++ )
                {
                     const int ip = (i + np0/2 ) % np0;
                     for ( int j = 0; j < np1; j++ )
                     {
                          const int jp = (j + np1/2 ) % np1;
                          for ( int k = 0; k < np2; k++ )
                          {
                                const int kp = (k + np2/2 ) % np2;
                                os << setw(13) << tmpr[ip+np0*(jp+np1*kp)];
                                if ( ( k % 6 ) == 5 ) os<<endl;
                          }
                          if ( ( np2 % 6 ) != 0 ) os<<endl;
                      }
                 }
                 os.close();
            }
     }
}  

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
#include "string.h"
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

TDNaturalOrbital::TDNaturalOrbital(const Sample& s):
s_(s),ref_((*s.proj_wf).sd(0,0)->c()),basis( (*s.proj_wf).sd(0,0)->basis()),ctxt_(ref_.context()),update_((*s.proj_wf).sd(0,0)->c()),
nto_coeff((*s.proj_wf).sd(0,0)->c()),hole_coeff((*s.proj_wf).sd(0,0)->c()),intermidiate(ctxt_,ref_.n(),ref_.n(),ref_.nb(),ref_.nb())
{ 
  m = ref_.m();
  n = ref_.n();
  mb= ref_.mb();
  nb= ref_.nb();
  nto_.resize(n);
  electron = new Wavefunction(*(s_.hamil_wf));
  
  //valarray<double> nto_(n);
  //ComplexMatrix  nto_coeff(ctxt_,m,n,mb,nb);
  //ComplexMatrix  elec_coeff(ctxt_,m,n,mb,nb);
  //ComplexMatrix  hole_coeff(ctxt_,m,n,mb,nb);
  electron ->clear();
  //hole_coeff = hole->sd(0,0)->c();
  np0 = basis.np(0);
  np1 = basis.np(1);
  np2 = basis.np(2); 
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
  intermidiate.gemm('c','n',1.0,ref_,update_,0.0);
}
void TDNaturalOrbital::update_NTO()
{
  ComplexMatrix z(intermidiate);
  ComplexMatrix ortho(intermidiate);
  
  //tmp.gemm('c','n',1.0,ref_,update_,0.0);
  ortho.gemm('c','n',1.0,intermidiate,intermidiate,0.0);
  ortho.heev('l',nto_,z);
  //(*nto_coeff).print(cout);
  nto_coeff.gemm('n','n',1.0,update_,z,0.0); 
}
void TDNaturalOrbital::update_hole()
{ 
  ComplexMatrix ortho(intermidiate);
  ortho.gemm('c','n',1.0,ref_,nto_coeff,0.0);
  hole_coeff.gemm('n','n',1.0,ref_,ortho,0.0);
}
void TDNaturalOrbital::update_effective_elec()
{
   // in current implementary we require projempty = no of occupied oribitals
   ComplexMatrix & elec_coeff = electron->sd(0,0)->c();
   int nvir = (s_.proj_wf_virtual)->sd(0,0)->c().n()-n;
   ComplexMatrix vir_coeff(ctxt_,m,nvir,mb,nb);
   vir_coeff.getsub(s_.proj_wf_virtual->sd(0,0)->c(),m,nvir,0,n);
   ComplexMatrix tmp(ctxt_,n,nvir,nb,nb);
   tmp.gemm('c','n',1.0,nto_coeff,vir_coeff,0.0);
   ComplexMatrix ortho(ctxt_,nvir,nvir,nb,nb);  
   ComplexMatrix phase(ortho);

   ortho.gemm('c','n',-1.0,tmp,tmp,0.0);
   tmp.clear();
   valarray<double> eig_;
   
   eig_.resize(nvir);
   ortho.heev('l',eig_,phase); 
   
   //ComplexMatrix tmp_result(ctxt_,m,nvir,mb,nb);
   elec_coeff.gemm('n','n',1.0,vir_coeff,phase,0.0);
   //vir_coeff.print(cout);
   //

}
void TDNaturalOrbital::update_elec()
{
   ComplexMatrix & elec_coeff = electron->sd(0,0)->c();
   int nloc = elec_coeff.nloc();
   int mloc = elec_coeff.mloc();
   for (int in = 0; in < nloc; in++)
   {
      int ist = nto_coeff.jglobal(in);
      double occ = pow(nto(ist),0.5);
      double elec = pow(abs(1-nto(ist)),0.5);
      std::complex<double> *ptrnto = nto_coeff.valptr(in*mb),*ptrhole = hole_coeff.valptr(in*mb),*ptrelec = elec_coeff.valptr(in*mb);
      for (int ii =0; ii<mloc; ii++)
      {
            //ptrhole[ii]= (ptrnto[ii]-occ*ptrelec[ii]);
            ptrelec[ii]= (ptrnto[ii]-ptrhole[ii])/elec;  
      }
   }
   //ComplexMatrix test(ctxt_,n,n,nb,nb);
   //test.gemm('c','n',1.0,hole_coeff,hole_coeff,0.0);
   //test.print(cout); 
}
void TDNaturalOrbital::print_elec_orbital(int m,string filename)
{
   ComplexMatrix & elec_coeff = electron->sd(0,0)->c();
   vector<complex<double>> wftmp(ft->np012());
   vector<double> wftmpr(ft->np012());
   int nloc = elec_coeff.y(m);
   ft->backward(elec_coeff.valptr(mb*nloc),&wftmp[0]);
   double *a = (double*) &wftmp[0];
   for ( int i = 0; i < ft->np012loc(); i++ ) 
          wftmpr[i] = sqrt(a[2*i]*a[2*i] + a[2*i+1]*a[2*i+1]);
   print_orbital(&wftmpr[0],filename);

}

void TDNaturalOrbital::print_hole_orbital(int m,string filename)
{
   vector<complex<double>> wftmp(ft->np012());
   vector<double> wftmpr(ft->np012());
   int nloc = hole_coeff.y(m);
   ft->backward(hole_coeff.cvalptr(mb*nloc),&wftmp[0]);
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
void TDNaturalOrbital::proj_sorted_hole_orbital()
{

   ComplexMatrix occ_coeff(ctxt_,m,n,mb,nb);
   occ_coeff.getsub(s_.proj_wf_virtual->sd(0,0)->c(),m,n,0,0);
   ComplexMatrix proj_hole(ctxt_,n,n,nb,nb);
   proj_hole.clear();
   proj_hole.gemm('c','n',1.0,occ_coeff,hole_coeff,0.0);
   int hole_index1 = s_.ctrl.holeindex1;
   int hole_index2 = s_.ctrl.holeindex2;
   if ( hole_index1 ==-1 &&  hole_index2 == -1)
   {
      hole_index1 = 0;
      hole_index2 = n;
   }
   Context ctxtl(1,1);

   ComplexMatrix t(ctxtl,n,n);
   t.getsub(proj_hole,n,n,0,0);
   if ( ctxt_.oncoutpe() )
   {
         //if (hole_index1== hole_index2)
         //{
         //    cout << " <hole_projection index=\"" << hole_index1 << "\">" << endl;
         //    for (int i =0 ; i <n ; i++)
         //    {
         //       cout.setf(ios::fixed, ios::floatfield);
         //       cout.setf(ios::right, ios::adjustfield);
         //       double pop = norm(t[hole_index1*n+i]) /nto_[hole_index1];
         //       cout << "   <population=\"" << setprecision(6) << setw(12) << pop << " \"/>"<<endl;
         //   }
         //   cout << " </hole_projection>" << endl;
        // }
         //else
         {        
             int size = hole_index2-hole_index1;
             int index[size];
             for (int i =0 ; i <size ; i++)
             {
                 double local_max = 0 ;
                 for  (int j =0 ; j <n ; j++)
                 {
                   int reverse = n-i-1-hole_index1;
                   double pop = norm(t[reverse*n+j])/nto_[i];
                   if ( pop>local_max)
                   {
                        index[i]=j;
                        local_max = pop;
                   }
                 }
             }
            cout << " <sorted_hole size=\"" << size << "\">" << endl;
            for ( int i = 0; i < size; i++ )
            {
                cout.setf(ios::fixed, ios::floatfield);
                cout.setf(ios::right, ios::adjustfield);
                cout << "   <population=\"" << setprecision(10)
                     << setw(12) << nto(index[i]) << " \"/>"<<endl;
            }
            cout << " </sorted_hole>" << endl;
            cout << " <hole_projection index=\"";
             for (int i =0 ; i <size ; i++)
             {
                  cout<<i<<"\t";
             }
             cout<<endl;
             for (int j =0 ; j <n ; j++)
             {

                cout.setf(ios::fixed, ios::floatfield);
                cout.setf(ios::right, ios::adjustfield);
                cout << "   <population=\"";
                for (int i = hole_index1 ; i < hole_index2 ; i++)
                {
                   double pop = norm(t[index[i]*n+j])/nto_[i] ;
                   cout<<setprecision(6) << setw(12) << pop ;
                }
                cout << " \"/>"<<endl;
            }


         }

   }
}
void TDNaturalOrbital::proj_hole_orbital()
{
  
   ComplexMatrix occ_coeff(ctxt_,m,n,mb,nb);
   occ_coeff.getsub(s_.proj_wf_virtual->sd(0,0)->c(),m,n,0,0);
   ComplexMatrix proj_hole(ctxt_,n,n,nb,nb);
   proj_hole.clear();
   proj_hole.gemm('c','n',1.0,occ_coeff,hole_coeff,0.0);
   int hole_index1 = s_.ctrl.holeindex1;
   int hole_index2 = s_.ctrl.holeindex2;
   if ( hole_index1 ==-1 &&  hole_index2 == -1)
   {
      hole_index1 = 0;
      hole_index2 = n;
   }
   Context ctxtl(1,1);

   ComplexMatrix t(ctxtl,n,n);
   t.getsub(proj_hole,n,n,0,0);
   if ( ctxt_.oncoutpe() )
   {
         if (hole_index1== hole_index2)
         {
             cout << " <hole_projection index=\"" << hole_index1 << "\">" << endl;
             for (int i =0 ; i <n ; i++)
             {
                cout.setf(ios::fixed, ios::floatfield);
                cout.setf(ios::right, ios::adjustfield);
                double pop = norm(t[hole_index1*n+i]) /nto_[hole_index1];
                cout << "   <population=\"" << setprecision(6) << setw(12) << pop << " \"/>"<<endl;
            }
            cout << " </hole_projection>" << endl;
         }
         else
         {
             int size = hole_index2-hole_index1;
             cout << " <hole_projection index=\"";
             for (int i =hole_index1 ; i <hole_index2 ; i++)
             {
                  cout<<i<<"\t";
             }
             cout<<endl;
             for (int j =0 ; j <n ; j++)
             {

                cout.setf(ios::fixed, ios::floatfield);
                cout.setf(ios::right, ios::adjustfield);
                cout << "   <population=\"";
                for (int i = hole_index1 ; i < hole_index2 ; i++)
                {
                   double pop = norm(t[i*n+j])/nto_[i] ;
                   cout<<setprecision(6) << setw(12) << pop ;
                }
                cout << " \"/>"<<endl;
            }


         }

   }
}
void TDNaturalOrbital::proj_elec_orbital()
{
   ComplexMatrix & elec_coeff = electron->sd(0,0)->c();
   int nvir = (s_.proj_wf_virtual)->sd(0,0)->c().n()-n;
   ComplexMatrix vir_coeff(ctxt_,m,nvir,mb,nb);
   vir_coeff.getsub(s_.proj_wf_virtual->sd(0,0)->c(),m,nvir,0,n);
   //vir_coeff.print(cout);
   //s_.proj_wf_virtual->sd(0,0)->c().print(cout);
   //s_.proj_wf_virtual->sd(0,0)->c().getsub( vir_coeff, 0,n,m,nvir);
  //vir_coeff.print(cout);
   ComplexMatrix proj_elec(ctxt_,nvir,n,nb,nb);
   proj_elec.clear();
   proj_elec.gemm('c','n',1.0,vir_coeff,elec_coeff,0.0);
   int elec_index1 = s_.ctrl.elecindex1;
   int elec_index2 = s_.ctrl.elecindex2;
   if ( elec_index1 ==-1 &&  elec_index2 == -1)
   {
      elec_index1 = 0;
      elec_index2 = nvir;
   }
   // To be done: add in memory buffer
   Context ctxtl(1,1);
   
   ComplexMatrix t(ctxtl,nvir,n);
   t.getsub(proj_elec,nvir,n,0,0);
   if ( ctxt_.oncoutpe() )
   {
         if (elec_index1== elec_index2)
         {
             cout << " <electron_projection index=\"" << elec_index1 << "\">" << endl;
             for (int i =0 ; i <nvir ; i++)
             {
                cout.setf(ios::fixed, ios::floatfield);
                cout.setf(ios::right, ios::adjustfield);
                double pop = norm(t[elec_index1*nvir+i]) ;
                cout << "   <population=\"" << setprecision(6) << setw(12) << pop << " \"/>"<<endl;
            }
            cout << " </elec_projection>" << endl;
         }
         else
         {
             int size = elec_index2-elec_index1;
             cout << " <elec_projection index=\"";
             for (int i =elec_index1 ; i <elec_index2 ; i++)
             {     
                  cout<<i<<"\t";
             }
             cout<<endl;
             for (int j =0 ; j <nvir ; j++)
             {
                
                cout.setf(ios::fixed, ios::floatfield);
                cout.setf(ios::right, ios::adjustfield);
                //double pop = pow(norm(t[hole_index1*nvir+i]),2) ;
                cout << "   <population=\"";
                for (int i =elec_index1 ; i <elec_index2 ; i++)
                {  
                   double pop = norm(t[i*nvir+j]) ;  
                   cout<<setprecision(6) << setw(12) << pop ;
                }
                cout << " \"/>"<<endl;
            }

             
         }
   } 
  // if ( t.active() )
   
}
/*
void TDNaturalOrbital::save_hole_orbital()
{
      string format = "binary";
      string filestr="hole/hole";
      if ( ctxt_.oncoutpe() ) {
      cout << "<!-- SaveCmd:  writing wf " << filestr << "... -->" << endl;
      string dirstr = filestr.substr(0, filestr.find_last_of('/'));
      //string distr = "./states"
      int mode = 0775;
      struct stat statbuf;
      int rc = stat(dirstr.c_str(), &statbuf);
      if (rc == -1) {
	cout << "<!-- Creating directory: " << dirstr << "/ -->" << endl;
	rc = mkdir(dirstr.c_str(), mode);
	rc = stat(dirstr.c_str(), &statbuf);
      }
      if (rc != 0 || !(statbuf.st_mode)) {
	cout << "<ERROR> Can't stat directory " << dirstr << " </ERROR> " << endl;
	MPI_Abort(MPI_COMM_WORLD,2);
	}
      
      }
      MPI_Barrier(MPI_COMM_WORLD);
      
      hole->write_states(filestr,format);
      hole->write_mditer(filestr,s_.ctrl.mditer);
}
*/
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
                //ostringstream oss;
                //oss.width(7);  
                //oss.fill('0');  
                //oss << endl;
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

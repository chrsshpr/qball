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
// jade_complex.C
//
////////////////////////////////////////////////////////////////////////////////
#include <config.h>
#include <cmath>
#include <cassert>
#include <vector>
#include <deque>
#include <algorithm>
#include <limits> // epsilon
#include <iostream>
#include <iomanip>
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef SCALAPACK
#include "blacs.h"
#endif

#include <complex>
#include "Context.h"
#include <math/matrix.h>
#include "jacobi.h"
#include "jacobi_eigenvalue.h"
#include <math/blas.h>
#include "blacs.h" 
#include "Timer.h"
using namespace std;

int jade_complex(int maxsweep, double tol, vector<ComplexMatrix*> a, ComplexMatrix& atmp1, ComplexMatrix& atmp2, ComplexMatrix& atmp3, ComplexMatrix& atmp4, ComplexMatrix& atmp5, ComplexMatrix& atmp6, ComplexMatrix& u, ComplexMatrix& tmpmat, vector<vector<complex<double> > >& adiag)
{
  Timer tm_comm;
  const bool debug_diag_sum = false;

  const double eps = numeric_limits<double>::epsilon();
  assert(tol>eps);
  const Context& ctxt = u.context();
  // The input matrices are *a[k]
  // the orthogonal transformation is returned in u
  // on exit, the matrices a[k] are maximally diagonal,
  // u contains the orthogonal transformation
  // adiag[k][i] contains the diagonal elements of the a[k]'s

  for ( int k = 0; k < a.size(); k++ )
  {
    assert(a[k]->context() == u.context());
    assert(a[k]->m()==a[k]->n());
    assert(a[k]->m()==u.n());
    assert(a[k]->m()==u.m());
    assert(a[k]->mb()==u.mb());
    assert(a[k]->nb()==u.nb());
  }

  int mloc = a[0]->mloc();
  int nloc = a[0]->nloc();

  // identify the last active process column
  // process columns beyond that column do not have any elements of a[k]
  // compute num_nblocks = total number of column blocks
  // if num_nblocks >= ctxt.npcol(), all process columns are active
  // otherwise, the last active process column has index num_nblocks-1
  const int num_nblocks = u.n() / u.nb() + ( u.n()%u.nb() == 0 ? 0 : 1 );
  const int last_active_process_col = min(ctxt.npcol()-1, num_nblocks-1);
  // initialize u with the identity
  u.identity();
  // eigenvalue array

  adiag.resize(a.size());
  for ( int k = 0; k < a.size(); k++ )
  {
    adiag[k].resize(a[k]->n());
  }

  // check if the local number of rows is odd
  const bool nloc_odd = ( a[0]->nloc()%2 != 0 );

  // if nloc is odd, auxiliary arrays are created to host an extra column
  // for both a[k] and u
  vector<vector<complex<double> > > a_aux(a.size());
  vector<complex<double> > u_aux;
  if ( nloc_odd )
  {
    for ( int k = 0; k < a.size(); k++ )
      a_aux[k].resize(mloc);
    u_aux.resize(mloc);
  }

  // compute local number of pairs nploc
  const int nploc = (a[0]->nloc()+1)/2;
  // dimension of top and bot arrays is nploc: local number of pairs
  deque<int> top(nploc), bot(nploc);

  // compute total number of pairs np
  int np = nploc;
  ctxt.isum('r',1,1,&np,1);
  // initialize top and bot arrays
  // the pair i is (top[i],bot[i])
  // top[i] is the local index of the top column of pair i
  // bot[i] is the local index of the bottom column of pair i
  for ( int i = 0; i < nploc; i++ )
    top[i] = i;
  for ( int i = 0; i < nploc; i++ )
    bot[nploc-i-1] = nploc+i;
  // if top[i] or bot[i] == nloc,the data resides in the array a_aux or u_aux

  // jglobal: global column index
  // jglobal[i] is the global column index of the column residing in
  // the local vector i. If nloc_odd and i==2*nploc-1, jglobal[i] == -1
  vector<int> jglobal(2*nploc,-1);
  for ( int jblock = 0; jblock < a[0]->nblocks(); jblock++ )
    for ( int y = 0; y < a[0]->nbs(jblock); y++ )
     {
      jglobal[y + jblock*a[0]->nb()] = a[0]->j(jblock,y);
     }
  // store addresses of columns of a and of u in acol and ucol
  vector<vector<complex<double>*> > acol(a.size());  //DCY
  vector<complex<double>*> ucol(2*nploc);	     //DCY

  //vector<vector<double*> > acol(a.size());
  //vector<double*> ucol(2*nploc);
  for ( int k = 0; k < a.size(); k++ )
  {
    acol[k].resize(2*nploc);
    for ( int i = 0; i < a[k]->nloc(); i++ )
     {
      acol[k][i] = a[k]->valptr(i*a[k]->mloc());
     }
    // if nloc is odd, store the address of vector 2*nploc-1
    if ( nloc_odd )
      acol[k][2*nploc-1] = &a_aux[k][0];
  }
  for ( int i = 0; i < u.nloc(); i++ )
   {
    ucol[i] = u.valptr(i*u.mloc());
   }
  // if nloc is odd, store the address of vector 2*nploc-1
  if ( nloc_odd )
    ucol[2*nploc-1] = &u_aux[0];

  // the vectors of the pair (top[i],bot[i]) are located at
  // addresses acol[top[i]] and acol[bot[i]]

  bool done = false;
  int nsweep = 0;
  // allocate matrix element packed array apq
  // apq[3*ipair   + k*3*nploc] = apq[k][ipair]
  // apq[3*ipair+1 + k*3*nploc] = app[k][ipair]
  // apq[3*ipair+2 + k*3*nploc] = aqq[k][ipair]
  vector<complex<double> > apq(a.size()*3*nploc); //DCY
  vector<double > tapq(a.size()*3*2*nploc); //DCY

  double diag_sum = 0.0, previous_diag_sum = 0.0;
  while ( !done )
  {
    // sweep: process local pairs and rotate 2*np-1 times
    nsweep++;
    double diag_change = 0.0;
    for ( int irot = 0; irot < 2*np-1; irot++ )
    {
/*
      cout << ctxt.mype() << ": top[i]: ";
      for ( int i = 0; i < nploc; i++ )
        cout << setw(3) << top[i];
      cout << endl;

      cout << ctxt.mype() << ": bot[i]: ";
      for ( int i = 0; i < nploc; i++ )
        cout << setw(3) << bot[i];
      cout << endl;

      cout << ctxt.mype() << ": jglobal[top[i]]: ";
      for ( int i = 0; i < nploc; i++ )
        cout << setw(3) << jglobal[top[i]];
      cout << endl;

      cout << ctxt.mype() << ": jglobal[bot[i]]: ";
      for ( int i = 0; i < nploc; i++ )
        cout << setw(3) << jglobal[bot[i]];
      cout << endl;
*/
      // perform Jacobi rotations for all local pairs

      // compute off-diagonal matrix elements apq for all pairs
      // skip the pair if one or both of the vectors is a dummy vector
      // i.e. a vector having jglobal==-1

      int mloc = a[0]->mloc();

#pragma omp parallel for 
      for ( int k = 0; k < a.size(); k++ )
      {
        for ( int ipair = 0; ipair < nploc; ipair++ )
        {
          const int iapq = 3*ipair + k*3*nploc;
          apq[iapq] = (0.0,0.0);
          apq[iapq+1] = (0.0,0.0);
          apq[iapq+2] = (0.0,0.0);
          if ( jglobal[top[ipair]] >= 0 && jglobal[bot[ipair]] >= 0 )
          {
            const complex<double> *ap = acol[k][top[ipair]];
            const complex<double> *aq = acol[k][bot[ipair]];
	    const complex<double> *up = ucol[top[ipair]];
            const complex<double> *uq = ucol[bot[ipair]];
	    int one = 1;
	    for (int ii=0; ii<mloc; ii++)
	    	{
		  apq[iapq]   += conj(ap[ii])*uq[ii];
		  apq[iapq+1] += conj(ap[ii])*up[ii];
		  apq[iapq+2] += conj(aq[ii])*uq[ii];
		}
          }
        }
      } // for k
      // apq now contains partial sums of matrix elements
      // create proxy array tapq with double the length of apq
      // so "dsum" will work
      // MPI_Allreduce works only for npcol = 1
      tm_comm.start();
      int len = apq.size()*2;
      double *tapq = (double*) &apq[0];
      ctxt.dsum('c',len,1,&tapq[0],len);
      tm_comm.stop();

      // apq now contains the matrix elements
#pragma omp parallel for
      for ( int ipair = 0; ipair < nploc; ipair++ )
      {
        if ( jglobal[top[ipair]] >= 0 && 
	     jglobal[bot[ipair]] >= 0 )
        {
          // compute rotation sine and cosine
          // Cardoso-Souloumiac expressions for the rotation angle

          // compute 3x3 matrix g
	  // ^ DCY. matrix g will no longer be 2x2 because of complex matrices. 
          double g11 = 0.0, g12 = 0.0, g13 = 0.0;
          double g21 = 0.0, g22 = 0.0, g23 = 0.0;
          double g31 = 0.0, g32 = 0.0, g33 = 0.0;

          for ( int k = 0; k < a.size(); k++ )
          {
            const int iapq = 3*ipair + k*3*nploc;
	    const complex<double> aij(apq[iapq]);
	    const complex<double> aii(apq[iapq+1]);
            const complex<double> ajj(apq[iapq+2]); 

	    const complex<double> i(0.0, 1.0);
            const complex<double> h1(aii - ajj); 
            const complex<double> h2(aij + conj(aij)); 
	    const complex<double> h3(i * (aij - conj(aij)));
	    //conjugate transposes of h1, h2, h3
	    const complex<double> h1ct(conj(h1));
	    const complex<double> h2ct(conj(h2));
	    const complex<double> h3ct(conj(h3));

            g11 += real(h1ct * h1);
            g12 += real(h1ct * h2);
            g13 += real(h1ct * h3);
	    g21 += real(h2ct * h1);
	    g22 += real(h2ct * h2);
	    g23 += real(h2ct * h3);
	    g31 += real(h3ct * h1);
	    g32 += real(h3ct * h2);
	    g33 += real(h3ct * h3);
          }

	int N = 3;
	double G[3*3] = { g11, g12, g13,   // matrix to be diagonalized
			  g21, g22, g23,
			  g31, g32, g33 };
	double D[3]; //array of eigenvalues, in ascending order
	int it_max = 10000;
	int it_num;
	int rot_num;
	double Q[3*3]; //matrix of eigenvectors
	jacobi_eigenvalue(3, G, it_max, Q, D, it_num, rot_num);

//DCY
/* 
	int Gmb = 3/ctxt.nprow();
	int Gnb = 3/ctxt.npcol();
	DoubleMatrix Gmat(ctxt,3,3,Gnb,Gnb);
	DoubleMatrix Qmat(ctxt,3,3,Gnb,Gnb);

	int Gind = 0;
	for ( int m = 0; m < Gmat.nblocks(); m++ )
	  for ( int l = 0; l < Gmat.mblocks(); l++ )
	    for ( int y = 0; y < Gmat.nbs(m); y++ )
	      for ( int x = 0; x < Gmat.mbs(l); x++ ) {
		int i = Gmat.i(l,x);
		int j = Gmat.j(m,y);
		int iii = x + l*Gmat.mb();
		int jjj = y + m*s.nb();
		int ival = iii + jjj * G.mloc();
		Gmat[ival] = G[Gind];
		Gind += 1;
		}
	Gind = 0; 
*/
	
////////////////////////////////////////////////////////////////////////////////

	// get eigenvector associated with largest eigenvalue of G
	double maxeig = 0.0;
	maxeig = D[2];
	double x, y, z;
	x = Q[6];
	y = Q[7];
	z = Q[8];


	// choose eigenvector with x positive to ensure small angle
	if ( x < 0.0 )
	{
	 x = -x; 
	 y = -y;
	 z = z;
	}
	

	double r = 1.0;//sqrt(x*x + y*y + z*z);
	complex<double> c(sqrt((x + r)/(2.0 * r)), 0.0);
	complex<double> s( y/sqrt(2.0 * r * (x + r)), -1.0*z/sqrt(2.0 * r * (x + r)));
	//complex<double> s( y / sqrt(2.0*(x + r)), 0.0);
          // apply the rotation R(p,q)
          //
          //           |  c      s |
          //  R(p,q) = |           |
          //           | -s      c |

          // U := U * R(p,q)^T
          // A := A * R(p,q)^T

          // apply rotation to columns of a and u

          // the drot function computes
          // c*x + s*y -> x
          //-s*x + c*y -> y
          // call drot with args c, conj(s)

	  complex<double> sconj(conj(s));
          int one = 1;
          for ( int k = 0; k < a.size(); k++ )
          {
            complex<double> *ap = acol[k][top[ipair]];
            complex<double> *aq = acol[k][bot[ipair]];
            int mloc = a[k]->mloc();
            zrot(&mloc,ap,&one,aq,&one,&c,&sconj);
          }
          complex<double> *up = ucol[top[ipair]];
          complex<double> *uq = ucol[bot[ipair]];
          int mloc = u.mloc();
	  zrot(&mloc,up,&one,uq,&one,&c,&sconj);

          // new value of off-diag element apq
          double diag_change_ipair = 0.0;
          for ( int k = 0; k < a.size(); k++ )
          {
            const int iapq = 3*ipair + k*3*nploc;
            const complex<double> aii(apq[iapq+1]);
            const complex<double> ajj(apq[iapq+2]);
	    const complex<double> v1(conj(c)*c - conj(s)*s);
	    const double apqnew = real(v1 * (aii - ajj) + 2.0 * c * s * apq[iapq] + 2.0 * conj(s) * c * apq[iapq]);
	    // accumulate change in sum of squares of diag elements
            // note negative sign: decrease in offdiag is increase in diag
            diag_change_ipair += 2.0 * abs(apqnew - real(aii - ajj));
	    //diag_change -= 2.0 * ( real(conj(apqnew)*apqnew) - real(conj(apq[iapq])*apq[iapq]) ); 
          }
#pragma omp critical 
	  diag_change += diag_change_ipair;
        }
      } // for ipair

      // all local pairs have been processed
      // rotate top and bot arrays
      if ( nploc > 0 )
      {
        bot.push_back(top.back());
        top.pop_back();
        top.push_front(bot.front());
        bot.pop_front();
        // make rotation skip element 0 on the first process column
        // if my process column is zero, swap top[0] and top[1]
        if ( ctxt.mycol() == 0 )
        {
          if ( nploc > 1 )
          {
            int tmp = top[0];
            top[0] = top[1];
            top[1] = tmp;
          }
          else
          {
            // if there is only one local pair, exchange top[0] and bot[0]
            int tmp = top[0];
            top[0] = bot[0];
            bot[0] = tmp;
          }
        }
        // exchange columns of a[k] and u

        int rbufi_left, rbufi_right, sbufi_left, sbufi_right;
        // send buffers contain k columns of a and one of u
	int bufsize = (a.size()+1)*a[0]->mloc();
	//DCY initiate proxy double arrays with twice the length
	//so that dsend and drecv blacs routines will work?
        vector< complex<double> > sbuf_left(bufsize), sbuf_right(bufsize);
        vector< complex<double> > rbuf_left(bufsize), rbuf_right(bufsize);
	
        // on each task except mycol==npcol-1
        // send jglobal[bot[nploc-1]] to the right
        // if jglobal != -1 send vector bot[nploc-1] to the right

        // on each task except mycol==npcol-1
        // recv jglobal from the right
        // if jglobal != -1 recv a vector from the right into bot[nploc-1]
        // set value of jglobal[bot[nploc-1]]

        // on each task except mycol==0
        // send jglobal[top[0]] to the left
        // if jglobal != -1 send vector top[0] to the left

        // on each task except mycol==0
        // recv jglobal from the left
        // if jglobal != -1 recv a vector from the left into top[0]
        // set value of jglobal[top[0]]
        // exchange jglobal values first
        tm_comm.start();
        if ( ctxt.mycol() < last_active_process_col )
        {
          sbufi_right = jglobal[bot[nploc-1]];
          ctxt.isend(1,1,&sbufi_right,1,ctxt.myrow(),ctxt.mycol()+1);
          ctxt.irecv(1,1,&rbufi_right,1,ctxt.myrow(),ctxt.mycol()+1);
          jglobal[bot[nploc-1]] = rbufi_right;
          //cout << ctxt.mype() << ": received jglobal="
          //               << jglobal[bot[nploc-1]] << " from right" << endl;
        }
        if ( ctxt.mycol() != 0 )
        {
          sbufi_left = jglobal[top[0]];
          ctxt.isend(1,1,&sbufi_left,1,ctxt.myrow(),ctxt.mycol()-1);
          ctxt.irecv(1,1,&rbufi_left,1,ctxt.myrow(),ctxt.mycol()-1);
          jglobal[top[0]] = rbufi_left;
          //cout << ctxt.mype() << ": received jglobal="
          //            << jglobal[top[0]] << " from left" << endl;
        }

        // exchange column vectors
        if ( ctxt.mycol() < last_active_process_col )
        {
          for ( int k = 0; k < a.size(); k++ )
          {
            memcpy(&sbuf_right[k*mloc],acol[k][bot[nploc-1]],
                   mloc*sizeof(complex<double>));
          }
          memcpy(&sbuf_right[a.size()*mloc], ucol[bot[nploc-1]],
                 mloc*sizeof(complex<double>) );
          ctxt.zsend(bufsize,1,&sbuf_right[0],bufsize,ctxt.myrow(),ctxt.mycol()+1);
          ctxt.zrecv(bufsize,1,&rbuf_right[0],bufsize,ctxt.myrow(),ctxt.mycol()+1);
          for ( int k = 0; k < a.size(); k++ )
          {
            memcpy(acol[k][bot[nploc-1]],&rbuf_right[k*mloc],
                   mloc*sizeof(complex<double>));
          }
          memcpy(ucol[bot[nploc-1]], &rbuf_right[a.size()*mloc],
                 mloc*sizeof(complex<double>) );
        }
        if ( ctxt.mycol() != 0 )
        { 
          for ( int k = 0; k < a.size(); k++ )
          {
            memcpy(&sbuf_left[k*mloc],acol[k][top[0]],mloc*sizeof(complex<double>));
          }
          memcpy(&sbuf_left[a.size()*mloc],ucol[top[0]],mloc*sizeof(complex<double>) );
	  ctxt.zsend(bufsize,1,&sbuf_left[0],bufsize,ctxt.myrow(),ctxt.mycol()-1);
	  ctxt.zrecv(bufsize,1,&rbuf_left[0],bufsize,ctxt.myrow(),ctxt.mycol()-1);
          for ( int k = 0; k < a.size(); k++ )
          {
            memcpy(acol[k][top[0]],&rbuf_left[k*mloc],mloc*sizeof(complex<double>) );
          }
          memcpy(ucol[top[0]],&rbuf_left[a.size()*mloc],mloc*sizeof(complex<double>) );
        }
        tm_comm.stop();
      } // if nploc > 0
      // end of step
	
    } // for irot
    // sweep is complete
    tm_comm.start();
    ctxt.dsum('r',1,1,&diag_change,1);
    tm_comm.stop();

    // compute sum of squares of diagonal elements

    if ( debug_diag_sum )
    {
      // compute sum of squares of diagonal elements using current values
      // (after rotation)
      previous_diag_sum = diag_sum;
      diag_sum = 0.0;
      for ( int k = 0; k < a.size(); k++ )
      {
        for ( int ipair = 0; ipair < nploc; ipair++ )
        {
          complex<double> tmp[2] = {(0.0, 0.0), (0.0, 0.0)};
          // compute the diagonal elements
          // skip dummy vectors
          int one = 1;
          //int mloc = a[k]->mloc();
          if ( jglobal[top[ipair]] >= 0 )
          {
            const complex<double> *ap = acol[k][top[ipair]];
            const complex<double> *up = ucol[top[ipair]];
	    //tmp[0] = zdotc(&mloc,ap,&one,up,&one);
            for ( int i = 0; i < mloc; i++){ tmp[0] += conj(ap[i]) * up[i];}
          }
          if ( jglobal[bot[ipair]] >= 0 )
          {
            const complex<double> *aq = acol[k][bot[ipair]];
            const complex<double> *uq = ucol[bot[ipair]];
	    //tmp[1] = zdotc(&mloc,aq,&one,uq,&one);
	    for ( int i = 0; i < mloc; i++){ tmp[1] += conj(aq[i]) * uq[i];}
          }
          // tmp now contains partial sums of app and aqq
          //ctxt.zsum('c',2,1,tmp,2);
          // tmp now contains the diagonal elements app and aqq
          //diag_sum += real(conj(tmp[0])*tmp[0]) + real(conj(tmp[1])*tmp[1]);
	  diag_sum += norm(conj(tmp[0])*tmp[0] + conj(tmp[1])*tmp[1]);
        }
      }
     ctxt.dsum('r',1,1,&diag_sum,1);
     const double diag_sum_increase = diag_sum - previous_diag_sum;
      if ( ctxt.onpe0() )
        cout << " jade: nsweep=" << nsweep
             << "zsum: "
             << setw(15) << setprecision(10) << diag_sum
             << " zsum_inc: "
             << setw(15) << setprecision(10) << diag_sum_increase << endl;
    }
     if ( ctxt.onpe0() )
       cout << " jade: nsweep=" << nsweep
           << " dchange: "
           << setw(15) << setprecision(10) << diag_change << endl;

    done = ( ( fabs(diag_change) < tol ) || ( nsweep >= maxsweep ) );

  } // while !done
  // if a dummy vector was used, (i.e. if nloc_odd), the dummy vector
  // may end up anywhere in the array after all rotations are completed.
  // The array a_aux may contain a (non-dummy) vector.

  // rotate columns of a and u to restore original order
  //complex<double> *tmpmat1 = new complex<double>[nloc*mloc];
  // rotate columns of a
  //for ( int k = 0; k < a.size(); k++ )
  //{
  if ( nloc_odd )
  {
    // find position of the dummy vector and copy a_aux onto it
    int idum = 0;
    while ( jglobal[idum] != -1 && idum < 2*nploc ) idum++;
    //cout << ctxt.mype() << ": idum=" << idum << endl;
    if ( idum != 2*nploc-1 )
    {
      for ( int k = 0; k < a.size(); k++ )
      {
        memcpy(acol[k][idum],&a_aux[k][0],mloc*sizeof(complex<double>));
      }
      memcpy(ucol[idum],&u_aux[0],mloc*sizeof(complex<double>));
    }
  }
/*
  // rotate columns of a and u to restore original order
  //complex<double> *tmpmat = new complex<double>[nloc*mloc];
  // rotate columns of a
  for ( int k = 1; k < a.size(); k++ )
  {
    for ( int ipair = 0; ipair < nploc; ipair++ )
    {
      // copy columns of a[k] to temporary array tmpmat in original order
      if ( jglobal[top[ipair]] >= 0 )
      {
        memmove(&tmpmat[ipair*mloc], acol[k][top[ipair]], mloc*sizeof(complex<double>));
      }
      if ( jglobal[bot[nploc-ipair-1]] >= 0 )
      {
        memmove(&tmpmat[(nploc+ipair)*mloc],acol[k][bot[nploc-ipair-1]],
               mloc*sizeof(complex<double>));
      }
    }
    // copy tmpmat back to a[k]
    memmove(acol[k][0],tmpmat,nloc*mloc*sizeof(complex<double>));
  }*/
    for ( int ipair = 0; ipair < nploc; ipair++ )
    {
      // copy columns of a[k] to temporary array tmpmat in original order
      if ( jglobal[top[ipair]] >= 0 )
      {
        memcpy(&atmp1[ipair*mloc], acol[0][top[ipair]], mloc*sizeof(complex<double>));
      }
      if ( jglobal[bot[nploc-ipair-1]] >= 0 )
      {
        memcpy(&atmp1[(nploc+ipair)*mloc],acol[0][bot[nploc-ipair-1]],
               mloc*sizeof(complex<double>));
      }
    }
    for ( int ipair = 0; ipair < nploc; ipair++ )
    {
      // copy columns of a[1] to temporary array tmpmat in original order
      if ( jglobal[top[ipair]] >= 0 )
      {
        memcpy(&atmp2[ipair*mloc], acol[1][top[ipair]], mloc*sizeof(complex<double>));
      }
      if ( jglobal[bot[nploc-ipair-1]] >= 0 )
      {
        memcpy(&atmp2[(nploc+ipair)*mloc],acol[1][bot[nploc-ipair-1]],
               mloc*sizeof(complex<double>));
      }
    }
    for ( int ipair = 0; ipair < nploc; ipair++ )
    {
      // copy columns of a[2] to temporary array tmpmat in original order
      if ( jglobal[top[ipair]] >= 0 )
      {
        memcpy(&atmp3[ipair*mloc], acol[2][top[ipair]], mloc*sizeof(complex<double>));
      }
      if ( jglobal[bot[nploc-ipair-1]] >= 0 )
      {
        memcpy(&atmp3[(nploc+ipair)*mloc],acol[2][bot[nploc-ipair-1]],
               mloc*sizeof(complex<double>));
      }
    }
    for ( int ipair = 0; ipair < nploc; ipair++ )
    {
      // copy columns of a[3] to temporary array tmpmat in original order
      if ( jglobal[top[ipair]] >= 0 )
      {
        memcpy(&atmp4[ipair*mloc], acol[3][top[ipair]], mloc*sizeof(complex<double>));
      }
      if ( jglobal[bot[nploc-ipair-1]] >= 0 )
      {
        memcpy(&atmp4[(nploc+ipair)*mloc],acol[3][bot[nploc-ipair-1]],
               mloc*sizeof(complex<double>));
      }
    }
    for ( int ipair = 0; ipair < nploc; ipair++ )
    {
      // copy columns of a[4] to temporary array tmpmat in original order
      if ( jglobal[top[ipair]] >= 0 )
      {
        memcpy(&atmp5[ipair*mloc], acol[4][top[ipair]], mloc*sizeof(complex<double>));
      }
      if ( jglobal[bot[nploc-ipair-1]] >= 0 )
      {
        memcpy(&atmp5[(nploc+ipair)*mloc],acol[4][bot[nploc-ipair-1]],
               mloc*sizeof(complex<double>));
      }
    }
    for ( int ipair = 0; ipair < nploc; ipair++ )
    {
      // copy columns of a[5] to temporary array tmpmat in original order
      if ( jglobal[top[ipair]] >= 0 )
      {
        memcpy(&atmp6[ipair*mloc], acol[5][top[ipair]], mloc*sizeof(complex<double>));
      }
      if ( jglobal[bot[nploc-ipair-1]] >= 0 )
      {
        memcpy(&atmp6[(nploc+ipair)*mloc],acol[5][bot[nploc-ipair-1]],
               mloc*sizeof(complex<double>));
      }
    }
    // copy tmpmat back to a[k]
    //memcpy(acol[k][0],tmpmat1,nloc*mloc*sizeof(complex<double>));
  //}
  
  // rotate columns of u
  for ( int ipair = 0; ipair < nploc; ipair++ )
  {
    // copy columns of u to temporary array tmpmat in original order
    if ( jglobal[top[ipair]] >= 0 )
    {
      memcpy(&tmpmat[ipair*mloc], ucol[top[ipair]], mloc*sizeof(complex<double>));
    }
    if ( jglobal[bot[nploc-ipair-1]] >= 0 )
    {
      memcpy(&tmpmat[(nploc+ipair)*mloc],ucol[bot[nploc-ipair-1]],
             mloc*sizeof(complex<double>));
    }
  }
  // copy tmpmat back to u
  //memcpy(ucol[0],&tmpmat,nloc*mloc*sizeof(complex<double>));
  //delete [] tmpmat1;
  {
    for ( int i = 0; i < atmp1.n(); i++ )
     {
      adiag[0][i] = (0.0,0.0);
     }
    for ( int jblock = 0; jblock < atmp1.nblocks(); jblock++ )
      for ( int y = 0; y < atmp1.nbs(jblock); y++ )
      {
        // j is the global column index
        int j = atmp1.j(jblock,y);
        int jjj = y + jblock*atmp1.nb();
        const complex<double> *ap = atmp1.valptr(jjj*atmp1.mloc());
        const complex<double> *up = tmpmat.valptr(jjj*tmpmat.mloc());
        //int mloc = a[k]->mloc();
        int one = 1;
        for (int ii = 0; ii < mloc; ii++)
        {
          adiag[0][j] += conj(ap[ii])*up[ii];
        }
      }

    // adiag[k][i] now contains the partial sums of the diagonal elements of a
    tm_comm.start();
    int len = 2*(atmp1.n());
    //double *tadiag = (double*) &adiag[k][0];
    //ctxt.dsum('c',len,1,&tadiag[0],len);
    MPI_Allreduce(MPI_IN_PLACE,&adiag[0][0],len,MPI_DOUBLE,MPI_SUM,ctxt.comm());
    tm_comm.stop();
    // adiag[k] contains the diagonal elements of a[k]
    // u contains the orthogonal transformation minimizing the spread
  }
  {
    for ( int i = 0; i < atmp2.n(); i++ )
     {
      adiag[1][i] = (0.0,0.0);
     }
    for ( int jblock = 0; jblock < atmp2.nblocks(); jblock++ )
      for ( int y = 0; y < atmp2.nbs(jblock); y++ )
      {
        // j is the global column index
        int j = atmp2.j(jblock,y);
        int jjj = y + jblock*atmp2.nb();
        const complex<double> *ap = atmp2.valptr(jjj*atmp2.mloc());
        const complex<double> *up = tmpmat.valptr(jjj*tmpmat.mloc());
        //int mloc = a[k]->mloc();
        int one = 1;
        for (int ii = 0; ii < mloc; ii++)
        {
          adiag[1][j] += conj(ap[ii])*up[ii];
        }
      }

    // adiag[k][i] now contains the partial sums of the diagonal elements of a
    tm_comm.start();
    int len = 2*(atmp2.n());
    //double *tadiag = (double*) &adiag[k][0];
    //ctxt.dsum('c',len,1,&tadiag[0],len);
    MPI_Allreduce(MPI_IN_PLACE,&adiag[1][0],len,MPI_DOUBLE,MPI_SUM,ctxt.comm());
    tm_comm.stop();
    // adiag[k] contains the diagonal elements of a[k]
    // u contains the orthogonal transformation minimizing the spread
  }
  {
    for ( int i = 0; i < atmp3.n(); i++ )
     {
      adiag[2][i] = (0.0,0.0);
     }
    for ( int jblock = 0; jblock < atmp3.nblocks(); jblock++ )
      for ( int y = 0; y < atmp3.nbs(jblock); y++ )
      {
        // j is the global column index
        int j = atmp3.j(jblock,y);
        int jjj = y + jblock*atmp3.nb();
        const complex<double> *ap = atmp3.valptr(jjj*atmp3.mloc());
        const complex<double> *up = tmpmat.valptr(jjj*tmpmat.mloc());
        //int mloc = a[k]->mloc();
        int one = 1;
        for (int ii = 0; ii < mloc; ii++)
        {
          adiag[2][j] += conj(ap[ii])*up[ii];
        }
      }

    // adiag[k][i] now contains the partial sums of the diagonal elements of a
    tm_comm.start();
    int len = 2*(atmp3.n());
    //double *tadiag = (double*) &adiag[k][0];
    //ctxt.dsum('c',len,1,&tadiag[0],len);
    MPI_Allreduce(MPI_IN_PLACE,&adiag[2][0],len,MPI_DOUBLE,MPI_SUM,ctxt.comm());
    tm_comm.stop();
    // adiag[k] contains the diagonal elements of a[k]
    // u contains the orthogonal transformation minimizing the spread
  }
  {
    for ( int i = 0; i < atmp4.n(); i++ )
     {
      adiag[3][i] = (0.0,0.0);
     }
    for ( int jblock = 0; jblock < atmp4.nblocks(); jblock++ )
      for ( int y = 0; y < atmp4.nbs(jblock); y++ )
      {
        // j is the global column index
        int j = atmp4.j(jblock,y);
        int jjj = y + jblock*atmp4.nb();
        const complex<double> *ap = atmp4.valptr(jjj*atmp4.mloc());
        const complex<double> *up = tmpmat.valptr(jjj*tmpmat.mloc());
        //int mloc = a[k]->mloc();
        int one = 1;
        for (int ii = 0; ii < mloc; ii++)
        {
          adiag[3][j] += conj(ap[ii])*up[ii];
        }
      }

    // adiag[k][i] now contains the partial sums of the diagonal elements of a
    tm_comm.start();
    int len = 2*(atmp4.n());
    //double *tadiag = (double*) &adiag[k][0];
    //ctxt.dsum('c',len,1,&tadiag[0],len);
    MPI_Allreduce(MPI_IN_PLACE,&adiag[3][0],len,MPI_DOUBLE,MPI_SUM,ctxt.comm());
    tm_comm.stop();
    // adiag[k] contains the diagonal elements of a[k]
    // u contains the orthogonal transformation minimizing the spread
  }
  {
    for ( int i = 0; i < atmp5.n(); i++ )
     {
      adiag[4][i] = (0.0,0.0);
     }
    for ( int jblock = 0; jblock < atmp5.nblocks(); jblock++ )
      for ( int y = 0; y < atmp5.nbs(jblock); y++ )
      {
        // j is the global column index
        int j = atmp5.j(jblock,y);
        int jjj = y + jblock*atmp5.nb();
        const complex<double> *ap = atmp5.valptr(jjj*atmp5.mloc());
        const complex<double> *up = tmpmat.valptr(jjj*tmpmat.mloc());
        //int mloc = a[k]->mloc();
        int one = 1;
        for (int ii = 0; ii < mloc; ii++)
        {
          adiag[4][j] += conj(ap[ii])*up[ii];
        }
      }

    // adiag[k][i] now contains the partial sums of the diagonal elements of a
    tm_comm.start();
    int len = 2*(atmp5.n());
    //double *tadiag = (double*) &adiag[k][0];
    //ctxt.dsum('c',len,1,&tadiag[0],len);
    MPI_Allreduce(MPI_IN_PLACE,&adiag[4][0],len,MPI_DOUBLE,MPI_SUM,ctxt.comm());
    tm_comm.stop();
    // adiag[k] contains the diagonal elements of a[k]
    // u contains the orthogonal transformation minimizing the spread
  }
  {
    for ( int i = 0; i < atmp6.n(); i++ )
     {
      adiag[5][i] = (0.0,0.0);
     }
    for ( int jblock = 0; jblock < atmp6.nblocks(); jblock++ )
      for ( int y = 0; y < atmp6.nbs(jblock); y++ )
      {
        // j is the global column index
        int j = atmp6.j(jblock,y);
        int jjj = y + jblock*atmp6.nb();
        const complex<double> *ap = atmp6.valptr(jjj*atmp6.mloc());
        const complex<double> *up = tmpmat.valptr(jjj*tmpmat.mloc());
        //int mloc = a[k]->mloc();
        int one = 1;
        for (int ii = 0; ii < mloc; ii++)
        {
          adiag[5][j] += conj(ap[ii])*up[ii];
        }
      }

    // adiag[k][i] now contains the partial sums of the diagonal elements of a
    tm_comm.start();
    int len = 2*(atmp6.n());
    //double *tadiag = (double*) &adiag[k][0];
    //ctxt.dsum('c',len,1,&tadiag[0],len);
    MPI_Allreduce(MPI_IN_PLACE,&adiag[5][0],len,MPI_DOUBLE,MPI_SUM,ctxt.comm());
    tm_comm.stop();
    // adiag[k] contains the diagonal elements of a[k]
    // u contains the orthogonal transformation minimizing the spread
  } /*
  // compute diagonal values
  for ( int k = 0; k < a.size(); k++ )
  {
    for ( int i = 0; i < a[k]->n(); i++ )
     {
      adiag[k][i] = (0.0,0.0);
     }
    for ( int jblock = 0; jblock < a[k]->nblocks(); jblock++ )
      for ( int y = 0; y < a[k]->nbs(jblock); y++ )
      {
        // j is the global column index
        int j = a[k]->j(jblock,y);
        int jjj = y + jblock*a[k]->nb();
        const complex<double> *ap = a[k]->valptr(jjj*a[k]->mloc());
        const complex<double> *up = u.valptr(jjj*u.mloc());
        //int mloc = a[k]->mloc();
        int one = 1;
	for (int ii = 0; ii < mloc; ii++) 
	{
	  adiag[k][j] += conj(ap[ii])*up[ii];
	}
      }

    // adiag[k][i] now contains the partial sums of the diagonal elements of a
    tm_comm.start();
    int len = 2*(a[k]->n());
    //double *tadiag = (double*) &adiag[k][0];
    //ctxt.dsum('c',len,1,&tadiag[0],len);
    MPI_Allreduce(MPI_IN_PLACE,&adiag[k][0],len,MPI_DOUBLE,MPI_SUM,ctxt.comm());
    tm_comm.stop();
    // adiag[k] contains the diagonal elements of a[k]
    // u contains the orthogonal transformation minimizing the spread
  }
*/
  if ( ctxt.onpe0() )
    cout << " jade_complex: comm time: " << tm_comm.real() << endl;
  return nsweep;
}

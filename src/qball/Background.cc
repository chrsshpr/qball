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
// Background.C
//
////////////////////////////////////////////////////////////////////////////////



#include <config.h>


#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
#include <vector>
#include <cmath>
#include <complex>
#include <math/d3vector.h>
#include "Background.h"
#include "UnitCell.h"
#include "FourierTransform.h"
using namespace std;


const double PI  =3.141592653589793238463;

Jellium::Jellium(string input)
{
      stringstream ss(input);
      string tempbuf;
      vector <string> info;
      int count=0;
      while ( ss >> tempbuf )
      {
           info.push_back(tempbuf);
           count++;
      }
      charge_= stod(info[1]);
      for(int i = 0; i < count; i++)
      {
          if(info[i] == "X" || info[i] == "x")
          {
              x=stod(info[i+1]);
          }
          if(info[i] == "y" || info[i] == "Y")
          {
              y=stod(info[i+1]);
          }
          if(info[i] == "Z" || info[i] == "z")
          {
              z=stod(info[i+1]);
          }
          if(info[i] == "sigma" || info[i] == "SIGMA")
          {
              sigma = true; 
              sigma_degree =stoi(info[i+1]);
          }
      }
}
void Jellium::update_rho_r(const UnitCell & cell,FourierTransform &vft, int ngloc)
{
      const int np012loc = vft.np012loc();
      rho_r.resize(np012loc);
      const int np0 = vft.np0();
      const int np1 = vft.np1();
      const int np2 = vft.np2();
      int idx0 = vft.np0() * vft.np1() * vft.np2_first();
      int idxx, i, j, k;
      double x_length,y_length,z_length;
      D3vector tmp;
      bool not_boundary;
      if (x==0)
      {       
             tmp = cell.a(0);
             x_length = tmp[0]/2;
      }else{  
             tmp = cell.a(0);
             x_length = x;
      }
      if (y==0)
      {
             tmp = cell.a(1);
             y_length = tmp[1]/2;
      }else{
             tmp = cell.a(1);
             y_length = y;
      }
      if (z==0)
      {
             tmp = cell.a(2);
             z_length = tmp[2]/2;
      }else{
             tmp = cell.a(2);
             z_length=z;
      }
      
      const double omega = 8*x_length*y_length*z_length;
      const double omega_inv = 1.0 / omega;
      const double local_charge = charge_* omega_inv;
      cout<< "length\t"<<x_length<<"\t"<<y_length<<"\t"<<z_length<<endl;
      for ( int ir = 0; ir < np012loc; ir++ ) 
      {
            D3vector r;
            idxx = idx0 + ir;
            k = idxx / ( np0 * np1 );
            idxx = idxx - ( np0 * np1 ) * k;
            j = idxx / np0;
            idxx = idxx - np0 * j;
            i = idxx;
            r = cell.a(0) / np0 * i + cell.a(1) / np1 * j + cell.a(2) / np2 * k;
            not_boundary = (i!=np0 && j!=np1 && k!=np2);
            cell.fold_in_ws(r);
            if (abs(r[0])<=x_length && abs(r[1])<=y_length && abs(r[2])<=z_length)
                     rho_r[ir] = complex<double>(local_charge,0.0);
            else
                     rho_r[ir] = complex<double>(0.0,0.0);
      }
      //vector <complex<double>> for_plot3,tmp_plot;
      //for_plot3.resize(np012loc);
      //for_plot3= rho_r;
      //tmp_plot.resize(np012loc);
      //for_plot2 = for_plot3;
      //vft.forward(&for_plot3[0],&tmp_plot[0]);
      //vft.backward(&tmp_plot[0],&for_plot3[0]);
      //vft.sigma_forward(&for_plot2[0],&tmp_plot[0]);
      //vft.backward(&tmp_plot[0],&for_plot2[0]);
      //for ( int ir = 0; ir < np012loc; ir++ )
      //{
      //      D3vector r;
      //      idxx = idx0 + ir;
      //      k = idxx / ( np0 * np1 );
      //      idxx = idxx - ( np0 * np1 ) * k;
      //      j = idxx / np0;
      //      idxx = idxx - np0 * j;
      //      i = idxx;
      //      r = cell.a(0) / np0 * i + cell.a(1) / np1 * j + cell.a(2) / np2 * k;
      //      not_boundary = (i!=np0 && j!=np1 && k!=np2);
      //      cell.fold_in_ws(r);
     //       if (r[1]==0 && r[2]==0 )
     //               cout<<"plot\t"<<r[0]<<"\t"<<real(for_plot3[ir])<<"\t"<<real(rho_r[ir])<<endl;
     // }
}

// Local Variables:
// mode: c++
// End:

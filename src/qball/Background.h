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
// Background.h
//
////////////////////////////////////////////////////////////////////////////////



#include <config.h>

#ifndef BACKGROUND_H
#define BACKGROUND_H


#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
#include <vector>

#include <complex>
#include "UnitCell.h"
#include "FourierTransform.h"
using namespace std;
class UnitCell;
class FourierTransform;
class Jellium
{
   private:
   
   double x=0;
   double y=0;
   double z=0;
   double charge_;
     
   public:
   
   bool sigma = false;
   int  sigma_degree;
   Jellium(string input) ;
   ~Jellium(){ rho_r.clear(); } ;
   vector< complex<double> > rho_r;
   double charge (void) {return charge_;}
   void update_rho_r(const UnitCell & cell,FourierTransform &vft, int ngloc);
};
#endif

// Local Variables:
// mode: c++
// End:

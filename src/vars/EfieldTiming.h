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
// EfieldTiming.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef EFIELDTIMING_H
#define EFIELDTIMING_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
//#include<math.h>

#include <qball/Sample.h>

class EfieldTiming : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "efield_timing"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 3 )
    {
      if ( ui->oncoutpe() )
      cout << " efield_timing takes two value: t0 and full width half max" << endl;
      return 1;
    }
    
    double ev2au = 0.0367493;
    //double sig = 2 * sqrt(2*log(2));
    double sig = 2.3548;
    double convert = atof(argv[2])*ev2au;
    double v0 = atof(argv[1]);
    double v1 = sig /convert; 

    if ( v0 < 0.0 || v1 < 0.0 )
    {
      if ( ui->oncoutpe() )
        cout << "timing values must be > 0" << endl;
      return 1;
    }

    s->ctrl.efield_timing[0] = v0;
    s->ctrl.efield_timing[1] = v1;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.efield_timing[0] << " " 
        << s->ctrl.efield_timing[1] << " ";
     return st.str();
  }

  EfieldTiming(Sample *sample) : s(sample)
  {
    s->ctrl.efield_timing[0] = 0;
    s->ctrl.efield_timing[1] = 0;
  }
};
#endif


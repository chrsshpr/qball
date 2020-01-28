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
// Polarization.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef POLARIZATION_H
#define POLARIZATION_H

#include<iostream>
#include<iomanip>
#include<sstream>

#include <qball/Sample.h>

class Polarization: public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "polarization"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      //if ( ui->onpe0() )
      cout << " polarization takes only one value" << endl;
      return 1;
    }

    string v = argv[1];

    if ( v == "OFF" || v == "MLWF" || v == "MLWF_REF" || v == "MLWF_REF_Q" ||
         v == "BERRY" || v == "TDMLWF" || v == "TDMLWF_REF")
      s->ctrl.polarization = v;
    else
    {
      //if ( ui->onpe0() )
      cout <<
      " polarization must be OFF, MLWF, TDMLWF, MLWF_REF, TDMLWF_REF, MLWF_REF_Q or BERRY" << endl;
      return 1;
    }
    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.polarization;
     return st.str();
  }

  Polarization(Sample *sample) : s(sample)
  {
    s->ctrl.polarization = "OFF";
  }
};
#endif

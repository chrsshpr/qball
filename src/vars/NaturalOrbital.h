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
// EfieldAnti.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NATURALORBITAL_H
#define NATURALORBITAL_H

#include<iostream>
#include<iomanip>
#include<sstream>

#include <qball/Sample.h>

class NaturalOrbital: public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "natural_orbital"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 &&  argc != 3)
    {
      //if ( ui->onpe0() )
      cout << "natural_orbital takes only one or two value" << endl;
      return 1;
    }

    string v = argv[1];

    if ( v == "ON")
      s->ctrl.natural_orbital = true;
    else if (v=="OFF")
       s->ctrl.natural_orbital = false;
    else
    {
      //if ( ui->onpe0() )
      cout <<
      " natural_orbital  must be OFF or ON" << endl;
      return 1;
    }
    if (argc == 3) 
    {   
       string w = argv[2];
       
       if (w == "-hole")
         s->ctrl.saveholestate = true;
       else
       {
         cout << "invalid commands" << endl;
         return 1;
        }
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

  NaturalOrbital(Sample *sample) : s(sample)
  {
    s->ctrl.natural_orbital = false;
    s->ctrl.saveholestate = false;
  }
};
#endif

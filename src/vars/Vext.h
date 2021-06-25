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
//
// Vext.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef VEXT_H
#define VEXT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>
#include <qball/ExternalPotential.h>

class Vext : public Var
{
  Sample *s;

  public:

  const char *name ( void ) const { return "vext"; };

  int set ( int argc, char **argv )
  {
    if ( argc > 2 )
    {
      if ( ui->oncoutpe() )
      cout << " vext takes only one value" << endl;
      return 1;
    }

    if ( !strcmp(argv[1],"NULL") )
    {
      // set vext NULL
      delete s->vext;
      s->vext = 0;
    }
    else
    {
      if ( s->vext )
        delete s->vext;
      s->vext = new ExternalPotential(*s,argv[1]);
    }

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     if ( s->vext )
     {
       st.setf(ios::right,ios::adjustfield);
       st << setw(10) << s->vext->filename();
     }
     return st.str();
  }

  Vext(Sample *sample) : s(sample) {}
};
#endif

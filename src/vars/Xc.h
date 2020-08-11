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
// Xc.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef XC_H
#define XC_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>

#ifdef HAVE_LIBXC //YY
#include "xc.h"
#endif //YY

class Xc : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "xc"; };

  int set ( int argc, char **argv )
  {
    string v; //YY
    v = argv[1];
#ifdef HAVE_LIBXC
    if ( (argc > 2) && (v == "LIBXC") ) {
      for ( int n = 2; n < argc ; n ++ ) {
        v = v + " " + argv[n];
        }
    }
#else
    if (false) {
    }
#endif
    else if ( argc != 2 )
    {
      if ( ui->oncoutpe() )
      cout << " <ERROR> xc takes only one value </ERROR>" << endl;
      return 1;
    }
    
    else if ( !( v == "LDA" || v == "PBE" || v == "PBEsol" || v == "PBErev" || v == "BLYP" || v == "HF" || v == "PBE0" || v == "RSH") )
    {
      if ( ui->oncoutpe() )
        cout << " <ERROR> xc must be LDA, PBE, PBEsol, PBErev , BLYP, HF, PBE0 or RSH </ERROR>" << endl;
      return 1;
    }
        
    s->ctrl.xc= v;
    if (v=="HF") s->ctrl.hf= 1.0;
    //if (v=="PBE0") s->ctrl.hf= 0.25;
    return 0;

  } //YY Allow LIBXC

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.xc;
     return st.str();
  }

  Xc(Sample *sample) : s(sample) { s->ctrl.xc = "LDA";s->ctrl.hf = 0.0; };
};
#endif

// Local Variables:
// mode: c++
// End:

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

#ifndef BACKGROUNDVAR_H
#define BACKGROUNDVAR_H


#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
#include <vector>

#include <complex>
#include <qball/Sample.h>

class BackgroundCharge : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "background_charge"; };

  int set ( int argc, char **argv )
  {
    string v;
    //if ( argc !=2 && argc !=4 && argc !=6 & argc !=8 )
    //{
    //  if ( ui->oncoutpe() )
    //  cout << " <ERROR> background_charge takes one value and x y z boundary </ERROR>" << endl;
    //  return 1;
    //}
    //else
    //{
    for (int i=0;i<argc;i++)
      v = v + " " + argv[i];
      s->ctrl.has_background_charge = true;
      s->ctrl.background_charge = v;
    //}

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << s->ctrl.background_charge;
     return st.str();
  }

  BackgroundCharge(Sample *sample) : s(sample) { s->ctrl.has_background_charge = false ; s->ctrl.background_charge="";};
};
#endif

// Local Variables:
// mode: c++
// End:

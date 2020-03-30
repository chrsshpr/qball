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
// SaveHoleFreq.h 
// Ruiyi
////////////////////////////////////////////////////////////////////////////////

#ifndef SAVEELECFREQ_H
#define SAVEELECFREQ_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/Sample.h>

class SaveElecFreq : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "saveelecfreq"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 3 && argc != 4)
    {
      if ( ui->oncoutpe() )
      cout << " saveholefreq takes only two or three values" << endl;
      return 1;
    }
    int v = atoi(argv[1]);
    s->ctrl.elecindex = v ;
    int w = atoi(argv[2]);
    s->ctrl.saveelecfreq = w ;

    if (argc == 4)
       s->ctrl.saveelecfilebase = argv[2];

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << name() << ":  ";
     st.setf(ios::right,ios::adjustfield);
     st << s->ctrl.saveelecfreq;
     return st.str();
  }

  SaveElecFreq(Sample *sample) : s(sample) {
     s->ctrl.saveelecfreq = -1 ;
     s->ctrl.elecindex = 0 ;
     s->ctrl.saveelecfilebase = "elec";
  }
};
#endif

//CS
// Local Variables:
// mode: c++
// End:

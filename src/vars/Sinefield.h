////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// Efield.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef SINEFIELD_H
#define SINEFIELD_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
#include <math/d3vector.h>

#include <qball/Sample.h>

class Sinefield : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "sine_field"; };

///input should be formatted direction amplitude period
  int set ( int argc, char **argv )
  {
    if ( argc != 4 )
    {
      //if ( ui->onpe0() )
      cout << " e_field takes 3 input:direction amplitude period" << endl;
      return 1;
    }
    if (argv[1] == "X" or argv[1] == "x" )
	     s->ctrl.sine_field[0] = 0;
    if (argv[1] == "Y" or argv[1] == "y" )
             s->ctrl.sine_field[0] = 1;
    if (argv[2] == "Z" or argv[2] == "z" )
             s->ctrl.sine_field[0] = 2;

    s->ctrl.sine_field[1] = atof(argv[2]);
    s->ctrl.sine_field[2] = atof(argv[3]);
    s->ctrl.compute_sine_field = true;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << "direction\t"  << s->ctrl.sine_field[0] << " "
        << "ampltitude\t" << s->ctrl.sine_field[1] << " "
        << "period\t"     << s->ctrl.sine_field[2] << " ";
     return st.str();
  }

  Sinefield(Sample *sample) : s(sample)
  {
    s->ctrl.compute_sine_field = false;	  
    s->ctrl.sine_field[0] = 0.0;
    s->ctrl.sine_field[1] = 0.0;
    s->ctrl.sine_field[2] = 0.0;
  }
};
#endif

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
// Cosinefield.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef COSINEFIELD_H
#define COSINEFIELD_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
#include <math/d3vector.h>

#include <qball/Sample.h>

class Cosinefield : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "cosine_field"; };

///input should be formatted direction amplitude period
  int set ( int argc, char **argv )
  {
    if ( argc != 4 )
    {
      ui->error ("cosine_field takes 4 input:direction strength frequency");
      return 1;
    }

    if (argv[1] == "X" or argv[1] == "x" )
	     s->ctrl.cosine_field[0] = 0;
    if (argv[1] == "Y" or argv[1] == "y" )
             s->ctrl.cosine_field[0] = 1;
    if (argv[1] == "Z" or argv[1] == "z" )
             s->ctrl.cosine_field[0] = 2;

    s->ctrl.cosine_field[1] = atof(argv[2]);
    s->ctrl.cosine_field[2] = atof(argv[3]);
//    s->ctrl.cosine_field[3] = atof(argv[4]);
//    s->ctrl.cosine_field[3] = atof(argv[5]);
    s->ctrl.compute_cosine_field = true;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << "direction\t"  << s->ctrl.cosine_field[0] << " "
        << "strength\t"   << s->ctrl.cosine_field[1] << " "
        << "frequency\t"  << s->ctrl.cosine_field[2] << " ";
//        << "gamma\t"      << s->ctrl.cosine_field[3] << " ";
//        << "t0\t"      << s->ctrl.cosine_field[4] << " ";
     return st.str();
  }

  Cosinefield(Sample *sample) : s(sample)
  {
    s->ctrl.compute_cosine_field = false;	  
    s->ctrl.cosine_field[0] = 0.0;
    s->ctrl.cosine_field[1] = 0.0;
    s->ctrl.cosine_field[2] = 0.0;
  }
};
#endif

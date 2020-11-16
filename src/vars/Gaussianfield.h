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
// Gaussianfield.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef GAUSSIANFIELD_H
#define GAUSSIANFIELD_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>
#include <math/d3vector.h>
//#include<math.h>

#include <qball/Sample.h>

class Gaussianfield : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "gaussian_field"; };

///input should be formatted direction amplitude period
  int set ( int argc, char **argv )
  {
    if ( argc != 4 )
    {
      ui->error ("gaussian_field takes 4 input:direction strength frequency");
      return 1;
    }

    if (argv[1] == "X" or argv[1] == "x" )
	     s->ctrl.gaussian_field[0] = 0;
    if (argv[1] == "Y" or argv[1] == "y" )
             s->ctrl.gaussian_field[0] = 1;
    if (argv[1] == "Z" or argv[1] == "z" )
             s->ctrl.gaussian_field[0] = 2;

    double ev2au = 0.0367493;
    s->ctrl.gaussian_field[1] = atof(argv[2]);
    double fwhm = atof(argv[3]);
    s->ctrl.gaussian_field[2] = fwhm*ev2au; //convert input in ev to au 
    s->ctrl.compute_gaussian_field = true;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << "direction\t"  << s->ctrl.gaussian_field[0] << " " 
        << "strength\t"   << s->ctrl.gaussian_field[1] << " "
        << "frequency\t"  << s->ctrl.gaussian_field[2] << " ";
     return st.str();
  }

  Gaussianfield(Sample *sample) : s(sample)
  {
    s->ctrl.compute_gaussian_field = false;	  
    s->ctrl.gaussian_field[0] = 0.0;
    s->ctrl.gaussian_field[1] = 0.0;
    s->ctrl.gaussian_field[2] = 0.0;
  }
};
#endif

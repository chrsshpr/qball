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
////
//// ApplyElectricFieldPulseCmd.cc:
////
//////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#include <ui/ApplyElectricFieldPulseCmd.h>
#include<iostream>
#include <qball/Context.h>
#include <qball/SlaterDet.h>
using namespace std;

int ApplyElectricFieldPulseCmd::action(int argc, char **argv)
{
  Wavefunction& wf = s->wf;
  SlaterDet& sd = *(wf.sd(0,0));

  if ( argc != 3 && argc != 4) 
  {
    if ( ui->oncoutpe() )
    cout << " <ERROR> apply_electric_field_pulse takes only two values </ERROR>" << endl;
    cout << " <!-- use: apply_electric_field_pulse e_direction e_strength -->" << endl;
    cout << " <!-- use: apply_electric_field_pulse e_direction e_strength state_index-->" << endl;
    return 1;
  }


  int e_direction = atoi(argv[1]);
  double e_strength = atof(argv[2]);
  int state_index;
  if ( argc == 4)
  {
    state_index = atoi(argv[3]);
  }
  else
  {
    state_index = -1;
  }

  if (wf.nkp() > 1 || wf.kpoint(0) != D3vector(0,0,0)) {
    if ( ui->oncoutpe() )
      cout << "<ERROR> apply_electric_field_pulse command only works for gamma-point calculations! </ERROR>" << endl;
    return 1;
  }

  if (not wf.force_complex_set()) {
    if ( ui->oncoutpe() )
      cout << "<ERROR> apply_electric_field_pulse command only works for complex wavefunction calculations! </ERROR>" << endl;
    return 1;
  }

  sd.apply_electric_field(e_direction, e_strength, state_index);

  return 0;
}

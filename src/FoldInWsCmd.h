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
// FoldInWsCmd.h:
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef FOLDINWSCMD_H
#define FOLDINWSCMD_H

#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

#include "UserInterface.h"
#include "Sample.h"

class FoldInWsCmd : public Cmd
{
  public:

  Sample *s;

  FoldInWsCmd(Sample *sample) : s(sample) {};

  char *name(void) const { return "fold_in_ws"; }
  char *help_msg(void) const
  {
    return
    "\n fold_in_ws\n\n"
    " syntax: fold_in_ws \n\n"
    "   The fold_in_ws command folds all atomic positions back in\n"
    "   the Wigner-Seitz cell of the current unit cell.\n\n";
  }

  int action(int argc, char **argv)
  {
    s->atoms.fold_in_ws();
    return 0;
  }
};
#endif

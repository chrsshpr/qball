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
//////////////////////////////////////////////////////////////////////////////////
//
// Bisection.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef BISECTION_H
#define BISECTION_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>

#include <qball/Context.h>
#include <qball/SlaterDet.h>
#include <math/matrix.h>

class FourierTransform;
class Bisection
{
  private:

    Context ctxt_;

    // bisection levels in each directions
    int nlevels_[3]; // bisection level
    int ndiv_[3];    // number of subdivisions
    int nlevelsmax_; // max level of bisection
    int nst_;        // number of states in SlaterDet
    int nstloc_;

    // real space grid
    int np_[3];            // grid size
    int np01_;             // size of one xy plane
    int np2loc_;           // local z size
    int np012loc_;         // local size
    FourierTransform *ft_;

    std::vector<long int> localization_; // localization indices

    // xy_proj[ir]: projector in xy plane associated with grid point ir
    std::vector<int> xy_proj_;

    // matrices of real space wave functions in subdomains
    std::vector<ComplexMatrix*> rmat_;

    // a matrices
    int nmat_;
    std::vector<ComplexMatrix*> amat_;
    std::vector<std::vector<std::complex<double> > > adiag_; 
    ComplexMatrix *u_;

    // test function
    bool check_amat(const ComplexMatrix &c);
    void trim_amat(const std::vector<double>& occ);

  public:

    Bisection(const SlaterDet& sd, const int nlevels[3]);
    void compute_transform(const SlaterDet& sd);
    void compute_localization(double epsilon);
    void forward(SlaterDet& sd);
    void forward(ComplexMatrix& u, SlaterDet& sd);
    void backward(SlaterDet& sd);
    void backward(ComplexMatrix& u, SlaterDet& sd);

    int nmat(void) const { return nmat_; }
    long int localization(int i) const { return localization_[i]; }
    const std::vector<long int>& localization(void) const
    { return localization_; }
    bool overlap(int i, int j) const;
    bool overlap(const std::vector<long int>& loc, int i, int j) const;
    const ComplexMatrix& u(void) const { return *u_; }
    double pair_fraction(void) const;
    double size(int i) const;
    double total_size(void) const;
    ~Bisection();
};
#endif

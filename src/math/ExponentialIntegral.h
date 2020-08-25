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
// ExponentialIntegral.h
//
////////////////////////////////////////////////////////////////////////////////
//
// Calculate exponential integral using the algorithm of
// Tseng, Lee, Journal of Hydrology, 205 (1998) 38-51
//
////////////////////////////////////////////////////////////////////////////////
#ifndef EXPONENTIAL_INTEGRAL_H
#define EXPONENTIAL_INTEGRAL_H

namespace util {

// constant at which we switch from series to gauss laguerre expansion
static const double series_cutoff = 4.0;

// Calculate the exponential integral E_1(x):
//
//        inf
//          /     -t
//         |     e
// E (x) = | dt -----
//  1      |      t
//        /
//         x
//
// Input:  x - position at which exponential integral is evaluated (x > 0)
// Return: E1(x)
double E1(const double x);

// Series expansion of the exponential integral
//
//                          n_cut
//                          -----     n  n
//                           \    (-1)  x
// E (x) = -gamma - ln(x) -   )   --------
//  1                        /     n * n!
//                          -----
//                          n = 1
//
// where gamma is the Euler constant.
// n_cut is set to 25
// Input:  x - position at which exponential integral is evaluated (x > 0)
// Return: approximation by series expansion for E_1(x)
double E1_series(const double x);

// The Gauss Laguerre expansion of the exponential integral can be written as
//
//             N
// E (x0)    -----     a
//  1         \         n
// ------ =    )   ---------
//   -x0      /     x  + x0
//  e        -----   n
//            n=1
//
// where the a_n and x_n are determined by least quadrature (see reference)
// Input: x0 - point at which Gaussian Laguerre quadrature is calculated
// Return: E_1(x0) / exp(-x0) in this approximation
double gauss_laguerre(const double x0);

}

#endif

//  Copyright (C) 2002 Zoltan Nagy
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include "bits/hhc-process.h"


namespace nlo {

  
  weight_hhc pdf_and_coupling_hhc::
  pdf(double x1, double x2, double mf2, unsigned int nu, unsigned int nd) 
  {
    int ia, iq;
    weight_hhc retval;
    static double __f1[13], __f2[13];
    static const int iu[3] = {2,4,6}, id[3] = {1,3,5};
    
    //----- calculat the pdfs -----
    double *f1 = __f1+6, *f2 = __f2+6;
    
    this -> hadronA(x1, mf2, nu, nd, f1);
    this -> hadronB(x2, mf2, nu, nd, f2);

    //----- gluon pdfs -----
    double A0 = f1[0];
    double B0 = f2[0];
    
    //---- up type quarks -----
    double q1, q2, a1, a2;
    double A = 0.0, B = 0.0, Ab = 0.0, Bb = 0.0, D = 0.0, Db = 0.0; 
    
    for(unsigned int u = 0; u < nu && u < 3; u++) {
      ia = -(iq = iu[u]);
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];
      
      A += q1; Ab += a1; B += q2; Bb += a2;
      D += q1*q2 + a1*a2; Db += q1*a2 + a1*q2;
    }
    
    //----- down type quarks -----
    for(unsigned int d = 0; d < nd && d < 3; d++) {
      ia = -(iq = id[d]);
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];
      
      A += q1; Ab += a1; B += q2; Bb += a2;
      D += q1*q2 + a1*a2; Db += q1*a2 + a1*q2;
    }
    
    retval[0] = A0*B0;
    retval[1] = (A + Ab)*B0;
    retval[2] = A0*(B + Bb);
    retval[3] = A*B + Ab*Bb - D;
    retval[4] = D;
    retval[5] = Db;
    retval[6] = A*Bb +Ab*B - Db;
    
    return retval;
  }
}

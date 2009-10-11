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
#include "bits/hhc2ph-process.h"


namespace nlo {
    
  
  weight_hhc2ph pdf_and_coupling_hhc2ph::
  pdf(double x1, double x2, double mf2, unsigned int nu, unsigned int nd)
  {
    int ia, iq;
    weight_hhc2ph retval;
    static double __f1[13], __f2[13];
    static const int iu[3] = {2,4,6}, id[3] = {1,3,5};
    
    //----- calculat the pdfs -----
    double *f1 = __f1+6, *f2 = __f2+6;
    
    this -> hadronA(x1, mf2, nu, nd, f1);
    this -> hadronB(x2, mf2, nu, nd, f2);
    
    //----- gluon pdfs -----
    double G1 = f1[0];
    double G2 = f2[0];
    
    //---- up type quarks -----
    double q1, q2, a1, a2;
    double U1 = 0.0, U2 = 0.0, Ub1 = 0.0, Ub2 = 0.0, X = 0.0, Xb = 0.0; 
    double D1 = 0.0, D2 = 0.0, Db1 = 0.0, Db2 = 0.0, Y = 0.0, Yb = 0.0; 
    
    for(unsigned int u = 0; u < nu && u < 3; u++) {
      ia = -(iq = iu[u]);
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];

      U1 += q1; Ub1 += a1; U2 += q2; Ub2 += a2;
      X += q1*q2 + a1*a2; Xb += q1*a2 + a1*q2;
    }
    
    //----- down type quarks -----
    for(unsigned int d = 0; d < nd && d < 3; d++) {
      ia = -(iq = id[d]);
      q1 = f1[iq]; q2 = f2[iq];
      a1 = f1[ia]; a2 = f2[ia];
      
      D1 += q1; Db1 += a1; D2 += q2; Db2 += a2;
      Y += q1*q2 + a1*a2; Yb += q1*a2 + a1*q2;
    }
        
    retval[0]  = G1*G2;
    retval[1]  = (U1 + Ub1)*G2;
    retval[2]  = (D1 + Db1)*G2;
    retval[3]  = G1*(U2 + Ub2);
    retval[4]  = G1*(D2 + Db2);
    retval[5]  = U1*U2 + Ub1*Ub2 - X;
    retval[6]  = X;
    retval[7]  = D1*D2 + Db1*Db2 - Y;
    retval[8]  = Y;
    retval[9]  = U1*D2 + Ub1*Db2;
    retval[10] = D1*U2 + Db1*Ub2;
    retval[11] = U1*Ub2 + Ub1*U2 - Xb;
    retval[12] = Xb;
    retval[13] = D1*Db2 + Db1*D2 - Yb;
    retval[14] = Yb;
    retval[15] = U1*Db2 + Ub1*D2;
    retval[16] = D1*Ub2 + Db1*U2;
    
    return retval;
  }
}

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

//   nlo includes
#include "bits/nlo-color.h"
#include "ampg6.h"
#include "defmacros.h"


namespace nlo {
  
#define ICPLX std::complex<double>(0,1)

  const short ampg6::perm[4][24] = 
  {{0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3},
   {1,1,2,2,3,3,2,2,3,3,0,0,3,3,0,0,1,1,0,0,1,1,2,2},
   {2,3,3,1,1,2,3,0,0,2,2,3,0,1,1,3,3,0,1,2,2,0,0,1},
   {3,2,1,3,2,1,0,3,2,0,3,2,1,0,3,1,0,3,2,1,0,2,1,0}};
  
  const short ampg6::hel[10][24] = 
  {{1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4},
   {2,2,3,3,4,4,3,3,4,4,1,1,4,4,1,1,2,2,1,1,2,2,3,3},
   {3,4,4,2,2,3,4,1,1,3,3,4,1,2,2,4,4,1,2,3,3,1,1,2},
   {4,3,2,4,3,2,1,4,3,1,4,3,2,1,4,2,1,4,3,2,1,3,2,1},
   {5,5,9,9,7,7,8,8,10,10,5,5,6,6,9,9,8,8,7,7,10,10,6,6},
   {6,6,10,10,8,8,7,7,9,9,6,6,5,5,10,10,7,7,8,8,9,9,5,5},
   {7,9,5,7,9,5,5,10,8,5,10,8,8,9,6,8,9,6,6,10,7,6,10,7},
   {8,10,6,8,10,6,6,9,7,6,9,7,7,10,5,7,10,5,5,9,8,5,9,8},
   {9,7,7,5,5,9,10,5,5,8,8,10,9,8,8,6,6,9,10,6,6,7,7,10},
   {10,8,8,6,6,10,9,6,6,7,7,9,10,7,7,5,5,10,9,5,5,8,8,9}};
  
  
  std::complex<double> 
  ampg6::APPPMMM(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    std::complex<double> BB, GG;
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s45 = S(4,5), 
      s56 = S(5,6), t123 = s12+s23+S(1,3), t234 = s23+s34+S(2,4),
      t345 = s34+s45+S(3,5);
  
    BB = B(2,3)*A(5,6)*(B(1,2)*A(2,4)+B(1,3)*A(3,4));
    GG = B(1,2)*A(4,5)*(B(3,1)*A(1,6)+B(3,2)*A(2,6));
    
    return ICPLX*(BB*BB*s12*s45/t234+GG*GG*s23*s56/t345+t123*BB*GG);
  }

  std::complex<double>
  ampg6::APPMPMM(int p1, int p2, int p3, int p4, int p5, int p6)  const
  {
    std::complex<double> AA, BB, GG;
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s45 = S(4,5), 
      s56 = S(5,6), s61 = S(6,1), t123 = s12+s23+S(1,3), 
      t234 = s23+s34+S(2,4), t345 = s34+s45+S(3,5);
    
    AA = B(2,1)*A(5,6)*(B(4,1)*A(1,3)+B(4,2)*A(2,3));
    BB = B(2,4)*A(5,6)*(B(1,2)*A(2,3)+B(1,4)*A(4,3));
    GG = B(1,2)*A(3,5)*(B(4,1)*A(1,6)+B(4,2)*A(2,6));
    
    return ICPLX*(AA*AA*s34*s61/t123+BB*BB*s12*s45/t234+GG*GG*s23*s56/t345
		  + t123*BB*GG+t234*GG*AA+t345*AA*BB);
  }

  std::complex<double>
  ampg6::APMPMPM(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    std::complex<double> AA, BB, GG;
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s45 = S(4,5), 
      s56 = S(5,6), s61 = S(6,1), t123 = s12+s23+S(1,3), 
      t234 = s23+s34+S(2,4), t345 = s34+s45+S(3,5);
    
    AA = B(1,3)*A(4,6)*(B(5,1)*A(1,2)+B(5,3)*A(3,2));
    BB = B(5,1)*A(2,4)*(B(3,1)*A(1,6)+B(3,5)*A(5,6));
    GG = B(3,5)*A(6,2)*(B(1,3)*A(3,4)+B(1,5)*A(5,4));
    
    return ICPLX*(AA*AA*s34*s61/t123+BB*BB*s12*s45/t234+GG*GG*s23*s56/t345
		  + t123*BB*GG+t234*GG*AA+t345*AA*BB);
  }
  
  double ampg6::Leading(std::complex<double> *cbase)
  {
    short i;
    double ret_val = 0.0;
    static std::complex<double> c[36];
    
    // Parameter adjustments
    --cbase;
    
    c[0]  = -cbase[24] - cbase[14] - cbase[17] - cbase[18];
    c[1]  = -cbase[14] - cbase[24] - cbase[21] - cbase[22];
    c[4]  = -cbase[10] - cbase[21] - cbase[24] - cbase[23];
    c[7]  = -cbase[21] - cbase[10] - cbase[ 7] - cbase[ 8];
    c[10] = -cbase[17] - cbase[ 7] - cbase[ 9] - cbase[10];
    c[11] = -cbase[ 7] - cbase[17] - cbase[14] - cbase[13];
    c[12] = -cbase[ 1] - cbase[11] - cbase[ 8] - cbase[ 7];
    c[13] = -cbase[ 2] - cbase[12] - cbase[ 9] - cbase[10];
    c[14] = -cbase[ 3] - cbase[16] - cbase[13] - cbase[14];
    c[15] = -cbase[ 4] - cbase[15] - cbase[17] - cbase[18];
    c[16] = -cbase[ 5] - cbase[19] - cbase[21] - cbase[22];
    c[17] = -cbase[ 6] - cbase[20] - cbase[23] - cbase[24];
    c[19] = -cbase[20] - cbase[ 6] - cbase[ 3] - cbase[ 4];
    c[21] = -cbase[16] - cbase[ 3] - cbase[ 5] - cbase[ 6];
    c[22] = -cbase[11] - cbase[ 1] - cbase[ 3] - cbase[ 4];
    c[23] = -cbase[12] - cbase[ 2] - cbase[ 5] - cbase[ 6];
    c[24] = -cbase[ 8] - cbase[18] - cbase[15] - cbase[16];
    c[25] = -cbase[ 9] - cbase[22] - cbase[19] - cbase[20];
    c[26] = -cbase[23] - cbase[13] - cbase[15] - cbase[16];
    c[27] = -cbase[13] - cbase[23] - cbase[19] - cbase[20];
    c[30] = -cbase[15] - cbase[ 4] - cbase[ 1] - cbase[ 2];
    c[31] = -cbase[18] - cbase[ 8] - cbase[11] - cbase[12];
    c[32] = -cbase[22] - cbase[ 9] - cbase[11] - cbase[12];
    c[35] = -cbase[19] - cbase[ 5] - cbase[ 2] - cbase[ 1];
    
    c[2]  = -cbase[ 1] - c[12] - c[ 0] - cbase[24];
    c[3]  = -cbase[ 2] - c[13] - c[ 1] - cbase[14];
    c[5]  = -cbase[ 4] - c[15] - c[ 7] - cbase[21];
    c[6]  = -cbase[ 3] - c[14] - c[ 4] - cbase[10];
    c[8]  = -cbase[ 5] - c[16] - c[10] - cbase[17];
    c[9]  = -cbase[ 6] - c[17] - c[11] - cbase[ 7];
    c[18] = -cbase[ 8] - c[24] - c[19] - cbase[20];
    c[20] = -cbase[ 9] - c[25] - c[21] - cbase[16];
    c[28] = -cbase[11] - c[22] - c[26] - cbase[23];
    c[29] = -cbase[12] - c[23] - c[27] - cbase[13];
    c[33] = -cbase[15] - c[30] - c[32] - cbase[22];
    c[34] = -cbase[18] - c[31] - c[35] - cbase[19];
    
    for(i = 1; i <= 24; i++) ret_val += real(cbase[i]*conj(cbase[i]));
    for(i = 0; i < 36; i++)  ret_val += real(c[i]*conj(c[i]));
    
    return ret_val;
  }
  
  double ampg6::SubLeading(std::complex<double> *cbase)
  {
    short i;
    double ret_val = 0.0;
    static std::complex<double> z[18];
    
    // Parameter adjustments
    --cbase;
    
    z[ 0] = cbase[13] + cbase[14] + cbase[21] + cbase[23] + cbase[24];
    z[ 1] = cbase[13] + cbase[14] + cbase[17] + cbase[23] + cbase[24];
    z[ 2] = cbase[ 7] + cbase[ 9] + cbase[10] + cbase[21] + cbase[22];
    z[ 3] = cbase[ 9] + cbase[10] + cbase[21] + cbase[22] + cbase[24];
    z[ 4] = cbase[ 7] + cbase[ 8] + cbase[14] + cbase[17] + cbase[18];
    z[ 5] = cbase[ 7] + cbase[ 8] + cbase[10] + cbase[17] + cbase[18];
    z[ 6] = cbase[19] + cbase[20];
    z[ 7] = cbase[19] + cbase[20] + cbase[23];
    z[ 8] = cbase[13] + cbase[15] + cbase[16];
    z[ 9] = cbase[15] + cbase[16];
    z[10] = cbase[13] + cbase[14] + cbase[20] + cbase[23] + cbase[24];
    z[11] = cbase[13] + cbase[14] + cbase[16] + cbase[23] + cbase[24];
    z[14] = cbase[19] + cbase[21] + cbase[22];
    z[15] = cbase[21] + cbase[22];
    z[16] = cbase[19] + cbase[20];
    z[17] = cbase[19] + cbase[20] + cbase[22];
    
    for (i = 1;  i <= 12; i++) ret_val += real(cbase[i]*conj(z[i-1]));
    for (i = 15; i <= 18; i++) ret_val += real(cbase[i]*conj(z[i-1]));
    
    return ret_val;
  } 
  
  double ampg6::su3_tree(int p1, int p2, int q3, int q4, int q5, int q6)
  {
    static const double coll = 4.0*Na*Nc2*Nc2; 
    static const double cols = 4.0*24.0*Na*Nc2;
    static std::complex<double> m[11][24];

    double ret_val, ptnum, tmp, prop;
    int i, p3, p4, p5, p6, p[4]={q3, q4, q5, q6};
    std::complex<double> temp;
    
    for(i = 0; i < 24; i++) {
      p3 = p[perm[0][i]];
      p4 = p[perm[1][i]];
      p5 = p[perm[2][i]];
      p6 = p[perm[3][i]];
      
      temp = 1.0/(A(1,2)*A(2,3)*A(3,4)*A(4,5)*A(5,6)*A(6,1));
      prop = real(temp*conj(temp));
      
      m[0][i] = temp;
      
      m[hel[0][i]][i] = prop*APPPMMM(p1,p2,p3,p4,p5,p6);
      m[hel[1][i]][i] = prop*APPMPMM(p1,p2,p3,p4,p5,p6);
      m[hel[2][i]][i] = prop*APPMMPM(p1,p2,p3,p4,p5,p6);
      m[hel[3][i]][i] = prop*APPMMMP(p1,p2,p3,p4,p5,p6);
      	
      m[hel[4][i]][i] = prop*APMPPMM(p1,p2,p3,p4,p5,p6);
      m[hel[5][i]][i] = prop*APMMMPP(p1,p2,p3,p4,p5,p6);
      m[hel[6][i]][i] = prop*APMPMMP(p1,p2,p3,p4,p5,p6);
      m[hel[7][i]][i] = prop*APMMPPM(p1,p2,p3,p4,p5,p6);
      	
      m[hel[8][i]][i] = prop*APMPMPM(p1,p2,p3,p4,p5,p6);
      m[hel[9][i]][i] = prop*APMMPMP(p1,p2,p3,p4,p5,p6);
    }
    
    tmp = Sij(p1,p2), tmp *= tmp; ptnum  = tmp*tmp;
    tmp = Sij(p1,q3), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(p1,q4), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(p1,q5), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(p1,q6), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(p2,q3), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(p2,q4), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(p2,q5), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(p2,q6), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(q3,q4), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(q3,q5), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(q3,q6), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(q4,q5), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(q4,q6), tmp *= tmp; ptnum += tmp*tmp;
    tmp = Sij(q5,q6), tmp *= tmp; ptnum += tmp*tmp;
    
    ret_val = ptnum*(coll*Leading(m[0]) + cols*SubLeading(m[0]));
    
    for(i = 1; i <= 10; i++)
      ret_val += coll*Leading(m[i]) + cols*SubLeading(m[i]);
    
    return ret_val;
  } 
}   //  namespace nlo

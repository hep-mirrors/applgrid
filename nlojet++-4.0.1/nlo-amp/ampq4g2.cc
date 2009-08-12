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
#include "ampq4g2.h"
#include "defmacros.h"


namespace nlo {

  const short ampq4g2::_S_itab[4][8] 
  = {{0,1,2,3,4,5,6,7},{0,1,2,3,7,6,5,4},
     {0,2,1,3,4,6,5,7},{0,2,1,3,7,5,6,4}};
  
  struct ampq4g2::_Colmat {
    _Colmat();
    double colmat[12][12];
  };
  
  ampq4g2::_Colmat::_Colmat()
  {
    double c1 = Na2, c3 = Na, c4 = -Na, c5 = Na2/Nc, c6 = -Na/Nc;
    double colfac[6] = {c1, 0.0, c3, c4, c5, c6};
    short i, j, ord[12][12] = 
    {{0,1,2,3,1,2,4,5,4,5,4,5}, {1,0,1,1,2,1,5,4,4,4,4,5},
     {2,1,0,2,1,3,4,4,4,5,5,5}, {3,1,2,0,1,2,5,4,5,4,5,4},
     {1,2,1,1,0,1,4,4,5,5,4,4}, {2,1,3,2,1,0,5,5,5,4,4,4},
     {4,5,4,5,4,5,0,1,2,3,1,2}, {5,4,4,4,4,5,1,0,1,1,2,1},
     {4,4,4,5,5,5,2,1,0,2,1,3}, {5,4,5,4,5,4,3,1,2,0,1,2},
     {4,4,5,5,4,4,1,2,1,1,0,1}, {5,5,5,4,4,4,2,1,3,2,1,0}};
    
    for(i = 0; i < 12; i++)
      for(j = 0; j < 12; j++)
	colmat[i][j] = colfac[ord[i][j]];
  }

  //    the color matrix
  const ampq4g2::_Colmat ampq4g2::_S_colmat;


  //  functions E20,E11,F20,F11 returns 8 helicity combinations
  //  p1,p2,p3,p4;p5,p6  :  0 = +-+-;++     1 =+--+;++       
  //                        2 = -++-;++     3 =-+-+;++       
  //                        4 = +-+-;+-     5 =+--+;+-       
  //                        6 = -++-;+-     7 =-+-+;+-       
  //  ICC determines whether to complex conjugate the last four amplitudes
  //  because the call is with G2,G1 instead of G1,G2. 
  void ampq4g2::E20(int p1, int p2, int p3, int p4, int p5, int p6,
		    const short *itab, _ComplexD *E, double fac, bool icc) const
  {
    std::complex<double> Z1, Z2, Z3, Z4;
    std::complex<double> I(0.0, fac), fac1, fac2, fac3, fac4;
    
    double s15 = S(1,5), s56 = S(5,6), s34 = S(3,4), s12 = S(1,2),  
      s46 = S(4,6), s125 = S3(1,2,5), s156 = S3(1,5,6), s564 = S3(5,6,4);
    
    double P1 = s15*s156*s56*s34, P2 = s12*s56*s564*s46,
      P3 = s12*s125*s15*s46*s34, P4 = s12*s15*s56*s46*s34;
    
    if(icc) swap();
    _ComplexD a13 = A(1,3), a14 = A(1,4), a16 = A(1,6), a23 = A(2,3), 
      a24 = A(2,4), a31 = A(3,1), a32 = A(3,2), a36 = A(3,6), a46 = A(4,6), 
      a56 = A(5,6), a61 = A(6,1), a63 = A(6,3), a64 = A(6,4), a65 = A(6,5),
      b15 = B(1,5), b23 = B(2,3), b25 = B(2,5), b31 = B(3,1), b32 = B(3,2), 
      b41 = B(4,1), b42 = B(4,2), b45 = B(4,5), b51 = B(5,1), b52 = B(5,2),
      b54 = B(5,4), b56 = B(5,6), b65 = B(6,5), c1253 = C(1,2,5,3), 
      c1254 = C(1,2,5,4), c1364 = C(1,3,6,4), c1465 = C(1,4,6,5), 
      c2364 = C(2,3,6,4), c2465 = C(2,4,6,5), c3165 = C(3,1,6,5), 
      c4165 = C(4,1,6,5), c6153 = C(6,1,5,3), c6154 = C(6,1,5,4),
      c6451 = C(6,4,5,1), c6452 = C(6,4,5,2);
    if(icc) swap();
    
    
    //  helicities: ****;++
    fac4 = I*A(1,4)/(A(1,5)*A(5,6)*A(6,4)*A(1,2)*A(3,4));
    
    E[itab[0]] = -fac4*A(2,4)*A(2,4);
    E[itab[1]] =  fac4*A(2,3)*A(2,3);
    E[itab[2]] =  fac4*A(1,4)*A(1,4);
    E[itab[3]] = -fac4*A(1,3)*A(1,3);
    
    // helicities: ****;+-
    fac1 = I/P1; fac2 = I/P2; fac3 = I/P3; fac4 = I/P4;
    
    Z1 = a61*a24*b15*b15*c6153;
    Z2 = a46*a46*b31*b54*c2465;
    Z3 = a46*b15*a46*b15*c1253*c2364;
    Z4 = s125*a46*a61*a24*b15*b31*b54 - a46*a46*a61*b15*b65*b31*c2364 
      -  a46*a56*a24*b15*b15*b54*c1253;
    
    E[itab[4]] = Z1*fac1 + Z2*fac2 + Z3*fac3 + Z4*fac4;
    
    Z1 = a32*a16*b15*b15*c6154;
    Z2 = a23*a46*b54*b54*c6451;
    Z3 = a63*a64*b15*b15*c1254*c2364;
    Z4 = a63*a64*a61*b65*b15*b41*c2364 + a64*a23*a65*b15*b15*b54*c1254
      +   s125*a61*a64*a23*b15*b41*b54;
    
    E[itab[5]] = Z1*fac1 + Z2*fac2 + Z3*fac3 - Z4*fac4;
    
    Z1 = a61*a61*b15*b32*c4165;
    Z2 = a46*a46*b23*b45*c1465;
    Z3 = a46*a46*b51*b52*c1364*c1253;
    Z4 = a61*a46*a46*b56*b32*b51*c1364 + a14*a46*a56*b51*b54*b52*c1253
      +  s125*a46*a14*a61*b51*b54*b32;
    
    E[itab[6]] = Z1*fac1 + Z2*fac2 + Z3*fac3 - Z4*fac4;
    
    Z1 = a61*a61*b42*b15*c3165;
    Z2 = a13*a46*b45*b45*c6452;
    Z3 = a46*a36*b25*b15*c1364*c1254;
    Z4 = a46*a36*a16*b56*b42*b15*c1364 + a46*a31*a56*b15*b45*b25*c1254
      +  s125*a16*a31*a46*b15*b45*b42;
    
    E[itab[7]] = Z1*fac1 + Z2*fac2 + Z3*fac3 - Z4*fac4;
  }

  void ampq4g2::E11(int p1, int p2, int p3, int p4, int p5, int p6,
		    const short *itab, _ComplexD *E, double fac, bool icc) const
  {
    std::complex<double> I(0.0, fac), fac1, fac2, fac3;
    
    double s12 = S(1,2), s15 = S(1,5), s26 = S(2,6), s34 = S(3,4), 
      s36 = S(3,6), s45 = S(4,5), s125 = S3(1,2,5), s345 = S3(3,4,5); 
    
    double P1 = s12*s15*s125*s34*s36, P2 = s12*s45*s345*s34*s26,
      P3 = s12*s15*s45*s34*s36*s26;
    
    if(icc) swap();
    _ComplexD a14 = A(1,4), a26 = A(2,6), a36 = A(3,6), a46 = A(4,6), 
      a61 = A(6,1), a62 = A(6,2), a63 = A(6,3), b15 = B(1,5), b32 = B(3,2),
      b45 = B(4,5), b51 = B(5,1), b52 = B(5,2), b53 = B(5,3), b54 = B(5,4),
      c1253 = C(1,2,5,3), c1254 = C(1,2,5,4), c1364 = C(1,3,6,4), 
      c2153 = C(2,1,5,3), c2154 = C(2,1,5,4), c3162 = C(3,1,6,2), 
      c3261 = C(3,2,6,1), c3451 = C(3,4,5,1), c4162 = C(4,1,6,2), 
      c4261 = C(4,2,6,1), c4351 = C(4,3,5,1), c4352 = C(4,3,5,2);
    if(icc) swap();
    
    //  helicities: ****;++
    fac3 = I*A(1,4)*A(3,2)/(A(1,5)*A(5,4)*A(3,6)*A(6,2)*A(1,2)*A(3,4));
    
    E[itab[0]] = -A(2,4)*A(2,4)*fac3;
    E[itab[1]] =  A(2,3)*A(2,3)*fac3;
    E[itab[2]] =  A(1,4)*A(1,4)*fac3;
    E[itab[3]] = -A(1,3)*A(1,3)*fac3;
    
    //  helicities: ****;+-
    fac1 = I/P1; fac2 = I/P2; fac3 = I*b51*a14*b45*a36*b32*a26/P3;
    
    E[itab[4]] = fac1*a46*a36*b51*b51*c1253*c2153
      +          fac2*a26*a26*b53*b54*c4162*c4261
      +          fac3*c2153*c4351;
    
    E[itab[5]] = fac1*a63*b51*a63*b51*c1253*c2154
      +          fac2*a26*b54*a26*b54*c3261*c4162
      +          fac3*c3451*c2154;
    
    E[itab[6]] = fac1*a46*a63*b52*b15*c1253*c1253
      +          fac2*a61*a62*b54*b53*c4162*c4162
      +          fac3*c1253*c4352;
    
    E[itab[7]] = fac1*a36*a36*b51*b52*c1253*c1254
      +          fac2*a62*a61*b54*b54*c3162*c4162
      +          fac3*c1364*c3162;
  }


  void ampq4g2::F20(int p1, int p2, int p3, int p4, int p5, int p6,
		    const short *itab, _ComplexD *F, double fac, bool icc) const
  {
    std::complex<double> I(0.0, fac), fac1, fac2, fac3;

    double s56 = S(5,6), s34 = S(3,4), s256 = S3(2,5,6), s156 = S3(1,5,6);
 
    if(icc) swap();
    _ComplexD a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), a31 = A(3,1), 
      a32 = A(3,2), a41 = A(4,1), a42 = A(4,2), a61 = A(6,1), a62 = A(6,2), 
      b13 = B(1,3), b14 = B(1,4), b15 = B(1,5), b23 = B(2,3), b24 = B(2,4), 
      b25 = B(2,5), b32 = B(3,2), b42 = B(4,2), b62 = B(6,2), 
      c3165 = C(3,1,6,5), c3265 = C(3,2,6,5), c4165 = C(4,1,6,5), 
      c4265 = C(4,2,6,5), c6153 = C(6,1,5,3), c6154 = C(6,1,5,4), 
      c6253 = C(6,2,5,3), c6254 = C(6,2,5,4);
    if(icc) swap();
        
    //  helicities: ****;++
    fac3 = I/(A(1,5)*A(5,6)*A(6,2)*A(3,4));
    F[itab[0]] = -fac3*A(2,4)*A(2,4);
    F[itab[1]] =  fac3*A(2,3)*A(2,3);
    F[itab[2]] =  fac3*A(1,4)*A(1,4);
    F[itab[3]] = -fac3*A(1,3)*A(1,3);
    
    //  helicities: ****;+-
    fac1 = I/(s34*s56*s156*a15);
    fac2 = I/(s34*s56*a15*b62);
    fac3 = I/(s34*s56*s256*b62);
     
    F[itab[4]] = fac1*a61*b15*c6153*a42 + fac2*c6153*c4265 
      +          fac3*b13*c4265*a62*b25;

    F[itab[5]] = fac1*a61*b15*c6154*a32 + fac2*c6154*c3265 
      +          fac3*b14*c3265*a62*b25;

    F[itab[6]] = fac1*a61*a61*b23*c4165 + fac2*a61*a14*b32*b25
      +          fac3*c6253*a41*b25*b25;

    F[itab[7]] = fac1*a61*a61*b24*c3165 + fac2*a61*a13*b42*b25
      +          fac3*c6254*a31*b25*b25;
  }

  void ampq4g2::F11(int p1, int p2, int p3, int p4, int p5, int p6,
		    const short *itab, _ComplexD *F, double fac, bool icc) const
  {
    std::complex<double> I(0.0, fac);

    //  helicities: ****;++
    _ComplexD fac2 = I/(A(1,5)*A(5,2)*A(3,6)*A(6,4));

    F[itab[0]] = -fac2*A(2,4)*A(2,4);
    F[itab[1]] =  fac2*A(2,3)*A(2,3);
    F[itab[2]] =  fac2*A(1,4)*A(1,4);
    F[itab[3]] = -fac2*A(1,3)*A(1,3);

    //  helicities: ****;+-
    if(icc) swap();
    _ComplexD fac1 = I/(A(1,5)*A(5,2)*B(3,6)*B(6,4)*S3(3,4,6));
    _ComplexD c1364 = C(1,3,6,4), c1463 = C(1,4,6,3),
      c2364 = C(2,3,6,4), c2463 = C(2,4,6,3); 

    F[itab[4]] = -fac1*c2463*c2463;
    F[itab[5]] =  fac1*c2364*c2364;
    F[itab[6]] =  fac1*c1463*c1463;
    F[itab[7]] = -fac1*c1364*c1364;
    if(icc) swap();
  }

  void ampq4g2::
  __su3_tree(int p1, int p2, int p3, int p4, int p5, int p6, double *out) const
  {
    static double coll = 2.0, cols = -2.0/Nc;
    static std::complex<double> m12[12][8], m14[12][8];

    F20(p3,p4,p1,p2,p5,p6, _S_itab[2], m12[0], cols, false);
    F11(p1,p2,p3,p4,p5,p6, _S_itab[0], m12[1], cols, false);
    F20(p1,p2,p3,p4,p5,p6, _S_itab[0], m12[2], cols, false);
    					  
    F20(p3,p4,p1,p2,p6,p5, _S_itab[3], m12[3], cols, true);
    F11(p1,p2,p3,p4,p6,p5, _S_itab[1], m12[4], cols, true);
    F20(p1,p2,p3,p4,p6,p5, _S_itab[1], m12[5], cols, true);
    					  
    E20(p3,p4,p1,p2,p5,p6, _S_itab[2], m12[6], coll, false);
    E11(p1,p2,p3,p4,p5,p6, _S_itab[0], m12[7], coll, false);
    E20(p1,p2,p3,p4,p5,p6, _S_itab[0], m12[8], coll, false);
    
    E20(p3,p4,p1,p2,p6,p5, _S_itab[3], m12[9],  coll, true);
    E11(p1,p2,p3,p4,p6,p5, _S_itab[1], m12[10], coll, true);
    E20(p1,p2,p3,p4,p6,p5, _S_itab[1], m12[11], coll, true);
    
    //----- same calls with quarks interchanged -----
    F20(p3,p2,p1,p4,p5,p6, _S_itab[2], m14[6], cols, false);
    F11(p1,p4,p3,p2,p5,p6, _S_itab[0], m14[7], cols, false);
    F20(p1,p4,p3,p2,p5,p6, _S_itab[0], m14[8], cols, false);
    
    F20(p3,p2,p1,p4,p6,p5, _S_itab[3], m14[9],  cols, true);
    F11(p1,p4,p3,p2,p6,p5, _S_itab[1], m14[10], cols, true);
    F20(p1,p4,p3,p2,p6,p5, _S_itab[1], m14[11], cols, true);
    	       	        	       	
    E20(p3,p2,p1,p4,p5,p6, _S_itab[2], m14[0], coll, false);
    E11(p1,p4,p3,p2,p5,p6, _S_itab[0], m14[1], coll, false);
    E20(p1,p4,p3,p2,p5,p6, _S_itab[0], m14[2], coll, false);
    	       	        	       		
    E20(p3,p2,p1,p4,p6,p5, _S_itab[3], m14[3], coll, true);
    E11(p1,p4,p3,p2,p6,p5, _S_itab[1], m14[4], coll, true);
    E20(p1,p4,p3,p2,p6,p5, _S_itab[1], m14[5], coll, true);

    //----- Color and helicity sums -----
    double colfac;
    unsigned int i, j;
    _ComplexD tmp12, tmp14, tmp0, tmp3, tmp4, tmp7;
    
    out[0] = out[1] = out[2] = 0.0; 
    for(int h = 0; h < 8; h++)
      for(i = 0; i < 12; i++) {
	colfac = 0.5*_S_colmat.colmat[i][i];
	tmp12 = colfac*m12[i][h];
	tmp14 = colfac*m14[i][h];
	for(j = i + 1; j < 12; j++) {
	  colfac = _S_colmat.colmat[i][j];
	  tmp12 += colfac*m12[j][h];
	  tmp14 += colfac*m14[j][h];
	}
	out[0] += real(tmp12*conj(m12[i][h]));
	out[1] += real(tmp14*conj(m14[i][h]));
      }
    
    for(i = 0; i < 12; i++) {
      tmp0 = tmp3 = tmp4 = tmp7 = 0.0;
      for(j = 0; j < 12; j++) {
	colfac = _S_colmat.colmat[i][j];
	tmp0 += colfac*m14[j][0];
	tmp3 += colfac*m14[j][3];
	tmp4 += colfac*m14[j][4];
	tmp7 += colfac*m14[j][7];
      }
      
      out[2] -= real(tmp0*conj(m12[i][0])) + real(tmp3*conj(m12[i][3]))
	+       real(tmp4*conj(m12[i][4])) + real(tmp7*conj(m12[i][7]));
    }
  }
}  //  namespace nlo

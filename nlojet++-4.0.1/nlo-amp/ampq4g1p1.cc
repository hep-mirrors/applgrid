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
#include "ampq4g1p1.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {

  
  std::complex<double>
  ampq4g1p1::g1(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    double s34 = S(3,4), s134 = S3(1,3,4), s156 = S3(1,5,6);
    _ComplexD a16 = A(1,6), a24 = A(2,4), a52 = A(5,2), a62 = A(6,2), 
      b15 = B(1,5), b16 = B(1,6), b31 = B(3,1), b52 = B(5,2), 
      d1364 = D(1,3,6,4), d3162 = D(3,1,6,2), d3165 = D(3,1,6,5),
      d6254 = D(6,2,5,4); 

    return ICPLX*(d1364*d3162/(a16*a62*b15*b52*s34)
		  + b16*a24*d3165/(a16*b15*s34*s156)
		  + a52*b31*d6254/(b52*a62*s34*s134));
  }
  
  std::complex<double>
  ampq4g1p1::g2(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    double s126 = S3(1,2,6);
    _ComplexD a16 = A(1,6), a62 = A(6,2), b35 = B(3,5), b54 = B(5,4), 
      d3162 = D(3,1,6,2); 
    
    return -ICPLX*d3162*d3162/(a16*a62*b35*b54*s126);
  }

  std::complex<double>
  ampq4g1p1::f1(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    double s34 = S(3,4), s126 = S3(1,2,6), s134 = S3(1,3,4);
    
    _ComplexD a16 = A(1,6), a45 = A(4,5), a52 = A(5,2),
      a62 = A(6,2), b31 = B(3,1), b35 = B(3,5), b52 = B(5,2), 
      d3162 = D(3,1,6,2), d3164 = D(3,1,6,4), d6254 = D(6,2,5,4); 
    
    return ICPLX*(d3162*d3164/(a16*a62*b35*b52*s34)
		  + a52*b31*d6254/(b52*a62*s34*s134)
		  + a45*d3162*d3162/(a16*a62*b35*s34*s126));
  }
  
  std::complex<double>
  ampq4g1p1::f2(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    double s34 = S(3,4), s126 = S3(1,2,6), s156 = S3(1,5,6);

    _ComplexD a16 = A(1,6), a24 = A(2,4), a53 = A(5,3), a62 = A(6,2), 
      b15 = B(1,5), b16 = B(1,6), b54 = B(5,4), d1452 = D(1,4,5,2), 
      d3162 = D(3,1,6,2), d3165 = D(3,1,6,5);
    
    return ICPLX*(d1452*d3162/(a16*a62*b15*b54*s34)
		  + b16*a24*d3165/(a16*b15*s34*s156)
		  + a53*d3162*d3162/(a16*a62*b54*s34*s126));
  }

  std::complex<double>
  ampq4g1p1::f3(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    double s34 = S(3,4), s126 = S3(1,2,6), s134 = S3(1,3,4), s146 = S3(1,4,6); 
    
    _ComplexD a16 = A(1,6), a52 = A(5,2), a54 = A(5,4), a62 = A(6,2), 
      b35 = B(3,5), b41 = B(4,1), b52 = B(5,2), d4162 = D(4,1,6,2), 
      d6253 = D(6,2,5,3); 

    return ICPLX*(-s146*d4162/(a16*a62*b35*b52*s34)
		  + a54*d4162*d4162/(a16*a62*b35*s34*s126)
		  + a52*b41*d6253/(b52*a62*s34*s134));
  }
  
  std::complex<double>
  ampq4g1p1::f4(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    double s34 = S(3,4), s126 = S3(1,2,6), s156 = S3(1,5,6);
    
    _ComplexD a16 = A(1,6), a23 = A(2,3), a35 = A(3,5), a62 = A(6,2), 
      b15 = B(1,5), b16 = B(1,6), b41 = B(4,1), b54 = B(5,4), 
      d4162 = D(4,1,6,2), d4165 = D(4,1,6,5);
    
    return ICPLX*(b41*a23*d4162/(a16*a62*b15*b54*s34)
		  + b16*a23*d4165/(a16*b15*s34*s156)
		  + a35*d4162*d4162/(a16*a62*b54*s34*s126));
  }


#define F(i,j) (A(i,j)/(A(i,5)*A(5,j)))

  void ampq4g1p1::
  Amhv(double Q1, double Q2, const _ComplexD& fhh, int p1, int p2, 
       int p3, int p4, int p5, int p6, _ComplexD *E) const
  {
    _ComplexD a21 = A(2,1), a43 = A(4,3);
    _ComplexD Fhh = ICPLX*fhh/(a21*a43)
      *(Q1*a21/(A(2,6)*A(6,1)) + Q2*a43/(A(4,6)*A(6,3)));
   
    E[0] = Fhh*F(2,3);  // A01
    E[1] = Fhh*F(4,1);  // A10
    E[2] = Fhh*F(2,1);  // B10
    E[3] = Fhh*F(4,3);  // B01
  }

  void ampq4g1p1::
  matrix_tree_pmpmmp(double Q1, double Q2, int p1, int p2, int p3, 
		     int p4, int p5, int p6, _ComplexD *E) const 
  {
    E[0] = Q1*f1(p1,p2, p3,p4, p5, p6) + Q2*f2(p3,p4, p1,p2, p5, p6);
    E[1] = Q1*f2(p1,p2, p3,p4, p5, p6) + Q2*f1(p3,p4, p1,p2, p5, p6);
    E[2] = Q1*g1(p1,p2, p3,p4, p5, p6) + Q2*g2(p3,p4, p1,p2, p5, p6);
    E[3] = Q1*g2(p1,p2, p3,p4, p5, p6) + Q2*g1(p3,p4, p1,p2, p5, p6);
  }

  void ampq4g1p1::
  matrix_tree_mpmpmp(double Q1, double Q2, int p1, int p2, int p3, 
		     int p4, int p5, int p6, _ComplexD *E) const 
  {
    E[1] = Q1*f1(p2,p1, p4,p3, p5, p6) + Q2*f2(p4,p3, p2,p1, p5, p6);
    E[0] = Q1*f2(p2,p1, p4,p3, p5, p6) + Q2*f1(p4,p3, p2,p1, p5, p6);
    E[2] = Q1*g1(p2,p1, p4,p3, p5, p6) + Q2*g2(p4,p3, p2,p1, p5, p6);
    E[3] = Q1*g2(p2,p1, p4,p3, p5, p6) + Q2*g1(p4,p3, p2,p1, p5, p6);
  }

  void ampq4g1p1::
  matrix_tree_pmmpmp(double Q1, double Q2, int p1, int p2, int p3, 
		     int p4, int p5, int p6, _ComplexD *E) const 
  {
    E[0] =  Q1*f3(p1,p2, p3,p4, p5, p6) + Q2*f3(p4,p3, p2,p1, p5, p6);
    E[1] =  Q1*f4(p1,p2, p3,p4, p5, p6) + Q2*f4(p4,p3, p2,p1, p5, p6);
    E[2] =  Q1*g1(p1,p2, p4,p3, p5, p6) - Q2*g2(p4,p3, p1,p2, p5, p6);
    E[3] = -Q1*g2(p1,p2, p4,p3, p5, p6) + Q2*g1(p4,p3, p1,p2, p5, p6);
  }

  void ampq4g1p1::
  matrix_tree_mppmmp(double Q1, double Q2, int p1, int p2, int p3, 
		     int p4, int p5, int p6, _ComplexD *E) const 
  {
    E[1] =  Q1*f3(p2,p1, p4,p3, p5, p6) + Q2*f3(p3,p4, p1,p2, p5, p6);
    E[0] =  Q1*f4(p2,p1, p4,p3, p5, p6) + Q2*f4(p3,p4, p1,p2, p5, p6);
    E[2] =  Q1*g1(p2,p1, p3,p4, p5, p6) - Q2*g2(p3,p4, p2,p1, p5, p6);
    E[3] = -Q1*g2(p2,p1, p3,p4, p5, p6) + Q2*g1(p3,p4, p2,p1, p5, p6);
  }

  double ampq4g1p1::amptree(_ComplexD *E) const
  {
    return 2.0*Na*(Nc2*(real(E[0]*conj(E[0])) + real(E[1]*conj(E[1])))
		   + real(E[2]*conj(E[2])) + real(E[3]*conj(E[3]))
		   - 2.0*real((E[0]+E[1])*conj(E[2]+E[3])))/Nc;
  }
  
  
  void ampq4g1p1::
  su3_tree(double Q, int p1, int p2, int p3,
	   int p4, int p5, int p6, double *out) const
  {
    double tmp;
    _ComplexD e14[4], E[8][4], F[8][4];
    out[0] = out[1] = 0.0;
    
    matrix_tree_pmpmpp(Q, Q, p1, p2, p3, p4, p5, p6, E[0]);
    matrix_tree_pmpmpp(Q, Q, p1, p4, p3, p2, p5, p6, F[0]);

    matrix_tree_mpmppp(Q, Q, p1, p2, p3, p4, p5, p6, E[1]);
    matrix_tree_mpmppp(Q, Q, p1, p4, p3, p2, p5, p6, F[1]);   

    matrix_tree_pmpmmp(Q, Q, p1, p2, p3, p4, p5, p6, E[2]);
    matrix_tree_pmpmmp(Q, Q, p1, p4, p3, p2, p5, p6, F[2]);

    matrix_tree_mpmpmp(Q, Q, p1, p2, p3, p4, p5, p6, E[3]);
    matrix_tree_mpmpmp(Q, Q, p1, p4, p3, p2, p5, p6, F[3]);


    for(unsigned i = 0; i < 4; i++) {
      e14[0] = E[i][0] + F[i][3]/Nc; 
      e14[1] = E[i][1] + F[i][2]/Nc; 
      e14[2] = E[i][2] + F[i][1]*Nc; 
      e14[3] = E[i][3] + F[i][0]*Nc;
     
      out[0] += amptree(E[i]); 
      out[1] += amptree(e14);
    }
    
    matrix_tree_pmmppp(Q, Q, p1, p2, p3, p4, p5, p6, E[4]);
    matrix_tree_pmmppp(Q, Q, p1, p4, p3, p2, p5, p6, F[4]);
    
    matrix_tree_mppmpp(Q, Q, p1, p2, p3, p4, p5, p6, E[5]);
    matrix_tree_mppmpp(Q, Q, p1, p4, p3, p2, p5, p6, F[5]);

    matrix_tree_pmmpmp(Q, Q, p1, p2, p3, p4, p5, p6, E[6]);
    matrix_tree_pmmpmp(Q, Q, p1, p4, p3, p2, p5, p6, F[6]);
    
    matrix_tree_mppmmp(Q, Q, p1, p2, p3, p4, p5, p6, E[7]);
    matrix_tree_mppmmp(Q, Q, p1, p4, p3, p2, p5, p6, F[7]);
    
    for(unsigned i = 4; i < 8; i++) {
      e14[0] = F[i][3]/Nc; 
      e14[1] = F[i][2]/Nc; 
      e14[2] = F[i][1]*Nc; 
      e14[3] = F[i][0]*Nc;
      
      out[0] += (tmp = amptree(E[i])); 
      out[1] += tmp + amptree(e14); 
    }
    
    out[0] *= 2.0; 
    out[1] *= 2.0;
  }
  
  
  double ampq4g1p1::
  su3_tree(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    double tmp = 0.0;
    _ComplexD E[8][4];
    
    matrix_tree_pmpmpp(Q1, Q2, p1, p2, p3, p4, p5, p6, E[0]);
    matrix_tree_mpmppp(Q1, Q2, p1, p2, p3, p4, p5, p6, E[1]);
    matrix_tree_pmmppp(Q1, Q2, p1, p2, p3, p4, p5, p6, E[2]);
    matrix_tree_mppmpp(Q1, Q2, p1, p2, p3, p4, p5, p6, E[3]);

    matrix_tree_pmpmmp(Q1, Q2, p1, p2, p3, p4, p5, p6, E[4]);
    matrix_tree_mpmpmp(Q1, Q2, p1, p2, p3, p4, p5, p6, E[5]);
    matrix_tree_pmmpmp(Q1, Q2, p1, p2, p3, p4, p5, p6, E[6]);
    matrix_tree_mppmmp(Q1, Q2, p1, p2, p3, p4, p5, p6, E[7]);
        
    for(unsigned i = 0; i < 8; i++)
      tmp += amptree(E[i]); 
    
    return 2.0*tmp;
  }
  
  void ampq4g1p1::
  su3_tree_mch(double Q, int p1, int p2, int p3,
	       int p4, int p5, int p6, double *out) const
  {
    double tmp;
    _ComplexD E[4], F[4], e14[4];
    unsigned int hel = (unsigned int) (8*_M_rng());
   
    switch(hel) {
    case 0:
      matrix_tree_pmpmpp(Q, Q, p1, p2, p3, p4, p5, p6, E);
      matrix_tree_pmpmpp(Q, Q, p1, p4, p3, p2, p5, p6, F);
      break;
    case 1:
      matrix_tree_mpmppp(Q, Q, p1, p2, p3, p4, p5, p6, E);
      matrix_tree_mpmppp(Q, Q, p1, p4, p3, p2, p5, p6, F);
      break;
    case 2:
      matrix_tree_pmpmmp(Q, Q, p1, p2, p3, p4, p5, p6, E);
      matrix_tree_pmpmmp(Q, Q, p1, p4, p3, p2, p5, p6, F);
      break;
    case 3:
      matrix_tree_mpmpmp(Q, Q, p1, p2, p3, p4, p5, p6, E);
      matrix_tree_mpmpmp(Q, Q, p1, p4, p3, p2, p5, p6, F);
      break;
    case 4:
      matrix_tree_pmmppp(Q, Q, p1, p2, p3, p4, p5, p6, E);
      matrix_tree_pmmppp(Q, Q, p1, p4, p3, p2, p5, p6, F);
      break;
    case 5:
      matrix_tree_mppmpp(Q, Q, p1, p2, p3, p4, p5, p6, E);
      matrix_tree_mppmpp(Q, Q, p1, p4, p3, p2, p5, p6, F);
      break;
    case 6:
      matrix_tree_pmmpmp(Q, Q, p1, p2, p3, p4, p5, p6, E);
      matrix_tree_pmmpmp(Q, Q, p1, p4, p3, p2, p5, p6, F);
      break;
    case 7:
      matrix_tree_mppmmp(Q, Q, p1, p2, p3, p4, p5, p6, E);
      matrix_tree_mppmmp(Q, Q, p1, p4, p3, p2, p5, p6, F);
      break;
    }
    
    if(hel < 4) {
      e14[0] = E[0] + F[3]/Nc; 
      e14[1] = E[1] + F[2]/Nc; 
      e14[2] = E[2] + F[1]*Nc; 
      e14[3] = E[3] + F[0]*Nc;
      
      out[0] = 16.0*amptree(E); 
      out[1] = 16.0*amptree(e14); 
    } else {
      e14[0] = F[3]/Nc; 
      e14[1] = F[2]/Nc; 
      e14[2] = F[1]*Nc; 
      e14[3] = F[0]*Nc;
      
      out[0] = 16.0*(tmp = amptree(E)); 
      out[1] = 16.0*(tmp + amptree(e14)); 
    }
  }
  
  double ampq4g1p1::
  su3_tree_mch(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    _ComplexD E[4];
    unsigned int hel = (unsigned int) (8*_M_rng());
   
    switch(hel) {
    case 0: matrix_tree_pmpmpp(Q1, Q2, p1, p2, p3, p4, p5, p6, E); break;
    case 1: matrix_tree_mpmppp(Q1, Q2, p1, p2, p3, p4, p5, p6, E); break;
    case 2: matrix_tree_pmpmmp(Q1, Q2, p1, p2, p3, p4, p5, p6, E); break;
    case 3: matrix_tree_mpmpmp(Q1, Q2, p1, p2, p3, p4, p5, p6, E); break;
    case 4: matrix_tree_pmmppp(Q1, Q2, p1, p2, p3, p4, p5, p6, E); break;
    case 5: matrix_tree_mppmpp(Q1, Q2, p1, p2, p3, p4, p5, p6, E); break;
    case 6: matrix_tree_pmmpmp(Q1, Q2, p1, p2, p3, p4, p5, p6, E); break;
    case 7: matrix_tree_mppmmp(Q1, Q2, p1, p2, p3, p4, p5, p6, E); break;
    }
    
    return 16.0*amptree(E);
  }

}

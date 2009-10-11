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
#include "ampq4p2.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {

#define Leading12 2.0*real(A1*conj(A1))
#define Leading14 2.0*real(A3*conj(A3))
#define SubLeading 4.0*real(A1*conj(A3))/Nc

  double ampq4p2::
  su3_tree(double Q1, double Q2, int p1, int p2,
	   int p3, int p4, int p5, int p6) const
  {
    double Amp;
    _ComplexD A1;
    
    A1 = Apmpmpp(Q1,Q2,p1,p2,p3,p4,p5,p6); Amp  = Leading12;
    A1 = Apmmppp(Q1,Q2,p1,p2,p3,p4,p5,p6); Amp += Leading12;
    A1 = Apmpmmp(Q1,Q2,p1,p2,p3,p4,p5,p6); Amp += Leading12;
    A1 = Apmmpmp(Q1,Q2,p1,p2,p3,p4,p5,p6); Amp += Leading12;
    A1 = Apmpmpm(Q1,Q2,p1,p2,p3,p4,p5,p6); Amp += Leading12;
    A1 = Apmmppm(Q1,Q2,p1,p2,p3,p4,p5,p6); Amp += Leading12;
    A1 = Apmpmmm(Q1,Q2,p1,p2,p3,p4,p5,p6); Amp += Leading12;
    A1 = Apmmpmm(Q1,Q2,p1,p2,p3,p4,p5,p6); Amp += Leading12;

    return 4.0*Na*Amp;
  }

  double ampq4p2::
  su3_tree_mch(double Q1, double Q2, int p1, int p2,
	       int p3, int p4, int p5, int p6) const
  {
    _ComplexD A1;
    unsigned int hel = (unsigned int) (8*_M_rng());
    
    switch(hel) {
    case 0: A1 = Apmpmpp(Q1,Q2,p1,p2,p3,p4,p5,p6); break;
    case 1: A1 = Apmmppp(Q1,Q2,p1,p2,p3,p4,p5,p6); break;
    case 2: A1 = Apmpmmp(Q1,Q2,p1,p2,p3,p4,p5,p6); break;
    case 3: A1 = Apmmpmp(Q1,Q2,p1,p2,p3,p4,p5,p6); break;
    case 4: A1 = Apmpmpm(Q1,Q2,p1,p2,p3,p4,p5,p6); break;
    case 5: A1 = Apmmppm(Q1,Q2,p1,p2,p3,p4,p5,p6); break;
    case 6: A1 = Apmpmmm(Q1,Q2,p1,p2,p3,p4,p5,p6); break;
    case 7: A1 = Apmmpmm(Q1,Q2,p1,p2,p3,p4,p5,p6); break;
    }
    
    return 32.0*Na*Leading12 ;
  }
  
  
  void ampq4p2::
  su3_tree(double Q, int p1, int p2, int p3,
	   int p4, int p5, int p6, double *Amp) const
  {
    _ComplexD A1, A3;
    
    A1 = Apmpmpp(Q,Q,p1,p2,p3,p4,p5,p6); 
    A3 = Apmpmpp(Q,Q,p1,p4,p3,p2,p5,p6);
    Amp[0] = Leading12; Amp[1] = Leading14 + SubLeading;

    A1 = Apmmppp(Q,Q,p1,p2,p3,p4,p5,p6);
    A3 = Apmmppp(Q,Q,p1,p4,p3,p2,p5,p6);
    Amp[0] += Leading12; Amp[1] += Leading14;
    
    A1 = Apmpmmp(Q,Q,p1,p2,p3,p4,p5,p6);
    A3 = Apmpmmp(Q,Q,p1,p4,p3,p2,p5,p6);
    Amp[0] += Leading12; Amp[1] += Leading14 + SubLeading;
    
    A1 = Apmmpmp(Q,Q,p1,p2,p3,p4,p5,p6);
    A3 = Apmmpmp(Q,Q,p1,p4,p3,p2,p5,p6);
    Amp[0] += Leading12; Amp[1] += Leading14;
    
    A1 = Apmpmpm(Q,Q,p1,p2,p3,p4,p5,p6);
    A3 = Apmpmpm(Q,Q,p1,p4,p3,p2,p5,p6);
    Amp[0] += Leading12; Amp[1] += Leading14 + SubLeading;
   
    A1 = Apmmppm(Q,Q,p1,p2,p3,p4,p5,p6);
    A3 = Apmmppm(Q,Q,p1,p4,p3,p2,p5,p6);
    Amp[0] += Leading12; Amp[1] += Leading14;
    
    A1 = Apmpmmm(Q,Q,p1,p2,p3,p4,p5,p6);
    A3 = Apmpmmm(Q,Q,p1,p4,p3,p2,p5,p6);
    Amp[0] += Leading12; Amp[1] += Leading14 + SubLeading;
    
    A1 = Apmmpmm(Q,Q,p1,p2,p3,p4,p5,p6);
    A3 = Apmmpmm(Q,Q,p1,p4,p3,p2,p5,p6);
    Amp[0] += Leading12; Amp[1] += Leading14;

    Amp[0] *= 4.0*Na; Amp[1] *= 4.0*Na; Amp[1] += Amp[0];
  }

  void ampq4p2::
  su3_tree_mch(double Q, int p1, int p2, int p3, 
	       int p4, int p5, int p6, double *Amp) const
  {
    _ComplexD A1, A3;
    unsigned int hel = (unsigned int) (8*_M_rng());
    Amp[0] = Amp[1] = 0.0;
    
    switch(hel) {
    case 0:
      A1 = Apmpmpp(Q,Q,p1,p2,p3,p4,p5,p6); 
      A3 = Apmpmpp(Q,Q,p1,p4,p3,p2,p5,p6);
      Amp[0] = Leading12; Amp[1] = Leading14 + SubLeading;
      break;
    case 1:
      A1 = Apmmppp(Q,Q,p1,p2,p3,p4,p5,p6);
      A3 = Apmmppp(Q,Q,p1,p4,p3,p2,p5,p6);
      Amp[0] = Leading12; Amp[1] = Leading14;
      break;
    case 2:
      A1 = Apmpmmp(Q,Q,p1,p2,p3,p4,p5,p6);
      A3 = Apmpmmp(Q,Q,p1,p4,p3,p2,p5,p6);
      Amp[0] = Leading12; Amp[1] = Leading14 + SubLeading;
      break;
    case 3:    
      A1 = Apmmpmp(Q,Q,p1,p2,p3,p4,p5,p6);
      A3 = Apmmpmp(Q,Q,p1,p4,p3,p2,p5,p6);
      Amp[0] = Leading12; Amp[1] = Leading14;
      break;
    case 4:    
      A1 = Apmpmpm(Q,Q,p1,p2,p3,p4,p5,p6);
      A3 = Apmpmpm(Q,Q,p1,p4,p3,p2,p5,p6);
      Amp[0] = Leading12; Amp[1] = Leading14 + SubLeading;
      break;
    case 5:  
      A1 = Apmmppm(Q,Q,p1,p2,p3,p4,p5,p6);
      A3 = Apmmppm(Q,Q,p1,p4,p3,p2,p5,p6);
      Amp[0] = Leading12; Amp[1] = Leading14;
      break;
    case 6:    
      A1 = Apmpmmm(Q,Q,p1,p2,p3,p4,p5,p6);
      A3 = Apmpmmm(Q,Q,p1,p4,p3,p2,p5,p6);
      Amp[0] = Leading12; Amp[1] = Leading14 + SubLeading;
      break;
    case 7:    
      A1 = Apmmpmm(Q,Q,p1,p2,p3,p4,p5,p6);
      A3 = Apmmpmm(Q,Q,p1,p4,p3,p2,p5,p6);
      Amp[0] = Leading12; Amp[1] = Leading14;
      break;
    }
    
    Amp[0] *= 32.0*Na; Amp[1] *= 32.0*Na; Amp[1] += Amp[0];
  }


  std::complex<double>
  ampq4p2::g1(int p1, int p2, int p3, int p4, int p5, int p6) const
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
  ampq4p2::g2(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    double s126 = S3(1,2,6);
    _ComplexD a16 = A(1,6), a62 = A(6,2), b35 = B(3,5), b54 = B(5,4), 
      d3162 = D(3,1,6,2); 
    
    return -ICPLX*d3162*d3162/(a16*a62*b35*b54*s126);
  }

  std::complex<double> ampq4p2::
  Amhv(double Q1, double Q2, int pi, int pj, int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    _ComplexD a21 = A(2,1), a25 = A(2,5), a26 = A(2,6), a43 = A(4,3), 
      a45 = A(4,5), a46 = A(4,6), a51 = A(5,1), a53 = A(5,3), a61 = A(6,1),
      a63 = A(6,3), aij = A(i, j); 
    
    return ICPLX*(Q1*a21/(a25*a51) + Q2*a43/(a45*a53))
      *(Q1*a21/(a26*a61) + Q2*a43/(a46*a63))*aij*aij/(a21*a43);
  }
}


//   double ampq4p2::
//   su3_tree(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const
//   {
//     double Amp;
//     _ComplexD A1, A2, A3, A4;
    
//     A1 = Apmpmpp(Q1,Q2,p1,p2,p3,p4,p5,p6); 
//     A2 = Apmpmpp(Q2,Q1,p3,p4,p1,p2,p5,p6);
//     A3 = Apmpmpp(Q1,Q2,p1,p4,p3,p2,p5,p6);
//     A4 = Apmpmpp(Q2,Q1,p3,p2,p1,p4,p5,p6);
//     Amp = Leading12 + Leading14 + SubLeading;

//     A1 = Apmmppp(Q1,Q2,p1,p2,p3,p4,p5,p6);
//     A2 = Amppmpp(Q2,Q1,p3,p4,p1,p2,p5,p6);
//     A3 = Apmmppp(Q1,Q2,p1,p4,p3,p2,p5,p6);
//     A4 = Amppmpp(Q2,Q1,p3,p2,p1,p4,p5,p6);
//     Amp += Leading12 + Leading14;
    
//     A1 = Apmpmmp(Q1,Q2,p1,p2,p3,p4,p5,p6);
//     A2 = Apmpmmp(Q2,Q1,p3,p4,p1,p2,p5,p6);
//     A3 = Apmpmmp(Q1,Q2,p1,p4,p3,p2,p5,p6);
//     A4 = Apmpmmp(Q2,Q1,p3,p2,p1,p4,p5,p6);
//     Amp += Leading12 + Leading14 + SubLeading;
    
//     A1 = Apmmpmp(Q1,Q2,p1,p2,p3,p4,p5,p6);
//     A2 = Amppmmp(Q2,Q1,p3,p4,p1,p2,p5,p6);
//     A3 = Apmmpmp(Q1,Q2,p1,p4,p3,p2,p5,p6);
//     A4 = Amppmmp(Q2,Q1,p3,p2,p1,p4,p5,p6);
//     Amp += Leading12 + Leading14;
    
//     A1 = Apmpmpm(Q1,Q2,p1,p2,p3,p4,p5,p6);
//     A2 = Apmpmpm(Q2,Q1,p3,p4,p1,p2,p5,p6);
//     A3 = Apmpmpm(Q1,Q2,p1,p4,p3,p2,p5,p6);
//     A4 = Apmpmpm(Q2,Q1,p3,p2,p1,p4,p5,p6);
//     Amp += Leading12 + Leading14 + SubLeading;
   
//     A1 = Apmmppm(Q1,Q2,p1,p2,p3,p4,p5,p6);
//     A2 = Amppmpm(Q2,Q1,p3,p4,p1,p2,p5,p6);
//     A3 = Apmmppm(Q1,Q2,p1,p4,p3,p2,p5,p6);
//     A4 = Amppmpm(Q2,Q1,p3,p2,p1,p4,p5,p6);
//     Amp += Leading12 + Leading14;
    
//     A1 = Apmpmmm(Q1,Q2,p1,p2,p3,p4,p5,p6);
//     A2 = Apmpmmm(Q2,Q1,p3,p4,p1,p2,p5,p6);
//     A3 = Apmpmmm(Q1,Q2,p1,p4,p3,p2,p5,p6);
//     A4 = Apmpmmm(Q2,Q1,p3,p2,p1,p4,p5,p6);
//     Amp += Leading12 + Leading14 + SubLeading;
    
//     A1 = Apmmpmm(Q1,Q2,p1,p2,p3,p4,p5,p6);
//     A2 = Amppmmm(Q2,Q1,p3,p4,p1,p2,p5,p6);
//     A3 = Apmmpmm(Q1,Q2,p1,p4,p3,p2,p5,p6);
//     A4 = Amppmmm(Q2,Q1,p3,p2,p1,p4,p5,p6);
//     Amp += Leading12 + Leading14;

//     return 4.0*Na*Amp;
//   }

//   double ampq4p2::
//   su3_tree_mch(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const
//   {
//     double Amp = 0.0;
//     _ComplexD A1, A2, A3, A4;
//     unsigned int hel = (unsigned int) (8*_M_rng());
   
//     switch(hel) {
//     case 0:
//       A1 = Apmpmpp(Q1,Q2,p1,p2,p3,p4,p5,p6); 
//       A2 = Apmpmpp(Q2,Q1,p3,p4,p1,p2,p5,p6);
//       A3 = Apmpmpp(Q1,Q2,p1,p4,p3,p2,p5,p6);
//       A4 = Apmpmpp(Q2,Q1,p3,p2,p1,p4,p5,p6);
//       Amp = Leading12 + Leading14 + SubLeading;
//       break;
//     case 1:
//       A1 = Apmmppp(Q1,Q2,p1,p2,p3,p4,p5,p6);
//       A2 = Amppmpp(Q2,Q1,p3,p4,p1,p2,p5,p6);
//       A3 = Apmmppp(Q1,Q2,p1,p4,p3,p2,p5,p6);
//       A4 = Amppmpp(Q2,Q1,p3,p2,p1,p4,p5,p6);
//       Amp = Leading12 + Leading14;
//       break;
//     case 2:
//       A1 = Apmpmmp(Q1,Q2,p1,p2,p3,p4,p5,p6);
//       A2 = Apmpmmp(Q2,Q1,p3,p4,p1,p2,p5,p6);
//       A3 = Apmpmmp(Q1,Q2,p1,p4,p3,p2,p5,p6);
//       A4 = Apmpmmp(Q2,Q1,p3,p2,p1,p4,p5,p6);
//       Amp = Leading12 + Leading14 + SubLeading;
//       break;
//     case 3:    
//       A1 = Apmmpmp(Q1,Q2,p1,p2,p3,p4,p5,p6);
//       A2 = Amppmmp(Q2,Q1,p3,p4,p1,p2,p5,p6);
//       A3 = Apmmpmp(Q1,Q2,p1,p4,p3,p2,p5,p6);
//       A4 = Amppmmp(Q2,Q1,p3,p2,p1,p4,p5,p6);
//       Amp = Leading12 + Leading14;
//       break;
//     case 4:    
//       A1 = Apmpmpm(Q1,Q2,p1,p2,p3,p4,p5,p6);
//       A2 = Apmpmpm(Q2,Q1,p3,p4,p1,p2,p5,p6);
//       A3 = Apmpmpm(Q1,Q2,p1,p4,p3,p2,p5,p6);
//       A4 = Apmpmpm(Q2,Q1,p3,p2,p1,p4,p5,p6);
//       Amp = Leading12 + Leading14 + SubLeading;
//       break;
//     case 5:  
//       A1 = Apmmppm(Q1,Q2,p1,p2,p3,p4,p5,p6);
//       A2 = Amppmpm(Q2,Q1,p3,p4,p1,p2,p5,p6);
//       A3 = Apmmppm(Q1,Q2,p1,p4,p3,p2,p5,p6);
//       A4 = Amppmpm(Q2,Q1,p3,p2,p1,p4,p5,p6);
//       Amp = Leading12 + Leading14;
//       break;
//     case 6:    
//       A1 = Apmpmmm(Q1,Q2,p1,p2,p3,p4,p5,p6);
//       A2 = Apmpmmm(Q2,Q1,p3,p4,p1,p2,p5,p6);
//       A3 = Apmpmmm(Q1,Q2,p1,p4,p3,p2,p5,p6);
//       A4 = Apmpmmm(Q2,Q1,p3,p2,p1,p4,p5,p6);
//       Amp = Leading12 + Leading14 + SubLeading;
//       break;
//     case 7:    
//       A1 = Apmmpmm(Q1,Q2,p1,p2,p3,p4,p5,p6);
//       A2 = Amppmmm(Q2,Q1,p3,p4,p1,p2,p5,p6);
//       A3 = Apmmpmm(Q1,Q2,p1,p4,p3,p2,p5,p6);
//       A4 = Amppmmm(Q2,Q1,p3,p2,p1,p4,p5,p6);
//       Amp = Leading12 + Leading14;
//       break;
//     }
    
//     return 32.0*Na*Amp;
//   }

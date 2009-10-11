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
#include "ampq2g2p2.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {

#define AmpTree (Na*real(A1*conj(A1) + A2*conj(A2)) - 2.0*real(A1*conj(A2)))

  double ampq2g2p2::
  su3_tree(int p1, int p2, int p3, int p4, int p5,int p6) const 
  {
    double Amp;
    _ComplexD A1,A2;
    
    A1 = Apmmppp(p1,p2,p3,p4,p5,p6); A2 = Apmpmpp(p1,p2,p4,p3,p5,p6); Amp  = AmpTree;
    A1 = Apmpmpp(p1,p2,p3,p4,p5,p6); A2 = Apmmppp(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmppmp(p1,p2,p3,p4,p5,p6); A2 = Apmppmp(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmpppm(p1,p2,p3,p4,p5,p6); A2 = Apmpppm(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmmmpp(p1,p2,p3,p4,p5,p6); A2 = Apmmmpp(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmmpmp(p1,p2,p3,p4,p5,p6); A2 = Apmpmmp(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmmppm(p1,p2,p3,p4,p5,p6); A2 = Apmpmpm(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmpmmp(p1,p2,p3,p4,p5,p6); A2 = Apmmpmp(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmpmpm(p1,p2,p3,p4,p5,p6); A2 = Apmmppm(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmppmm(p1,p2,p3,p4,p5,p6); A2 = Apmppmm(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmpmmm(p1,p2,p3,p4,p5,p6); A2 = Apmmpmm(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmmpmm(p1,p2,p3,p4,p5,p6); A2 = Apmpmmm(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmmmpm(p1,p2,p3,p4,p5,p6); A2 = Apmmmpm(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    A1 = Apmmmmp(p1,p2,p3,p4,p5,p6); A2 = Apmmmmp(p1,p2,p4,p3,p5,p6); Amp += AmpTree;
    
    return 8.0*Na*Amp/Nc;
  }

  double ampq2g2p2::
  su3_tree_mch(int p1, int p2, int p3, int p4, int p5,int p6) const 
  {
    _ComplexD A1,A2;
    unsigned int hel = (unsigned int) (14*_M_rng());

    switch(hel) {
    case 0:  A1 = Apmmppp(p1,p2,p3,p4,p5,p6); A2 = Apmpmpp(p1,p2,p4,p3,p5,p6); break;
    case 1:  A1 = Apmpmpp(p1,p2,p3,p4,p5,p6); A2 = Apmmppp(p1,p2,p4,p3,p5,p6); break;
    case 2:  A1 = Apmppmp(p1,p2,p3,p4,p5,p6); A2 = Apmppmp(p1,p2,p4,p3,p5,p6); break;
    case 3:  A1 = Apmpppm(p1,p2,p3,p4,p5,p6); A2 = Apmpppm(p1,p2,p4,p3,p5,p6); break;
    case 4:  A1 = Apmmmpp(p1,p2,p3,p4,p5,p6); A2 = Apmmmpp(p1,p2,p4,p3,p5,p6); break;
    case 5:  A1 = Apmmpmp(p1,p2,p3,p4,p5,p6); A2 = Apmpmmp(p1,p2,p4,p3,p5,p6); break;
    case 6:  A1 = Apmmppm(p1,p2,p3,p4,p5,p6); A2 = Apmpmpm(p1,p2,p4,p3,p5,p6); break;
    case 7:  A1 = Apmpmmp(p1,p2,p3,p4,p5,p6); A2 = Apmmpmp(p1,p2,p4,p3,p5,p6); break;
    case 8:  A1 = Apmpmpm(p1,p2,p3,p4,p5,p6); A2 = Apmmppm(p1,p2,p4,p3,p5,p6); break;
    case 9:  A1 = Apmppmm(p1,p2,p3,p4,p5,p6); A2 = Apmppmm(p1,p2,p4,p3,p5,p6); break;
    case 10: A1 = Apmpmmm(p1,p2,p3,p4,p5,p6); A2 = Apmmpmm(p1,p2,p4,p3,p5,p6); break;
    case 11: A1 = Apmmpmm(p1,p2,p3,p4,p5,p6); A2 = Apmpmmm(p1,p2,p4,p3,p5,p6); break;
    case 12: A1 = Apmmmpm(p1,p2,p3,p4,p5,p6); A2 = Apmmmpm(p1,p2,p4,p3,p5,p6); break;
    case 13: A1 = Apmmmmp(p1,p2,p3,p4,p5,p6); A2 = Apmmmmp(p1,p2,p4,p3,p5,p6); break;
    }

    return 112.0*Na*AmpTree/Nc;
  }
  
  std::complex<double> 
  ampq2g2p2::Amhv(int pi, int p1, int p2, int p3, int p4, int p5,int p6) const 
  {
    _ComplexD a12 = A(1,2), a15 = A(1,5), a16 = A(1,6), a23 = A(2,3), 
      a34 = A(3,4), a41 = A(4,1), a52 = A(5,2), a62 = A(6,2),
      a1i = Aij(p1,pi), a2i = Aij(p2,pi);

    return ICPLX*a1i*std::pow(a2i, 3)/(a12*a23*a34*a41)
      *a12/(a15*a52)*a12/(a16*a62);
  }
 
  std::complex<double>
  ampq2g2p2::Apmmmpp(int p5, int p6, int p1, int p2, int p3, int p4) const 
  {
    double s126 = S3(1,2,6), s136 = S3(1,3,6), s146 = S3(1,4,6);
   
    _ComplexD a16 = A(1,6), a35 = A(3,5), a45 = A(4,5), a53 = A(5,3), 
      a54 = A(5,4), a56 = A(5,6), a63 = A(6,3), a64 = A(6,4), b16 = B(1,6),
      b21 = B(2,1), b35 = B(3,5), b45 = B(4,5), b52 = B(5,2);
    
    _ComplexD d3162 = D(3,1,6,2), d4162 = D(4,1,6,2), d5126 = D(5,1,2,6);

    return ICPLX*(d5126*a56*s126/(b52*b16*b21*a63*a64*a53*a54)
		  - a16*b45*d3162/(b16*a45*b52*a63*s136)
		  - a16*b35*d4162/(b16*a35*b52*a64*s146));
  }
  
  std::complex<double>
  ampq2g2p2::Apmpmmp(int p5, int p6, int p1, int p2, int p3, int p4) const
  {
    double s12 = S(1,2), s125 = S3(1,2,5), s126 = S3(1,2,6), s136 = S3(1,3,6),
      s145 = S3(1,4,5); 
    
    _ComplexD a16 = A(1,6), a26 = A(2,6), a36 = A(3,6), a45 = A(4,5), 
      a46 = A(4,6), a54 = A(5,4), a61 = A(6,1), a62 = A(6,2), a64 = A(6,4),
      b15 = B(1,5), b25 = B(2,5), b36 = B(3,6), b45 = B(4,5), b51 = B(5,1),
      b52 = B(5,2), b53 = B(5,3), b54 = B(5,4);
    
    _ComplexD d1263 = D(1,2,6,3), d1362 = D(1,3,6,2), d4362 = D(4,3,6,2),
      d5146 = D(5,1,4,6), d5236 = D(5,2,3,6);
    
    return ICPLX*(b45*a26*a26*d1263/(a45*b53*a16*s12*s126)
		  + a36*b54*d1362/(b36*a54*b52*a16*s136)
		  + a36*b15*b15*d4362/(b36*b25*a64*s12*s125)
		  + b51*a62*d5146*s145/(b36*b52*b53*a61*a64*a45*s12)
		  + d1362*(b51*a16*d5236+b54*a46*b52*a26)/
		  (b36*b52*b53*a61*a64*a45*s12));
  }
   
  std::complex<double>
  ampq2g2p2::Apmmpmp(int p5,int p6,int p1,int p2,int p3,int p4) const
  {
    double s12 = S(1,2), s126 = S3(1,2,6), s136 = S3(1,3,6), s146 = S3(1,4,6),
      s346 = S3(3,4,6);
    
    _ComplexD a16 = A(1,6), a25 = A(2,5), a36 = A(3,6), a45 = A(4,5), 
      a46 = A(4,6), a51 = A(5,1), a54 = A(5,4), a64 = A(6,4), b12 = B(1,2),
      b16 = B(1,6), b25 = B(2,5), b26 = B(2,6), b35 = B(3,5), b36 = B(3,6),
      b45 = B(4,5), b53 = B(5,3), b62 = B(6,2), b65 = B(6,5);
    
    _ComplexD  d2453 = D(2,4,5,3), d4251 = D(4,2,5,1), d4253 = D(4,2,5,3);

    return ICPLX*(a16*b25*d4253/(b16*a25*b53*a64*s146)
		  + a16*b45*b62*d2453/(b16*a45*b53*s12*s126)
		  + a36*b25*a51*d4251/(b36*a25*a64*s12*s346)
		  + a16*b25*d2453/(b16*b35*a54*a46*s12)
		  + s136*(a51*b26*b53*a36 + b12*b65*a16*a51)/
		  (b16*b53*b36*a25*a54*a46*s12));
  }

}  //  namespace nlo

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

//   Standard includes
#include <cstring>

//   nlo includes
#include "bits/nlo-color.h"
#include "ampq4g1l2.h"
#include "defmacros.h"


namespace nlo {


  std::complex<double> 
  ampq4g1l2::A1ppp(int p1, int p2, int p3, int p4, int p5, int p6, int p7)
  {
    double s15 = S(1,5), s34 = S(3,4), t267 = S3(2,6,7), t234 = S3(2,3,4),
      t167 = S3(1,6,7), t345 = S3(3,4,5);
    
    _ComplexD b15 = B(1,5), a62 = A(6,2), a45 = A(4,5), b17 = B(1,7), 
      a42 = A(4,2), b23 = B(2,3), a15 = A(1,5), a54 = A(5,4), b53 = B(5,3),
      a34 = A(3,4), b35 = B(3,5);
    
    _ComplexD c4153 = C(4,1,5,3), c4267 = C(4,2,6,7), c6175 = C(6,1,7,5),
      c4157 = C(4,1,5,7), c6243 = C(6,2,4,3), c4351 = C(4,3,5,1),
      c6173 = C(6,1,7,3);
    
    _ComplexD M = - b15*c4153*c4267*a62/(a45*s15*s34*t267)        //  --> +H1
      - b17*c6175*a42*a42*b23/(a45*s34*t234*t167)                //  --> +H2  
      - c4157*c6243*a42/(a15*a54*s34*t234)                       //  --> +H3
      + b53*c4351*c4267*a62/(a45*s34*t345*t267)                  //  --> -H4 
      + b17*(c6173*a34 + c6175*a54)*b35*a42/(a45*s34*t345*t167); //  --> -H5
    
    return M/S(6,7);
  }
  
  std::complex<double> 
  ampq4g1l2::A1ppm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s15 = S(1,5), s34 = S(3,4), t267 = S3(2,6,7), t234 = S3(2,3,4),
      t167 = S3(1,6,7), t345 = S3(3,4,5);
    
    _ComplexD b13 = B(1,3), a51 = A(5,1), a62 = A(6,2), b35 = B(3,5), 
      b17 = B(1,7), a42 = A(4,2), b15 = B(1,5), b53 = B(5,3), a54 = A(5,4), 
      b34 = B(3,4);
    
    _ComplexD c4267 = C(4,2,6,7), c6243 = C(6,2,4,3), c6173 = C(6,1,7,3), 
      c5243 = C(5,2,4,3), c5267 = C(5,2,6,7), c2453 = C(2,4,5,3);  
    
    _ComplexD M = b13*b13*a51*c4267*a62/(b35*s15*s34*t267)
      + b17*c6173*c5243*a42/(b35*s34*t234*t167)
      - b13*b17*c6243*a42/(b15*b53*s34*t234)
      + b13*a54*(b34*c4267 + b35*c5267)*a62/(b35*s34*t345*t267)
      - b17*c6173*a54*c2453/(b35*s34*t345*t167);
    
    return M/S(6,7);
  }

  std::complex<double> 
  ampq4g1l2::A1pmp(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s15 = S(1,5), s34 = S(3,4), t267 = S3(2,6,7), t243 = S3(2,4,3),
      t167 = S3(1,6,7), t345 = S3(3,4,5);
    
    _ComplexD b15 = B(1,5), a62 = A(6,2), a35 = A(3,5), b17 = B(1,7), 
      a32 = A(3,2), b24 = B(2,4), a15 = A(1,5), a53 = A(5,3), b54 = B(5,4), 
      a43 = A(4,3), b45 = B(4,5);
    
    _ComplexD c3451 = C(3,4,5,1), c6174 = C(6,1,7,4), c3157 = C(3,1,5,7),
      c6234 = C(6,2,3,4), c6175 = C(6,1,7,5), c3154 = C(3,1,5,4), 
      c3267 = C(3,2,6,7); 
    
    _ComplexD M = - b15*c3154*c3267*a62/(a35*s15*s34*t267)
      - b17*c6175*a32*a32*b24/(a35*s34*t243*t167)
      - c3157*c6234*a32/(a15*a53*s34*t243)
      + b54*c3451*c3267*a62/(a35*s34*t345*t267)
      + b17*(c6174*a43 + c6175*a53)*b45*a32/(a35*s34*t345*t167);
    
    return M/S(6,7);
  }

  std::complex<double> 
  ampq4g1l2::A1pmm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s15 = S(1,5), s34 = S(3,4), t267 = S3(2,6,7), t243 = S3(2,4,3),
      t167 = S3(1,6,7), t435 = S3(4,3,5);
    
    _ComplexD b14 = B(1,4), a51 = A(5,1), a62 = A(6,2), b45 = B(4,5), 
      b17 = B(1,7), a32 = A(3,2), b15 = B(1,5), b54 = B(5,4), a53 = A(5,3),
      b43 = B(4,3);
    
    _ComplexD c5267 = C(5,2,6,7), c3267 = C(3,2,6,7), c6234 = C(6,2,3,4),
      c6174 = C(6,1,7,4), c5234 = C(5,2,3,4), c2354 = C(2,3,5,4);
    
    _ComplexD M = b14*b14*a51*c3267*a62/(b45*s15*s34*t267)
      + b17*c6174*c5234*a32/(b45*s34*t243*t167)
      - b14*b17*c6234*a32/(b15*b54*s34*t243)
      + b14*a53*(b43*c3267 + b45*c5267)*a62/(b45*s34*t435*t267)
      - b17*c6174*a53*c2354/(b45*s34*t435*t167);
    
    return M/S(6,7);
  }

  std::complex<double> 
  ampq4g1l2::A2ppp(int p1, int p2, int p3, int p4, int p5, int p6, int p7)
  {
    double s34 = S(3,4), s25 = S(2,5), t345 = S3(3,4,5), t267 = S3(2,6,7),
      t167 = S3(1,6,7), t134 = S3(1,3,4);
    
    _ComplexD b53 = B(5,3), a62 = A(6,2), a45 = A(4,5), b17 = B(1,7),
      a34 = A(3,4), a54 = A(5,4), b35 = B(3,5), a42 = A(4,2), b13 = B(1,3),
      a52 = A(5,2), b25 = B(2,5);
    
    _ComplexD c4351 = C(4,3,5,1), c4267 = C(4,2,6,7), c6173 = C(6,1,7,3), 
      c6175 = C(6,1,7,5), c4135 = C(4,1,3,5), c4137 = C(4,1,3,7);
    
    _ComplexD M =- b53*c4351*c4267*a62/(a45*s34*t345*t267)         // --> H4 
      - b17*(c6173*a34 + c6175*a54)*b35*a42/(a45*s34*t345*t167)    // --> H5
      - b13*c4135*c4267*a62/(a45*s34*t134*t267)                    // --> H6
      - b13*c4137*a62*a42/(a45*a52*s34*t134)                       // --> H7
      - b17*c6173*b25*a42*a42/(a45*s25*s34*t167);
    
    return M/S(6,7);
  }

  std::complex<double> 
  ampq4g1l2::A2ppm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s34 = S(3,4), s25 = S(2,5), t345 = S3(3,4,5), t267 = S3(2,6,7),
      t167 = S3(1,6,7), t134 = S3(1,3,4);
    
    _ComplexD b13 = B(1,3), a54 = A(5,4), b34 = B(3,4), b35 = B(3,5), 
      a62 = A(6,2), b17 = B(1,7), a41 = A(4,1), b52 = B(5,2), a52 = A(5,2);
    
    _ComplexD c4267 = C(4,2,6,7), c5267 = C(5,2,6,7), c6173 = C(6,1,7,3),
      c2453 = C(2,4,5,3), c4137 = C(4,1,3,7), c6253 = C(6,2,5,3),
      c4253 = C(4,2,5,3);
    
    _ComplexD M = - b13*a54*(b34*c4267 + b35*c5267)*a62/(b35*s34*t345*t267)
      + b17*c6173*a54*c2453/(b35*s34*t345*t167)
      + b13*b13*a41*c5267*a62/(b35*s34*t134*t267)
      - b13*c4137*c6253/(b35*b52*s34*t134)
      + b17*c6173*a52*c4253/(b35*s25*s34*t167);
    
    return M/S(6,7);
  }
  
  std::complex<double> 
  ampq4g1l2::A2pmp(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s34 = S(3,4), s25 = S(2,5), t345 = S3(3,4,5), t267 = S3(2,6,7),
      t167 = S3(1,6,7), t143 = S3(1,4,3);
    
    _ComplexD b54 = B(5,4), a62 = A(6,2), a35 = A(3,5), b17 = B(1,7), 
      a43 = A(4,3), a53 = A(5,3), b45 = B(4,5), a32 = A(3,2), b14 = B(1,4), 
      a52 = A(5,2), b25 = B(2,5); 
    
    _ComplexD c6174 = C(6,1,7,4), c6175 = C(6,1,7,5), c3145 = C(3,1,4,5), 
      c3147 = C(3,1,4,7), c3451 = C(3,4,5,1), c3267 = C(3,2,6,7); 
    
    _ComplexD M = - b54*c3451*c3267*a62/(a35*s34*t345*t267)
      - b17*(c6174*a43 + c6175*a53)*b45*a32/(a35*s34*t345*t167)
      - b14*c3145*c3267*a62/(a35*s34*t143*t267)
      - b14*c3147*a62*a32/(a35*a52*s34*t143)
      - b17*c6174*b25*a32*a32/(a35*s25*s34*t167);
    
    return M/S(6,7);
  }

  std::complex<double> 
  ampq4g1l2::A2pmm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s25 = S(2,5), s34 = S(3,4), t435 = S3(4,3,5), t267 = S3(2,6,7), 
      t167 = S3(1,6,7), t143 = S3(1,4,3);
    
    _ComplexD b14 = B(1,4), a53 = A(5,3), b43 = B(4,3), b45 = B(4,5), 
      a62 = A(6,2), b17 = B(1,7), a31 = A(3,1), b52 = B(5,2), a52 = A(5,2);
    
    _ComplexD c5267 = C(5,2,6,7), c3267 = C(3,2,6,7), c6174 = C(6,1,7,4), 
      c2354 = C(2,3,5,4), c3254 = C(3,2,5,4), c3147 = C(3,1,4,7),
      c6254 = C(6,2,5,4);
    
    _ComplexD M = - b14*a53*(b43*c3267 + b45*c5267)*a62/(b45*s34*t435*t267)
      + b17*c6174*a53*c2354/(b45*s34*t435*t167)
      + b14*b14*a31*c5267*a62/(b45*s34*t143*t267)
      - b14*c3147*c6254/(b45*b52*s34*t143)
      + b17*c6174*a52*c3254/(b45*s25*s34*t167);
    
    return M/S(6,7);
  }

  std::complex<double>
  ampq4g1l2::A3ppp(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s35 = S(3,5), t345 = S3(3,4,5), t267 = S3(2,6,7), t167 = S3(1,6,7);
    
    _ComplexD b53 = B(5,3), a62 = A(6,2), a45 = A(4,5), b17 = B(1,7), 
      a34 = A(3,4), a54 = A(5,4), b35 = B(3,5), a42 = A(4,2); 
    
    _ComplexD c6175 = C(6,1,7,5), c6173 = C(6,1,7,3), c4351 = C(4,3,5,1),
      c4267 = C(4,2,6,7);
    
    _ComplexD M =- b53*c4351*c4267*a62/(a45*s35*t345*t267)         // --> H9 
      - b17*(c6173*a34 + c6175*a54)*b35*a42/(a45*s35*t345*t167);   // --> H10
    
    return M/S(6,7);
  }

  std::complex<double> 
  ampq4g1l2::A3pmm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s35 = S(3,5), t435 = S3(4,3,5), t267 = S3(2,6,7), t167 = S3(1,6,7);
    
    _ComplexD b14 = B(1,4), a53 = A(5,3), b43 = B(4,3), b45 = B(4,5),
      a62 = A(6,2), b17 = B(1,7);
    
    _ComplexD c3267 = C(3,2,6,7), c5267 = C(5,2,6,7), c6174 = C(6,1,7,4), 
      c2354 = C(2,3,5,4);
    
    _ComplexD M = - b14*a53*(b43*c3267 + b45*c5267)*a62/(b45*s35*t435*t267)
      + b17*c6174*a53*c2354/(b45*s35*t435*t167);
    
    return M/S(6,7);
  }

  std::complex<double> 
  ampq4g1l2::A4ppm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s45 = S(4,5), t345 = S3(3,4,5), t267 = S3(2,6,7), t167 = S3(1,6,7);
    
    _ComplexD b13 = B(1,3), a54 = A(5,4), b34 = B(3,4), b35 = B(3,5), 
      a62 = A(6,2), b17 = B(1,7); 
    
    _ComplexD c4267 = C(4,2,6,7), c5267 = C(5,2,6,7), c6173 = C(6,1,7,3), 
      c2453 = C(2,4,5,3);
    
    _ComplexD M = b13*a54*(b34*c4267 + b35*c5267)*a62/(b35*s45*t345*t267)// H11
      - b17*c6173*a54*c2453/(b35*s45*t345*t167);                         // H12
    
    return M/S(6,7);
  }

  std::complex<double> 
  ampq4g1l2::A4pmp(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s45 = S(4,5), t345 = S3(3,4,5), t267 = S3(2,6,7), t167 = S3(1,6,7);
    
    _ComplexD b54 = B(5,4), a62 = A(6,2), a35 = A(3,5), b17 = B(1,7), 
      a43 = A(4,3), a53 = A(5,3), b45 = B(4,5), a32 = A(3,2); 
    
    _ComplexD c6175 = C(6,1,7,5), c6174 = C(6,1,7,4), c3451 = C(3,4,5,1), 
      c3267 = C(3,2,6,7);
    
    _ComplexD M = b54*c3451*c3267*a62/(a35*s45*t345*t267)
      + b17*(c6174*a43 + c6175*a53)*b45*a32/(a35*s45*t345*t167);
    
    return M/S(6,7);
  }
  
  void ampq4g1l2::
  matrix_tree_pmpmp(int p1, int p4, int p3, int p2,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = A1ppp(p1, p2, p3, p4, p5, p6, p7);
    M[1] = A2ppp(p1, p2, p3, p4, p5, p6, p7);
    M[2] = A3ppp(p1, p2, p3, p4, p5, p6, p7);
    M[3] = 0.0;
    
    M[4] = A3ppp(p3, p4, p1, p2, p5, p6, p7);
    M[5] = 0.0;
    M[6] = A1ppp(p3, p4, p1, p2, p5, p6, p7);
    M[7] = A2ppp(p3, p4, p1, p2, p5, p6, p7);
    
    
    M[8]  = A1ppp(p1, p4, p3, p2, p5, p6, p7);
    M[9]  = A2ppp(p1, p4, p3, p2, p5, p6, p7);
    M[10] = A3ppp(p1, p4, p3, p2, p5, p6, p7);
    M[11] = 0.0;
    
    M[12] = A3ppp(p3, p2, p1, p4, p5, p6, p7);
    M[13] = 0.0;
    M[14] = A1ppp(p3, p2, p1, p4, p5, p6, p7);
    M[15] = A2ppp(p3, p2, p1, p4, p5, p6, p7);
  }
  
  void ampq4g1l2::
  matrix_tree_pmpmm(int p1, int p4, int p3, int p2,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = A1ppm(p1, p2, p3, p4, p5, p6, p7);
    M[1] = A2ppm(p1, p2, p3, p4, p5, p6, p7);
    M[2] = 0;
    M[3] = A4ppm(p1, p2, p3, p4, p5, p6, p7);
    
    M[4] = 0.0;
    M[5] = A4ppm(p3, p4, p1, p2, p5, p6, p7);
    M[6] = A1ppm(p3, p4, p1, p2, p5, p6, p7);
    M[7] = A2ppm(p3, p4, p1, p2, p5, p6, p7);
    
    M[8]  = A1ppm(p1, p4, p3, p2, p5, p6, p7);
    M[9]  = A2ppm(p1, p4, p3, p2, p5, p6, p7);
    M[10] = 0.0;
    M[11] = A4ppm(p1, p4, p3, p2, p5, p6, p7);
    
    M[12] = 0.0;
    M[13] = A4ppm(p3, p2, p1, p4, p5, p6, p7);
    M[14] = A1ppm(p3, p2, p1, p4, p5, p6, p7);
    M[15] = A2ppm(p3, p2, p1, p4, p5, p6, p7);
  }
  
  void ampq4g1l2::
  matrix_tree_ppmmp(int p1, int p4, int p3, int p2,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = A1pmp(p1, p2, p3, p4, p5, p6, p7);
    M[1] = A2pmp(p1, p2, p3, p4, p5, p6, p7);
    M[2] = 0;
    M[3] = A4pmp(p1, p2, p3, p4, p5, p6, p7);
    
    swap();
    M[4] = -A3pmm(p3, p4, p1, p2, p5, p7, p6);
    M[5] = 0.0;
    M[6] = -A1pmm(p3, p4, p1, p2, p5, p7, p6);
    M[7] = -A2pmm(p3, p4, p1, p2, p5, p7, p6);
    swap();
    
    M[8] = M[9] = M[10] = M[11] = 0.0;
    M[12] = M[13] = M[14] = M[15] = 0.0;
  }
  
  void ampq4g1l2::
  matrix_tree_ppmmm(int p1, int p4, int p3, int p2,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = A1pmm(p1, p2, p3, p4, p5, p6, p7);
    M[1] = A2pmm(p1, p2, p3, p4, p5, p6, p7);
    M[2] = A3pmm(p1, p2, p3, p4, p5, p6, p7);
    M[3] = 0.0;
    
    swap();
    M[4] = 0.0;
    M[5] = -A4pmp(p3, p4, p1, p2, p5, p7, p6);
    M[6] = -A1pmp(p3, p4, p1, p2, p5, p7, p6);
    M[7] = -A2pmp(p3, p4, p1, p2, p5, p7, p6);
    swap();
    
    M[8] = M[9] = M[10] = M[11] = 0.0;
    M[12] = M[13] = M[14] = M[15] = 0.0;
  }
  
  void ampq4g1l2::
  matrix_tree_pmmpp(int p1, int p4, int p3, int p2,
		    int p5, int p6, int p7, _ComplexD *M)
  {
    M[0] = M[1] = M[2] = M[3] = 0;  
    M[4] = M[5] = M[6] = M[7] = 0;
    
    M[8]  = A1pmp(p1, p4, p3, p2, p5, p6, p7);
    M[9]  = A2pmp(p1, p4, p3, p2, p5, p6, p7);
    M[10] = 0.0;
    M[11] = A4pmp(p1, p4, p3, p2, p5, p6, p7);
    
    swap();
    M[12] = -A3pmm(p3, p2, p1, p4, p5, p7, p6);
    M[13] = 0.0;
    M[14] = -A1pmm(p3, p2, p1, p4, p5, p7, p6);
    M[15] = -A2pmm(p3, p2, p1, p4, p5, p7, p6);
    swap();
  }
  
  void ampq4g1l2::
  matrix_tree_pmmpm(int p1, int p4, int p3, int p2,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = M[1] = M[2] = M[3] = 0;  
    M[4] = M[5] = M[6] = M[7] = 0;
    
    M[8]  = A1pmm(p1, p4, p3, p2, p5, p6, p7);
    M[9]  = A2pmm(p1, p4, p3, p2, p5, p6, p7);
    M[10] = A3pmm(p1, p4, p3, p2, p5, p6, p7);
    M[11] = 0;
    
    swap();
    M[12] = 0;
    M[13] = -A4pmp(p3, p2, p1, p4, p5, p7, p6);
    M[14] = -A1pmp(p3, p2, p1, p4, p5, p7, p6);
    M[15] = -A2pmp(p3, p2, p1, p4, p5, p7, p6);
    swap();
  }
  
  double ampq4g1l2::su3_tree_aa(_ComplexD *A, _ComplexD *B) 
  {
    double D = real((A[0]+A[1])*conj(B[0]+B[1])+(A[2]+A[3])*conj(B[2]+B[3]));
    double E = real((A[0]+A[1])*conj(B[2]+B[3])+(A[2]+A[3])*conj(B[0]+B[1]));
    double F = real((A[0]+A[3])*conj(B[1]+B[2])+(A[1]+A[2])*conj(B[0]+B[3]));
    
    return 2.0*Na*(Na*D + (Nc2-2.0)*E - Nc2*F)/Nc;
  }
  
  double ampq4g1l2::su3_tree_ac(_ComplexD *A, _ComplexD *B) 
  {
    double D = real(A[0]*conj(B[0])+A[2]*conj(B[2]) + A[1]*conj(B[3])+A[3]*conj(B[1]));
    double E = real(A[1]*conj(B[1])+A[3]*conj(B[3]) + A[0]*conj(B[2])+A[2]*conj(B[0]));
    double F = real((A[0]+A[2])*conj(B[1]+B[3])     + (A[1]+A[3])*conj(B[0]+B[2]));
    
    return 2.0*Na*(-Na*D + (Nc2+1.0)*E + F)/Nc2;
  }
  
  void ampq4g1l2::su3_amptree(_ComplexD *M, double *Amp) 
  {
    Amp[0] +=     su3_tree_aa(M,   M);               //   A*A
    Amp[1] +=     su3_tree_aa(M+4, M+4);             //   B*B
    Amp[2] += 2.0*su3_tree_aa(M,   M+4);             // 2*A*B 
    
    Amp[3] +=     su3_tree_aa(M+8,  M+8);            //   C*C
    Amp[4] +=     su3_tree_aa(M+12, M+12);           //   D*D
    Amp[5] += 2.0*su3_tree_aa(M+8,  M+12);           // 2*C*D
    
    Amp[6] -= 2.0*su3_tree_ac(M,   M+8);             // 2*A*C
    Amp[7] -= 2.0*su3_tree_ac(M,   M+12);            // 2*A*D
    Amp[8] -= 2.0*su3_tree_ac(M+4, M+8);             // 2*B*C
    Amp[9] -= 2.0*su3_tree_ac(M+4, M+12);            // 2*B*D 
  }

  void ampq4g1l2::
  su3_tree(int p1, int p2, int p3, int p4, int p5, int p6, int p7, double *Amp) 
  {
    std::complex<double>  M[16];
    std::memset(Amp, 0, 10*sizeof(double));
    
    matrix_tree_pmpmp(p1,p2,p3,p4,p5,p6,p7, M); su3_amptree(M, Amp);
    matrix_tree_pmpmm(p1,p2,p3,p4,p5,p6,p7, M); su3_amptree(M, Amp);
    matrix_tree_ppmmp(p1,p2,p3,p4,p5,p6,p7, M); su3_amptree(M, Amp);
    matrix_tree_ppmmm(p1,p2,p3,p4,p5,p6,p7, M); su3_amptree(M, Amp);
    matrix_tree_pmmpp(p1,p2,p3,p4,p5,p6,p7, M); su3_amptree(M, Amp);
    matrix_tree_pmmpm(p1,p2,p3,p4,p5,p6,p7, M); su3_amptree(M, Amp);
    
    matrix_tree_pmpmp(p1,p2,p3,p4,p5,p7,p6, M); su3_amptree(M, Amp);
    matrix_tree_pmpmm(p1,p2,p3,p4,p5,p7,p6, M); su3_amptree(M, Amp);
    matrix_tree_ppmmp(p1,p2,p3,p4,p5,p7,p6, M); su3_amptree(M, Amp);
    matrix_tree_ppmmm(p1,p2,p3,p4,p5,p7,p6, M); su3_amptree(M, Amp);
    matrix_tree_pmmpp(p1,p2,p3,p4,p5,p7,p6, M); su3_amptree(M, Amp);
    matrix_tree_pmmpm(p1,p2,p3,p4,p5,p7,p6, M); su3_amptree(M, Amp);
  }
  void ampq4g1l2::
  su3_tree_mch(int p1, int p2, int p3, int p4, int p5, int p6, int p7, double *Amp) 
  {
    std::complex<double>  M[16];
    std::memset(Amp, 0, 10*sizeof(double));
    int hel = (int) (12*_M_rng());
 
    switch(hel){
    case 0:  matrix_tree_pmpmp(p1,p2,p3,p4,p5,p6,p7, M); break;
    case 1:  matrix_tree_pmpmm(p1,p2,p3,p4,p5,p6,p7, M); break;
    case 2:  matrix_tree_ppmmp(p1,p2,p3,p4,p5,p6,p7, M); break;
    case 3:  matrix_tree_ppmmm(p1,p2,p3,p4,p5,p6,p7, M); break;
    case 4:  matrix_tree_pmmpp(p1,p2,p3,p4,p5,p6,p7, M); break;
    case 5:  matrix_tree_pmmpm(p1,p2,p3,p4,p5,p6,p7, M); break;
    case 6:  matrix_tree_pmpmp(p1,p2,p3,p4,p5,p7,p6, M); break;
    case 7:  matrix_tree_pmpmm(p1,p2,p3,p4,p5,p7,p6, M); break;
    case 8:  matrix_tree_ppmmp(p1,p2,p3,p4,p5,p7,p6, M); break;
    case 9:  matrix_tree_ppmmm(p1,p2,p3,p4,p5,p7,p6, M); break;
    case 10: matrix_tree_pmmpp(p1,p2,p3,p4,p5,p7,p6, M); break;
    case 11: matrix_tree_pmmpm(p1,p2,p3,p4,p5,p7,p6, M); break;
    }

    su3_amptree(M, Amp); 
    for(unsigned int i = 0; i < 10; i++) Amp[i] *= 12.0;
  }
}  //  namespace nlo

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
#include "ampq2g3l2.h"
#include "defmacros.h"


namespace nlo {

  std::complex<double> 
  ampq2g3l2::Appp(int p1, int p2, int p3, int p4, int p5, int p6, int p7) {
    return A(6,5)*A3(5,6,7)/(A(1,2)*A(2,3)*A(3,4)*A(4,5)*S(6,7));
  }
  
  
  std::complex<double> 
  ampq2g3l2::Ammm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) {
    return -B(1,7)*A3(6,7,1)/(B(1,2)*B(2,3)*B(3,4)*B(4,5)*S(6,7));
  }

  std::complex<double> 
  ampq2g3l2::Appm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s34 = S(3,4), t534 = S3(5,3,4), t567 = S3(5,6,7), t234 = S3(2,3,4),
      t167 = S3(1,6,7); 
    
    _ComplexD a65 = A(6,5), a23 = A(2,3), a34 = A(3,4), b34 = B(3,4), 
      b23 = B(2,3), a12 = A(1,2), b45 = B(4,5), a45 = A(4,5), b53 = B(5,3),
      a24 = A(2,4), b17 = B(1,7), a42 = A(4,2);
    
    _ComplexD c4231 = C(4,2,3,1), c4123 = C(4,1,2,3), c6543 = C(6,5,4,3),
      c4127 = C(4,1,2,7), c6172 = C(6,1,7,2), c6173 = C(6,1,7,3),
      c4567 = C(4,5,6,7);
    
    _ComplexD m = a65*c4567/(a23*a34*b34*t567)*(b23*c4231/t234 + c4123/a12)
      + c4567*c6543/(a12*a23*a34*b34*b45)
      - c4127*c6543*a45*b53/(a12*a24*b45*s34*t534)
      - b17*c6172*a45*a45*b53*b53/(b45*a42*s34*t534*t167)
      + b17*a45*b53/(a23*b34*b45*t167)*(c6172/a34 + c6173/a24)
      - b17*a45*b23/(a23*b34*t234*t167)*(c6172*a24/a34 + c6173);
    
    return m/S(6,7);
  }
  
  std::complex<double> 
  ampq2g3l2::Apmp(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s23 = S(2,3), s34 = S(3,4), t234 = S3(2,3,4), t567 = S3(5,6,7),
      t123 = S3(1,2,3), t345 = S3(3,4,5), t167 = S3(1,6,7);
    
    _ComplexD b42 = B(4,2), a65 = A(6,5), b23 = B(2,3), a12 = A(1,2), 
      a23 = A(2,3), a31 = A(3,1), b12 = B(1,2), a34 = A(3,4), a35 = A(3,5),
      a45 = A(4,5), b17 = B(1,7), b54 = B(5,4), a43 = A(4,3);
    
    _ComplexD c3124 = C(3,1,2,4), c3127 = C(3,1,2,7), c6172 = C(6,1,7,2), 
      c6534 = C(6,5,3,4), c3542 = C(3,5,4,2), c6174 = C(6,1,7,4),
      c3241 = C(3,2,4,1), c3567 = C(3,5,6,7); 
    
    _ComplexD m = b42*b42*c3241*c3567*a65/(s23*s34*t234*t567)
      + b42*c3124*c3567*a65/(b23*a12*a23*s34*t567)
      + a31*b12*c3124*c3567*a65/(a12*a34*s23*t123*t567)
      + a31*b12*c3127*a65*a35/(a12*a34*a45*s23*t123)
      - c3127*c6172*a35*a35/(a12*a23*a34*a45*b23*t345)
      + c3127*c6534*a35*b42/(b23*a12*a23*s34*t345)
      + b17*c6172*a35*a35/(s23*t345*t167)*(c3542/(a34*a45) + b42*b54/s34)
      - b42*b42*b17/(s23*s34*t234*t167)*(c6172*a23*a35 + c6174*a43*a35);
    
    return m/S(6,7);
  }
  
  
  std::complex<double> 
  ampq2g3l2::Ampp(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s23 = S(2,3), s34 = S(3,4), t234 = S3(2,3,4), t567 = S3(5,6,7),
      t345 = S3(3,4,5), t167 = S3(1,6,7), t123 = S3(1,2,3); 
    
    _ComplexD b43 = B(4,3), a65 = A(6,5), b13 = B(1,3), b12 = B(1,2), 
      a34 = A(3,4), a42 = A(4,2), a24 = A(2,4), a25 = A(2,5), a45 = A(4,5),
      b17 = B(1,7), a32 = A(3,2);
    
    _ComplexD c2134 = C(2,1,3,4), c6534 = C(6,5,3,4), c2137 = C(2,1,3,7), 
      c6543 = C(6,5,4,3), c6173 = C(6,1,7,3), c2543 = C(2,5,4,3),
      c2341 = C(2,3,4,1), c2567 = C(2,5,6,7), c2534 = C(2,5,3,4),
      c6174 = C(6,1,7,4);
    
    _ComplexD m = b43*b43*c2341*c2567*a65/(s23*s34*t234*t567)
      + b13*c2341*c2567*a65/(b12*a34*a42*s23*t567)
      - b13*b13*c2134*c2567*a65/(b12*a24*s23*t123*t567)
      - b13*b13*c2137*a65*a25/(b12*a24*a45*s23*t123)
      + b13*b17*a25/(b12*s23*a24*t345)*(c6543*(a25/a45-a32/a34)-c6534*a42/a34)
      + b17*c6173*a25/(a24*s23*t345*t167)
      *(c2543*(a25/a45-a32/a34)-c2534*a42/a34)
      - b17*a25*b43*b43/(s23*s34*t234*t167)*(c6173*a32 + c6174*a42);
    
    return m/S(6,7);
  }

  std::complex<double> 
  ampq2g3l2::Apmm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s23 = S(2,3), s34 = S(3,4), t234 = S3(2,3,4), t167 = S3(1,6,7), 
      t123 = S3(1,2,3), t567 = S3(5,6,7);
    
    _ComplexD b12 = B(1,2), b23 = B(2,3), b24 = B(2,4), a65 = A(6,5), 
      a43 = A(4,3), a31 = A(3,1), b34 = B(3,4), b42 = B(4,2), a12 = A(1,2),
      b45 = B(4,5), b17 = B(1,7); 
    
    _ComplexD c3127 = C(3,1,2,7), c4567 = C(4,5,6,7), c3567 = C(3,5,6,7), 
      c6172 = C(6,1,7,2), c6542 = C(6,5,4,2), c3542 = C(3,5,4,2), 
      c5342 = C(5,3,4,2);
    
    _ComplexD m = b12*(b23*c3567 + b24*c4567)*a65/(s23*t567)
      *(a43*a43/(s34*t234) - a31/(b34*b42*a12))
      - a31*a31*b12*b12*c4567*a65/(a12*b24*s23*t123*t567)
      + a31*b12*c3127*c6542/(a12*b24*b45*s23*t123)  
      + c3127*c6172/(a12*b34*b45*s23)
      + b17*c6172*c3542/(b34*b45*s23*t167)
      - b17*c6172*c5342*a43*a43/(s23*s34*t234*t167);
    
    return m/S(6,7);
  }
  

  std::complex<double>
  ampq2g3l2::Ampm(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s23 = S(2,3), s34 = S(3,4), t234 = S3(2,3,4), t123 = S3(1,2,3),
      t345 = S3(3,4,5), t167 = S3(1,6,7), t567 = S3(5,6,7);
    
    _ComplexD a42 = A(4,2), b13 = B(1,3), b32 = B(3,2), b34 = B(3,4), 
      a65 = A(6,5), b12 = B(1,2), a21 = A(2,1), b45 = B(4,5), b17 = B(1,7),
      a45 = A(4,5);
    
    _ComplexD c2567 = C(2,5,6,7), c4567 = C(4,5,6,7), c2137 = C(2,1,3,7), 
      c6543 = C(6,5,4,3), c2543 = C(2,5,4,3), c6173 = C(6,1,7,3),
      c5243 = C(5,2,4,3);
    
    _ComplexD m = a42*a42*b13*(b32*c2567 + b34*c4567)*a65/(s23*s34*t234*t567)
      - a42*b13*b13*c4567*a65/(b12*s23*s34*t567)
      + b13*b13*b13*a21*c4567*a65/(b12*b34*s23*t123*t567)
      - b13*b13*c2137*c6543/(b12*b34*b45*s23*t123)
      - b13*b17*c6543/(b12*s23*t345)*(a45*a42/s34 - c2543/(b34*b45))
      - b17*c6173*c2543/(s23*t345*t167)*(a45*a42/s34 - c2543/(b34*b45))
      - b17*c6173*a42*a42*c5243/(s23*s34*t234*t167);
    
    return m/S(6,7);
  }

  std::complex<double> 
  ampq2g3l2::Ammp(int p1, int p2, int p3, int p4, int p5, int p6, int p7) 
  {
    double s23 = S(2,3), s34 = S(3,4), t234 = S3(2,3,4), t567 = S3(5,6,7), 
      t123 = S3(1,2,3), t345 = S3(3,4,5), t167 = S3(1,6,7); 
    
    _ComplexD a23 = A(2,3), b14 = B(1,4), b42 = B(4,2), b43 = B(4,3),
      a65 = A(6,5), b23 = B(2,3), b12 = B(1,2), a35 = A(3,5), a45 = A(4,5), 
      b17 = B(1,7), b24 = B(2,4); 
    
    _ComplexD c2567 = C(2,5,6,7), c3567 = C(3,5,6,7), c2134 = C(2,1,3,4),
      c3124 = C(3,1,2,4), c2137 = C(2,1,3,7), c3127 = C(3,1,2,7),
      c6534 = C(6,5,3,4), c6174 = C(6,1,7,4), c2534 = C(2,5,3,4),
      c5234 = C(5,2,3,4);
    
    _ComplexD m = a23*a23*b14*(b42*c2567 + b43*c3567)*a65/(s23*s34*t234*t567)
      - b14*c3567*a65/(b42*s34*t123*t567)
      *((b42*c2134 + b43*c3124)/b23 - b14*c3124/b12)
      - b14*a65*a35/(b42*a45*s34*t123)
      *((b42*c2137 + b43*c3127)/b23 - b14*c3127/b12)
      + b14*b17*c6534*a35*a35/(b12*b24*a45*s34*t345)
      - b17*c6174*c2534*a35*a35/(b42*a45*s34*t345*t167)
      - b17*c6174*c5234*a35/(a45*b42*b23*s34*t167)
      - b17*c6174*c5234*a23*a23/(s23*s34*t234*t167);
    
    return m/S(6,7);
  }

  void ampq2g3l2::
  matrix_tree_ppppm(int p1, int p2, int p3, int p4, 
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = Appp(p1, p2, p3, p4, p5, p6, p7);
    M[1] = Appp(p1, p3, p4, p2, p5, p6, p7);
    M[2] = Appp(p1, p4, p2, p3, p5, p6, p7);
    M[3] = Appp(p1, p2, p4, p3, p5, p6, p7);
    M[4] = Appp(p1, p4, p3, p2, p5, p6, p7);
    M[5] = Appp(p1, p3, p2, p4, p5, p6, p7);
  }
  
  void ampq2g3l2::
  matrix_tree_pmmmm(int p1, int p2, int p3, int p4,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = Ammm(p1, p2, p3, p4, p5, p6, p7);
    M[1] = Ammm(p1, p3, p4, p2, p5, p6, p7);
    M[2] = Ammm(p1, p4, p2, p3, p5, p6, p7);
    M[3] = Ammm(p1, p2, p4, p3, p5, p6, p7);
    M[4] = Ammm(p1, p4, p3, p2, p5, p6, p7);
    M[5] = Ammm(p1, p3, p2, p4, p5, p6, p7);
  }
  
  void ampq2g3l2::
  matrix_tree_pppmm(int p1, int p2, int p3, int p4,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = Appm(p1, p2, p3, p4, p5, p6, p7);
    M[1] = Apmp(p1, p3, p4, p2, p5, p6, p7);
    M[2] = Ampp(p1, p4, p2, p3, p5, p6, p7);
    M[3] = Apmp(p1, p2, p4, p3, p5, p6, p7);
    M[4] = Ampp(p1, p4, p3, p2, p5, p6, p7);
    M[5] = Appm(p1, p3, p2, p4, p5, p6, p7);
  }
  
  void ampq2g3l2::
  matrix_tree_ppmpm(int p1, int p2, int p3, int p4,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = Apmp(p1, p2, p3, p4, p5, p6, p7);
    M[1] = Ampp(p1, p3, p4, p2, p5, p6, p7);
    M[2] = Appm(p1, p4, p2, p3, p5, p6, p7);
    M[3] = Appm(p1, p2, p4, p3, p5, p6, p7);
    M[4] = Apmp(p1, p4, p3, p2, p5, p6, p7);
    M[5] = Ampp(p1, p3, p2, p4, p5, p6, p7);
  }
  
  void ampq2g3l2::
  matrix_tree_pmppm(int p1, int p2, int p3, int p4,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = Ampp(p1, p2, p3, p4, p5, p6, p7);
    M[1] = Appm(p1, p3, p4, p2, p5, p6, p7);
    M[2] = Apmp(p1, p4, p2, p3, p5, p6, p7);
    M[3] = Ampp(p1, p2, p4, p3, p5, p6, p7);
    M[4] = Appm(p1, p4, p3, p2, p5, p6, p7);
    M[5] = Apmp(p1, p3, p2, p4, p5, p6, p7);
  }
  
  void ampq2g3l2::
  matrix_tree_ppmmm(int p1, int p2, int p3, int p4,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = Apmm(p1, p2, p3, p4, p5, p6, p7);
    M[1] = Ammp(p1, p3, p4, p2, p5, p6, p7);
    M[2] = Ampm(p1, p4, p2, p3, p5, p6, p7);
    M[3] = Apmm(p1, p2, p4, p3, p5, p6, p7);
    M[4] = Ammp(p1, p4, p3, p2, p5, p6, p7);
    M[5] = Ampm(p1, p3, p2, p4, p5, p6, p7);
  }
  
  void ampq2g3l2::
  matrix_tree_pmpmm(int p1, int p2, int p3, int p4,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = Ampm(p1, p2, p3, p4, p5, p6, p7);
    M[1] = Apmm(p1, p3, p4, p2, p5, p6, p7);
    M[2] = Ammp(p1, p4, p2, p3, p5, p6, p7);
    M[3] = Ammp(p1, p2, p4, p3, p5, p6, p7);
    M[4] = Ampm(p1, p4, p3, p2, p5, p6, p7);
    M[5] = Apmm(p1, p3, p2, p4, p5, p6, p7);
  }
  
  void ampq2g3l2::
  matrix_tree_pmmpm(int p1, int p2, int p3, int p4,
		    int p5, int p6, int p7, _ComplexD *M) 
  {
    M[0] = Ammp(p1, p2, p3, p4, p5, p6, p7);
    M[1] = Ampm(p1, p3, p4, p2, p5, p6, p7);
    M[2] = Apmm(p1, p4, p2, p3, p5, p6, p7);
    M[3] = Ampm(p1, p2, p4, p3, p5, p6, p7);
    M[4] = Apmm(p1, p4, p3, p2, p5, p6, p7);
    M[5] = Ammp(p1, p3, p2, p4, p5, p6, p7);
  }
 
  double ampq2g3l2::su3_amptree(_ComplexD *m) 
  {
    _ComplexD temp = m[0]+m[1]+m[2]+m[3]+m[4]+m[5];
    
    double A0 = real(temp*conj(temp));
    double A2 = real(m[0]*conj(m[0])+m[1]*conj(m[1])+m[2]*conj(m[2])
		     + m[3]*conj(m[3])+m[4]*conj(m[4])+m[5]*conj(m[5]));
    double A1 = -2.0*A2 - 2.0*real(m[0]*conj(m[3]+m[5]-m[4])
				   + m[1]*conj(m[5]+m[4]-m[3])
				   + m[2]*conj(m[4]+m[3]-m[5]));
  
    return Na*(Nc2*A2 + A1 + A0/Nc2);
  }
  
  double ampq2g3l2::
  su3_tree(int p1,int p2,int p3,int p4,int p5,int p6,int p7)
  {
    double A;
    _ComplexD M[6];
        
    matrix_tree_ppppm(p1, p2, p3, p4, p5, p6, p7, M); A  = su3_amptree(M);
    matrix_tree_pppmm(p1, p2, p3, p4, p5, p6, p7, M); A += su3_amptree(M);
    matrix_tree_ppmpm(p1, p2, p3, p4, p5, p6, p7, M); A += su3_amptree(M);
    matrix_tree_pmppm(p1, p2, p3, p4, p5, p6, p7, M); A += su3_amptree(M);
    matrix_tree_ppmmm(p1, p2, p3, p4, p5, p6, p7, M); A += su3_amptree(M);
    matrix_tree_pmpmm(p1, p2, p3, p4, p5, p6, p7, M); A += su3_amptree(M);
    matrix_tree_pmmpm(p1, p2, p3, p4, p5, p6, p7, M); A += su3_amptree(M);
    matrix_tree_pmmmm(p1, p2, p3, p4, p5, p6, p7, M); A += su3_amptree(M);
    
    matrix_tree_ppppm(p1, p2, p3, p4, p5, p7, p6, M); A += su3_amptree(M);
    matrix_tree_pppmm(p1, p2, p3, p4, p5, p7, p6, M); A += su3_amptree(M);
    matrix_tree_ppmpm(p1, p2, p3, p4, p5, p7, p6, M); A += su3_amptree(M);
    matrix_tree_pmppm(p1, p2, p3, p4, p5, p7, p6, M); A += su3_amptree(M);
    matrix_tree_ppmmm(p1, p2, p3, p4, p5, p7, p6, M); A += su3_amptree(M);
    matrix_tree_pmpmm(p1, p2, p3, p4, p5, p7, p6, M); A += su3_amptree(M);
    matrix_tree_pmmpm(p1, p2, p3, p4, p5, p7, p6, M); A += su3_amptree(M);
    matrix_tree_pmmmm(p1, p2, p3, p4, p5, p7, p6, M); A += su3_amptree(M);
    
    return 2.0*A;
  }

  double ampq2g3l2::
  su3_tree_mch(int p1,int p2,int p3,int p4,int p5,int p6,int p7)
  {
    _ComplexD M[6];
    int hel = (int) (16*_M_rng());
      
    switch(hel){
    case 0:  matrix_tree_ppppm(p1, p2, p3, p4, p5, p6, p7, M); break;
    case 1:  matrix_tree_pppmm(p1, p2, p3, p4, p5, p6, p7, M); break;
    case 2:  matrix_tree_ppmpm(p1, p2, p3, p4, p5, p6, p7, M); break;
    case 3:  matrix_tree_pmppm(p1, p2, p3, p4, p5, p6, p7, M); break;
    case 4:  matrix_tree_ppmmm(p1, p2, p3, p4, p5, p6, p7, M); break;
    case 5:  matrix_tree_pmpmm(p1, p2, p3, p4, p5, p6, p7, M); break;
    case 6:  matrix_tree_pmmpm(p1, p2, p3, p4, p5, p6, p7, M); break;
    case 7:  matrix_tree_pmmmm(p1, p2, p3, p4, p5, p6, p7, M); break;
    case 8:  matrix_tree_ppppm(p1, p2, p3, p4, p5, p7, p6, M); break;
    case 9:  matrix_tree_pppmm(p1, p2, p3, p4, p5, p7, p6, M); break;
    case 10: matrix_tree_ppmpm(p1, p2, p3, p4, p5, p7, p6, M); break;
    case 11: matrix_tree_pmppm(p1, p2, p3, p4, p5, p7, p6, M); break;
    case 12: matrix_tree_ppmmm(p1, p2, p3, p4, p5, p7, p6, M); break;
    case 13: matrix_tree_pmpmm(p1, p2, p3, p4, p5, p7, p6, M); break;
    case 14: matrix_tree_pmmpm(p1, p2, p3, p4, p5, p7, p6, M); break;
    case 15: matrix_tree_pmmmm(p1, p2, p3, p4, p5, p7, p6, M); break;
    }
    return 32.0*su3_amptree(M);
  }
}   //   namespace nlo

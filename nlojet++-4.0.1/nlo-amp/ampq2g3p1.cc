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
#include "ampq2g3p1.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {


  std::complex<double> 
  ampq2g3p1::Amhv(int pi, int p1, int p2, int p3, int p4, int p5,int p6) const 
  {
    _ComplexD a51 = A(5,1), a16 = A(1,6), a23 = A(2,3), 
      a34 = A(3,4), a45 = A(4,5), a62 = A(6,2), 
      a1i = Aij(p1,pi), a2i = Aij(p2,pi);

    return -ICPLX*a1i*std::pow(a2i,3)/(a23*a34*a45*a51*a16*a62);
  }
 

  std::complex<double>
  ampq2g3p1::Apmmmpp(int p5, int p6, int p1, int p2, int p3, int p4) const 
  {
    double s23 = S(2,3), s46 = S(4,6), s123 = S3(1,2,3), s126 = S3(1,2,6), 
      s146 = S3(1,4,6);
    
    _ComplexD a12 = A(1,2), a35 = A(3,5), a45 = A(4,5), b12 = B(1,2), 
      b16 = B(1,6), b45 = B(4,5), d3125 = D(3,1,2,5), d3126 = D(3,1,2,6), 
      d3251 = D(3,2,5,1), d4352 = D(4,3,5,2); 

    return ICPLX*(d3125*d4352*s126/(b12*a45*b16*a35*s23*s46)
		  + d3251*d4352*d4352/(a35*b16*s23*s46*s146)
		  - a12*b45*d3126*d3125/(b12*a45*s23*s46*s123));
  }
  
  std::complex<double>
  ampq2g3p1::Apmpmmp(int p5, int p6, int p1, int p2, int p3, int p4) const 
  {
    double s12 = S(1,2), s16 = S(1,6), s35 = S(3,5), s123 = S3(1,2,3),
      s126 = S3(1,2,6), s146 = S3(1,4,6);
  
    _ComplexD a12 = A(1,2), a23 = A(2,3), a26 = A(2,6), a46 = A(4,6), 
      a54 = A(5,4), b23 = B(2,3), b41 = B(4,1), d1236 = D(1,2,3,6), 
      d1263 = D(1,2,6,3), d4352 = D(4,3,5,2), d5146 = D(5,1,4,6); 
    
    return ICPLX*(a23*d1236*d1236/(b23*a54*a46*s12*s123)
		  + a26*d1236*d1263*s146/(b23*a54*a46*s12*s16*s35)
		  + a23*a26*b41*d5146/(b23*a12*a46*s16*s35)
		  + a26*d1263*d1263*d4352/(a54*s12*s16*s35*s126));
  }
  
  std::complex<double>
  ampq2g3p1::Apmmpmp(int p5, int p6, int p1, int p2, int p3, int p4) const 
  {
    double s12 = S(1,2), s23 = S(2,3), s35 = S(3,5), s123 = S3(1,2,3), 
      s126 = S3(1,2,6), s146 = S3(1,4,6); 
    
    _ComplexD a13 = A(1,3), a16 = A(1,6), a21 = A(2,1), a46 = A(4,6),
      a53 = A(5,3), a54 = A(5,4), a61 = A(6,1), b16 = B(1,6), b23 = B(2,3), 
      b25 = B(2,5), b26 = B(2,6), b45 = B(4,5), b53 = B(5,3), 
      d2135 = D(2,1,3,5), d2136 = D(2,1,3,6), d2163 = D(2,1,6,3),
      d2453 = D(2,4,5,3), d4253 = D(4,2,5,3);
    
    return ICPLX*(-a13*a13*d2136*d2136/(a54*a46*s12*s23*s123)
		  + b26*b45*a61*a53*d2453/(b16*a54*s12*s35*s126)
		  + b25*b25*a61*d4253/(b16*b53*a46*s23*s146)
		  + a16*b25*d2135*d2163/(b16*b23*a54*a46*s12*s35)
		  + a13*d2136*d2163/(b16*a54*a46*s12*s23)
		  + a16*b25*a13*d4253/(b16*a46*a21*s23*s35));
  }

#define MatrixElement(h3,h4,h5,h6)				\
  void ampq2g3p1::						\
  matrix_tree_pm##h3##h4##h5##h6(int p1, int p2, int p3, 	\
				 int p4, int p5, int p6, 	\
				 _ComplexD *M) const 		\
  {								\
    M[0] = Apm##h3##h4##h5##h6(p1,p2, p3,p4,p5, p6);		\
    M[1] = Apm##h4##h5##h3##h6(p1,p2, p4,p5,p3, p6);		\
    M[2] = Apm##h5##h3##h4##h6(p1,p2, p5,p3,p4, p6);		\
    M[3] = Apm##h3##h5##h4##h6(p1,p2, p3,p5,p4, p6);		\
    M[4] = Apm##h5##h4##h3##h6(p1,p2, p5,p4,p3, p6);		\
    M[5] = Apm##h4##h3##h5##h6(p1,p2, p4,p3,p5, p6);		\
  } 

MatrixElement(m,m,p,p)
MatrixElement(p,m,m,p)
MatrixElement(m,p,m,p) 
MatrixElement(m,p,p,p)
MatrixElement(p,m,p,p)
MatrixElement(p,p,m,p)
MatrixElement(m,m,m,p)
MatrixElement(p,p,p,m)
MatrixElement(m,m,p,m)
MatrixElement(m,p,m,m)
MatrixElement(p,m,m,m)
MatrixElement(p,p,m,m)
MatrixElement(p,m,p,m)
MatrixElement(m,p,p,m)
  
  double ampq2g3p1::amptree(_ComplexD *m) 
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
  
  double ampq2g3p1::su3_tree(int p1,int p2,int p3,int p4,int p5,int p6) const
  {
    _ComplexD M[6];
    double A;
   
    matrix_tree_pmpmmm(p1, p2, p3, p4, p5, p6, M); A  = amptree(M);
    matrix_tree_pmmpmm(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmmmpm(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmmmmp(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
 
    matrix_tree_pmpppm(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmppmp(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmpmpp(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmmppp(p1, p2, p3, p4, p5, p6, M); A += amptree(M);

    matrix_tree_pmppmm(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmpmpm(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmpmmp(p1, p2, p3, p4, p5, p6, M); A += amptree(M);

    matrix_tree_pmmppm(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmmpmp(p1, p2, p3, p4, p5, p6, M); A += amptree(M);
    matrix_tree_pmmmpp(p1, p2, p3, p4, p5, p6, M); A += amptree(M);

    return 4.0*A;
  }

  double ampq2g3p1::
  su3_tree_mch(int p1,int p2, int p3, int p4, int p5, int p6) const
  {
    _ComplexD M[6];
    unsigned int hel = (unsigned int) (14*_M_rng());

    switch(hel) {
    case 0:  matrix_tree_pmmmpp(p1, p2, p3, p4, p5, p6, M); break;
    case 1:  matrix_tree_pmpmmp(p1, p2, p3, p4, p5, p6, M); break;
    case 2:  matrix_tree_pmmpmp(p1, p2, p3, p4, p5, p6, M); break;
    case 3:  matrix_tree_pmmppp(p1, p2, p3, p4, p5, p6, M); break;
    case 4:  matrix_tree_pmpmpp(p1, p2, p3, p4, p5, p6, M); break;
    case 5:  matrix_tree_pmppmp(p1, p2, p3, p4, p5, p6, M); break;
    case 6:  matrix_tree_pmmmmp(p1, p2, p3, p4, p5, p6, M); break;
    case 7:  matrix_tree_pmpppm(p1, p2, p3, p4, p5, p6, M); break;
    case 8:  matrix_tree_pmmmpm(p1, p2, p3, p4, p5, p6, M); break;
    case 9:  matrix_tree_pmmpmm(p1, p2, p3, p4, p5, p6, M); break;
    case 10: matrix_tree_pmpmmm(p1, p2, p3, p4, p5, p6, M); break;
    case 11: matrix_tree_pmppmm(p1, p2, p3, p4, p5, p6, M); break;
    case 12: matrix_tree_pmpmpm(p1, p2, p3, p4, p5, p6, M); break;
    case 13: matrix_tree_pmmppm(p1, p2, p3, p4, p5, p6, M); break;
    }
    
    return 56.0*amptree(M);
  }


}















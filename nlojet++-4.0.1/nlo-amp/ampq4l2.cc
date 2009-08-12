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
#include "ampq4l2.h"
#include "defmacros.h"



namespace nlo {

  std::complex<double> 
  ampq4l2::App(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
    _ComplexD Atree = _ComplexD(0.0,1.0)
      *(B(1,2)*A(5,4)*C(3,1,2,6)/S3(1,2,3)
	+ A(3,4)*B(6,1)*C(5,3,4,2)/S3(2,3,4))/(S(2,3)*S(5,6));

    return Atree;
  }
  
  std::complex<double> 
  ampq4l2::Apm(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
    _ComplexD Atree = _ComplexD(0.0,-1.0)
      *(B(1,3)*A(5,4)*C(2,1,3,6)/S3(1,2,3)
	+ A(2,4)*B(6,1)*C(5,2,4,3)/S3(2,3,4))/(S(2,3)*S(5,6));

    return Atree;
  }

  std::complex<double> 
  ampq4l2::FApp(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s56 = S(5,6),
      t123 = S3(1,2,3); 
    
    _ComplexD a23 = A(2,3), b56 = B(5,6), b12 = B(1,2), a45 = A(4,5), 
      b23 = B(2,3), a56 = A(5,6);
    
    _ComplexD a346 = A3(3,4,6), a321 = A3(3,2,1), c3126 = C(3,1,2,6), 
      c1234 = C(1,2,3,4), c4231 = C(4,2,3,1), c1236 = C(1,2,3,6); 
    
    _ComplexD F = (c3126*c3126/(a23*b56*t123*c1234)
		   - b12*b12*a45*a45/(b23*a56*t123*c4231))
      *(Ls_1(-s12, -t123, -s23, -t123) + tLs_1_2mh(s34, t123, s12, s56))
      
      - 2.0*c3126/(b56*c1234)*(c1236*b12/t123*L0(-s23, -t123)/t123
			       + a346/a23*L0(-s56, -t123)/t123)
      
    - 0.5/(a23*b56*t123*c1234)
      *(a321*c1236*a321*c1236*L1(-t123, -s23)/(s23*s23)
	+ a346*t123*a346*t123*L1(-s56, -t123)/(t123*t123));
    
    return F;
  }

  std::complex<double> 
  ampq4l2::FApm(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
    double g3 = G3, s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s56 = S(5,6),
      t123 = S3(1,2,3), t124 = S3(1,2,4);
    
    _ComplexD b13 = B(1,3), a45 = A(4,5), b23 = B(2,3), a56 = A(5,6), 
      a12 = A(1,2), a23 = A(2,3), b56 = B(5,6), a13 = A(1,3), b64 = B(6,4),
      a42 = A(4,2), b16 = B(1,6), b26 = B(2,6), b12 = B(1,2), a25 = A(2,5), 
      a15 = A(1,5);
    
    _ComplexD a214 = A3(2,1,4), a536 = A3(5,3,6), a546 = A3(5,4,6), 
      a234 = A3(2,3,4), a213 = A3(2,1,3);
    
    _ComplexD c1236 = C(1,2,3,6), c4231 = C(4,2,3,1), c1234 = C(1,2,3,4), 
      c2341 = C(2,3,4,1), c1342 = C(1,3,4,2), c5126 = C(5,1,2,6), 
      c5346 = C(5,3,4,6), c2136 = C(2,1,3,6), c3126 = C(3,1,2,6), 
      c2134 = C(2,1,3,4), c3124 = C(3,1,2,4);
    
    _ComplexD F = (- b13*a45*b13*a45/(b23*a56*t123*c4231)
		   + a12*c3126*a12*c3126/(a23*b56*t123*a13*a13*c1234)
		   )*Ls_1(-s12, -t123, -s23, -t123)
      
      + (- b13*a45*b13*a45/(b23*a56*t123*c4231)
	 + c3126*c2134*c3126*c2134/(a23*b56*t123*c1234*c3124*c3124)
	 )*Ls_1_2mh(s34, t123, s12, s56)
      
      + (0.5*D56*c5126*(s56*a214 - t123*c2134)/(c1234*c3124*g3)
	 - 0.5*(a214*(a536 - a546) + a234*c5346)/(c1234*c3124)
	 + 0.5*(t123*D34 + 2.0*s12*s56)/(t123*t123)
	 *(b13*a45*b13*a45/(b23*a56*c4231) - c2136*c2136/(a23*b56*c1234))
       - 2.0*a213*a546/(t123*c1234))*I3_3m(s12, s34, s56)

      + b13*a12*c3126*c3126/(b56*t123*a13*c3124)*L0(-t123, -s12)/s12
      - 0.5*b13*b13*a23*c1236*c1236/(b56*t123*c1234)*L1(-t123, -s23)/(s23*s23)
      + b13*c1236*c2136/(b56*t123*c1234)*L0(-t123, -s23)/s23
      + b13*a12*c1236*c3126/(b56*t123*a13*c1234)*L0(-t123, -s23)/s23
      - 0.5*b64*a42*b64*a42*t123/(a23*b56*c1234)*L1(-s56, -t123)/(t123*t123)
      - b64*b64*a42*t123/(b56*c1234*c3124)*L0(-t123, -s56)/s56
      - 2.0*b64*a42*c2136/(a23*b56*c1234)*L0(-t123, -s56)/s56
      
      + 1.0/(c3124*g3)
      *(a12*t123/b56*(b16*b16 + b26*b26*c2341/c1342)
	+ b12*t124/a56*(a25*a25 + a15*a15*c2341/c1342)
	+ 2.0*s12*(-a25*b16 + a15*b26*c2341/c1342))*(Log(-s12) - Log(-s56));
  
    return F;
  }
  
  std::complex<double> 
  ampq4l2::FCsl(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
    double g3 = G3, d34 = D34, d56 = D56, s34 = S(3,4), s12 = S(1,2),
      s56 = S(5,6), t123 = S3(1,2,3), t124 = S3(1,2,4);
    
    _ComplexD b64 = B(6,4), a42 = A(4,2), a12 = A(1,2), b56 = B(5,6), 
      a23 = A(2,3), a24 = A(2,4), b12 = B(1,2), b46 = B(4,6), a25 = A(2,5),
      b16 = B(1,6);
    
    _ComplexD a236 = A3(2,3,6), a246 = A3(2,4,6), c2341 = C(2,3,4,1), 
      c5346 = C(5,3,4,6), c3124 = C(3,1,2,4), c2136 = C(2,1,3,6); 
    
    _ComplexD F = 0.5*b64*a42*b64*a42*t123/(a12*b56*c3124)
      *L1(-s56, -t123)/(t123*t123)
      + 2.0*b64*a42*c2136/(a12*b56*c3124)*L0(-t123, -s56)/s56
      - a23*a24*b64*b64*t123/(a12*b56*c3124*c3124)*L0(- t123, -s56)/s56
      - (0.5*a23*b64*c2136/(a12*b56*c3124*c3124)
	 + 0.75*c2136*c2136/(a12*b56*t123*c3124))
      *(Log(-t123) + Log(-s34) - 2.0*Log(-s56))
      
      + (1.5*d56*(t123 - t124)*c2341*c5346/(c3124*g3*g3)
	 - b12*a23*b46/(c3124*c3124*g3)*(a25*(t123 - t124) - 2.0*A4(2,1,6,5))
	 + b12*a25*(a236 - a246)/(c3124*g3)
	 + b16/(b56*c3124*g3)*(a236*t123 - a246*t124 + a23*b46*d34*t123/c3124) 
	 )*(Log(-s12) - Log(-s34))
      
      - 0.25*b16*(t123 - t124)*(b16*d34 - 2.0*B4(1,2,5,6))/(b12*b56*c3124*g3);
    
    return F;
  }

  std::complex<double> 
  ampq4l2::F1sl(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
    double g3 = G3, d56 = D56, d12 = D12, d34 = D34, s12 = S(1,2),
      s56 = S(5,6), s34 = S(3,4), t123 = S3(1,2,3), t124 = S3(1,2,4);
    
    _ComplexD a25 = A(2,5), b16 = B(1,6), b12 = B(1,2), a56 = A(5,6),
      a23 = A(2,3), b46 = B(4,6), b13 = B(1,3), a45 = A(4,5), a12 = A(1,2), 
      b56 = B(5,6);
    
    _ComplexD a521 = A3(5,2,1), a561 = A3(5,6,1), c2341 = C(2,3,4,1),
      c2134 = C(2,1,3,4), c3126 = C(3,1,2,6), c4123 = C(4,1,2,3),
      c5231 = C(5,2,3,1), c2346 = C(2,3,4,6), c5346 = C(5,3,4,6),
      c3124 = C(3,1,2,4);
    
    _ComplexD T = 0.5*(3*s34*d34 - g3)*(t123 - t124)*c2341*c5346/(c3124*g3*g3)
      + 0.5*s34*(t123 - t124)*a25*b16/(c3124*g3)
      - b12*a56*c2346*c2346/(c3124*g3)
      + a23*b46*d34*(a521*d12 - a561*d56)/(c3124*c3124*g3)
      + a23*b46*c5231/(c3124*c3124)
      - 2.0*b13*a45*c2134*c3126/(t123*t123*c3124);
  
    _ComplexD F = (b13*b13*a45*a45/(b12*a56*t123*c4123)
		   - c3126*c3126*c2134*c2134/(a12*b56*t123*c3124*c3124*c3124))
      *Ls_1_2mh(s34, t123, s12, s56)
      
      + T*I3_3m(s12, s34, s56) 
      + FCsl(p1, p2, p3, p4, p5, p6);
    
    swap();
    F -= FCsl(p5, p6, p4, p3, p1, p2);
    swap();
    
    return F;
  }
  
  std::complex<double>
  ampq4l2::Aax(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
    double s12 = S(1,2), s34 = S(3,4), s56 = S(5,6);  
    return _ComplexD(0.0, -1.233700550136169)
      *(f(MASS_OF_TOP, s12,s34,s56)-f(0, s12,s34,s56))/s56
      *(B(6,3)*A(4,2)*A(2,5)/A(1,2)-B(6,1)*B(1,3)*A(4,5)/B(1,2));
  }

  void ampq4l2::
  su3_m1_ppmm(unsigned int nf, int p1, int p2, int p3, 
	      int p4, int p5, int p6, _ComplexD *M) const 
  {
    _ComplexD I(0.0,1.0), atree, app, asl, apm, vpp, vsl, asf;
    double __nf = (double) nf;
    
    atree = App(p1, p2, p3, p4, p5, p6);
    vpp   = -13.0*(Log(-S(2,3))-Log(-S(5,6)))/6.0 + 10.0/9.0;
    vsl   = -1.5*(2.0*Log(-S(5,6))-Log(-S(1,4))-Log(-S(2,3))) - 7.5;
    
    asl = -vsl*atree + I*Fsl(p2, p3, p1, p4, p5, p6);
    app =  vpp*atree + I*Fpp(p1, p2, p3, p4, p5, p6);
    apm = -vpp*atree + I*Fpm(p1, p3, p2, p4, p5, p6);
    asf = -__nf*(2.0*(Log(-S(5,6))-Log(-S(2,3)))/3.0+10.0/9.0)*atree
      -    2.0*S(2,3)*atree/(15.0*mtsq);
    
    M[0] = atree;
    M[1] = app - 2.0*(app+apm)/Nc2 + asl/Nc2 + asf/Nc;
    M[2] = apm + (app+apm)/Nc2 - asl/Nc2 -asf/Nc;
  }
   
  void ampq4l2::
  su3_m1_pmpm(unsigned int nf, int p1, int p2, int p3, 
	      int p4, int p5, int p6, _ComplexD *M) const
  {
    _ComplexD I(0.0, 1.0), atree, app, asl, apm, vpp, vsl, asf;
    double __nf = (double) nf;
    
    atree = -App(p1, p3, p2, p4, p5, p6);
    vpp   = -13.0*(Log(-S(2,3))-Log(-S(5,6)))/6.0 + 10.0/9.0;
    vsl   = -1.5*(2.0*Log(-S(5,6))-Log(-S(1,4))-Log(-S(3,2))) - 7.5;
    
    asl =  vsl*atree + I*Fsl(p3, p2, p1, p4, p5, p6);
    app = -vpp*atree + I*Fpp(p1, p3, p2, p4, p5, p6);
    apm =  vpp*atree + I*Fpm(p1, p2, p3, p4, p5, p6);
    asf = -__nf*(2.0*(Log(-S(5,6))-Log(-S(2,3)))/3.0+10.0/9.0)*atree
      -   2.0*S(2,3)*atree/(15.0*mtsq);
    
    M[0] = atree;
    M[1] = apm - 2.0*(apm+app)/Nc2 - asl/Nc2 + asf/Nc;
    M[2] = app + (app+apm)/Nc2 + asl/Nc2 - asf/Nc;
  }
  
  void ampq4l2::
  color_m1_ppmm(unsigned int nf, int p1, int p2, int p3, 
		int p4, int p5, int p6, _ComplexD *M) const 
  {
    _ComplexD I(0.0,1.0), atree, app, asl, apm, vpp, vsl;
    double __nf = (double) nf;
    
    atree = App(p1, p2, p3, p4, p5, p6);
    vpp   = -13.0*(Log(-S(2,3))-Log(-S(5,6)))/6.0 + 10.0/9.0;
    vsl   = -1.5*(2.0*Log(-S(5,6))-Log(-S(1,4))-Log(-S(2,3))) - 7.5;
    
    asl  = -vsl*atree + I*Fsl(p2, p3, p1, p4, p5, p6);
    app  =  vpp*atree + I*Fpp(p1, p2, p3, p4, p5, p6);
    apm  = -vpp*atree + I*Fpm(p1, p3, p2, p4, p5, p6);
    M[3] = -__nf*(2.0*(Log(-S(5,6))-Log(-S(2,3)))/3.0+10.0/9.0)*atree
      -    2.0*S(2,3)*atree/(15.0*mtsq);
    
    M[0] =  atree;
    M[1] = -asl;
    M[2] =  asl - apm;
    M[4] =  app + apm;
    M[5] =  Aax(p1, p4, p2, p3, p5, p6); 
  }

  void ampq4l2::
  color_m1_pmpm(unsigned int nf, int p1, int p2, int p3, 
		int p4, int p5, int p6, _ComplexD *M) const
  {
    _ComplexD I(0.0, 1.0), atree, app, asl, apm, vpp, vsl;
    double __nf = (double) nf;
    
    atree = -App(p1, p3, p2, p4, p5, p6);
    vpp   = -13.0*(Log(-S(2,3))-Log(-S(5,6)))/6.0 + 10.0/9.0;
    vsl   = -1.5*(2.0*Log(-S(5,6))-Log(-S(1,4))-Log(-S(3,2))) - 7.5;
    
    asl  =  vsl*atree + I*Fsl(p3, p2, p1, p4, p5, p6);
    app  = -vpp*atree + I*Fpp(p1, p3, p2, p4, p5, p6);
    apm  =  vpp*atree + I*Fpm(p1, p2, p3, p4, p5, p6);
    M[3] = -__nf*(2.0*(Log(-S(5,6))-Log(-S(2,3)))/3.0+10.0/9.0)*atree
      -    2.0*S(2,3)*atree/(15.0*mtsq);
    
    M[0] =  atree;
    M[1] =  asl;
    M[2] = -asl - app;
    M[4] =  app + apm;
    M[5] = -Aax(p1, p4, p3, p2, p5, p6);
  }
  
  void ampq4l2::
  matrix_1loop_ppmm(unsigned int nf, int p1, int p2, int p3, 
		    int p4, int p5, int p6, _ComplexD *M) const 
  {
    su3_m1_ppmm(nf, p1, p2, p3, p4, p5, p6, M);
    swap();
    su3_m1_ppmm(nf, p3, p4, p1, p2, p6, p5, M+3);
    swap();
    M[6] = M[7] = M[8] = M[9] = M[10] = M[11] = 0.0;
  }
  
  void ampq4l2::
  matrix_1loop_pmmp(unsigned int nf, int p1, int p2, int p3, 
		    int p4, int p5, int p6, _ComplexD *M) const 
  {
    su3_m1_ppmm(nf, p1, p4, p3, p2, p5, p6, M+6);
    swap();
    su3_m1_ppmm(nf, p3, p2, p1, p4, p6, p5, M+9);
    swap();
    M[0] = M[1] = M[2] = M[3] = M[4] = M[5] = 0.0;
  }
  
  void ampq4l2::
  matrix_1loop_pmpm(unsigned int nf, int p1, int p2, int p3, 
		    int p4, int p5, int p6, _ComplexD *M) const 
  {
    su3_m1_pmpm(nf, p1, p2, p3, p4, p5, p6, M);
    su3_m1_pmpm(nf, p3, p4, p1, p2, p5, p6, M+3);
    su3_m1_pmpm(nf, p1, p4, p3, p2, p5, p6, M+6);
    su3_m1_pmpm(nf, p3, p2, p1, p4, p5, p6, M+9);
  }
  
  void ampq4l2::
  color_1loop_ppmm(unsigned int nf, int p1, int p2, int p3,
		   int p4, int p5, int p6, _ComplexD *M) const
  {
    color_m1_ppmm(nf, p1, p2, p3, p4, p5, p6, M);
    swap();
    color_m1_ppmm(nf, p3, p4, p1, p2, p6, p5, M+6);
    swap();
    
    for(short i = 12; i <= 23; i++) M[i] = 0.0;
  }

  void ampq4l2::
  color_1loop_pmmp(unsigned int nf, int p1, int p2, int p3,
		   int p4, int p5, int p6, _ComplexD *M) const
  {
    color_m1_ppmm(nf, p1, p4, p3, p2, p5, p6, M+12);
    swap();
    color_m1_ppmm(nf, p3, p2, p1, p4, p6, p5, M+18);
    swap();
    
    for(short i = 0; i <= 11; i++) M[i] = 0.0;
  }

  void ampq4l2::
  color_1loop_pmpm(unsigned int nf, int p1, int p2, int p3,
		   int p4, int p5, int p6, _ComplexD *M) const
  {
    color_m1_pmpm(nf, p1, p2, p3, p4, p5, p6, M);
    color_m1_pmpm(nf, p3, p4, p1, p2, p5, p6, M+6);
    color_m1_pmpm(nf, p1, p4, p3, p2, p5, p6, M+12);
    color_m1_pmpm(nf, p3, p2, p1, p4, p5, p6, M+18);
  }

  //--------       AMPLITUDES      -------
#define TreeAA(x,y) 2.0*Na*real(x*conj(y))
#define TreeAC(x,y) Na*4.0*real(x*conj(y))/Nc
  
  
  void ampq4l2::su3_amptree(_ComplexD *M, double *A) 
  {
    A[0] +=   TreeAA(M[0], M[0]);                      //   A*A
    A[1] +=   TreeAA(M[1], M[1]);                      //   B*B
    A[2] += 2.0*TreeAA(M[0], M[1]);                    // 2*A*B 
    
    A[3] +=   TreeAA(M[2], M[2]);                      //   C*C
    A[4] +=   TreeAA(M[3], M[3]);                      //   D*D
    A[5] += 2.0*TreeAA(M[2], M[3]);                    // 2*C*D
    
    A[6] += TreeAC(M[0], M[2]);                        // A*C
    A[7] += TreeAC(M[0], M[3]);                        // A*D
    A[8] += TreeAC(M[1], M[2]);                        // B*C
    A[9] += TreeAC(M[1], M[3]);                        // B*D 
  }
  
  void ampq4l2::
  su3_tree(int p1, int p2, int p3, int p4, int p5, int p6, double *Amp) const
  {
    _ComplexD M[4];
    std::memset(Amp, 0, 10*sizeof(double));
    
    matrix_tree_pmpm(p1, p2, p3, p4, p5, p6, M); su3_amptree(M, Amp);
    matrix_tree_ppmm(p1, p2, p3, p4, p5, p6, M); su3_amptree(M, Amp);
    matrix_tree_pmmp(p1, p2, p3, p4, p5, p6, M); su3_amptree(M, Amp);
    
    matrix_tree_pmpm(p1, p2, p3, p4, p6, p5, M); su3_amptree(M, Amp);
    matrix_tree_ppmm(p1, p2, p3, p4, p6, p5, M); su3_amptree(M, Amp);
    matrix_tree_pmmp(p1, p2, p3, p4, p6, p5, M); su3_amptree(M, Amp);
  }

  void ampq4l2::
  su3_tree_mch(int p1, int p2, int p3, int p4, int p5, int p6, double *Amp) const
  {
    _ComplexD M[4];
    std::memset(Amp, 0, 10*sizeof(double));
    int hel = (int) (6*_M_rng());
    
    switch(hel) {
    case 0: matrix_tree_pmpm(p1, p2, p3, p4, p5, p6, M); break;
    case 1: matrix_tree_ppmm(p1, p2, p3, p4, p5, p6, M); break;
    case 2: matrix_tree_pmmp(p1, p2, p3, p4, p5, p6, M); break;
    case 3: matrix_tree_pmpm(p1, p2, p3, p4, p6, p5, M); break;
    case 4: matrix_tree_ppmm(p1, p2, p3, p4, p6, p5, M); break;
    case 5: matrix_tree_pmmp(p1, p2, p3, p4, p6, p5, M); break;
    }

    su3_amptree(M, Amp);
    for(unsigned int i = 0; i < 10; i++) 
      Amp[i] *= 6.0;
  }
  
#undef TreeAA
#undef TreeAC


#define LoopAA(x,y) Nc*Na*real(x*conj(y))
#define LoopAC(x,y) (-Na*real(x*conj(y)))

  void ampq4l2::su3_amploop(_ComplexD *L, double *Amp) 
  {
    _ComplexD *A = L, *B = L+3, *C= L+6, *D = L+9;
    
    Amp[0] += LoopAA(A[1], A[0]);                            // A1*A0
    Amp[1] += LoopAA(B[1], B[0]);                            // B1*B0
    Amp[2] += LoopAA(A[1], B[0]) + LoopAA(B[1], A[0]);       // A1*B0+A0*B1 
    
    Amp[3] += LoopAA(C[1], C[0]);                            // C1*C0
    Amp[4] += LoopAA(D[1], D[0]);                            // D1*D0
    Amp[5] += LoopAA(C[1], D[0]) + LoopAA(D[1], C[0]);       // C1*D0+C0*D1 
    
    Amp[6] += LoopAC(A[2], C[0]) + LoopAC(C[2], A[0]);       // A2*C0+C2*A0  
    Amp[7] += LoopAC(A[2], D[0]) + LoopAC(D[2], A[0]);       // A2*D0+D2*A0  
    Amp[8] += LoopAC(B[2], C[0]) + LoopAC(C[2], B[0]);       // B2*C0+C2*B0
    Amp[9] += LoopAC(B[2], D[0]) + LoopAC(D[2], B[0]);       // B2*D0+D2*B0
  }
  
#undef LoopAA
#undef LoopAC

  void ampq4l2::
  su3_1loop(unsigned int nf, int p1, int p2, int p3, 
	    int p4, int p5, int p6, double *Amp) const 
  {
    _ComplexD M[12];
    std::memset(Amp, 0, 10*sizeof(double));
    
    matrix_1loop_pmpm(nf, p1,p2,p3,p4,p5,p6, M); su3_amploop(M, Amp);  
    matrix_1loop_ppmm(nf, p1,p2,p3,p4,p5,p6, M); su3_amploop(M, Amp);
    matrix_1loop_pmmp(nf, p1,p2,p3,p4,p5,p6, M); su3_amploop(M, Amp);
    matrix_1loop_pmpm(nf, p1,p2,p3,p4,p6,p5, M); su3_amploop(M, Amp);
    matrix_1loop_ppmm(nf, p1,p2,p3,p4,p6,p5, M); su3_amploop(M, Amp);
    matrix_1loop_pmmp(nf, p1,p2,p3,p4,p6,p5, M); su3_amploop(M, Amp);

    //-----the opposit helicity states (e+e- : means factor 2) -----
    swap();
    matrix_1loop_pmpm(nf, p1,p2,p3,p4,p5,p6, M); su3_amploop(M, Amp);
    matrix_1loop_ppmm(nf, p1,p2,p3,p4,p5,p6, M); su3_amploop(M, Amp);
    matrix_1loop_pmmp(nf, p1,p2,p3,p4,p5,p6, M); su3_amploop(M, Amp);
    matrix_1loop_pmpm(nf, p1,p2,p3,p4,p6,p5, M); su3_amploop(M, Amp);
    matrix_1loop_ppmm(nf, p1,p2,p3,p4,p6,p5, M); su3_amploop(M, Amp);
    matrix_1loop_pmmp(nf, p1,p2,p3,p4,p6,p5, M); su3_amploop(M, Amp);
    swap();
  }

  void ampq4l2::
  su3_1loop_mch(unsigned int nf, int p1, int p2, int p3, 
		int p4, int p5, int p6, double *ret_val) const 
  {
    _ComplexD M[12];
    int hpm = (int) (2*_M_rng());
    int hel = (int) (6*_M_rng());
    std::memset(ret_val, 0, 10*sizeof(double));
    
    if(hpm == 1) swap();
    switch(hel){
    case 0: matrix_1loop_pmpm(nf, p1,p2,p3,p4,p5,p6, M); break;
    case 1: matrix_1loop_ppmm(nf, p1,p2,p3,p4,p5,p6, M); break;
    case 2: matrix_1loop_pmmp(nf, p1,p2,p3,p4,p5,p6, M); break;
    case 3: matrix_1loop_pmpm(nf, p1,p2,p3,p4,p6,p5, M); break;
    case 4: matrix_1loop_ppmm(nf, p1,p2,p3,p4,p6,p5, M); break;
    case 5: matrix_1loop_pmmp(nf, p1,p2,p3,p4,p6,p5, M); break;
    }
    if(hpm == 1) swap();
    
    su3_amploop(M, ret_val);
    for(unsigned int i = 0; i < 10; i++) 
      ret_val[i] *= 12.0;
  }


  //
  //            COLOR CORRALETED AMPLITUDES
  //
#define Prod(x,y) (-Na*real(x*conj(y))/Nc)
  
  void ampq4l2::
  su3_ampcc(_ComplexD *M, double s1, double s2, double s3, double *A) 
  {
    A[0] +=     s1*Prod(M[0], M[0]);        //   A*A
    A[1] +=     s1*Prod(M[1], M[1]);        //   B*B
    A[2] += 2.0*s1*Prod(M[0], M[1]);        // 2*A*B
    
    A[3] +=     s2*Prod(M[2], M[2]);        //   C*C
    A[4] +=     s2*Prod(M[3], M[3]);        //   D*D
    A[5] += 2.0*s2*Prod(M[2], M[3]);        // 2*C*D
    
    A[6] += 2.0*s3*Prod(M[0], M[2]);        // A*C
    A[7] += 2.0*s3*Prod(M[0], M[3]);        // A*D
    A[8] += 2.0*s3*Prod(M[1], M[2]);        // B*C
    A[9] += 2.0*s3*Prod(M[1], M[3]);        // B*D
  }
#undef Prod
  
#define Cond14(i,j)       			\
(i==p1 && j==p4) || (i==p4 && j==p1) ||		\
(i==p3 && j==p2) || (i==p2 && j==p3)
  
#define Cond12(i,j)				\
(i==p1 && j==p2) || (i==p2 && j==p1) ||		\
(i==p3 && j==p4) || (i==p4 && j==p3)

#define Cond13(i,j)				\
(i==p1 && j==p3) || (i==p3 && j==p1) ||		\
(i==p2 && j==p4) || (i==p4 && j==p2)

 
  void ampq4l2::
  su3_cc(int i, int j, int p1, int p2, int p3, 
	 int p4, int p5, int p6, double *cc) const 
  {
    _ComplexD M[4];
    double s1, s2, s3;
    std::memset(cc, 0, 10*sizeof(double));
    
    if(Cond14(i,j))      s1 = -1.0, s2 = Nc2-2.0, s3 = -1.0/Nc;
    else if(Cond12(i,j)) s1 = Nc2-2.0, s2 = -1.0, s3 = -1.0/Nc;
    else if(Cond13(i,j)) s1 = s2 = 2.0, s3 = Nc+1.0/Nc; 
    else throw("Error in ampq4l2::su3_cc");

    matrix_tree_pmpm(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, s1, s2, s3, cc);
    matrix_tree_ppmm(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, s1, s2, s3, cc);
    matrix_tree_pmmp(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, s1, s2, s3, cc);
    
    matrix_tree_pmpm(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, s1, s2, s3, cc);
    matrix_tree_ppmm(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, s1, s2, s3, cc);
    matrix_tree_pmmp(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, s1, s2, s3, cc);
  }

  void ampq4l2::su3_ampcc(_ComplexD *M, double *cc) 
  {
    double s1 = -1.0, s2 = Nc2-2.0, s3 = -1.0/Nc;
    su3_ampcc(M, s1, s2, s3, cc+20);
    
    s1 = Nc2-2.0; s2 = -1.0;
    su3_ampcc(M, s1, s2, s3, cc);
    
    s1 = s2 = 2.0; s3 = Nc+1.0/Nc; 
    su3_ampcc(M, s1, s2, s3, cc+10);
  }
  
#define X(i,j) std::log(std::abs(s/S(i,j)))
#define l(a,b,c,d)  2.0*(Xq(S(a,b),s) + Xq(S(c,d), s))
  
  void ampq4l2::su3_kp(int pa, int p1, int p2, int p3, int p4,
		       int p5, int p6, su3_kp_i1 *res, double al) const 
  {
    _ComplexD M[4];
    double l12,l13, l14, cc[3][10], s = S(5,6);
    std::memset((double *) cc, 0, 30*sizeof(double));
    
    matrix_tree_pmpm(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, (double *) cc);
    matrix_tree_ppmm(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, (double *) cc);
    matrix_tree_pmmp(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, (double *) cc);
    
    matrix_tree_pmpm(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, (double *) cc);
    matrix_tree_ppmm(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, (double *) cc);
    matrix_tree_pmmp(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, (double *) cc);
    
    if(pa==p1)      l12 = X(1,2), l13 = X(1,3), l14 = X(1,4);
    else if(pa==p2) l12 = X(1,2), l13 = X(2,4), l14 = X(2,3);
    else if(pa==p3) l12 = X(3,4), l13 = X(1,3), l14 = X(2,3);
    else if(pa==p4) l12 = X(3,4), l13 = X(2,4), l14 = X(1,4);
    else throw("Error in ampq4l2::su3_kp(...)");
 
    double x12 = l(1,2,3,4), x13 = l(1,3,2,4), x14 = l(1,4,2,3), 
      gk = 4.0*(Gq + Kq(al)) - 2.0*Cf + Ca/3.0;
    
    for(unsigned int i = 0; i < 10; i++) {
      res[i].tree = -(cc[0][i]+cc[1][i]+cc[2][i])/Cf;
      res[i].loop = x12*cc[0][i]+x13*cc[1][i]+x14*cc[2][i] + (res[i].tree)*gk;
      
      res[i].ga = -Gq*(res[i].tree);
      res[i].pa = (l12*cc[0][i]+l13*cc[1][i]+l14*cc[2][i])/Cf;
    }
  }

  void ampq4l2::su3_kp_mch(int pa, int p1, int p2, int p3, int p4, 
			   int p5, int p6, su3_kp_i1 *res, double al) const 
  {
    _ComplexD M[4];
    double l12,l13, l14, cc[3][10], s = S(5,6);
    int hel = (int) (6*_M_rng());
    std::memset((double *) cc, 0, 30*sizeof(double));
    
    switch(hel) {
    case 0: matrix_tree_pmpm(p1, p2, p3, p4, p5, p6, M); break;
    case 1: matrix_tree_ppmm(p1, p2, p3, p4, p5, p6, M); break;
    case 2: matrix_tree_pmmp(p1, p2, p3, p4, p5, p6, M); break;
    case 3: matrix_tree_pmpm(p1, p2, p3, p4, p6, p5, M); break;
    case 4: matrix_tree_ppmm(p1, p2, p3, p4, p6, p5, M); break;
    case 5: matrix_tree_pmmp(p1, p2, p3, p4, p6, p5, M); break;
    }

    su3_ampcc(M, (double *) cc);
    if(pa==p1)      l12 = X(1,2), l13 = X(1,3), l14 = X(1,4);
    else if(pa==p2) l12 = X(1,2), l13 = X(2,4), l14 = X(2,3);
    else if(pa==p3) l12 = X(3,4), l13 = X(1,3), l14 = X(2,3);
    else if(pa==p4) l12 = X(3,4), l13 = X(2,4), l14 = X(1,4);
    else throw("Error in ampq4l2::su3_kp_mch(...)");
    
    double x12 = l(1,2,3,4), x13 = l(1,3,2,4), x14 = l(1,4,2,3), 
      gk = 4.0*(Gq + Kq(al)) - 2.0*Cf + Ca/3.0;

    for(unsigned int i = 0; i < 10; i++) {
      res[i].tree = -6.0*(cc[0][i]+cc[1][i]+cc[2][i])/Cf;
      res[i].loop = 6.0*(x12*cc[0][i]+x13*cc[1][i]+x14*cc[2][i]) + (res[i].tree)*gk;
      
      res[i].ga = -Gq*(res[i].tree);
      res[i].pa =  6.0*(l12*cc[0][i]+l13*cc[1][i]+l14*cc[2][i])/Cf;
    }
  }

  void ampq4l2::su3_ins(int p1, int p2, int p3, int p4, int p5,
			int p6, double *res, double al) const 
  {
    _ComplexD M[4];
    double tree, cc[3][10], s = S(5,6);
    std::memset((double *) cc, 0, 30*sizeof(double));
    
    matrix_tree_pmpm(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, (double *) cc);
    matrix_tree_ppmm(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, (double *) cc);
    matrix_tree_pmmp(p1, p2, p3, p4, p5, p6, M); su3_ampcc(M, (double *) cc);
    
    matrix_tree_pmpm(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, (double *) cc);
    matrix_tree_ppmm(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, (double *) cc);
    matrix_tree_pmmp(p1, p2, p3, p4, p6, p5, M); su3_ampcc(M, (double *) cc);
     
    double x12 = l(1,2,3,4), x13 = l(1,3,2,4), x14 = l(1,4,2,3), 
      gk = 4.0*(Gq + Kq(al)) - 2.0*Cf + Ca/3.0;
    
    for(unsigned int i = 0; i < 10; i++) {
      tree = -(cc[0][i]+cc[1][i]+cc[2][i])/Cf;
      res[i] = x12*cc[0][i]+x13*cc[1][i]+x14*cc[2][i] + tree*gk;
    }
  }

  void ampq4l2::su3_ins_mch(int p1, int p2, int p3, int p4, int p5,
			    int p6, double *res, double al) const 
  {
    _ComplexD M[4];
    double tree, cc[3][10], s = S(5,6);
    int hel = (int) (6*_M_rng());
    std::memset((double *) cc, 0, 30*sizeof(double));
    
    switch(hel) {
    case 0: matrix_tree_pmpm(p1, p2, p3, p4, p5, p6, M); break;
    case 1: matrix_tree_ppmm(p1, p2, p3, p4, p5, p6, M); break;
    case 2: matrix_tree_pmmp(p1, p2, p3, p4, p5, p6, M); break;
    case 3: matrix_tree_pmpm(p1, p2, p3, p4, p6, p5, M); break;
    case 4: matrix_tree_ppmm(p1, p2, p3, p4, p6, p5, M); break;
    case 5: matrix_tree_pmmp(p1, p2, p3, p4, p6, p5, M); break;
    }

    su3_ampcc(M, (double *) cc);

    double x12 = l(1,2,3,4), x13 = l(1,3,2,4), x14 = l(1,4,2,3), 
      gk = 4.0*(Gq + Kq(al)) - 2.0*Cf + Ca/3.0;

    for(unsigned int i = 0; i < 10; i++) {
      tree = -6.0*(cc[0][i]+cc[1][i]+cc[2][i])/Cf;
      res[i] = 6.0*(x12*cc[0][i]+x13*cc[1][i]+x14*cc[2][i]) + tree*gk;
    }
  }
}   //   namespace nlo

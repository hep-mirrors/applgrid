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
#include "ampq2g1p2.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {

  double ampq2g1p2::su3_tree(int p1, int p2, int p3,int p4,int p5) const
  {
    double s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), s12 = S(1,2),
      s24 = S(2,4), s25 = S(2,5), s23 = S(2,3);

    return 8.0*Na*s12/(s23*s13*s14*s24*s15*s25)
      *(s13*s23*(s13*s13+s23*s23) + s14*s24*(s14*s14+s24*s24)
	+ s15*s25*(s15*s15+s25*s25));
  }

  std::complex<double> 
  ampq2g1p2::amphtree(int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5),
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5); 

    return -8.0*Na*a12*a12*(a15*a24/(a14*a25)+a14*a25/(a15*a24))/(a13*a13*a23*a23);
  }
  
  //  COLOR CORRALETED AMPLITUDES
#define Cond12(i,j)                                     \
(i == p1 && j == p2) || (i == p2 && j == p1)

#define Cond13(i,j)                                     \
(i == p1 && j == p3) || (i == p3 && j == p1) ||         \
(i == p3 && j == p2) || (i == p2 && j == p3)

  std::pair<double, std::complex<double> > ampq2g1p2::
  su3_cc(int pi,int pj, int p1, int p2, int p3, int p4, int p5) const
  {
    double s = 0.0, cc = su3_tree(p1, p2, p3, p4, p5);
    _ComplexD hcc = 0.0;

    if(Cond12(pi,pj))      s = 0.5/Nc;
    else if(Cond13(pi,pj)) s = -0.5*Nc;
    else throw("Error in ampq2g1p2::su3_cc(...)");
    
    if(pj == p3) hcc = amphtree(p1, p2, p3, p4, p5);
    
    return _Pair(s*cc, s*hcc);
  }

#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampq2g1p2::su3_kp(unsigned int nf, int pa, int pb, int p1, 
			 int p2, int p3, int p4, int p5, su3_kp_i2 *res,
			 double al) const
  {
    double cc, xi, xj, xq = Gq/Cf, xg = Gg(nf)/Nc;
    int p[3] = {p1,p2,p3};

    double Ta = (pa == p1 || pa == p2 ? Cf : Ca);
    double Tb = (pb == p1 || pb == p2 ? Cf : Ca);
    double x[3] = {xq, xq, xg};

    res->tree = su3_tree(p1,p2, p3, p4,p5);
    res->pa = res->pb = res->ga = res->gb = res->loop = 0.0;
 
    for(unsigned i = 0; i < 3; i++)
      for(unsigned j = i+1; j < 3; j++) {
	if(Cond12(p[i], p[j]))      cc = 0.5*(res->tree)/Nc;
	else if(Cond13(p[i], p[j])) cc = -0.5*Nc*(res->tree);
	else throw("Error in ampq2g1p2::su3_kp(...)");
	
	//----- loop contribution -----
	if(i < 2) {
	  xi = Xq(Sij(p[i],p[j]), 1.0);
	  xj = j < 2 ? xi : Xg(Sij(p[i],p[j]), 1.0, nf);
	} else {
	  xi = Xg(Sij(p[i],p[j]), 1.0, nf);
	  xj = j < 2 ? Xq(Sij(p[i],p[j]), 1.0) : xi;
	}

	res->loop += (xi+xj)*cc;

	//----- finite contributions -----
	if(p[i] == pa || p[j] == pa) {
	  int iq = p[i] == pa ? j : i;
	  res->pa += X(pa, p[iq])*cc;
	  if(p[iq] != pb) res->ga += x[iq]*cc;
	  else res->cca = -cc/Ta;
	}
	
	if(p[i] == pb || p[j] == pb) {
	  int iq = p[i] == pb ? j : i;
	  res->pb += X(pb, p[iq])*cc;
	  if(p[iq] != pa) res->gb += x[iq]*cc;
	  else res->ccb = -cc/Tb;
	}
      }
    
    res->pa /= Ta; res->pb /= Tb;
    res->loop += (Gg(nf)+Kg(nf,al) + 2.0*(Gq+Kq(al)) - Cf)*(res->tree);
  }

#define Qu2 0.44444444444444444444
#define Qd2 0.11111111111111111111

  void ampq2g1p2::
  matrix_1loop(unsigned int nu, unsigned int nd, const _AmpPrim *A, const _AmpPrim *B, amp_1loop *res)
  { 
    double TrQ2 = (nu*Qu2 + nd*Qd2);
    _ComplexD Asum(0.0), Bsum, Csum(0.0), Tsum(0.0);
    
    
    for(unsigned i = 0; i < 6; i++) {
      Tsum += A[i].Atree;
      Asum += A[i].Asusy - A[i].AL;
      Csum += A[i].Asf;
    }
    
    Bsum = B[0].Asusy-B[0].AL+B[0].Asf + B[1].Asusy-B[1].AL+B[1].Asf;
    
    res -> U0 = Qu2*Tsum;
    res -> D0 = Qd2*Tsum;
    
    res -> U1 = -Qu2*Bsum - Qu2*(Asum+Csum)/Nc2 - TrQ2*Csum/Nc;
    res -> D1 = -Qd2*Bsum - Qd2*(Asum+Csum)/Nc2 - TrQ2*Csum/Nc;
  } 
  
#define __PrimitiveAmplitudes(h3,h4,h5)		\
  A1##h3##h4##h5(p1,p2, p3,p4,p5, A[0]);	\
  A1##h4##h5##h3(p1,p2, p4,p5,p3, A[1]);	\
  A1##h5##h3##h4(p1,p2, p5,p3,p4, A[2]);	\
  			       	                \
  A1##h3##h5##h4(p1,p2, p3,p5,p4, A[3]);	\
  A1##h5##h4##h3(p1,p2, p5,p4,p3, A[4]);	\
  A1##h4##h3##h5(p1,p2, p4,p3,p5, A[5]);	\
						\
						\
  A2##h3##h4##h5(p1, p3, p2, p4,p5, B[0]);	\
  A2##h3##h5##h4(p1, p3, p2, p5,p4, B[1])
  
  void ampq2g1p2::
  matrix_1loop_mpmpp(unsigned int nu, unsigned int nd, int p1, int p2,
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(m,p,p); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mppmp(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(p,m,p); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mpppm(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(p,p,m); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mppmm(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(p,m,m); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mpmpm(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(m,p,m); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g1p2::
  matrix_1loop_mpmmp(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[2];
    __PrimitiveAmplitudes(m,m,p); 
    matrix_1loop(nu, nd, A, B, res);
  }

#define Amp1Loop				\
  out[0] += real(M.U1*conj(M.U0));		\
  out[1] += real(M.D1*conj(M.D0))
  
  
  void ampq2g1p2::
  su3_1loop(unsigned int nu, unsigned int nd, int p1, int p2,
	    int p3, int p4, int p5, double *out) const
  {
    static amp_1loop M;
    out[0] = out[1] = 0.0;
    
    matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    		  
    swap();	  
    matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); Amp1Loop;
    swap();

    out[0] *= 4.0*Na*Nc;
    out[1] *= 4.0*Na*Nc;
  }

  void ampq2g1p2::
  su3_1loop_mch(unsigned int nu, unsigned int nd, int p1, int p2,
		int p3, int p4, int p5, double *out) const
  {
    static amp_1loop M;
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (6*_M_rng());
    
    if(hpm == 1) swap();
    switch(hel){
    case 0: matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 1: matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 2: matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 3: matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 4: matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); break;
    case 5: matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); break;
    }
    if(hpm == 1) swap();

    out[0] = 48.0*Na*Nc*real(M.U1*conj(M.U0));
    out[1] = 48.0*Na*Nc*real(M.D1*conj(M.D0));
  }

}  //  namespace nlo

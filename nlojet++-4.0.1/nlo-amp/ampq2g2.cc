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
#include "ampq2g2.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {


  double ampq2g2::su3_tree(int p1, int p2, int p3, int p4) const
  {
    static double cfac = 2.0*(Nc2-1.0)/Nc;
    double s = S(1,2), t = S(1,3), u = S(1,4);
    return cfac*(t/u + u/t)*(Nc2*(1.0-2.0*u*t/(s*s))-1.0);
  }

  double ampq2g2::ampcc12(int p1, int p2, int p3, int p4) const
  {
    double s = S(1,2), t = S(1,3), u = S(1,4);
    return -(Nc2-1.0)*(u*u + t*t)*(2.0/(s*s) + 1.0/(u*t*Nc2));
  }

  double ampq2g2::ampcc13(int p1, int p2, int p3, int p4) const
  {
    double s = S(1,2), t = S(1,3), u = S(1,4);
    return -(Nc2-1.0)*(u*u + t*t)*(Nc2*u*u - s*s)/(u*t*s*s);
  }
  
  double ampq2g2::ampcc14(int p1, int p2, int p3, int p4) const
  {
    double s = S(1,2), t = S(1,3), u = S(1,4);
    return -(Nc2-1.0)*(u*u + t*t)*(Nc2*t*t - s*s)/(u*t*s*s);
  }

  double ampq2g2::ampcc34(int p1, int p2, int p3, int p4) const
  {
    double s = S(1,2), t = S(1,3), u = S(1,4);
    double ut = u*u + t*t;
    return -(Nc2-1.0)*Nc2*ut*ut/(u*t*s*s);
  }

 
#define __cc_switch12 (pi == p1 && pj == p2) || (pi == p2 && pj == p1)
#define __cc_switch34 (pi == p3 && pj == p4) || (pi == p4 && pj == p3)

#define __cc_switch13 					\
  (pi == p1 && pj == p3) || (pi == p2 && pj == p4) ||	\
  (pi == p3 && pj == p1) || (pi == p4 && pj == p2)

#define __cc_switch14 					\
  (pi == p1 && pj == p4) || (pi == p2 && pj == p3) ||	\
  (pi == p4 && pj == p1) || (pi == p3 && pj == p2)
  
  
  double ampq2g2::su3_cc(int pi,int pj, int p1, int p2, int p3, int p4) const
  {
    if(__cc_switch12) return ampcc12(p1, p2, p3, p4);
    else if(__cc_switch13) return ampcc13(p1, p2, p3, p4);
    else if(__cc_switch14) return ampcc14(p1, p2, p3, p4);
    else if(__cc_switch34) return ampcc34(p1, p2, p3, p4);
    else throw("Error in ampq2g2::amp_cc");	
  }

#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampq2g2::su3_kp(unsigned int nf, int pa, int pb, int p1, 
		     int p2, int p3, int p4, su3_kp_i2 *res, double al) const
  {
    double cc, xi, xj, xq = Gq/Cf, xg = Gg(nf)/Nc;
    int p[4] = {p1,p2,p3,p4};

    double Ta = (pa == p1 || pa == p2 ? Cf : Ca);
    double Tb = (pb == p1 || pb == p2 ? Cf : Ca);
    double x[4] = {xq, xq, xg, xg};

    res->tree = res->pa = res->pb = res->ga = res->gb = res->loop = 0.0;
    for(unsigned i = 0; i < 4; i++)
      for(unsigned j = i+1; j < 4; j++) {
	cc = su3_cc(p[i], p[j], p1, p2, p3, p4);
	
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
	  res->tree -= cc;
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
    
    res->tree /= Ta; res->pa /= Ta; res->pb /= Tb;
    res->loop += (2.0*(Gg(nf)+Kg(nf,al) + Gq+Kq(al)) - Cf)*(res->tree);
  }

  //  
  //   Loop contribution
  // 
#define NLO_PI2 9.86960440108935861881

  void ampq2g2::A1mp(int p1, int p2, int p3, int p4, _AmpPrim& res) const
  {
    double s12 = S(1,2), s13 = S(1,3), s23 = S(2,3);
    _ComplexD a12 = A(1,2), a13 = A(1,3), a34 = A(3,4), a41 = A(4,1);
    
    _ComplexD l23 = Log(-s23), l12 = Log(-s12);
    _ComplexD Vg = std::pow(l12 - l23,2) + NLO_PI2;
    _ComplexD Vsusy = Vg + 1.5*(l12+l23) - 6.0;
    _ComplexD VL = Vg + 1.5*l12 - 2.5;
    
    _ComplexD Fsusy = 1.5*((s12-s23)/s13*(l12-l23) - s12*s23*Vg/(s13*s13));
    _ComplexD FL = s12/(2.0*s13)*((3.0-2.0*s23/s13)*(l12-l23)
				  + (s12*s12/(s13*s13)-3.0*s23/s13)*Vg + 1.0);
    
    res.Atree = ICPLX*a13*a13*a13/(a12*a34*a41);
    res.Asusy = res.Atree*(Vsusy + Fsusy); res.AL = res.Atree*(VL + FL);
  }
  
  void ampq2g2::A1pm(int p1, int p2, int p3, int p4, _AmpPrim& res) const
  {
    double s12 = S(1,2), s13 = S(1,3), s23 = S(2,3);
    _ComplexD a12 = A(1,2), a24 = A(2,4), a23 = A(2,3), a34 = A(3,4),
      a41 = A(4,1);
    
    _ComplexD l23 = Log(-s23), l12 = Log(-s12);
    _ComplexD Vg = std::pow(l12 - l23,2) + NLO_PI2;
    _ComplexD Vsusy = Vg + 3.0*l12 - 6.0;
    _ComplexD VL = Vg + 1.5*l12 - 2.5;
    
    _ComplexD FL = s12/(2.0*s13)*Vg;
    
    res.Atree = -ICPLX*a41*a41*a24/(a12*a23*a34);
    res.Asusy = res.Atree*Vsusy; res.AL = res.Atree*(VL + FL);
  }
  
  void ampq2g2::A2mp(int p1, int p2, int p3, int p4, _AmpPrim& res) const
  {
    double s12 = S(1,2), s23 = S(2,3);
    _ComplexD a12 = A(1,2), a34 = A(3,4), a14 = A(1,4);
    
    _ComplexD l23 = Log(-s23), l12 = Log(-s12);
    _ComplexD Vg = std::pow(l12 - l23,2) + NLO_PI2;
    _ComplexD Vsusy = Vg + 3.0*l23 - 6.0;
    
    res.Atree = ICPLX*a12*a12/(a34*a14);
    res.AL = 0.5*(res.Asusy = res.Atree*Vsusy);
  }

#define __PrimitiveAmplitudes(h3,h4)		\
  A1##h3##h4(p1,p2, p3,p4, A[0]);		\
  A1##h4##h3(p1,p2, p4,p3, A[1]);		\
  						\
  A2##h3##h4(p1, p3, p2, p4, B[0]);		\
  A2##h4##h3(p1, p4, p2, p3, B[1])

  void ampq2g2::
  matrix_1loop_pmpm(int p2, int p3, int p4, int p1, amp_1loop *res) const
  {
    static _AmpPrim A[2], B[2];
    __PrimitiveAmplitudes(m,p); 
    matrix_1loop(A, B, res);
  }


  void ampq2g2::
  matrix_1loop_ppmm(int p2, int p3, int p4, int p1, amp_1loop *res) const
  {
    static _AmpPrim A[2], B[2];
    __PrimitiveAmplitudes(p,m); 
    matrix_1loop(A, B, res);
  }

  void ampq2g2::
  matrix_1loop(const _AmpPrim *A, const _AmpPrim *B, amp_1loop *res) const
  { 
    double onc2 = 1.0/Nc2;
    _ComplexD Asum = A[0].Asusy + A[1].Asusy;
    
    for(unsigned i = 0; i < 2; i++) {
      res[i].A0 = A[i].Atree;
      res[i].A1 = (1.0+onc2)*A[i].AL - onc2*A[i].Asusy;
      res[i].A3 = Asum + B[i].Asusy;
    }
  } 

  double ampq2g2::amp1loop(amp_1loop *M) const
  {
    _ComplexD X0 = Na*M[0].A1 - M[1].A1 + M[0].A3;
    _ComplexD X1 = Na*M[1].A1 - M[0].A1 + M[1].A3;
    return Na*real(X0*conj(M[0].A0) + X1*conj(M[1].A0));
  }			      

  double ampq2g2::su3_1loop(int p1, int p2, int p3, int p4) const
  {
    static amp_1loop M[2];
    double ret_val;

    matrix_1loop_ppmm(p1,p3,p4,p2, M); ret_val  = amp1loop(M);
    matrix_1loop_pmpm(p1,p3,p4,p2, M); ret_val += amp1loop(M);

    swap();
    matrix_1loop_ppmm(p1,p3,p4,p2, M); ret_val += amp1loop(M);
    matrix_1loop_pmpm(p1,p3,p4,p2, M); ret_val += amp1loop(M);
    swap();

    return ret_val;
  }

  double ampq2g2::su3_1loop_mch(int p1, int p2, int p3, int p4) const
  {
    static amp_1loop M[2];
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (2*_M_rng());
    
    if(hpm == 1) swap();
    switch(hel){
    case 0: matrix_1loop_ppmm(p1,p3,p4,p2, M); break;
    case 1: matrix_1loop_pmpm(p1,p3,p4,p2, M); break;
    }
    if(hpm == 1) swap();

    return 4.0*amp1loop(M);
  }

  
  
}


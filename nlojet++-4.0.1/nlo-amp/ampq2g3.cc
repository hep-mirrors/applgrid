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
#include "ampq2g3.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {

  std::complex<double> 
  ampq2g3::amp2p1m(int I, int p1, int p2, int p3, int p4, int p5) const
  {
    std::complex<double> res(0.0, 1.0);
    
    res *= Aij(p1,I)*pow(Aij(p5,I), 3);
    res /= A(1,2)*A(2,3)*A(3,4)*A(4,5)*A(5,1);
    return res;
  }

  std::complex<double> 
  ampq2g3::amp1p2m(int I, int p1, int p2, int p3, int p4, int p5) const
  {
    std::complex<double> res(0.0, 1.0);
    
    res *= Bij(I,p5)*pow(Bij(I,p1), 3);
    res /= B(1,5)*B(5,4)*B(4,3)*B(3,2)*B(2,1);
    return res;
  }
    
  void ampq2g3::
  matrix_tree_pppmm(int p1, int p2, int p3, int p4, int p5, _ComplexD *M) const
  {
    M[0] = amp2p1m(p4, p1, p2, p3, p4, p5);
    M[1] = amp2p1m(p4, p1, p3, p4, p2, p5);
    M[2] = amp2p1m(p4, p1, p4, p2, p3, p5);
    M[3] = amp2p1m(p4, p1, p2, p4, p3, p5);
    M[4] = amp2p1m(p4, p1, p4, p3, p2, p5);
    M[5] = amp2p1m(p4, p1, p3, p2, p4, p5);
  }
  
  void ampq2g3::
  matrix_tree_ppmpm(int p1, int p2, int p3, int p4, int p5, _ComplexD *M) const
  {
    M[0] = amp2p1m(p3, p1, p2, p3, p4, p5);
    M[1] = amp2p1m(p3, p1, p3, p4, p2, p5);
    M[2] = amp2p1m(p3, p1, p4, p2, p3, p5);
    M[3] = amp2p1m(p3, p1, p2, p4, p3, p5);
    M[4] = amp2p1m(p3, p1, p4, p3, p2, p5);
    M[5] = amp2p1m(p3, p1, p3, p2, p4, p5);
  }
  
  void ampq2g3::
  matrix_tree_pmppm(int p1, int p2, int p3, int p4, int p5, _ComplexD *M) const
  {
    M[0] = amp2p1m(p2, p1, p2, p3, p4, p5);
    M[1] = amp2p1m(p2, p1, p3, p4, p2, p5);
    M[2] = amp2p1m(p2, p1, p4, p2, p3, p5);
    M[3] = amp2p1m(p2, p1, p2, p4, p3, p5);
    M[4] = amp2p1m(p2, p1, p4, p3, p2, p5);
    M[5] = amp2p1m(p2, p1, p3, p2, p4, p5);
  }
  
  void ampq2g3::
  matrix_tree_ppmmm(int p1, int p2, int p3, int p4, int p5, _ComplexD *M) const
  {
    M[0] = amp1p2m(p2, p1, p2, p3, p4, p5);
    M[1] = amp1p2m(p2, p1, p3, p4, p2, p5);
    M[2] = amp1p2m(p2, p1, p4, p2, p3, p5);
    M[3] = amp1p2m(p2, p1, p2, p4, p3, p5);
    M[4] = amp1p2m(p2, p1, p4, p3, p2, p5);
    M[5] = amp1p2m(p2, p1, p3, p2, p4, p5);
  }
  
  void ampq2g3::
  matrix_tree_pmpmm(int p1, int p2, int p3, int p4, int p5, _ComplexD *M) const
  {
    M[0] = amp1p2m(p3, p1, p2, p3, p4, p5);
    M[1] = amp1p2m(p3, p1, p3, p4, p2, p5);
    M[2] = amp1p2m(p3, p1, p4, p2, p3, p5);
    M[3] = amp1p2m(p3, p1, p2, p4, p3, p5);
    M[4] = amp1p2m(p3, p1, p4, p3, p2, p5);
    M[5] = amp1p2m(p3, p1, p3, p2, p4, p5);
  }
  
  void ampq2g3::
  matrix_tree_pmmpm(int p1, int p2, int p3, int p4, int p5, _ComplexD *M) const
  {
    M[0] = amp1p2m(p4, p1, p2, p3, p4, p5);
    M[1] = amp1p2m(p4, p1, p3, p4, p2, p5);
    M[2] = amp1p2m(p4, p1, p4, p2, p3, p5);
    M[3] = amp1p2m(p4, p1, p2, p4, p3, p5);
    M[4] = amp1p2m(p4, p1, p4, p3, p2, p5);
    M[5] = amp1p2m(p4, p1, p3, p2, p4, p5);
  }

  double ampq2g3::su3_tree(int p1, int p2, int p3,int p4,int p5) const
  {
    double res, s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), s12 = S(1,2),
      s34 = S(3,4), s35 = S(3,5), s32 = S(3,2), s45 = S(4,5), 
      s42 = S(4,2), s52 = S(5,2);

    double prop = s13*s14*s15*s32*s42*s52;
    double num  = s13*s32*(s13*s13 + s32*s32) + s14*s42*(s14*s14 + s42*s42)
      + s15*s52*(s15*s15 + s52*s52);

    double sum1 = (s13*s42+s14*s32)/s34 + (s13*s52+s15*s32)/s35 
      + (s14*s52+s15*s42)/s45;
    
    double sum2 = s13*s32*(s14*s52 + s15*s42)/(s34*s35)
      + s14*s42*(s15*s32 + s13*s52)/(s45*s34) 
      + s15*s52*(s13*s42 + s14*s32)/(s45*s35);
    
    res  = 2.0*Na*num*(s12 + Nc2*(s12-sum1) + Nc2*Nc2*sum2/s12)/(prop*Nc2);
  
    return res;
  }

  double ampq2g3::ampcc13(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5),
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), 
      s35 = S(3,5), s45 = S(4,5);

    double num  = s13*s23*(s13*s13 + s23*s23) + s14*s24*(s14*s14 + s24*s24)
      + s15*s25*(s15*s15 + s25*s25);

    double si = (Nc2-1.0)*s12/(Nc*s13*s23*s14*s24*s15*s25)
      + Nc*(Nc2*Nc2-1.0)/(s45*s13*s12)*(1.0/(s34*s25)+1.0/(s35*s24))
      - Nc*(Nc2-1.0)/(s34*s35*s45*s13*s23*s14*s24*s15*s25*s12)
      *(-3.0*s45*s13*s13*s24*s25-s45*s45*s13*s23*s12+s45*s23*s23*s14*s15

	+ (s35*s13*s14*s24*s25 + s45*s13*s23*s24*s15 + s35*s45*s13*s24*s12
	   +s35*s35*s24*s14*s12 + s35*s23*s14*s14*s25 - s35*s13*s15*s24*s24)

	+ (s34*s13*s15*s25*s24 + s45*s13*s23*s25*s14 + s34*s45*s13*s25*s12
	   +s34*s34*s25*s15*s12 + s34*s23*s15*s15*s24 - s34*s13*s14*s25*s25)
	);
    
    return -num*si;
  }

  double ampq2g3::ampcc12(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4),
      s35 = S(3,5), s45 = S(4,5);
    double D = s34*s35*s45*s13*s23*s14*s24*s15*s25*s12;
    
    double num  = s13*s23*(s13*s13 + s23*s23) + s14*s24*(s14*s14 + s24*s24)
      + s15*s25*(s15*s15 + s25*s25);

    double si = (1.0-Nc2*Nc2)/(Nc2*Nc)*s12/(s13*s23*s14*s24*s15*s25)
      + (Nc2-1.0)
      *(Nc*(-2.0*s34*s15*s15*s23*s24-2.0*s34*s25*s25*s13*s14
	    +(s34*s15*s25+s35*s45*s12)*(s13*s24+s23*s14))
	-(s34*s35*s14*s25*s12-2.0*s45*s13*s23*s24*s25-2.0*s13*s13*s45*s24*s25
	  +s34*s35*s45*s13*s23+s34*s35*s24*s15*s12-2.0*s45*s23*s13*s14*s15
	  -2.0*s23*s23*s45*s14*s15+s34*s35*s45*s23*s13)/Nc 

	+ Nc*(-2.0*s45*s13*s13*s24*s25-2.0*s45*s23*s23*s14*s15
	      +(s45*s13*s23+s34*s35*s12)*(s14*s25+s24*s15))
	-(s45*s34*s15*s23*s12-2.0*s35*s14*s24*s25*s23-2.0*s14*s14*s35*s25*s23
	  +s45*s34*s35*s14*s24+s45*s34*s25*s13*s12-2.0*s35*s24*s14*s15*s13
	  -2.0*s24*s24*s35*s15*s13+s45*s34*s35*s24*s14)/Nc 
	
	+ Nc*(-2.0*s35*s14*s14*s25*s23-2.0*s35*s24*s24*s15*s13
	      +(s35*s14*s24+s45*s34*s12)*(s15*s23+s25*s13))
	-(s35*s45*s13*s24*s12-2.0*s34*s15*s25*s23*s24-2.0*s15*s15*s34*s23*s24
	  +s35*s45*s34*s15*s25+s35*s45*s23*s14*s12-2.0*s34*s25*s15*s13*s14
	  -2.0*s25*s25*s34*s13*s14+s35*s45*s34*s25*s15)/Nc)/D; 

    return -num*si;
  }

  double ampq2g3::ampcc34(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), 
      s35 = S(3,5), s45 = S(4,5);

    double num  = s13*s23*(s13*s13 + s23*s23) + s14*s24*(s14*s14 + s24*s24)
      + s15*s25*(s15*s15 + s25*s25);

    double si = Nc2*Nc*(Nc2-1.0)/s12
      *(1.0/(s34*s35*s14*s25) + 1.0/(s34*s45*s13*s25)
	+1.0/(s34*s35*s24*s15) + 1.0/(s34*s45*s23*s15))
      -Nc*(Nc2-1.0)/(s34*s35*s45*s13*s23*s14*s24*s15*s25*s12)
      *(-2.0*s34*s35*s45*(s13*s23+s14*s24)-2.0*s34*s34*s15*s25*s12
	-s35*s45*s12*(s13*s24+s23*s14)
	+2.0*s45*(s13+s23)*(s23*s14*s15+s13*s24*s25)
	+2.0*s35*(s14+s24)*(s23*s14*s25+s13*s24*s15));

    return -num*si;
  }

#define __hamp_denom(i,j,k) (a1##i*a##i##j*a##j##k*a##k##2)

  std::complex<double> 
  ampq2g3::amphtree(int p, int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), a32 = A(3,2),
      a34 = A(3,4), a35 = A(3,5), a42 = A(4,2), a43 = A(4,3), a45 = A(4,5),
      a52 = A(5,2), a53 = A(5,3), a54 = A(5,4), a12 = A(1,2);
    
    int pi, pj;
    if(p == p3) pi = p4, pj = p5;
    else if(p == p4) pi = p3, pj = p5;
    else if(p == p5) pi = p3, pj = p4;
    else throw("Error in ampq2g3::amphtree");

    _ComplexD a1i = Aij(p1,pi), a2i = Aij(p2,pi), a1j = Aij(p1,pj), 
      a2j = Aij(p2,pj), num;

    num = 2.0*a1i*a2i*a1j*a2j*(a1i*a1i*a2j*a2j+a1j*a1j*a2i*a2i)/(a12*a12);
   
    _ComplexD m[6];
    m[0] = 1.0/__hamp_denom(3,4,5);
    m[1] = 1.0/__hamp_denom(4,5,3);
    m[2] = 1.0/__hamp_denom(5,3,4);
    m[3] = 1.0/__hamp_denom(3,5,4);
    m[4] = 1.0/__hamp_denom(5,4,3);
    m[5] = 1.0/__hamp_denom(4,3,5);
    
    _ComplexD A0, A1, A2;
    A0 = m[0]+m[1]+m[2]+m[3]+m[4]+m[5];
    A1 = m[0]*(m[3]-m[4]+m[5])+m[1]*(m[5]-m[3]+m[4])+m[2]*(m[4]-m[5]+m[3]);
    A2 = m[0]*m[0]+m[1]*m[1]+m[2]*m[2]+m[3]*m[3]+m[4]*m[4]+m[5]*m[5];
    
    return Na*num*((Nc2-2.0)*A2 - 2.0*A1 + A0*A0/Nc2);
  }
  

#define __cc_switch(i,j)					\
  (pi == p##i && pj == p##j) || (pi == p##j && pj == p##i)


  std::pair<double, std::complex<double> > ampq2g3::
  su3_cc(int pi,int pj, int p1, int p2, int p3, int p4, int p5) const
  {
    double cc;
    _ComplexD hcc = 0.0;

    if(__cc_switch(1,2)) cc = ampcc12(p1, p2, p3, p4, p5);
    else if(__cc_switch(1,3)) cc = ampcc13(p1, p2, p3, p4, p5);
    else if(__cc_switch(1,4)) cc = ampcc13(p1, p2, p4, p5, p3);
    else if(__cc_switch(1,5)) cc = ampcc13(p1, p2, p5, p3, p4);
    else if(__cc_switch(2,3)) cc = ampcc13(p2, p1, p3, p4, p5);
    else if(__cc_switch(2,4)) cc = ampcc13(p2, p1, p4, p5, p3);
    else if(__cc_switch(2,5)) cc = ampcc13(p2, p1, p5, p3, p4);
    else if(__cc_switch(3,4)) cc = ampcc34(p1, p2, p3, p4, p5);
    else if(__cc_switch(3,5)) cc = ampcc34(p1, p2, p5, p3, p4);
    else if(__cc_switch(4,5)) cc = ampcc34(p1, p2, p4, p5, p3);
    else throw("Error in ampq2g3::su3_cc");
   
    if(pj == p3 || pj == p4 || pj == p5)
      hcc = Nc*amphtree(pj, p1, p2, p3, p4, p5)/4.0;

    return _Pair(cc, hcc);
  }

  double ampq2g3::
  ampcc(int pi,int pj, int p1, int p2, int p3, int p4, int p5) const
  {
    if(__cc_switch(1,2)) return ampcc12(p1, p2, p3, p4, p5);
    else if(__cc_switch(1,3)) return ampcc13(p1, p2, p3, p4, p5);
    else if(__cc_switch(1,4)) return ampcc13(p1, p2, p4, p5, p3);
    else if(__cc_switch(1,5)) return ampcc13(p1, p2, p5, p3, p4);
    else if(__cc_switch(2,3)) return ampcc13(p2, p1, p3, p4, p5);
    else if(__cc_switch(2,4)) return ampcc13(p2, p1, p4, p5, p3);
    else if(__cc_switch(2,5)) return ampcc13(p2, p1, p5, p3, p4);
    else if(__cc_switch(3,4)) return ampcc34(p1, p2, p3, p4, p5);
    else if(__cc_switch(3,5)) return ampcc34(p1, p2, p5, p3, p4);
    else if(__cc_switch(4,5)) return ampcc34(p1, p2, p4, p5, p3);
    else throw("Error in ampq2g3::su3_cc");
  }


#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampq2g3::su3_kp(unsigned int nf, int pa, int pb, int p1, int p2, 
		       int p3, int p4, int p5, su3_kp_i2 *res, double al) const
  {
    double cc, xi, xj, xq = Gq/Cf, xg = Gg(nf)/Nc;
    int p[5] = {p1,p2,p3,p4,p5};

    double Ta = (pa == p1 || pa == p2 ? Cf : Ca);
    double Tb = (pb == p1 || pb == p2 ? Cf : Ca);
    double x[5] = {xq, xq, xg, xg, xg};

    res->tree = res->pa = res->pb = res->ga = res->gb = res->loop = 0.0;
    for(unsigned i = 0; i < 5; i++)
      for(unsigned j = i+1; j < 5; j++) {
	cc = ampcc(p[i], p[j], p1, p2, p3, p4, p5);
	
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
    res->loop += (3.0*(Gg(nf)+Kg(nf,al)) + 2.0*(Gq+Kq(al)) - Cf)*(res->tree);
  }


  //
  //    1-loop amplitudes
  //
#define LL(x,y,z,w) (x-y)*(z-w)

  std::complex<double> 
  ampq2g3::Vg(int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD  l12 = Log(-S(1,2)), l23 = Log(-S(2,3)),
      l34 = Log(-S(3,4)), l45 = Log(-S(4,5)), l51 = Log(-S(5,1));
    
    return LL(l12,l23,l34,l45) + LL(l23,l34,l45,l51) + LL(l34,l45,l51,l12)
      + LL(l45,l51,l12,l23) + LL(l51,l12,l23,l34) + 8.224670334241132; 
  }  

  void ampq2g3::A1mpp(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const
  {
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s51 = S(5,1); 
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), 
      a23 = A(2,3), a34 = A(3,4), a35 = A(3,5), a41 = A(4,1), a43 = A(4,3),
      a45 = A(4,5), a51 = A(5,1), b12 = B(1,2), b15 = B(1,5), b23 = B(2,3),
      b24 = B(2,4), b25 = B(2,5), b34 = B(3,4), b45 = B(4,5), b52 = B(5,2);

    _ComplexD Vsusy, VL, Fsusy, Ff, Fs, FL;
    _ComplexD l34 = Log(-s34), l51 = Log(-s51), vg = Vg(p1, p2, p3, p4, p5);
    
    res.Atree = ICPLX*a13*a13*a13*a23/(a12*a23*a34*a45*a51);
    Vsusy = vg + 1.5*(l34 + l51) - 6.0;
    VL    = vg + 1.5*l34 + l51 - Log(-s12) - 3.0;
    
    _ComplexD tmp1 = a13*a23*a41*b24*b24/(a45*a51)*Ls1(-s23,-s51,-s34,-s51)/(s51*s51)
      +             a13*a23*a35*b25*b25/(a34*a45)*Ls1(-s12,-s34,-s51,-s34)/(s34*s34);
    
    _ComplexD tmp2 = L0(-s34, -s51)/s51;
    
    Fsusy = - 3.0*tmp1
      - 3.0*a13*a13/(2.0*a12*a34*a45*a51)*(a15*b52*a23 + a12*b24*a43)*tmp2;
    
    Ff = - a13*a13*b45/(a12*a45)*L0(-s12, -s34)/s34;
    Fs = (2.0*a15*b25*a34*a35*b45*b45/a45*L2(-s12, -s34)/(s34*s34*s34)
          - a13*a15*a34*b45*b45/(a12*a45)*L1(-s12, -s34)/(s34*s34)
          - a13*b24*b45/(a12*b12*b34*a45) + b24*b24*b25/(b12*b23*b34*a45))/3.0;
    
    FL = Fs - 2.0*tmp1
      - a12*a23*a34*a41*b24*b24*b24/(a45*a51)*Ls2(-s23,-s51,-s34,-s51)/(s51*s51*s51)
      - a12*a23*a35*a35*b25*b25*b25/(a34*a45)*Ls2(-s12,-s34,-s51,-s34)/(s34*s34*s34)
      
      - a13*a13*b24/(a45*a51)*Ls0(-s23,-s51,-s34,-s51)/s51
      - a13*a13*a35*b25/(a34*a45*a51)*Ls0(-s12,-s34,-s51,-s34)/s34

      - (a13*a23*a23*b25*b25*a15/(a12*a34*a45)
	 + 0.5*a13*a13*b12*a23*b25/(a34*a45))*L1(-s34,-s51)/(s51*s51)

      + 0.5*a13*a14*a23*b24*b24/(a45*a51)*L1(-s23,-s51)/(s51*s51)
      - 0.5*a13*a15*a34*b45*b45/(a12*a45)*L1(-s12,-s34)/(s34*s34)
      
      + 0.5*a13*a13*b24/(a45*a51)*tmp2
      - (2.0*a13*a13*b45/(a12*a45) + a13*a13*a35*b25/(a34*a45*a51))*L0(-s12,-s34)/s34
      
      + 0.5*a14*b24*b24*b45/(a45*b23*b34*s51)
      - a13*a23*b25*b45/(a12*a34*b34*a45*b15)
      - 0.5*a13*a13*b12*a23*b25/(a34*a34*b34*a45*a15*b15);
    
    res.Asusy = res.Atree*Vsusy + ICPLX*Fsusy;
    res.Asf = ICPLX*(Fs+Ff);
    res.AL = res.Atree*VL + ICPLX*FL;
  }
  
  void ampq2g3::A1pmp(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const
  {
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s45 = S(4,5),
      s51 = S(5,1);
    
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5),
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a34 = A(3,4), a45 = A(4,5),
      a51 = A(5,1), a54 = A(5,4), b12 = B(1,2), b13 = B(1,3), b15 = B(1,5),
      b23 = B(2,3), b25 = B(2,5), b34 = B(3,4), b35 = B(3,5), b45 = B(4,5),
      b53 = B(5,3); 

    _ComplexD Vsusy, VL, Fsusy, Ff, Fs, FL;
    _ComplexD l34 = Log(-s34), vg = Vg(p1, p2, p3, p4, p5);
    res.Atree = ICPLX*a14*a14*a14*a24/(a12*a23*a34*a45*a51);
    Vsusy = vg + 1.5*(Log(-s12) + l34) - 6.0;
    VL    = vg + 1.5*l34 - 2.0;
    
    Fsusy = -3.0*a13*a14*a24*b35*b35/(a12*a23)
      *Ls1(-s34,-s12,-s45,-s12)/(s12*s12)
      
      -3.0*a14*a24*a24*b25*b25/(a23*a34)
      *Ls1(-s51,-s34,-s12,-s34)/(s34*s34)
      
      -3.0*a14*a14*a24/(2.0*a12*a23*a34*a45*a51)
      *(a12*b25*a54 + a15*b53*a34)*L0(-s12,-s34)/s34;
    
    Fs = 2.0*a13*a34*a45*a15*b35*b35*b35*b35/a12
      *Ls3(-s34,-s12,-s45,-s12)/(s12*s12*s12*s12)
      + 2.0*a14*a15*a34*b35*b35*b35/(3.0*a12)*L2(-s12,-s34)/(s34*s34*s34)
      + 2.0*a13*a14*a45*b35*b35*b35/(3.0*a12)*L2(-s12,-s45)/(s45*s45*s45)
      - a13*a15*b35*b35*b35*b35/(3.0*a12*b34*b45*s12*s12);

    Ff = -std::pow(a14*b35,2)/a12*Ls1(-s34,-s12,-s45,-s12)/(s12*s12);
    
    FL = Fs - a13*a34*a45*b35*b35*b35/a23
      *Ls2(-s34,-s12,-s45,-s12)/(s12*s12*s12)
      
      - a12*a24*a24*a45*b25*b25*b25/(a23*a34)
      *Ls2(-s51,-s34,-s12,-s34)/(s34*s34*s34)

      + a14*b35*b35*(a14*a23*a15+2.0*a13*a15*a24+a13*a14*a25)/(a12*a23*a51)
      *Ls1(-s34,-s12,-s45,-s12)/(s12*s12)

      + a14*a24*b25*b25*(2.0*a15*a24+a14*a25)/(a23*a34*a51)
      *Ls1(-s51,-s34,-s12,-s34)/(s34*s34)
      
      + a14*a14*a14*a13*b13*b13/(a23*a45*a51)
      *Ls1(-s12,-s45,-s23,-s45)/(s45*s45)
      
      - (3.0*a24/(2.0*a12*a23*a34*a45*a51)*(a12*b25*a54+a15*b53*a34)
	 + b35/(a23*a15))*a14*a14*L0(-s12,-s34)/s34

      + a14*a14*a24*b23/(a45*a15*a23)*L0(-s12,-s45)/s45
      - a14*a14*a14*b13/(a23*a45*a15)*L0(-s23,-s45)/s45
      - a14*a14*a24*b25/(a23*a34*a15)*L0(-s34,-s51)/s51

      + 0.5*(a14*a34*b35*b35/a23*L1(-s12,-s34)/(s34*s34)
	     + a12*a14*a24*a45*b25*b25/(a23*a34*a15)*L1(-s34,-s12)/(s12*s12)
	     + a14*a24*a24*b25*b25/(a23*a34)*L1(-s34,-s51)/(s51*s51)
	     + a24*b25*b25*b35/(b12*a23*s34*b45)
	     + a13*a14*b15*b35*b35/(s12*a23*s34*b45)
	     + a14*a24*a24*b25*b25/(a23*a34*s34*s51)
	     - a14*a14*a24*b23/(s12*a23*a45*a15));

    res.Asusy = res.Atree*Vsusy + ICPLX*Fsusy;
    res.Asf = ICPLX*(Fs+Ff);
    res.AL = res.Atree*VL + ICPLX*FL;
  }

  void ampq2g3::A1ppm(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const
  {
    double s12 = S(1,2), s23 = S(2,3), s45 = S(4,5), s51 = S(5,1); 
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), 
      a21 = A(2,1), a23 = A(2,3), a25 = A(2,5), a31 = A(3,1), a34 = A(3,4),
      a35 = A(3,5), a45 = A(4,5), a51 = A(5,1), a53 = A(5,3), a54 = A(5,4),
      b12 = B(1,2), b13 = B(1,3), b14 = B(1,4), b15 = B(1,5), b24 = B(2,4),
      b32 = B(3,2), b34 = B(3,4), b43 = B(4,3), b45 = B(4,5);

    _ComplexD Vsusy, VL, Fsusy, Ff, Fs, FL;
    _ComplexD l45 = Log(-s45), vg = Vg(p1, p2, p3, p4, p5);
    res.Atree = ICPLX*a15*a15*a15*a25/(a12*a23*a34*a45*a51);
    Vsusy = vg + 1.5*(Log(-s12) + l45) - 6.0;
    VL    = vg + 1.5*l45 - 2.5;
    
    Ff = -a15*a15*b34/(a12*a34)*L0(-s12,-s45)/s45;
    Fsusy = 3.0*a15*a25/(2.0*a12*a23*a34*a45)
      *(a54*b43*a31 + a53*b32*a21)*L0(-s12,-s45)/s45;
    
    Fs = (- std::pow(a13*a45*b34,2)*b34/(a12*a34)*2.0
	      *L2(-s45,-s12)/(s12*s12*s12)
	      + a13*a45*a15*b34*b34/(a12*a34)*3.0*L1(-s45,-s12)/(s12*s12)
	      + 2.0*a15*a15*b34/(a12*a12*b12*a34)
	      + a15*b13*b24/(a12*b12*a34*b15)
	      + a14*b14*b24*b34/(a12*b12*a34*b45*b15))/3.0;
    
    FL = Fs - (a14*b14*b14/(a23*a34)*Ls1(-s45,-s23,-s51,-s23)/(s23*s23)
	 + a35*a13*b13*b13/(a23*a34*a45)*Ls1(-s12,-s45,-s23,-s45)/(s45*s45)
	 + b14/(a23*a34)*L0(-s23,-s51)/s51
	 - b12*a25/(a23*a34*a45)*L0(-s23,-s45)/s45)*a15*a15
      
      - (a15*a15*b34/(a12*a34) - 0.5*a15*a25
	 *(3.0*a54*b43*a31+a53*b32*a21)/(a12*a23*a34*a45))*L0(-s12,-s45)/s45

      + 0.5*(-a13*a35*a45*b34*b34/(a23*a34)*L1(-s12,-s45)/(s45*s45)
	     + b24*b34/(a34*b45*b15));

    res.Asusy = res.Atree*Vsusy + ICPLX*Fsusy;
    res.Asf = ICPLX*(Fs+Ff);
    res.AL = res.Atree*VL + ICPLX*FL;
  }

  void ampq2g3::A2mpp(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const
  {
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s51 = S(5,1); 
    
    _ComplexD a12 = A(1,2), a15 = A(1,5), a23 = A(2,3), a24 = A(2,4), 
      a32 = A(3,2), a34 = A(3,4), a41 = A(4,1), a42 = A(4,2), a45 = A(4,5),
      a51 = A(5,1), b12 = B(1,2), b13 = B(1,3), b15 = B(1,5), b23 = B(2,3),
      b34 = B(3,4), b35 = B(3,5), b45 = B(4,5), b54 = B(5,4);

    _ComplexD Vsusy, VL, Fsusy, Fs, FL;
    _ComplexD l23 = Log(-s23), l51 = Log(-s51), vg = Vg(p1, p2, p3, p4, p5);
    res.Atree = ICPLX*a12*a12*a12*a32/(a12*a23*a34*a45*a51);
    Vsusy = vg + 1.5*(l23 + l51) - 6.0;
    VL = vg + 1.5*l23 + Ls_1(-s34,-s51,-s23,-s51) 
      + Ls_1(-s12,-s34,-s51,-s34) - 3.0;
    
    _ComplexD tmp = L0(-s23,-s51)/s51;
    Fsusy = 1.5*a12*(a23*b34*a41+a24*b45*a51)/(a34*a45*a51)*tmp;
    Fs = b34*b35/(3.0*b12*b23*a45);
    FL = Fs + 0.5*a15*b54*b54*a42*a42/(a34*a45)*L1(-s23,-s51)/(s51*s51)
      + 2.0*a12*a24*b45*a51/(a34*a45*a51)*tmp
      + 0.5*b45*b35*b13/(b23*b12*b15*a45);
    
    res.Asusy = res.Atree*Vsusy + ICPLX*Fsusy;
    res.Asf = ICPLX*Fs;
    res.AL = res.Atree*VL + ICPLX*FL;
  }
  
  void ampq2g3::A2pmp(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const
  {
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s45 = S(4,5), 
      s51 = S(5,1); 

    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5),
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a34 = A(3,4), a45 = A(4,5),
      a51 = A(5,1), a54 = A(5,4), b15 = B(1,5), b25 = B(2,5), b34 = B(3,4),
      b35 = B(3,5), b45 = B(4,5), b53 = B(5,3); 

    _ComplexD Vsusy, VL, Fsusy, FL;
    _ComplexD l12 = Log(-s12), l34 = Log(-s34), vg = Vg(p1, p2, p3, p4, p5);
    res.Atree = ICPLX*a14*a14*a14*a34/(a12*a23*a34*a45*a51);
    Vsusy = vg + 1.5*(l12 + l34) - 6.0;
    VL = vg + 2.0*l34 - 0.5*l12 + Ls_1(-s23,-s45,-s12,-s45) - 3.0;
    
    Fsusy = -3.0*a13*a14*a34*b35*b35/(a12*a23)*Ls1(-s45,-s12,-s34,-s12)/(s12*s12)
      - 3.0*a14*a24*b25*b25/a23*Ls1(-s51,-s34,-s12,-s34)/(s34*s34)
      + 1.5*a14*a14*a34*(a12*b25*a54 + a15*b53*a34)/(a12*a23*a34*a45*a15)*L0(-s12,-s34)/s34;
  
    FL = -std::pow(a13*b35,2)*a34*a45*b35/(a12*a23)*Ls2(-s34,-s12,-s45,-s12)/(s12*s12*s12)
      - a12*a45*a24*b25*b25*b25/a23*Ls2(-s12,-s34,-s51,-s34)/(s34*s34*s34)
      - 2.0*a14*a34*a13*b35*b35/(a12*a23)*Ls1(-s34,-s12,-s45,-s12)/(s12*s12)
      - 2.0*a14*a24*b25*b25/a23*Ls1(-s12,-s34,-s51,-s34)/(s34*s34)
      - a14*a14*a13*b35/(a12*a23*a51)*Ls0(-s34,-s12,-s45,-s12)/s12
      - a14*a14*b25/(a23*a51)*Ls0(-s12,-s34,-s51,-s34)/s34
      - 0.5*a14*a24*b25*b25/a23*L1(-s51,-s34)/(s34*s34)
      + 0.5*a12*a14*a45*b15*b25/(a23*a25)*L1(-s12,-s34)/(s34*s34)
      + 0.5*a13*a14*a34*b35*b35/(a12*a23)*L1(-s12,-s34)/(s34*s34)
      - (2.0*a14*a14*a34*b35/(a12*a23*a45)-0.5*a12*a14*a45*b25/(a15*a23*a25))
      *L0(-s12,-s34)/s34 + 0.5*a14*a24*b25/(a23*a25*a34*b34)
      + 0.5*a13*a13*b35*b35*b35/(s12*a12*a23*b34*b45);
    
    res.Asusy = res.Atree*Vsusy + ICPLX*Fsusy;
    res.Asf = 0.0;
    res.AL = res.Atree*VL + ICPLX*FL;
  }
  
  void ampq2g3::A2ppm(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const
  {
    double s12 = S(1,2), s23 = S(2,3), s34 = S(3,4), s45 = S(4,5),
      s51 = S(5,1); 
    
    _ComplexD a12 = A(1,2), a13 = A(1,3), a15 = A(1,5), a21 = A(2,1), 
      a23 = A(2,3), a31 = A(3,1), a34 = A(3,4), a35 = A(3,5), a45 = A(4,5),
      a51 = A(5,1), a53 = A(5,3), a54 = A(5,4), b14 = B(1,4), b15 = B(1,5),
      b24 = B(2,4), b32 = B(3,2), b34 = B(3,4), b43 = B(4,3), b45 = B(4,5);
    
    _ComplexD Vsusy, VL, Fsusy, FL;
    _ComplexD l12 = Log(-s12), vg = Vg(p1, p2, p3, p4, p5);
    res.Atree = ICPLX*a15*a15*a15*a35/(a12*a23*a34*a45*a51);
    Vsusy = vg + 1.5*(l12 + Log(-s45)) - 6.0;
    VL = vg + 1.5*l12 + Ls_1(-s23,-s45,-s12,-s45) -3.0;

    Fsusy = 1.5*a15*a35*(a54*b43*a31+a53*b32*a21)/(a12*a23*a34*a45)*L0(-s12,-s45)/s45;
    FL = a15*a15*b24/(a12*a34)*Ls0(-s34,-s51,-s23,-s51)/s51
      + a15*a15*a13*b14/(a12*a23*a34)*Ls0(-s51,-s23,-s45,-s23)/s23
      - 0.5*a35*a45*std::pow(b34*a13,2)/(a12*a23*a34)*L1(-s12,-s45)/(s45*s45)
      - 2.0*a15*a35*a13*b34/(a12*a23*a34)*L0(-s12,-s45)/s45
      - 0.5*b24*b24/(a34*b15*b45);
    
    res.Asusy = res.Atree*Vsusy + ICPLX*Fsusy;
    res.Asf = 0.0;
    res.AL = res.Atree*VL + ICPLX*FL;
  }

  void ampq2g3::
  matrix_1loop(unsigned int nf, const _AmpPrim *A, const _AmpPrim *B, amp_1loop *res) const
  { 
    double onc2 = 1.0/Nc2, nf1 = nf/Nc+1.0/Nc2, nf2 = nf/Nc+1.0;
    _ComplexD Lsum(0.0), Csum(0.0), Asum(0.0);
    
    for(unsigned i = 0; i < 6; i++) {
      res[i].A0 = A[i].Atree;
      res[i].A1 = (1.0+onc2)*A[i].AL - onc2*A[i].Asusy - nf1*A[i].Asf;
      
      Lsum += A[i].AL + B[i].AL;
      if(i < 3) Csum += A[i].Asusy + B[i].Asusy;
      else Asum += A[i].Asusy + B[i].Asusy;
    }
    
    for(unsigned i = 0; i < 6; i++) {
      res[i].A3 =  Lsum - B[i].Asusy - B[(i+3)%6].Asusy;
      res[i].A4 = -Lsum + (i < 3 ? Csum : Asum) + nf2*(A[i].Asf + B[i].Asf); 
    }
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
  A2##h4##h5##h3(p1, p4, p2, p5,p3, B[1]);	\
  A2##h5##h3##h4(p1, p5, p2, p3,p4, B[2]);	\
  		               	                \
  A2##h3##h5##h4(p1, p3, p2, p5,p4, B[3]);	\
  A2##h5##h4##h3(p1, p5, p2, p4,p3, B[4]);	\
  A2##h4##h3##h5(p1, p4, p2, p3,p5, B[5])
  
  
  void ampq2g3::
  matrix_1loop_pmppm(unsigned int nf, int p2, int p3, int p4,
		     int p5, int p1, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(m,p,p); 
    matrix_1loop(nf, A, B, res);
  }

  void ampq2g3::
  matrix_1loop_ppmpm(unsigned int nf, int p2, int p3, int p4,
		     int p5, int p1, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(p,m,p); 
    matrix_1loop(nf, A, B, res);
  }

  void ampq2g3::
  matrix_1loop_pppmm(unsigned int nf, int p2, int p3, int p4,
		     int p5, int p1, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(p,p,m); 
    matrix_1loop(nf, A, B, res);
  }

  void ampq2g3::
  matrix_1loop_ppmmm(unsigned int nf, int p2, int p3, int p4,
		     int p5, int p1, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(p,m,m); 
    matrix_1loop(nf, A, B, res);
  }

  void ampq2g3::
  matrix_1loop_pmpmm(unsigned int nf, int p2, int p3, int p4,
		     int p5, int p1, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(m,p,m); 
    matrix_1loop(nf, A, B, res);
  }

  void ampq2g3::
  matrix_1loop_pmmpm(unsigned int nf, int p2, int p3, int p4,
		     int p5, int p1, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(m,m,p); 
    matrix_1loop(nf, A, B, res);
  }

  const unsigned int ampq2g3::s[5][6] = 
    {{5,4,3,2,1,0}, {3,5,4,0,2,1}, {1,2,0,4,5,3},
     {4,3,5,1,0,2}, {2,0,1,5,3,4}};
  
  
  double ampq2g3::amp1loop(amp_1loop *M) const
  {
    _ComplexD tmp;
    double ret_val = 0.0;
    
    for(unsigned int i = 0; i < 6; i++) {
      tmp = Na2*M[i].A1
	- Na*(M[s[0][i]].A1 + M[s[1][i]].A1
	      - M[s[2][i]].A3 - M[i].A3)
	+ (Nc2+1.0)*M[s[3][i]].A1 + (Nc2-2.0)*M[i].A4
	+ M[s[2][i]].A1 + M[s[4][i]].A1
	- M[s[4][i]].A3 - 2.0*M[s[3][i]].A4;
      
      ret_val += real(conj(M[i].A0)*tmp);
    }
    
    return Na*ret_val/Nc;
  }			      
							
  double ampq2g3::
  su3_1loop(unsigned int nf, int p1, int p2, int p3, int p4, int p5) const
  {
    static amp_1loop M[6];
    double ret_val;

    matrix_1loop_pppmm(nf, p1,p3,p4,p5,p2, M); ret_val  = amp1loop(M);
    matrix_1loop_ppmpm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_pmppm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_ppmmm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_pmpmm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_pmmpm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    
    swap();
    matrix_1loop_pppmm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_ppmpm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_pmppm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_ppmmm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_pmpmm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    matrix_1loop_pmmpm(nf, p1,p3,p4,p5,p2, M); ret_val += amp1loop(M);
    swap();

    return ret_val;
  }

  double ampq2g3::
  su3_1loop_mch(unsigned int nf, int p1, int p2, int p3, int p4, int p5) const
  {
    static amp_1loop M[6];
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (6*_M_rng());
    
    if(hpm == 1) swap();
    switch(hel){
    case 0: matrix_1loop_pppmm(nf, p1,p3,p4,p5,p2, M); break;
    case 1: matrix_1loop_ppmpm(nf, p1,p3,p4,p5,p2, M); break;
    case 2: matrix_1loop_pmppm(nf, p1,p3,p4,p5,p2, M); break;
    case 3: matrix_1loop_ppmmm(nf, p1,p3,p4,p5,p2, M); break;
    case 4: matrix_1loop_pmpmm(nf, p1,p3,p4,p5,p2, M); break;
    case 5: matrix_1loop_pmmpm(nf, p1,p3,p4,p5,p2, M); break;
    }
    if(hpm == 1) swap();

    return 12.0*amp1loop(M);
  }
}  //  namespace nlo

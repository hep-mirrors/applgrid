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
#include "ampq4g1.h"
#include "defmacros.h"

#define E(i,j) (s##i##j/(s##i##5*s##j##5))
#define F(i,j) (A(i,j)/(A(i,5)*A(5,j)))
#define G(i,j) conj(f##i##j)


namespace nlo {


  void ampq4g1::
  su3_tree(int p1, int p2, int p3, int p4, int p5, double *out) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s32 = S(3,2), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5),
      s45 = S(4,5);

    double e12 = E(1,2), e13 = E(1,3), e14 = E(1,4), e32 = E(3,2),
      e24 = E(2,4), e34 = E(3,4), AA, BB, CC, pa, pb, pc; 

    AA = (Nc2-2.0)*(e14+e32)-e12-e34+2.0*(e13+e24);
    BB = (Nc2-2.0)*(e12+e34)-e14-e32+2.0*(e13+e24);
    CC = ((Nc2+1.0)*(e13+e24)-e12-e14-e32-e34)/Nc;
    
    pa = (s13*s13 + s24*s24 + s14*s14 + s32*s32)/(s34*s12);
    pb = (s13*s13 + s24*s24 + s12*s12 + s34*s34)/(s14*s32);
    pc = (s24*s24 + s13*s13)/(s12*s34*s14*s32)*(s12*s34 + s14*s32 - s13*s24);

    out[0] = 2.0*(Nc2-1.0)/Nc*pa*AA;
    out[1] = 2.0*(Nc2-1.0)/Nc*(pa*AA + pb*BB + pc*CC);
  }

  void ampq4g1::
  ampcc12(int p1, int p2, int p3, int p4, int p5, double *cc) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5),
      s45 = S(4,5);

    double e12 = E(1,2), e13 = E(1,3), e14 = E(1,4), e23 = E(2,3),
      e24 = E(2,4), e34 = E(3,4), A12, B12, C12, pa, pb, pc; 

    _ComplexD f12 = F(1,2), f14 = F(1,4), f32 = F(3,2), f34 = F(3,4); 
    _ComplexD g12 = G(1,2), g32 = G(3,2), g34 = G(3,4); 
    _ComplexD g = A(1,2)*A(3,4)/(A(1,4)*A(3,2));
    _ComplexD G = 0.5*(g/Nc+Nc*conj(g))/real(g);

    A12 = ((Nc2-2.0)*(e23+e14) - (Nc2+1.0)*e12 + 2.0*(e24+e13) - e34)/Nc;
    B12 = ((3.0-Nc2)*(e24+e13-e12) + (3.0*Nc2-3.0-Nc2*Nc2)*e34 
	   + (Nc2-2.0)*(e23 + e14))/Nc;
    
    C12 = 2.0*real(G*(f14+f32)*(g12-Na*g34)/Nc+f14*g32)+Na*e34-e12;
    
    pa = (s13*s13 + s24*s24 + s14*s14 + s23*s23)/(s34*s12);
    pb = (s13*s13 + s24*s24 + s12*s12 + s34*s34)/(s14*s23);
    pc = (s24*s24 + s13*s13)/(s12*s34)*2.0*real(g);

    cc[0] = (Nc2-1.0)/Nc*pa*A12;
    cc[1] = (Nc2-1.0)/Nc*(pa*A12 + pb*B12 - pc*C12);
  }

  void ampq4g1::
  ampcc23(int p1, int p2, int p3, int p4, int p5, double *cc) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5),
      s45 = S(4,5); 
    
    double e12 = E(1,2), e13 = E(1,3), e14 = E(1,4), e23 = E(2,3),
      e24 = E(2,4), e34 = E(3,4), A23, B23, C23, pa, pb, pc; 
    
    _ComplexD f12 = F(1,2), f14 = F(1,4), f32 = F(3,2), f34 = F(3,4); 
    _ComplexD g12 = G(1,2), g34 = G(3,4), g = A(1,2)*A(3,4)/(A(1,4)*A(3,2));
    _ComplexD G = 0.5*(g/Nc+Nc*conj(g))/real(g);
    
    A23 = ((3.0-Nc2)*(e24+e13-e23) + (3.0*Nc2-3.0-Nc2*Nc2)*e14 
	      + (Nc2-2.0)*(e12+e34))/Nc;
    B23 = ((Nc2-2.0)*(e12+e34) - (Nc2+1.0)*e23 + 2.0*(e24+e13) - e14)/Nc;
    C23 = 2.0*real(G*(g12+g34)*(f32-Na*f14)/Nc+f12*g34)+Na*e14-e23;


    pa = (s13*s13 + s24*s24 + s14*s14 + s23*s23)/(s34*s12);
    pb = (s13*s13 + s24*s24 + s12*s12 + s34*s34)/(s14*s23);
    pc = (s24*s24 + s13*s13)/(s12*s34)*2.0*real(g);

    cc[0] = (Nc2-1.0)/Nc*pa*A23;
    cc[1] = (Nc2-1.0)/Nc*(pa*A23 + pb*B23 - pc*C23);
  }

  void ampq4g1::
  ampcc24(int p1, int p2, int p3, int p4, int p5, double *cc) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5),
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5),
      s45 = S(4,5); 

    double e12 = E(1,2), e13 = E(1,3), e14 = E(1,4), e23 = E(2,3), 
      e24 = E(2,4), e34 = E(3,4), A24, B24, C24, pa, pb, pc; 
    
    _ComplexD f12 = F(1,2), f14 = F(1,4), f32 = F(3,2), f34 = F(3,4);
    _ComplexD g12 = G(1,2), g32 = G(3,2), g34 = G(3,4); 
    _ComplexD g = A(1,2)*A(3,4)/(A(1,4)*A(3,2));
    _ComplexD G = 0.5*(g/Nc+Nc*conj(g))/real(g);

    A24 = ((3.0-Nc2)*(e14 + e23 - e24) + 2.0*e12 - (3.0+Nc2)*e13 + 2.0*e34)/Nc;
    B24 = ((3.0-Nc2)*(e12 + e34 - e24) + 2.0*e14 - (3.0+Nc2)*e13 + 2.0*e23)/Nc;
    C24 = real(G*((f14*g34+f32*g12)-Na*(f14*g12+f32*g34))/Nc+f14*g32+f12*g34);
    
    pa = (s13*s13 + s24*s24 + s14*s14 + s23*s23)/(s34*s12);
    pb = (s13*s13 + s24*s24 + s12*s12 + s34*s34)/(s14*s23);
    pc = (s24*s24 + s13*s13)/(s12*s34)*2.0*real(g);

    cc[0] = (Nc2-1.0)/Nc*pa*A24;
    cc[1] = (Nc2-1.0)/Nc*(pa*A24 + pb*B24 + 2.0*pc*C24);
  }

  void ampq4g1::
  ampcc25(int p1, int p2, int p3, int p4, int p5, double *cc) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5),
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), 
      s35 = S(3,5), s45 = S(4,5);

    double e12 = E(1,2), e23 = E(2,3), e24 = E(2,4),
      A25, B25, C25, pa, pb, pc; 
    
    _ComplexD f12 = F(1,2), f14 = F(1,4), f32 = F(3,2), f34 = F(3,4), 
      g12 = G(1,2), g34 = G(3,4), g = A(1,2)*A(3,4)/(A(1,4)*A(3,2));
    _ComplexD G = 0.5*(g/Nc+Nc*conj(g))/real(g);
    
    A25 = Nc*(e12 - 2.0*e24 + (2.0-Nc2)*e23);
    B25 = Nc*(e23 - 2.0*e24 + (2.0-Nc2)*e12);
    C25 = Nc*(-Nc*(e12+e23) + 2.0*real(G*((f14+f32)*g12+(f32-f14)*g34)));

    pa = (s13*s13 + s24*s24 + s14*s14 + s23*s23)/(s34*s12);
    pb = (s13*s13 + s24*s24 + s12*s12 + s34*s34)/(s14*s23);
    pc = (s24*s24 + s13*s13)/(s12*s34)*2.0*real(g);

    cc[0] = (Nc2-1.0)/Nc*pa*A25;
    cc[1] = (Nc2-1.0)/Nc*(pa*A25 + pb*B25 + pc*C25);
  }

  void ampq4g1::
  amphtree(int p1, int p2, int p3, int p4, int p5, _ComplexD *cc) const
  {
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a23 = A(2,3),
      a24 = A(2,4), a34 = A(3,4), f12 = F(1,2), f14 = F(1,4), f32 = F(3,2),
      f34 = F(3,4);
    
    _ComplexD fa = f14*f14+f32*f32, fb = f12*f12+f34*f34, 
      fc = (f14+f32)*(f12+f34), ha, hb, hc;
    
    ha = 2.0*(a13*a13*a24*a24+a14*a14*a23*a23)/(a12*a12*a34*a34);
    hb = 2.0*(a13*a13*a24*a24+a12*a12*a34*a34)/(a14*a14*a23*a23);
    hc = 4.0*a13*a13*a24*a24/(a12*a34*a14*a23);
    
    cc[0] = ha*(Nc2*fa+fb-2.0*fc)/Nc;
    cc[1] = cc[0] + hb*(Nc2*fb+fa-2.0*fc)/Nc - hc*(fa+fb-(1.0+1.0/Nc2)*fc);
    cc[0] *= Nc2-1.0; cc[1] *= Nc2-1.0;
  }
  
  
#define __cc_switch(i,j)					\
  (pi == p##i && pj == p##j) || (pi == p##j && pj == p##i)

  void ampq4g1::
  su3_cc(int pi, int pj, int p1, int p2, int p3, int p4, int p5, _Pair *out) const
  {
    double cc[2];
    _ComplexD hcc[2];
    hcc[0] = hcc[1] = 0.0;
    
    if(__cc_switch(1,2)) ampcc12(p1, p2, p3, p4, p5, cc);
    else if(__cc_switch(1,3)) ampcc24(p2, p1, p4, p3, p5, cc);
    else if(__cc_switch(1,4)) ampcc23(p3, p4, p1, p2, p5, cc);
    else if(__cc_switch(1,5)) ampcc25(p2, p1, p4, p3, p5, cc);
    else if(__cc_switch(2,3)) ampcc23(p1, p2, p3, p4, p5, cc);
    else if(__cc_switch(2,4)) ampcc24(p1, p2, p3, p4, p5, cc);
    else if(__cc_switch(2,5)) ampcc25(p1, p2, p3, p4, p5, cc);
    else if(__cc_switch(3,4)) ampcc12(p3, p4, p1, p2, p5, cc);
    else if(__cc_switch(3,5)) ampcc25(p4, p3, p2, p1, p5, cc);
    else if(__cc_switch(4,5)) ampcc25(p3, p4, p1, p2, p5, cc);
    else throw("Error in ampq2g3::su3_cc");

    if(pj == p5) {
      amphtree(p1, p2, p3, p4, p5, hcc);
      hcc[0] *= -Nc/4.0; hcc[1] *= -Nc/4.0;
    }
    
    out[0] = _Pair(cc[0], hcc[0]);
    out[1] = _Pair(cc[1], hcc[1]);
  }
  
  void ampq4g1::
  ampcc(int pi, int pj, int p1, int p2, int p3, int p4, int p5, double *cc) const
  {
    if(__cc_switch(1,2)) ampcc12(p1, p2, p3, p4, p5, cc);
    else if(__cc_switch(1,3)) ampcc24(p2, p1, p4, p3, p5, cc);
    else if(__cc_switch(1,4)) ampcc23(p3, p4, p1, p2, p5, cc);
    else if(__cc_switch(1,5)) ampcc25(p2, p1, p4, p3, p5, cc);
    else if(__cc_switch(2,3)) ampcc23(p1, p2, p3, p4, p5, cc);
    else if(__cc_switch(2,4)) ampcc24(p1, p2, p3, p4, p5, cc);
    else if(__cc_switch(2,5)) ampcc25(p1, p2, p3, p4, p5, cc);
    else if(__cc_switch(3,4)) ampcc12(p3, p4, p1, p2, p5, cc);
    else if(__cc_switch(3,5)) ampcc25(p4, p3, p2, p1, p5, cc);
    else if(__cc_switch(4,5)) ampcc25(p3, p4, p1, p2, p5, cc);
    else throw("Error in ampq2g3::su3_cc");
  }
  

#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampq4g1::su3_kp(unsigned int nf, int pa, int pb, int p1, int p2, 
		       int p3, int p4, int p5, su3_kp_i2 *res, double al) const
  {
    double cc[2], xi, xj, xq = Gq/Cf, xg = Gg(nf)/Nc;
    int p[5] = {p1,p2,p3,p4,p5};

    double Ta = (pa == p5 ? Ca : Cf);
    double Tb = (pb == p5 ? Ca : Cf);
    double x[5] = {xq, xq, xq, xq, xg};

    res[0].tree=res[0].pa=res[0].pb=res[0].ga=res[0].gb=res[0].loop = 0.0;
    res[1].tree=res[1].pa=res[1].pb=res[1].ga=res[1].gb=res[1].loop = 0.0;
    for(unsigned i = 0; i < 5; i++)
      for(unsigned j = i+1; j < 5; j++) {
	ampcc(p[i], p[j], p1, p2, p3, p4, p5, cc);
	
	//----- loop contribution -----
	if(i == 4) {
	  xi = Xg(Sij(p[i],p[j]), 1.0, nf);
	  xj = Xq(Sij(p[i],p[j]), 1.0);
	} else {
	  xi = Xq(Sij(p[i],p[j]), 1.0);
	  xj = j == 4 ? Xg(Sij(p[i],p[j]), 1.0, nf) : xi;
	}

	res[0].loop += (xi+xj)*cc[0];
	res[1].loop += (xi+xj)*cc[1];

	//----- finite contributions -----
	if(p[i] == pa || p[j] == pa) {
	  int iq = p[i] == pa ? j : i;
	  res[0].tree -= cc[0];
	  res[1].tree -= cc[1];

	  xi = X(pa, p[iq]);
	  res[0].pa += xi*cc[0];
	  res[1].pa += xi*cc[1];

	  if(p[iq] != pb){
	    res[0].ga += x[iq]*cc[0];
	    res[1].ga += x[iq]*cc[1];
	  } else {
	    res[0].cca = -cc[0]/Ta;
	    res[1].cca = -cc[1]/Ta;
	  }
	}
	
	if(p[i] == pb || p[j] == pb) {
	  int iq = p[i] == pb ? j : i;
	  
	  xj = X(pb, p[iq]);
	  res[0].pb += xj*cc[0];
	  res[1].pb += xj*cc[1];

	  if(p[iq] != pa) {
	    res[0].gb += x[iq]*cc[0];
	    res[1].gb += x[iq]*cc[1];
	  } else {
	    res[0].ccb = -cc[0]/Tb;
	    res[1].ccb = -cc[1]/Tb;
	  }
	}
      }
    
    res[0].tree /= Ta; res[0].pa /= Ta; res[0].pb /= Tb;
    res[1].tree /= Ta; res[1].pa /= Ta; res[1].pb /= Tb;
    
    xi = Gg(nf)+Kg(nf,al) + 4.0*(Gq+Kq(al)) - 2.0*Cf + Ca/3.0;
    res[0].loop += xi*res[0].tree;
    res[1].loop += xi*res[1].tree;
  }

#undef F

  //   1-loop amplitude
#define ICPLX std::complex<double>(0.0,1.0)

#define F0(sij, skl) (Log(-(sij))-Log(-(skl)))
#define F1(sij, skl) (L0(-(sij), -(skl))/(skl))
#define F2(sij, skl) (L1(-(sij), -(skl))/((skl)*(skl)))
#define F3(sij, skl) (L2(-(sij), -(skl))/((skl)*(skl)*(skl)))
#define F(smn, sij, sjk) Ls_1(-(sij), -(smn), -(sjk), -(smn))
#define H(i,j,k,l) std::pow(a##i##k*a##j##l/(a##i##l*a##j##k), 2)
  


  std::complex<double> 
  ampq4g1::D34ppp(unsigned int, int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s14 = S(1,4), s23 = S(2,3), s34 = S(3,4), 
      s45 = S(4,5);

    _ComplexD a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), a21 = A(2,1),
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a35 = A(3,5), a43 = A(4,3),
      a45 = A(4,5), b13 = B(1,3), b15 = B(1,5), b35 = B(3,5);

    _ComplexD D = Nc*(a23*a23*a45*b35/(a21*a35*a35) 
		      + 3.0*a24*a23*b35/(2.0*a21*a35))*F1(s45, s12)
      - (Nc2+1.0)*(a24*a43*b13/(2.0*a45*a35*s12) + a24*b15/(a35*s12)
		   + 0.5*a23*a43*b13*b35/a35*F2(s45, s12)
		   - a24*a23/(a21*a35*a35)*F0(s34, s12) 
		   + a23*a45*b15/(a35*a35)*F1(s34, s12)
		   - a25*a43*b15*b35/a35*F2(s34, s12))/Nc
      - (a24/(a15*a35)*F0(s45, s23) + a24*a23/(a21*a35*a35)*F0(s45, s12)
	 + 3.0*a24*a24/(2.0*a21*a45*a35)*F0(s34, s45)
	 - (a23*a15*a43*a43*b13/(a13*a45*a35*a35)
	    - a24*a43*b13/(2.0*a45*a35))*F1(s45, s12)
	 + a25*a14*b15/(a15*a35)*F1(s14, s23) 
	 - a23*a14*a14*b13/(a13*a15*a45)*F1(s45, s23))/Nc;

    return ICPLX*D;
  }

  std::complex<double> 
  ampq4g1::D34mpp(unsigned int, int p1, int p2, int p3, int p4, int p5) const
  {
    double s13 = S(1,3), s21 = S(2,1), s24 = S(2,4), s43 = S(4,3), 
      s45 = S(4,5); 

    _ComplexD a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), a21 = A(2,1),
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a35 = A(3,5), a43 = A(4,3),
      a45 = A(4,5), b23 = B(2,3), b25 = B(2,5), b35 = B(3,5);

    _ComplexD D = Nc*(3.0*a14*a14/(2.0*a21*a45*a35)*F0(s45, s21)
		      + 3.0*a14*a43*b23/(a45*a35)*F1(s45, s21))
      + 3.0*a14*a14/(2.0*Nc*a21*a45*a35)*F0(s43, s45)
      - (Nc2+1.0)*(a14*b25/(a35*s21) + a14/(a25*a35)*F0(s13, s45)
		   + a14*a13/(a21*a35*a35)*F0(s43, s45)
		   - a24*a15*b25/(a25*a35)*F1(s13, s24)
		   + a13*a45*b25/(a35*a35)*F1(s43, s21)
		   - 2.0*a15*a43*a43*b23/(a45*a35*a35)*F1(s45, s21)
		   + a25*a13*a43*a43*b23/(a23*a45*a35*a35)*F1(s45, s21)
		   + a24*a24*a13*b23/(a23*a25*a45)*F1(s45, s13)
		   - a15*a43*b25*b35/a35*F2(s43, s21)
		   - a21*a43*a43*b23*b23/(2.0*a45*a35)*F2(s45, s21))/Nc;
    
    return ICPLX*D;
  }

  std::complex<double> 
  ampq4g1::D32ppp(unsigned int Nf, int p1, int p2, int p3, int p4, int p5) const
  {
    double s14 = S(1,4), s21 = S(2,1), s25 = S(2,5), s43 = S(4,3); 

    _ComplexD a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), a21 = A(2,1), 
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a35 = A(3,5), a41 = A(4,1), 
      a43 = A(4,3), a45 = A(4,5), b13 = B(1,3), b15 = B(1,5), b35 = B(3,5), 
      b45 = B(4,5);

    _ComplexD D = 2.0*(Nc-Nf)*(a24*a24*b45/(a21*a35)*F1(s43, s21)
			       - 0.5*a24*b15*b35*F2(s43, s21)
			       + a25*a43*b15*b35*b35*F3(s43, s21))/3.0
      + (Nc2+1.0)*(-a24*a24*a13/(2.0*a21*a15*a43*a35) - a24*b15/(2.0*a35*s21) 
		   - a21*a45*a24/(a15*a15*a43*a25)*F0(s43, s25) 
		   + a21*a21*b13*(a41*a35+a43*a15)/(a25*a13*a15*a15)*F1(s43, s25)
		   - 2.0*a23*a23*a14*b13/(a25*a13*a35)*F1(s14, s25) 
		   + (-1.5*a24*b15/a35 + a25*a14*b35/(a15*a15)
		      - 2.0*a23*a45*b35/(a15*a35))*F1(s43, s21)
		   + 0.5*a43*(a25*b15*b35/a35*F2(s43, s21)
			      - a25*b35*b35/a15*F2(s21, s43)
			      - a21*a21*b13*b13/(a25*a15)*F2(s25, s43)))/Nc
      + Nc*(a24*a23*a45/(a25*a15*a43*a35)*F0(s43, s14)
	    + (a25*a14/(a15*a24) + 0.5)*a24*a24*a13/(a21*a15*a43*a35)*F0(s43, s21)
	    + 1.5*a24*a24/(a25*a15*a43)*F0(s43, s25) 
	    + (a14*a35/(a15*a43) - 2.0)*a24*a23/(a25*a13*a35)*F0(s14, s25)
	    - 3.0*a23*(a43*b35/(a13*a35)*F1(s14, s25) 
		       + a25*a14*b35/(a21*a15*a35)*F1(s43, s21) 
		       - a21*a14*b13/(a25*a13*a15)*F1(s43, s25)))
      + (a21*a24*a45/(a25*a15*a15*a43)*F0(s43, s21) 
	 + a24*a24/(2.0*a25*a15*a43)*F0(s25, s21) 
	 + 2.0*a24/(a15*a35)*F0(s25, s21))/Nc;

    return ICPLX*D;
  }

  std::complex<double> 
  ampq4g1::D32mpp(unsigned int __nf, int p1, int p2, int p3, int p4, int p5) const
  {
    double nf = (double) __nf, s13 = S(1,3), s21 = S(2,1), s24 = S(2,4), s25 = S(2,5), 
      s35 = S(3,5), s43 = S(4,3); 

    _ComplexD a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), a21 = A(2,1), 
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a35 = A(3,5), a41 = A(4,1),
      a43 = A(4,3), a45 = A(4,5), b25 = B(2,5), b35 = B(3,5), b43 = B(4,3),
      b45 = B(4,5); 

    _ComplexD D = (Nc2+1.0)*(a23*a14*a14/(2.0*a21*a25*a43*a35)
			     - a14*b25/(2.0*a35*s21)
			     + a15*b35*b35/(2.0*a25*b43*s21)
			     - a14*b35/a25*F1(s43, s21) 
			     + a15*a43*b35*(s35-s25)/(2.0*a25*a35)*F2(s43, s21))/Nc
      + nf*(-2.0*a14*a14*a13/(a21*a15*a43*a35)*F0(s43, s21)
	    - 2.0*a14*a14*b25/(a15*a43)*F1(s43, s21)
	    + a14*b25*b35*F2(s21, s43) 
	    - 2.0*a21*a45*b25*b25*b35*F3(s43, s21))/3.0
      - (2.0*a14/(a25*a35)*F0(s24, s43)
	 - 2.0*a13*a45*b35/(a25*a35)*F1(s13, s24) 
	 -(1.5*a41/a45 + 2.0*a21/a25)*a45*b25/a35*F1(s43, s21))/Nc
      + Nc*(-3.0*a14*a14*a24/(2.0*a21*a43*a25*a45)*F0(s43, s21)
	    + (1.5*a14*a15*b35*(a23*a45+a24*a35)/(a21*a25*a45*a35)
	       - 2.0*a14*a14*b45/(3.0*a21*a35))*F1(s43, s21)
	    - a14*b25*b35/3.0*F2(s43, s21) 
	    + 2.0*a15*a43*b25*b35*b35/3.0*F3(s43, s21));
    
    return ICPLX*D;
  }
  
  std::complex<double> 
  ampq4g1::D32pmp(unsigned int Nf, int p1, int p2, int p3, int p4, int p5) const
  {
    double s13 = S(1,3), s21 = S(2,1), s24 = S(2,4), s25 = S(2,5),
      s35 = S(3,5), s43 = S(4,3); 

    _ComplexD a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), a21 = A(2,1), 
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a35 = A(3,5), a43 = A(4,3), 
      a45 = A(4,5), b14 = B(1,4), b15 = B(1,5), b25 = B(2,5), b35 = B(3,5),
      b43 = B(4,3), b45 = B(4,5); 

    _ComplexD D = 2.0*(Nf-Nc)*(a23*a23*b45/(a21*a35)*F1(s43, s21)
			       + 0.5*a23*b15*b45*F2(s43, s21)
			       + a25*a43*b15*b45*b45*F3(s43, s21))/3.0
      + Nc*(a23*a23*a23/(2.0*a21*a25*a43*a35)
	    -(a25*a13*a13*b15/(a15*a15*a43) + 1.5*a23*a13*b15/(a15*a43))*F1(s43, s25)
	    -(a25*a13*a13*b15/(a15*a15*a43) 
	      + a23*a14*a35*b15/(2.0*a15*a43*a45) - a23*a23*b35/(2.0*a21*a45)
	      + a24*a24*a35*b45/(a21*a45*a45))*F1(s43, s21)
	    - (1.5*a24*a23*b45/(a21*a45) + a24*a24*a35*b45/(a21*a45*a45))*F1(s35, s21)
	    - a24*a43*b14*b45/(2.0*a45)*F2(s35, s21)
	    + a25*a43*b15*b45/(2.0*a45)*F2(s43, s21) 
	    + a21*a35*b15*b45/(2.0*a15)*F2(s21, s43)
	    - a21*a13*b14*b15/(2.0*a15)*F2(s25, s43))
      + (a23*a23*a13/(2.0*a21*a15*a43*a35) + a23*a43*b14/(2.0*a45*a35*s21)
	 + a25*a43*b45/(2.0*a15*a45*s21) + a25*b45*b45/(2.0*a15*b43*s21)
	 - 2.0*a23*a13/(a14*a15*a35)*F0(s24, s21) 
	 - a23*a23/(2.0*a25*a15*a43)*F0(s25, s21)
	 - 2.0*a21*a23/(a25*a14*a15)*F0(s13, s25)
	 + (a25*a43*a23/(a45*a45*a35*a21) - a23*a21*a35/(a25*a15*a15*a43))*F0(s43, s21)
	 + a21*a35*b45/(a15*a15*s43)*F0(s43, s25) 
	 + (2.0*a23*a43/(a14*a45*a35) - a25*a43*b15/(a45*a45*s21)
	    - a23*a23/(2.0*a21*a45*a35))*F0(s35, s21)
	 - (a25*a43*b15/(a45*a45) + 2.0*(a23*a45+a25*a43)*b45/(a15*a45)
	    + a25*a25*a13*a43*b45/(2.0*a21*a15*a45*a35)
	    - a25*a13*b45/(a15*a15))*F1(s43, s21)
	 - 2.0*a24*a35*b45/(a15*a45)*F1(s13, s24)
	 + 2.0*a24*a43*b45/(a14*a45)*F1(s13, s25)
	 + (2.0*a21*a21*a43*b14/(a25*a14*a15)
	    + a21*a21*a35*b25*b14/(a15*a15*s43))*F1(s43, s25)
	 + 2.0*a21*a13*b15/(a14*a15)*F1(s35, s24)
	 + (2.0*a21*a43*a43*b14/(a14*a45*a35)
	    + a25*a43*a43*b14*b35/(a45*a45*s21))*F1(s35, s21)
	 + a21*a21*a43*b14*b14/(2.0*a25*a15)*F2(s25, s43)
	 + a25*a43*b45*(s13-s24)/(2.0*a15*a45)*F2(s43, s21) 
	 + a21*a43*a43*b14*b14/(2.0*a45*a35)*F2(s35, s21))/Nc;
    
    return ICPLX*D;
  }

  std::complex<double> 
  ampq4g1::E32ppp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5), 
      s45 = S(4,5); 

    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), 
      a21 = A(2,1), a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a32 = A(3,2), 
      a34 = A(3,4), a35 = A(3,5), a51 = A(5,1), a52 = A(5,2), a53 = A(5,3), 
      a54 = A(5,4); 

    _ComplexD h2134 = H(2,1,3,4), h2144 = H(2,1,4,4), h2154 = H(2,1,5,4), 
      h2234 = H(2,2,3,4), h2244 = H(2,2,4,4), h2254 = H(2,2,5,4), 
      h2514 = H(2,5,1,4), h2534 = H(2,5,3,4), h2544 = H(2,5,4,4); 

    _ComplexD f123 = F(s45,s12,s23), f124 = F(s35,s12,s24), 
      f125 = F(s34,s12,s25), f134 = F(s25,s13,s34), f135 = F(s24,s13,s35), 
      f143 = F(s25,s14,s34), f145 = F(s23,s14,s45), f153 = F(s24,s15,s35), 
      f154 = F(s23,s15,s45), f213 = F(s45,s12,s13), f214 = F(s35,s12,s14), 
      f215 = F(s34,s12,s15), f234 = F(s15,s23,s34), f235 = F(s14,s23,s35), 
      f243 = F(s15,s24,s34), f245 = F(s13,s24,s45), f253 = F(s14,s25,s35), 
      f254 = F(s13,s25,s45), f351 = F(s24,s35,s15), f435 = F(s12,s34,s35), 
      f451 = F(s23,s45,s15), f513 = F(s24,s15,s13), f514 = F(s23,s15,s14), 
      f521 = F(s34,s25,s12), f523 = F(s14,s25,s23), f524 = F(s13,s25,s24), 
      f534 = F(s12,s35,s34), f543 = F(s12,s45,s34);


    _ComplexD E = -a32*(Nc*(f253 + h2244*f214 + h2154*f125
			    + h2544*f534 + h2134*f143)
			-(-h2144*(f134 + f124) + f154
			  -h2234*(f243 + f213) + f253
			  +h2244*(f234 + f214) - f254
			  +h2134*(f123 + f143) - f153)/Nc)/(a35*a52)
      
      + a34*(Nc*h2544*f534
	     -(-f254 + h2254*f215 + f451 - h2154*f125 - f435 
	       - h2144*(f134 - f124) + h2244*(f234 - f214)
	       + h2254*f245 - f235 + f135 - h2154*f145)/Nc)/(a35*a54)
      
      + a12*(Nc*h2154*f125
	     -(f451 - h2544*f534 - f351 + h2534*f543 - h2514*f521
	       - h2144*(f124 - f134) + h2134*(f123 - f143)
	       + h2544*f524 - f514 + f513 - h2534*f523)/Nc)/(a15*a52);
    
    return -ICPLX*E*a24*a24/(a12*a34);
  }

  std::complex<double> 
  ampq4g1::E32mpp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5), 
      s45 = S(4,5); 

    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), 
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a32 = A(3,2), a34 = A(3,4), 
      a35 = A(3,5), a52 = A(5,2), a53 = A(5,3), a54 = A(5,4);
    
    _ComplexD h1134 = H(1,1,3,4), h1144 = H(1,1,4,4), h1154 = H(1,1,5,4),
      h1234 = H(1,2,3,4), h1244 = H(1,2,4,4), h1254 = H(1,2,5,4), 
      h1534 = H(1,5,3,4), h1544 = H(1,5,4,4); 

    _ComplexD f123 = F(s45,s12,s23), f124 = F(s35,s12,s24), 
      f125 = F(s34,s12,s25), f134 = F(s25,s13,s34), f135 = F(s24,s13,s35), 
      f143 = F(s25,s14,s34), f145 = F(s23,s14,s45), f153 = F(s24,s15,s35), 
      f154 = F(s23,s15,s45), f213 = F(s45,s12,s13), f214 = F(s35,s12,s14), 
      f215 = F(s34,s12,s15), f234 = F(s15,s23,s34), f235 = F(s14,s23,s35), 
      f243 = F(s15,s24,s34), f245 = F(s13,s24,s45), f253 = F(s14,s25,s35), 
      f254 = F(s13,s25,s45), f351 = F(s24,s35,s15), f435 = F(s12,s34,s35), 
      f451 = F(s23,s45,s15), f513 = F(s24,s15,s13), f514 = F(s23,s15,s14),
      f521 = F(s34,s25,s12), f523 = F(s14,s25,s23), f524 = F(s13,s25,s24), 
      f534 = F(s12,s35,s34), f543 = F(s12,s45,s34); 

    _ComplexD E = -a32*(Nc*(f253 + h1244*f214 + h1154*f125
			    + h1544*f534 + h1134*f143)
			-(-h1144*(f134 + f124) + f154
			  -h1234*(f243 + f213) + f253
			  +h1244*(f234 + f214) - f254
			  +h1134*(f123 + f143) - f153)/Nc)/(a35*a52)
      
      + a34*(Nc*h1544*f534
	     -(-f254 + h1254*f215 + f451 - h1154*f125 - f435 
	       - h1144*(f134 - f124) + h1244*(f234 - f214)
	       + h1254*f245 - f235 + f135 - h1154*f145)/Nc)/(a35*a54)
      
      + a12*(Nc*h1154*f125
	     -(f451 - h1544*f534 - f351 + h1534*f543 - f521
	       - h1144*(f124 - f134) + h1134*(f123 - f143)
	       + f524 - h1544*f514 + h1534*f513 - f523)/Nc)/(a15*a52);
    
    return ICPLX*E*a14*a14/(a12*a34);
  }

  std::complex<double> 
  ampq4g1::E32pmp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5), 
      s45 = S(4,5); 

    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), 
      a21 = A(2,1), a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a32 = A(3,2),
      a34 = A(3,4), a35 = A(3,5), a43 = A(4,3), a45 = A(4,5), a51 = A(5,1), 
      a52 = A(5,2), a53 = A(5,3), a54 = A(5,4);
    
    _ComplexD h2133 = H(2,1,3,3), h2143 = H(2,1,4,3), h2153 = H(2,1,5,3),
      h2233 = H(2,2,3,3), h2243 = H(2,2,4,3), h2253 = H(2,2,5,3), 
      h2453 = H(2,4,5,3), h2513 = H(2,5,1,3), h2533 = H(2,5,3,3),
      h2543 = H(2,5,4,3);
    
    _ComplexD f123 = F(s45,s12,s23), f124 = F(s35,s12,s24), 
      f125 = F(s34,s12,s25), f134 = F(s25,s13,s34), f135 = F(s24,s13,s35), 
      f143 = F(s25,s14,s34), f145 = F(s23,s14,s45), f153 = F(s24,s15,s35), 
      f154 = F(s23,s15,s45), f213 = F(s45,s12,s13), f214 = F(s35,s12,s14), 
      f215 = F(s34,s12,s15), f234 = F(s15,s23,s34), f235 = F(s14,s23,s35), 
      f243 = F(s15,s24,s34), f245 = F(s13,s24,s45), f253 = F(s14,s25,s35), 
      f254 = F(s13,s25,s45), f351 = F(s24,s35,s15), f435 = F(s12,s34,s35), 
      f451 = F(s23,s45,s15), f513 = F(s24,s15,s13), f514 = F(s23,s15,s14),
      f521 = F(s34,s25,s12), f523 = F(s14,s25,s23), f524 = F(s13,s25,s24), 
      f534 = F(s12,s35,s34), f543 = F(s12,s45,s34);

    _ComplexD E = -a32*(Nc*(f253 + h2243*f214 + h2153*f125
			    + h2543*f534 + h2133*f143)
			-(-h2143*(f134 + f124) + f154
			  -h2233*(f243 + f213) + f253
			  +h2243*(f234 + f214) - f254
			  +h2133*(f123 + f143) - f153)/Nc)/(a35*a52)
      
      + a34*(Nc*h2543*f534
	     -(-f254 + h2253*f215 + f451 - h2153*f125 - h2453*f435 
	       - h2143*(f134 - f124) + h2243*(f234 - f214)
	       + f245 - h2253*f235 + h2153*f135 - f145)/Nc)/(a35*a54)
      
      + a12*(Nc*h2153*f125
	     -(f451 - h2543*f534 - f351 + h2533*f543 - h2513*f521
	       - h2143*(f124 - f134) + h2133*(f123 - f143)
	       + h2543*f524 - f514 + f513 - h2533*f523)/Nc)/(a15*a52);
    
    return ICPLX*E*a23*a23/(a12*a34);
  }

  std::complex<double> 
  ampq4g1::E34ppp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5), 
      s45 = S(4,5); 

    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), 
      a21 = A(2,1), a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a32 = A(3,2), 
      a34 = A(3,4), a35 = A(3,5), a51 = A(5,1), a52 = A(5,2), a53 = A(5,3), 
      a54 = A(5,4); 

    _ComplexD h2134 = H(2,1,3,4), h2144 = H(2,1,4,4), h2154 = H(2,1,5,4), 
      h2354 = H(2,3,5,4), h2514 = H(2,5,1,4), h2534 = H(2,5,3,4), 
      h2544 = H(2,5,4,4); 

    _ComplexD f124 = F(s35,s12,s24), f125 = F(s34,s12,s25), 
      f134 = F(s25,s13,s34), f135 = F(s24,s13,s35), f145 = F(s23,s14,s45), 
      f152 = F(s34,s15,s25), f153 = F(s24,s15,s35), f213 = F(s45,s12,s13), 
      f214 = F(s35,s12,s14), f215 = F(s34,s12,s15), f235 = F(s14,s23,s35), 
      f245 = F(s13,s24,s45), f251 = F(s34,s25,s15), f254 = F(s13,s25,s45), 
      f312 = F(s45,s13,s12), f315 = F(s24,s13,s15), f321 = F(s45,s23,s12),
      f342 = F(s15,s34,s24), f345 = F(s12,s34,s45), f351 = F(s24,s35,s15),
      f421 = F(s35,s24,s12), f435 = F(s12,s34,s35), f452 = F(s13,s45,s25),
      f512 = F(s34,s15,s12), f521 = F(s34,s25,s12), f524 = F(s13,s25,s24),
      f534 = F(s12,s35,s34), f543 = F(s12,s45,s34), H34, H32, H14;
    
    H34 = Nc*(f251+f315-f452-f351-f312-f512+f524-f124-h2514*f521)
      - (f213-f214+f421+f245+f135-f235-f435-h2134*f321-h2154*f145-h2354*f345)/Nc;
    
    H32 = Nc*(f254-f251+f312-f315-f342-h2544*f524+h2514*f521+h2534*f543);
    H14 = Nc*(f153-f152-f315+f215+h2144*f124-h2544*f524+h2544*f534
	      -h2144*f134+h2514*f521-h2154*f125);
        
    return -ICPLX*a24*a24*(a34*H34/(a35*a54) + a32*H32/(a35*a52) 
			   + a14*H14/(a15*a54))/(a12*a34);
  }
  
  std::complex<double> 
  ampq4g1::E34pmp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), s35 = S(3,5), 
      s45 = S(4,5);
 
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), 
      a21 = A(2,1), a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a32 = A(3,2), 
      a34 = A(3,4), a35 = A(3,5), a43 = A(4,3), a45 = A(4,5), a51 = A(5,1), 
      a52 = A(5,2), a53 = A(5,3), a54 = A(5,4);
    
    _ComplexD h2143 = H(2,1,4,3), h2153 = H(2,1,5,3), h2453 = H(2,4,5,3), 
      h2513 = H(2,5,1,3), h2533 = H(2,5,3,3), h2543 = H(2,5,4,3);

    _ComplexD f124 = F(s35,s12,s24), f125 = F(s34,s12,s25),
      f134 = F(s25,s13,s34), f135 = F(s24,s13,s35), f145 = F(s23,s14,s45),
      f152 = F(s34,s15,s25), f153 = F(s24,s15,s35), f213 = F(s45,s12,s13),
      f214 = F(s35,s12,s14), f215 = F(s34,s12,s15), f235 = F(s14,s23,s35),
      f245 = F(s13,s24,s45), f251 = F(s34,s25,s15), f254 = F(s13,s25,s45),
      f312 = F(s45,s13,s12), f315 = F(s24,s13,s15), f321 = F(s45,s23,s12),
      f342 = F(s15,s34,s24), f345 = F(s12,s34,s45), f351 = F(s24,s35,s15), 
      f435 = F(s12,s34,s35), f452 = F(s13,s45,s25), f512 = F(s34,s15,s12), 
      f521 = F(s34,s25,s12), f524 = F(s13,s25,s24), f534 = F(s12,s35,s34), 
      f543 = F(s12,s45,s34), H34, H32, H14;
    
    H34 = Nc*(f251+f315-f452-f351-f312-f512+h2543*f524-h2143*f124-h2513*f521)
      - (f213-f214-f321+f245-f345-f235-f145+h2143*f124+h2153*f135-h2453*f435)/Nc;
    
    H32 = Nc*(f254-f251+f312-f315-f342-h2543*f524+h2513*f521+h2533*f543);
    H14 = Nc*(f153-f152-f315+f215+h2143*f124-h2543*f524+h2543*f534
	      -h2143*f134+h2513*f521-h2153*f125);
    
    return ICPLX*a23*a23*(a34*H34/(a35*a54) + a32*H32/(a35*a52) 
			  + a14*H14/(a15*a54))/(a12*a34);
  }

  void ampq4g1::A34ppp(unsigned int nf, int p1, int p2, int p3,
		       int p4, int p5, _ComplexD *A) const
  {
    _ComplexD l12 = -Log(-S(1,2)), l34 = -Log(-S(3,4));
    _ComplexD V = 2.0*(Nc-nf)*l12/3.0 + 1.5*(l12+l34)/Nc
      + 29.0*Nc/18.0 + 6.5/Nc - 10.0*nf/9.0;
    
    A[0] = -ICPLX*A(2,4)*A(2,4)/(A(1,2)*A(3,5)*A(5,4));
    A[1] = V*A[0] + D34ppp(nf, p1,p2,p3,p4,p5) + E34ppp(p1,p2,p3,p4,p5);
  }

  void ampq4g1::A34mpp(unsigned int nf, int p1, int p2, int p3,
		       int p4, int p5, _ComplexD *A) const
  {
    _ComplexD l12 = -Log(-S(1,2)), l34 = -Log(-S(3,4));
    _ComplexD V = 2.0*(Nc-nf)*l12/3.0 + 1.5*(l12+l34)/Nc
      + 29.0*Nc/18.0 + 6.5/Nc - 10.0*nf/9.0;
    
    A[0] = ICPLX*A(1,4)*A(1,4)/(A(1,2)*A(3,5)*A(5,4));
    A[1] = V*A[0] + D34mpp(nf, p1,p2,p3,p4,p5) - E34pmp(p2,p1,p4,p3,p5);
  }

  void ampq4g1::A32ppp(unsigned int nf, int p1, int p2, int p3,
		       int p4, int p5, _ComplexD *A) const
  {
    _ComplexD l12 = -Log(-S(1,2)), l34 = -Log(-S(3,4));
    _ComplexD V = 2.0*(Nc-nf)*l34/3.0 + 1.5*(l12+l34)/Nc
      + 29.0*Nc/18.0 + 6.5/Nc - 10.0*nf/9.0;

    A[0] = -ICPLX*A(2,4)*A(2,4)/(A(1,2)*A(3,4))*A(3,2)/(A(3,5)*A(5,2));
    A[1] = V*A[0] + D32ppp(nf, p1,p2,p3,p4,p5) + E32ppp(p1,p2,p3,p4,p5);
  }

  void ampq4g1::A32mpp(unsigned int nf, int p1, int p2, int p3,
		       int p4, int p5, _ComplexD *A) const
  {
    _ComplexD l12 = -Log(-S(1,2)), l34 = -Log(-S(3,4));
    _ComplexD V = 2.0*(Nc-nf)*l34/3.0 + 1.5*(l12+l34)/Nc
      + 29.0*Nc/18.0 + 6.5/Nc - 10.0*nf/9.0;

    A[0] = ICPLX*A(1,4)*A(1,4)/(A(1,2)*A(3,4))*A(3,2)/(A(3,5)*A(5,2));
    A[1] = V*A[0] + D32mpp(nf, p1,p2,p3,p4,p5) + E32mpp(p1,p2,p3,p4,p5);
  }

  void ampq4g1::A32pmp(unsigned int nf, int p1, int p2, int p3,
		       int p4, int p5, _ComplexD *A) const
  {
    _ComplexD l12 = -Log(-S(1,2)), l34 = -Log(-S(3,4));
    _ComplexD V = 2.0*(Nc-nf)*l34/3.0 + 1.5*(l12+l34)/Nc
      + 29.0*Nc/18.0 + 6.5/Nc - 10.0*nf/9.0;

    A[0] = ICPLX*A(2,3)*A(2,3)/(A(1,2)*A(3,4))*A(3,2)/(A(3,5)*A(5,2));
    A[1] = V*A[0] + D32pmp(nf, p1,p2,p3,p4,p5) + E32pmp(p1,p2,p3,p4,p5);
  }

#undef F
#define F(i,j) (A(i,j)/(A(i,5)*A(5,j)))

  void ampq4g1::
  matrix_tree_pmpmp(int p1, int p2, int p3, int p4, int p5, amp_tree *M) const
  {
    _ComplexD H = -ICPLX*A(2,4)*A(2,4)/(A(1,2)*A(3,4));
    _ComplexD g = A(1,2)*A(3,4)/(A(1,4)*A(3,2));

    M[0].A12 = F(1,2)*H; M[0].A34 = F(3,4)*H;
    M[0].A32 = F(3,2)*H; M[0].A14 = F(1,4)*H;

    M[1].A12 = M[0].A12*(1.0+g*Nc); M[1].A34 = M[0].A34*(1.0+g*Nc);
    M[1].A32 = M[0].A32*(1.0+g/Nc); M[1].A14 = M[0].A14*(1.0+g/Nc);
  }

  void ampq4g1::
  matrix_tree_mpmpp(int p1, int p2, int p3, int p4, int p5, amp_tree *M) const
  {
    _ComplexD H = -ICPLX*A(1,3)*A(1,3)/(A(1,2)*A(3,4));
    _ComplexD g = A(1,2)*A(3,4)/(A(1,4)*A(3,2));

    M[0].A12 = F(1,2)*H; M[0].A34 = F(3,4)*H;
    M[0].A32 = F(3,2)*H; M[0].A14 = F(1,4)*H;

    M[1].A12 = M[0].A12*(1.0+g*Nc); M[1].A34 = M[0].A34*(1.0+g*Nc);
    M[1].A32 = M[0].A32*(1.0+g/Nc); M[1].A14 = M[0].A14*(1.0+g/Nc);
  }

  void ampq4g1::
  matrix_tree_pmmpp(int p1, int p2, int p3, int p4, int p5, amp_tree *M) const
  {
    _ComplexD H = ICPLX*A(2,3)*A(2,3)/(A(1,2)*A(3,4));

    M[1].A12 = M[0].A12 = F(1,2)*H; M[1].A34 = M[0].A34 = F(3,4)*H;
    M[1].A32 = M[0].A32 = F(3,2)*H; M[1].A14 = M[0].A14 = F(1,4)*H;
  }

  void ampq4g1::
  matrix_tree_mppmp(int p1, int p2, int p3, int p4, int p5, amp_tree *M) const
  {
    _ComplexD H = ICPLX*A(1,4)*A(1,4)/(A(1,2)*A(3,4));

    M[1].A12 = M[0].A12 = F(1,2)*H; M[1].A34 = M[0].A34 = F(3,4)*H;
    M[1].A32 = M[0].A32 = F(3,2)*H; M[1].A14 = M[0].A14 = F(1,4)*H;
  }

  void ampq4g1::
  matrix_tree_ppmmp(int p1, int p2, int p3, int p4, int p5, amp_tree *M) const
  {
    _ComplexD H = ICPLX*A(3,4)*A(3,4)/(A(1,4)*A(3,2));
    
    M[0].A12 = M[0].A34 = M[0].A32 = M[0].A14 = 0.0;
    M[1].A12 = F(1,2)*H*Nc; M[1].A34 = F(3,4)*H*Nc;
    M[1].A32 = F(3,2)*H/Nc; M[1].A14 = F(1,4)*H/Nc;
  }

  void ampq4g1::
  matrix_tree_mmppp(int p1, int p2, int p3, int p4, int p5, amp_tree *M) const
  {
    _ComplexD H = ICPLX*A(1,2)*A(1,2)/(A(1,4)*A(3,2));
    
    M[0].A12 = M[0].A34 = M[0].A32 = M[0].A14 = 0.0;
    M[1].A12 = F(1,2)*H*Nc; M[1].A34 = F(3,4)*H*Nc;
    M[1].A32 = F(3,2)*H/Nc; M[1].A14 = F(1,4)*H/Nc;
  }






  void ampq4g1::
  matrix_1loop_pmpmp(unsigned int nf, int p1, int p2, int p3,
		     int p4, int p5, amp_1loop *M) const
  {
    _ComplexD A[2];

    //   a12
    A12ppp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A12 = M[1].tree.A12 = A[0]; 
    M[0].loop.A12 = M[1].loop.A12 = A[1];
  
    A14ppp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A12 += A[0]*Nc; 
    M[1].loop.A12 += A[1]*Nc;
    
    //   a14
    A14ppp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A14 = M[1].tree.A14 = A[0]; 
    M[0].loop.A14 = M[1].loop.A14 = A[1];
  
    A12ppp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A14 += A[0]/Nc; 
    M[1].loop.A14 += A[1]/Nc;
    
    //   a32
    A32ppp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A32 = M[1].tree.A32 = A[0]; 
    M[0].loop.A32 = M[1].loop.A32 = A[1];
  
    A34ppp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A32 += A[0]/Nc; 
    M[1].loop.A32 += A[1]/Nc;
    
    //   a34
    A34ppp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A34 = M[1].tree.A34 = A[0]; 
    M[0].loop.A34 = M[1].loop.A34 = A[1];
  
    A32ppp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A34 += A[0]*Nc; 
    M[1].loop.A34 += A[1]*Nc;
  }

  void ampq4g1::
  matrix_1loop_mpmpp(unsigned int nf, int p1, int p2, int p3,
		     int p4, int p5, amp_1loop *M) const
  {
    _ComplexD A[2];

    //   a12
    A12mmp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A12 = M[1].tree.A12 = A[0]; 
    M[0].loop.A12 = M[1].loop.A12 = A[1];
  
    A14mmp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A12 += A[0]*Nc; 
    M[1].loop.A12 += A[1]*Nc;
    
    //   a14
    A14mmp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A14 = M[1].tree.A14 = A[0]; 
    M[0].loop.A14 = M[1].loop.A14 = A[1];
  
    A12mmp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A14 += A[0]/Nc; 
    M[1].loop.A14 += A[1]/Nc;
    
    //   a32
    A32mmp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A32 = M[1].tree.A32 = A[0]; 
    M[0].loop.A32 = M[1].loop.A32 = A[1];
  
    A34mmp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A32 += A[0]/Nc; 
    M[1].loop.A32 += A[1]/Nc;
    
    //   a34
    A34mmp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A34 = M[1].tree.A34 = A[0]; 
    M[0].loop.A34 = M[1].loop.A34 = A[1];
  
    A32mmp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A34 += A[0]*Nc; 
    M[1].loop.A34 += A[1]*Nc;
  }

  void ampq4g1::
  matrix_1loop_pmmpp(unsigned int nf, int p1, int p2, int p3,
		     int p4, int p5, amp_1loop *M) const
  {
    _ComplexD A[2];

    //   a12
    A12pmp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A12 = M[1].tree.A12 = A[0]; 
    M[0].loop.A12 = M[1].loop.A12 = A[1];
      
    //   a14
    A14pmp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A14 = M[1].tree.A14 = A[0]; 
    M[0].loop.A14 = M[1].loop.A14 = A[1];
      
    //   a32
    A32pmp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A32 = M[1].tree.A32 = A[0]; 
    M[0].loop.A32 = M[1].loop.A32 = A[1];
      
    //   a34
    A34pmp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A34 = M[1].tree.A34 = A[0]; 
    M[0].loop.A34 = M[1].loop.A34 = A[1];
  }

  void ampq4g1::
  matrix_1loop_mppmp(unsigned int nf, int p1, int p2, int p3,
		     int p4, int p5, amp_1loop *M) const
  {
    _ComplexD A[2];

    //   a12
    A12mpp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A12 = M[1].tree.A12 = A[0]; 
    M[0].loop.A12 = M[1].loop.A12 = A[1];
      
    //   a14
    A14mpp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A14 = M[1].tree.A14 = A[0]; 
    M[0].loop.A14 = M[1].loop.A14 = A[1];
      
    //   a32
    A32mpp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A32 = M[1].tree.A32 = A[0]; 
    M[0].loop.A32 = M[1].loop.A32 = A[1];
      
    //   a34
    A34mpp(nf, p1, p2, p3, p4, p5, A); 
    M[0].tree.A34 = M[1].tree.A34 = A[0]; 
    M[0].loop.A34 = M[1].loop.A34 = A[1];
  }

  void ampq4g1::
  matrix_1loop_ppmmp(unsigned int nf, int p1, int p2, int p3,
		     int p4, int p5, amp_1loop *M) const
  {
    _ComplexD A[2];

    //   a12
    M[0].tree.A12 = M[0].loop.A12 = 0.0;
    A14pmp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A12 = A[0]*Nc; 
    M[1].loop.A12 = A[1]*Nc;
    
    //   a14
    M[0].tree.A14 = M[0].loop.A14 = 0.0;
    A12pmp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A14 = A[0]/Nc; 
    M[1].loop.A14 = A[1]/Nc;
    
    //   a32
    M[0].tree.A32 = M[0].loop.A32 = 0.0;  
    A34pmp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A32 = A[0]/Nc; 
    M[1].loop.A32 = A[1]/Nc;
    
    //   a34
    M[0].tree.A34 = M[0].loop.A34 = 0.0;
    A32pmp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A34 = A[0]*Nc; 
    M[1].loop.A34 = A[1]*Nc;
  }

  void ampq4g1::
  matrix_1loop_mmppp(unsigned int nf, int p1, int p2, int p3,
		     int p4, int p5, amp_1loop *M) const
  {
    _ComplexD A[2];

    //   a12
    M[0].tree.A12 = M[0].loop.A12 = 0.0;
    A14mpp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A12 = A[0]*Nc; 
    M[1].loop.A12 = A[1]*Nc;
    
    //   a14
    M[0].tree.A14 = M[0].loop.A14 = 0.0;
    A12mpp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A14 = A[0]/Nc; 
    M[1].loop.A14 = A[1]/Nc;
    
    //   a32
    M[0].tree.A32 = M[0].loop.A32 = 0.0;  
    A34mpp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A32 = A[0]/Nc; 
    M[1].loop.A32 = A[1]/Nc;
    
    //   a34
    M[0].tree.A34 = M[0].loop.A34 = 0.0;
    A32mpp(nf, p1, p4, p3, p2, p5, A);
    M[1].tree.A34 = A[0]*Nc; 
    M[1].loop.A34 = A[1]*Nc;
  }

  double ampq4g1::amp1loop(const amp_1loop&m) const
  {
    return Nc2*real(m.loop.A14*conj(m.tree.A14))
      +    Nc2*real(m.loop.A32*conj(m.tree.A32))
      + real(m.loop.A12*conj(m.tree.A12)) 
      + real(m.loop.A34*conj(m.tree.A34))
      - real((m.loop.A14+m.loop.A32)*conj(m.tree.A12+m.tree.A34))
      - real((m.loop.A12+m.loop.A34)*conj(m.tree.A14+m.tree.A32));
  }
 
  void ampq4g1::
  su3_1loop(unsigned int nf, int p1, int p2, int p3, int p4, int p5, double *out) const
  {
    amp_1loop m[2];
  
    matrix_1loop_pmpmp(nf, p1,p2,p3,p4,p5, m);
    out[0] = amp1loop(m[0]); out[1] = amp1loop(m[1]);

    matrix_1loop_pmmpp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_mpmpp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_mppmp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_ppmmp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_mmppp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    swap();
    matrix_1loop_pmpmp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_pmmpp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_mpmpp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_mppmp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_ppmmp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);

    matrix_1loop_mmppp(nf, p1,p2,p3,p4,p5, m);
    out[0] += amp1loop(m[0]); out[1] += amp1loop(m[1]);
    swap();

    out[0] *= Na/Nc; out[1] *= Na/Nc;
  }

  void ampq4g1::
  su3_1loop_mch(unsigned int nf, int p1, int p2, int p3, int p4, int p5, double *out) const
  {
    amp_1loop m[2];
  
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (6*_M_rng());
    
    if(hpm == 1) swap();
    switch(hel){
    case 0: matrix_1loop_pmpmp(nf, p1,p2,p3,p4,p5, m); break;
    case 1: matrix_1loop_pmmpp(nf, p1,p2,p3,p4,p5, m); break;
    case 2: matrix_1loop_mpmpp(nf, p1,p2,p3,p4,p5, m); break;
    case 3: matrix_1loop_mppmp(nf, p1,p2,p3,p4,p5, m); break;
    case 4: matrix_1loop_ppmmp(nf, p1,p2,p3,p4,p5, m); break;
    case 5: matrix_1loop_mmppp(nf, p1,p2,p3,p4,p5, m); break;
    }
    if(hpm == 1) swap();

    out[0] = 12.0*Na*amp1loop(m[0])/Nc; 
    out[1] = 12.0*Na*amp1loop(m[1])/Nc;
  }

}  //  namespace nlo

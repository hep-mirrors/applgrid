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
#include "ampg4.h"
#include "defmacros.h"


namespace nlo {
  
#define ICPLX std::complex<double>(0.0,1.0)


  std::complex<double> 
  ampg4::Atree(int p1, int p2, int p3, int p4, int pi, int pj) const
  {
    _ComplexD aij = Aij(pi,pj);
    return ICPLX*aij*aij*aij*aij/(A(1,2)*A(2,3)*A(3,4)*A(4,1));
  }

  double ampg4::su3_tree(int p1, int p2, int p3, int p4) const
  {
    double w12 = S(1,2), w13 = S(1,3), w14 = S(1,4);
    w12 *= w12; w13 *= w13; w14 *= w14;
    
    return Nc2*(Nc2-1.0)*(w12*w12 + w13*w13 + w14*w14)
      *(4.0/(w12*w13)+4.0/(w12*w14)+4.0/(w13*w14));
  }

  double ampg4::ampcc12(int p1, int p2, int p3, int p4) const
  {
    static double cfac = 2.0*Nc2*Nc*(Nc2-1.0);
    double w12 = S(1,2), w13 = S(1,3), w14 = S(1,4);
    w12 *= w12; w13 *= w13; w14 *= w14;

    return -cfac*(w13+w14)*(w12*w12 + w13*w13 + w14*w14)/(w12*w13*w14);
    
  }
  
  double ampg4::su3_cc(int pi,int pj, int p1, int p2, int p3, int p4) const
  {
    int q[4] = {p1,p2,p3,p4};
    
    if(pi == pj) throw("Error in ampg4::amp_cc");
    
    for(int k = 0; k < 4; k++){
      if(q[k] == pi && k != 0) std::swap(q[0], q[k]);
      if(q[k] == pj && k != 1) std::swap(q[1], q[k]);
    }
    
    return ampcc12(q[0], q[1], q[2], q[3]);
  }

#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampg4::su3_kp(unsigned int nf, int pa, int pb, int p1, 
		     int p2, int p3, int p4, su3_kp_i2 *res, double al) const
  {
    double cc;
    int p[4] = {p1,p2,p3,p4};
    
    res->tree = res->pa = res->pb = res->ga = res->gb = res->loop = 0.0;
    for(unsigned i = 0; i < 4; i++)
      for(unsigned j = i+1; j < 4; j++) {
	cc = su3_cc(p[i], p[j], p1, p2, p3, p4);
	res->loop += 2.0*Xg(Sij(p[i],p[j]), 1.0, nf)*cc;
	
	if(p[i] == pa || p[j] == pa) {
	  int q = (p[i] == pa ? p[j] : p[i]);
	  
	  res->tree -= cc;
	  res->pa += X(pa, q)*cc;
	  
	  if(q == pb) 
	    res->cca = res->ccb = -cc/Ca;
	}
	
	if(p[i] == pb || p[j] == pb) 
	  res->pb += X(pb, p[i] == pb ? p[j] : p[i])*cc;
      }
    
    res->tree /= Ca; res->pa /= Ca; res->pb /= Ca;
    res->ga = res->gb = -Gg(nf)*(res->tree - res->cca);
    res->loop += (4.0*(Gg(nf) + Kg(nf,al)) - Ca/3.0)*(res->tree);
  }


#define NLO_PI2 9.86960440108935861881

  void ampg4::
  A1mmpp(unsigned nf, int p1, int p2, int p3, int p4, _ComplexD *A) const
  {
    _ComplexD l23 = Log(-S(2,3));
    _ComplexD Vg = std::pow(Log(-S(1,2)) - l23,2) + NLO_PI2;
    _ComplexD Vf = l23 - 2.0;
    _ComplexD Vs = -Vf/3.0 + 2.0/9.0;
    
    A[0] = Atree(p1,p2,p3,p4, p1,p2);
    A[1] = A[0]*(Nc*Vg + (4.0*Nc-nf)*Vf + (Nc-nf)*Vs)/Nc;
  }
  
  void ampg4::
  A1mpmp(unsigned nf, int p1, int p2, int p3, int p4, _ComplexD *A) const
  {
    double s12 = S(1,2), s13 = S(1,3), s23 = S(2,3); 

    _ComplexD l12 = Log(-s12), l23 = Log(-s23);
    _ComplexD Vg = std::pow(l12-l23,2) + NLO_PI2;
    _ComplexD Vf = 0.5*(l12+l23) - 2.0;
    _ComplexD Vs = -Vf/3.0 + 2.0/9.0;
    
    double tmp1 = s13*s13/(s12*s23);
    _ComplexD Ff, Fs;
    
    Ff = -0.5*Vg/tmp1 + 0.5*(s12-s23)/s13*(l12-l23);
    Fs = -1.0/tmp1*(-Vg/tmp1 + (s12-s23)/s13*(1.0+tmp1/6.0)*(l12-l23) + 1.0);
    
    A[0] = Atree(p1,p2,p3,p4, p1,p3);
    A[1] = A[0]*(Nc*Vg + (4.0*Nc-nf)*(Vf+Ff) + (Nc-nf)*(Vs+Fs))/Nc; 
  }
  
#define __hel_switch(i, j)						\
  ((p##i == pi) && (p##j == pj)) || ((p##j == pi) && (p##i == pj))
  
  void ampg4::matrix_1loop(unsigned nf, int pi, int pj, int p1,
			   int p2, int p3, int p4, _ComplexD *d) const
  {
    if(__hel_switch(1,2))      A1mmpp(nf, p1,p2,p3,p4, d);
    else if(__hel_switch(1,4)) A1mmpp(nf, p4,p1,p2,p3, d);
    else if(__hel_switch(2,3)) A1mmpp(nf, p2,p3,p4,p1, d);
    else if(__hel_switch(3,4)) A1mmpp(nf, p3,p4,p1,p2, d);
    
    else if(__hel_switch(1,3)) A1mpmp(nf, p1,p2,p3,p4, d);
    else if(__hel_switch(2,4)) A1mpmp(nf, p2,p3,p4,p1, d);
    else throw "ampg4::matrix_1loop(...) : bad helicity configuration";
  }

  double ampg4::
  su3_1loop(unsigned nf, int p1, int p2, int p3, int p4) const
  {
    int pi, pj, p[4] = {p1, p2, p3, p4};
    double leading = 0.0;
    
    static _ComplexD d[2];
    static const struct {unsigned int pi, pj;}
    helconf[6]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    
    for(unsigned int hel = 0; hel < 6; hel++) {
      pi = p[helconf[hel].pi];
      pj = p[helconf[hel].pj];
      
      matrix_1loop(nf, pi, pj, p1, p2, p3, p4, d);
      leading += real(d[1]*conj(d[0]));
      
      matrix_1loop(nf, pi, pj, p1, p2, p4, p3, d);
      leading += real(d[1]*conj(d[0]));
      
      matrix_1loop(nf, pi, pj, p1, p3, p2, p4, d);
      leading += real(d[1]*conj(d[0]));
    }
    
    return 2.0*Nc*Nc2*Na*leading;
  }

  double ampg4::
  su3_1loop_mch(unsigned nf, int p1, int p2, int p3, int p4) const
  {
    int pi, pj, p[4] = {p1, p2, p3, p4};
    double leading;
    unsigned int hel = (unsigned int) (6*_M_rng());

    static _ComplexD d[2];
    static const struct {unsigned int pi, pj;}
    helconf[6]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    
    pi = p[helconf[hel].pi];
    pj = p[helconf[hel].pj];
    
    matrix_1loop(nf, pi, pj, p1, p2, p3, p4, d);
    leading  = real(d[1]*conj(d[0]));
    
    matrix_1loop(nf, pi, pj, p1, p2, p4, p3, d);
    leading += real(d[1]*conj(d[0]));
    
    matrix_1loop(nf, pi, pj, p1, p3, p2, p4, d);
    leading += real(d[1]*conj(d[0]));
    
    return 12.0*Nc*Nc2*Na*leading;
  }
}

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
#include "ampg5.h"
#include "defmacros.h"


namespace nlo {
  
#define ICPLX std::complex<double>(0.0,1.0)

  const int ampg5::perm1[12][4] = 
  {{0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1},
   {0,3,1,2}, {0,3,2,1}, {1,0,2,3}, {1,0,3,2}, 
   {1,2,0,3}, {1,3,0,2}, {2,0,1,3}, {2,1,0,3}}; 
  
  const int ampg5::perm2[6][3] = 
  {{0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}};
  
  
  std::complex<double> 
  ampg5::ap3m2(int p1, int p2, int p3, int p4, int p5) const {
    return ICPLX/(A(1,2)*A(2,3)*A(3,4)*A(4,5)*A(5,1));
  }
  
  double ampg5::su3_tree(int p1, int q2, int q3, int q4, int q5) const
  {
    int p2, p3, p4, p5, p[4] = {q2,q3,q4,q5};
    double sum = 0.0;
  
    for(int i = 0; i < 12; i++) {
      p2 = p[perm1[i][0]];
      p3 = p[perm1[i][1]];
      p4 = p[perm1[i][2]];
      p5 = p[perm1[i][3]];
      sum += 1.0/(S(1,2)*S(2,3)*S(3,4)*S(4,5)*S(5,1));
    }
    
    sum *= std::pow(S(1,2),4) + std::pow(S(1,3),4)
      +    std::pow(S(1,4),4) + std::pow(S(1,5),4)
      +    std::pow(S(2,3),4) + std::pow(S(2,4),4)
      +    std::pow(S(2,5),4) + std::pow(S(3,4),4)
      +    std::pow(S(3,5),4) + std::pow(S(4,5),4);
    
    return 4.0*Nc2*Nc*(Nc2-1.0)*sum;
  }
  
  std::pair<double, std::complex<double> >
  ampg5::ampcc12(int p1, int p2, int p3, int p4, int p5) const
  {
    int q3, q4, q5, p[3] = {p3,p4,p5};
    std::complex<double> A1, A2, sum2 = 0.0;
    static const double cfac = 2.0*Nc2*Nc2*(Nc2-1.0);
    double sum1 = 0.0;
  
    for(int i = 0; i < 6; i++) {
      q3 = p[perm2[i][0]];
      q4 = p[perm2[i][1]];
      q5 = p[perm2[i][2]];
      
      A1 = ap3m2(p1,p2,q3,q4,q5);
      A2 = ap3m2(p1,p2,q4,q5,q3)+ap3m2(p1,p2,q5,q3,q4)+ap3m2(p1,p2,q5,q4,q3);
    
      sum1 -= real(conj(A1)*(A1 + 12.0*A2/Nc2));
      sum2 += A1*(A1 + 12.0*A2/Nc2);    
    }
  
    sum1 *= pow(S(1,2),4)+pow(S(1,3),4)+pow(S(1,4),4)+pow(S(1,5),4)
      + pow(S(2,3),4)+pow(S(2,4),4)+pow(S(2,5),4)
      + pow(S(3,4),4)+pow(S(3,5),4)
      + pow(S(4,5),4);
    
    sum2 *= pow(A(1,3)*A(4,5),4)+pow(A(1,4)*A(3,5),4)+pow(A(1,5)*A(3,4),4);
    
    return std::pair<double, std::complex<double> >(cfac*sum1, cfac*sum2);
  }

  std::pair<double, std::complex<double> >
  ampg5::su3_cc(int i,int j, int p1, int p2, int p3, int p4, int p5) const
  {
    int q[5] = {p1,p2,p3,p4,p5}, temp;
    
    if(i == j) throw("Error in ampg5::amp_cc");
    
    for(int k = 0; k < 5; k++){
      if(q[k] == i && k != 0){ temp=q[0]; q[0]=q[k]; q[k]=temp;}
      if(q[k] == j && k != 1){ temp=q[1]; q[1]=q[k]; q[k]=temp;}
    }
    
    return ampcc12(q[0], q[1], q[2], q[3], q[4]);
  }

#define LL(x,y,z,w) (x-y)*(z-w)
#define F0(sij, skl) (Log(-(sij))-Log(-(skl)))
#define F1(sij, skl) (L0(-(sij), -(skl))/(skl))
#define F2(sij, skl) (L1(-(sij), -(skl))/((skl)*(skl)))
#define F3(sij, skl) (L2(-(sij), -(skl))/((skl)*(skl)*(skl)))
#define F(smn, sij, sjk) Ls_1(-(sij), -(smn), -(sjk), -(smn))


  std::complex<double> 
  ampg5::Vg(int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD  l12 = Log(-S(1,2)), l23 = Log(-S(2,3)),
      l34 = Log(-S(3,4)), l45 = Log(-S(4,5)), l51 = Log(-S(5,1));
    
    return LL(l12,l23,l34,l45) + LL(l23,l34,l45,l51) + LL(l34,l45,l51,l12)
      + LL(l45,l51,l12,l23) + LL(l51,l12,l23,l34) + 8.224670334241132; 
  }  

  std::complex<double> 
  ampg5::Vf(bool hel, int p1, int p2, int p3, int p4, int p5) const
  {              
    double s23 = (hel ? S(2,3) : S(3,4));
    return 0.5*(Log(-s23) + Log(-S(5,1))) - 2.0;
  }
  
  std::complex<double>
  ampg5::Ff(bool hel, int p1, int p2, int p3, int p4, int p5) const
  {              
    if(hel) 
      return 0.5*(R(1,3,4,2)*S(3,4) + R(1,4,5,2)*S(4,5))*F1(S(2,3), S(5,1));  
    else {
      double s34 = S(3,4), s51 = S(5,1), s24 = S(2,4), s23 = S(2,3), 
	s12 = S(1,2), s25 = S(2,5);
      _ComplexD r1253 = R(1,2,5,3), r1423 = R(1,4,2,3), f13451 = F1(s34, s51);
      
      return R(1,2,4,3)*r1423*(F(s51, s23, s34) + s24*(F1(s23, s51) + f13451))
	+    R(1,5,2,3)*r1253*(F(s34, s12, s51) + s25*(F1(s12, s34) + f13451))
	-    0.5*(r1253*s25 + r1423*s24)*f13451;
    }
  }

  std::complex<double>
  ampg5::Fs(bool hel, int p1, int p2, int p3, int p4, int p5) const
  {    
    if(hel) { 
      _ComplexD r1342 = R(1,3,4,2), r1452 = R(1,4,5,2), r1352 = R(1,3,5,2);
      double s34 = S(3,4), s45 = S(4,5), s23 = S(2,3), s51 = S(5,1), 
	s12 = S(1,2), s35 = S(3,5);
    
      return ((r1342*s34 + r1452*s45)*
	      (r1342*r1452*s34*s45*F3(s23, s51)-0.5*F1(s23, s51))
	      -r1352/conj(r1352)*(s35*s35/(s12*s12)-s35/s12)
	      +0.5*r1342*r1452*s34*s45/(s23*s51))/3.0;
    } else {

      _ComplexD r1243 = R(1,2,4,3), r1423 = R(1,4,2,3), r1253 = R(1,2,5,3),
	r1523 = R(1,5,2,3);
      double s51 = S(5,1), s23 = S(2,3), s34 = S(3,4), s24 = S(2,4),
	s12 = S(1,2), s25 = S(2,5);


      return -r1243*r1423*r1243*r1423*
	(2.0*F(s51, s23, s34)
	 + s24*s24*(F2(s23, s51) + F2(s34, s51))
	 + 2.0*s24*(F1(s23, s51) + F1(s34, s51)))
	   
	-r1253*r1523*r1253*r1523*
	(2.0*F(s34, s12, s51) 
	 + s25*s25*(F2(s12, s34) + F2(s51, s34))
	 + 2.0*s25*(F1(s12, s34) + F1(s51, s34)))
	      
	+ 2.0*(r1423*r1243*r1243*r1243*s24*s24*s24*F3(s23, s51)
	       +r1253*r1523*s25*r1523*s25*r1523*s25*F3(s12, s34))/3.0
	      
	+ (2.0*(r1243*r1423*r1423*r1423*s24*s24*s24
		+ r1523*r1253*s25*r1253*s25*r1253*s25)
	   - r1423*r1253*s24*s25*(r1423*s24 + r1253*s25))*F3(s34, s51)/3.0
	      
	+ (r1423*s24 + r1253*s25)*F1(s34, s51)/6.0
	+ (r1423*r1253*s24*s25*r1423*r1253*s24*s25/(s12*s23*s34*s51)
	   - 0.5*r1423*r1253*s24*s25/(s34*s51)
	   + r1243*r1423*s24*r1243*r1423*s24*s24/(s23*s34*s51)
	   + r1253*r1523*s25*r1253*r1523*s25*s25/(s12*s34*s51))/3.0;
    }
  } 
  
#define __hel_switch(i, j)						\
  ((p##i == pi) && (p##j == pj)) || ((p##j == pi) && (p##i == pj))
  
  void ampg5::matrix_1loop(unsigned nf, int pi, int pj, int p1, int p2,
			   int p3, int p4, int p5, _ComplexD *d) const
  {
    if(__hel_switch(1,2))      A1mmppp(nf, p1,p2,p3,p4,p5, d);
    else if(__hel_switch(5,1)) A1mmppp(nf, p5,p1,p2,p3,p4, d);
    else if(__hel_switch(4,5)) A1mmppp(nf, p4,p5,p1,p2,p3, d);
    else if(__hel_switch(3,4)) A1mmppp(nf, p3,p4,p5,p1,p2, d);
    else if(__hel_switch(2,3)) A1mmppp(nf, p2,p3,p4,p5,p1, d);
    				      			 
    else if(__hel_switch(1,3)) A1mpmpp(nf, p1,p2,p3,p4,p5, d);
    else if(__hel_switch(5,2)) A1mpmpp(nf, p5,p1,p2,p3,p4, d);
    else if(__hel_switch(4,1)) A1mpmpp(nf, p4,p5,p1,p2,p3, d);
    else if(__hel_switch(3,5)) A1mpmpp(nf, p3,p4,p5,p1,p2, d);
    else if(__hel_switch(2,4)) A1mpmpp(nf, p2,p3,p4,p5,p1, d);
    else throw "ampg5::matrix_1loop(...) : bad helicity configuration";
  }

  void ampg5::
  A1mmppp(unsigned nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *d) const
  {
    _ComplexD vg = Vg(p1, p2, p3, p4, p5);
    _ComplexD vf = Vf(true, p1, p2, p3, p4, p5);
    _ComplexD ff = Ff(true, p1, p2, p3, p4, p5);
    _ComplexD fs = Fs(true, p1, p2, p3, p4, p5);
    
    d[0] = Atree(p1,p2,p3,p4,p5, p1, p2);
    d[1] = Nc*vg + B0(nf)*vf + 2.0*(Nc-nf)/9.0 + (4.0*Nc-nf)*ff + (Nc-nf)*fs;
    d[2] = d[1] + (nf*10.0)*vf/18.0 + 5.0*nf/27.0 + (5.0*nf)*(ff + fs)/6.0;
    d[1] *= d[0]; d[2] *= d[0];
  }
  
  void ampg5::
  A1mpmpp(unsigned nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *d) const
  {
    _ComplexD vg = Vg(p1, p2, p3, p4, p5);
    _ComplexD vf = Vf(false, p1, p2, p3, p4, p5);
    _ComplexD ff = Ff(false, p1, p2, p3, p4, p5);
    _ComplexD fs = Fs(false, p1, p2, p3, p4, p5);
    
    d[0] = Atree(p1,p2,p3,p4,p5, p1, p3);
    d[1] = Nc*vg + B0(nf)*vf + 2.0*(Nc-nf)/9.0 + (4.0*Nc-nf)*ff + (Nc-nf)*fs;
    d[2] = d[1] + (nf*10.0)*vf/18.0 + 5.0*nf/27.0 + (5.0*nf)*(ff + fs)/6.0;
    d[1] *= d[0]; d[2] *= d[0];
  }

  std::complex<double> 
  ampg5::Atree(int p1, int p2, int p3, int p4, int p5, int pi, int pj) const
  {
    _ComplexD aij = Aij(pi,pj);
    return ICPLX*aij*aij*aij*aij/(A(1,2)*A(2,3)*A(3,4)*A(4,5)*A(5,1));
  }

  double ampg5::
  su3_1loop_mmppp(unsigned int nf, int p1, int q2, int q3, int q4, int q5) const
  {
    int pi, pj, p2, p3, p4, p5;
    int p[4] = {q2, q3, q4, q5};
    int q[5] = {p1, q2, q3, q4, q5};
    
    double leading = 0.0, subleading = 0.0;

    static _ComplexD d[12][3];
    static const double cfac = 2.0*Nc*(Nc2-1.0);
    static const struct {unsigned int pi, pj;}
    helconf[10]={{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
    
    for(unsigned int hel = 0; hel < 10; hel++) {
      pi = q[helconf[hel].pi];
      pj = q[helconf[hel].pj];
      
      for(unsigned int i = 0; i < 12; i++) {
	p2 = p[perm1[i][0]];
	p3 = p[perm1[i][1]];
	p4 = p[perm1[i][2]];
	p5 = p[perm1[i][3]];
	
	matrix_1loop(nf, pi, pj, p1, p2, p3, p4, p5, d[i]);
	leading += real(d[i][1]*conj(d[i][0]));
      }
      
      subleading += real(d[0][2]*conj(d[ 9][0]) - d[ 9][2]*conj(d[0][0]))
	+           real(d[1][2]*conj(d[ 8][0]) - d[ 8][2]*conj(d[1][0]))
	-           real(d[2][2]*conj(d[ 7][0]) - d[ 7][2]*conj(d[2][0]))
	+           real(d[3][2]*conj(d[11][0]) - d[11][2]*conj(d[3][0]))
	-           real(d[4][2]*conj(d[ 6][0]) - d[ 6][2]*conj(d[4][0]))
	-           real(d[5][2]*conj(d[10][0]) - d[10][2]*conj(d[5][0])); 
    }
    
    return cfac*(Nc2*leading + 12.0*subleading);
  }
  
  double ampg5::
  su3_1loop_mch(unsigned int nf, int p1, int q2, int q3, int q4, int q5) const
  {
    int pi, pj, p2, p3, p4, p5;
    int p[4] = {q2, q3, q4, q5};
    int q[5] = {p1, q2, q3, q4, q5};
    
    double leading = 0.0, subleading = 0.0;
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (10*_M_rng());
    
    static _ComplexD d[12][3];
    static const double cfac = 40.0*Nc*(Nc2-1.0);
    static const struct {unsigned int pi, pj;}
    helconf[10]={{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
    
    pi = q[helconf[hel].pi];
    pj = q[helconf[hel].pj];
    
    if(hpm == 1) swap();
    for(unsigned int i = 0; i < 12; i++) {
      p2 = p[perm1[i][0]];
      p3 = p[perm1[i][1]];
      p4 = p[perm1[i][2]];
      p5 = p[perm1[i][3]];
      
      matrix_1loop(nf, pi, pj, p1, p2, p3, p4, p5, d[i]);
      leading += real(d[i][1]*conj(d[i][0]));
    }
    if(hpm == 1) swap();

    subleading += real(d[0][2]*conj(d[ 9][0]) - d[ 9][2]*conj(d[0][0]))
      +           real(d[1][2]*conj(d[ 8][0]) - d[ 8][2]*conj(d[1][0]))
      -           real(d[2][2]*conj(d[ 7][0]) - d[ 7][2]*conj(d[2][0]))
      +           real(d[3][2]*conj(d[11][0]) - d[11][2]*conj(d[3][0]))
      -           real(d[4][2]*conj(d[ 6][0]) - d[ 6][2]*conj(d[4][0]))
      -           real(d[5][2]*conj(d[10][0]) - d[10][2]*conj(d[5][0])); 
    
    return cfac*(Nc2*leading + 12.0*subleading);
  }

#define X(i,j) (-std::log(std::abs(Sij(i,j))))


  void ampg5::su3_kp(unsigned int nf, int pa, int pb, int p1, int p2,
		     int p3, int p4, int p5, su3_kp_i2 *res, double al) const
  {
    double cc;
    int p[5] = {p1,p2,p3,p4,p5};
    
    res->tree = res->pa = res->pb = res->ga = res->gb = res->loop = 0.0;
    for(unsigned i = 0; i < 5; i++) 
      for(unsigned j = i+1; j < 5; j++) {
	cc = ampcc(p[i], p[j], p1, p2, p3, p4, p5);
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
    res->loop += (5.0*(Gg(nf) + Kg(nf,al)) - Ca/3.0)*(res->tree);
  }
  
  double ampg5::cc12(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), 
      s23 = S(2,3), s24 = S(2,4), s25 = S(2,5), s34 = S(3,4), 
      s35 = S(3,5), s45 = S(4,5); 
    
    //----- helicity sum -----
    double pp = pow(s12,4)+pow(s13,4)+pow(s14,4)+pow(s15,4)+pow(s23,4)
      +  pow(s24,4)+pow(s25,4)+pow(s34,4)+pow(s35,4)+pow(s45,4);

    //----- leading color -----
    double tmp = Nc2/s12;
    double ccl = tmp/(s23*s34*s45*s15) + tmp/(s24*s45*s35*s13) 
      +          tmp/(s25*s35*s34*s14) + tmp/(s23*s35*s45*s14) 
      +          tmp/(s25*s45*s34*s13) + tmp/(s24*s35*s34*s15);
    
    //----- subleading color -----
    double ccsl = 1.0/(s13*s14*s15*s34*s35) + 1.0/(s14*s15*s23*s24*s35)
      +           1.0/(s14*s15*s23*s25*s34) - 1.0/(s13*s14*s15*s24*s25)
      -           1.0/(s13*s24*s25*s34*s35) + 1.0/(s13*s15*s23*s24*s45)
      +           1.0/(s13*s14*s23*s25*s45) - 2.0*s12/(s13*s14*s15*s23*s24*s25)
      -           s34/(s13*s14*s23*s24*s35*s45) - s35/(s13*s15*s23*s25*s34*s45);
    
    return 2.0*Nc2*Na*pp*(12.0*ccsl-ccl);
  }


  double 
  ampg5::ampcc(int pi, int pj, int p1, int p2, int p3, int p4, int p5) const
  {
    int p[5] = {p1,p2,p3,p4,p5};
    
    if(pi == pj) throw("Error in ampg5::ampcc");
    
    for(int k = 0; k < 5; k++){
      if(p[k] == pi && k != 0) std::swap(p[0], p[k]);
      if(p[k] == pj && k != 1) std::swap(p[1], p[k]);
    }
    
    return cc12(p[0], p[1], p[2], p[3], p[4]);
  }



}   //   namespace nlo






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
#include "bits/nlo-color.h"
#include "ampq4p1.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {


  std::complex<double> ampq4p1::
  Amhv(double Q1, double Q2, const _ComplexD& fhh, int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD a21 = A(2,1), a43 = A(4,3);
    return ICPLX*fhh/(a21*a43)
      *(Q1*a21/(A(2,5)*A(5,1)) + Q2*a43/(A(4,5)*A(5,3)));
  }
  
  void ampq4p1::su3_tree(double Q, int p1, int p2, int p3, 
			 int p4, int p5, double *out) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s23 = S(2,3), 
      s24 = S(2,4), s34 = S(3,4);
    
    _ComplexD F1 = Q/(A(1,5)*A(5,2)*A(3,4)) + Q/(A(3,5)*A(5,4)*A(1,2));
    _ComplexD F2 = Q/(A(1,5)*A(5,4)*A(3,2)) + Q/(A(3,5)*A(5,2)*A(1,4));
    
    double S1 = 2.0*(s13*s13 + s14*s14 + s23*s23 + s24*s24)*real(F1*conj(F1));
    double S2 = 2.0*(s13*s13 + s12*s12 + s34*s34 + s24*s24)*real(F2*conj(F2));
    double S3 = 4.0*(s13*s13 + s24*s24)*real(F1*conj(F2));
    
    out[0] = 2.0*Na*S1;
    out[1] = 2.0*Na*(S1 + S2 + S3/Nc);
  }

  double ampq4p1::su3_tree(double Q1, double Q2, int p1, 
			   int p2, int p3, int p4, int p5) const
  {
    double s13 = S(1,3), s14 = S(1,4), s23 = S(2,3), s24 = S(2,4);
    _ComplexD F1 = Q1/(A(1,5)*A(5,2)*A(3,4)) + Q2/(A(3,5)*A(5,4)*A(1,2));
   
    return 4.0*Na*(s13*s13 + s14*s14 + s23*s23 + s24*s24)*real(F1*conj(F1));
  }

  void ampq4p1::
  ampcc(double c1, double c2, double c3, double Q,
	int p1, int p2, int p3, int p4, int p5, double *out) const
  {
    double s12 = S(1,2), s13 = S(1,3), s14 = S(1,4), s23 = S(2,3), 
      s24 = S(2,4), s34 = S(3,4);
    
    _ComplexD F1 = Q/(A(1,5)*A(5,2)*A(3,4)) + Q/(A(3,5)*A(5,4)*A(1,2));
    _ComplexD F2 = Q/(A(1,5)*A(5,4)*A(3,2)) + Q/(A(3,5)*A(5,2)*A(1,4));
    
    double S1 = 2.0*(s13*s13 + s14*s14 + s23*s23 + s24*s24)*real(F1*conj(F1));
    double S2 = 2.0*(s13*s13 + s12*s12 + s34*s34 + s24*s24)*real(F2*conj(F2));
    double S3 = 4.0*(s13*s13 + s24*s24)*real(F1*conj(F2));
    
    out[0] = -Na*c1*S1/Nc;
    out[1] = -Na*(c1*S1 + c2*S2 + c3*S3)/Nc;
  }
  
  double ampq4p1::
  ampcc(double c, double Q1, double Q2, int p1, int p2, int p3, int p4, int p5) const
  {
    double s13 = S(1,3), s14 = S(1,4), s23 = S(2,3), s24 = S(2,4);
    _ComplexD F1 = Q1/(A(1,5)*A(5,2)*A(3,4)) + Q2/(A(3,5)*A(5,4)*A(1,2));
    
    return -2.0*Na*c*(s13*s13+s14*s14+s23*s23+s24*s24)*real(F1*conj(F1))/Nc;
  }
  
#define __cc_switch12                                   \
  (pi == p1 && pj == p2) || (pi == p2 && pj == p1) ||   \
  (pi == p3 && pj == p4) || (pi == p4 && pj == p3)
  
#define __cc_switch13                                   \
  (pi == p1 && pj == p3) || (pi == p2 && pj == p4) ||   \
  (pi == p3 && pj == p1) || (pi == p4 && pj == p2)

#define __cc_switch14                                   \
  (pi == p1 && pj == p4) || (pi == p2 && pj == p3) ||   \
  (pi == p4 && pj == p1) || (pi == p3 && pj == p2)

  void ampq4p1::
  su3_cc(double Q, int pi,int pj, int p1, 
	 int p2, int p3, int p4, int p5, double *cc) const
  {
    double onc = 1.0/Nc, nc2 = Nc2-2.0; 
    if(__cc_switch12)      ampcc(-1.0, nc2, -onc,   Q, p1,p2, p3,p4, p5, cc);
    else if(__cc_switch13) ampcc( 2.0, 2.0, Nc+onc, Q, p1,p2, p3,p4, p5, cc);
    else if(__cc_switch14) ampcc(nc2, -1.0, -onc,   Q, p1,p2, p3,p4, p5, cc);
    else throw("Error in ampq4p1::su3_cc");     
  }
  
  double ampq4p1::
  su3_cc(double Q1, double Q2, int pi,int pj, 
	 int p1, int p2, int p3, int p4, int p5) const
  {
    if(__cc_switch12)      return ampcc(-1.0,    Q1,Q2, p1,p2, p3,p4, p5);
    else if(__cc_switch13) return ampcc( 2.0,    Q1,Q2, p1,p2, p3,p4, p5);
    else if(__cc_switch14) return ampcc(Nc2-2.0, Q1,Q2, p1,p2, p3,p4, p5);
    else throw("Error in ampq4p1::su3_cc");     
  }
  

#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampq4p1::su3_kp(double Q, int pa, int p1, int p2, int p3, 
		     int p4, int p5, su3_kp_i1 *res, double al) const
  {
    double cc[2], xi;
    int p[4] = {p1,p2,p3,p4};
    
    res[0].tree = res[0].pa = res[0].loop = 0.0;
    res[1].tree = res[1].pa = res[1].loop = 0.0;
    for(unsigned i = 0; i < 4; i++)
      for(unsigned j = i+1; j < 4; j++) {
	su3_cc(Q, p[i], p[j], p1, p2, p3, p4, p5, cc);
	
	//----- loop contribution -----
	xi = 2.0*Xq(Sij(p[i],p[j]), 1.0);
	res[0].loop += xi*cc[0];
	res[1].loop += xi*cc[1];
	
	//----- finite contributions -----
	if(p[i] == pa || p[j] == pa) {
	  int q = p[i] == pa ? p[j] : p[i];
	  res[0].tree -= cc[0];
	  res[1].tree -= cc[1];

	  xi = X(pa, q);
	  res[0].pa += xi*cc[0];
	  res[1].pa += xi*cc[1];
	}
      }
    
    res[0].tree /= Cf; res[0].pa /= Cf;
    res[1].tree /= Cf; res[1].pa /= Cf;
    
    res[0].ga = -Gq*res[0].tree;
    res[1].ga = -Gq*res[1].tree;
    
    xi = 4.0*(Gq+Kq(al)) - 2.0*Cf + Ca/3.0; 
    res[0].loop += xi*res[0].tree;
    res[1].loop += xi*res[1].tree;
  }
  
  void ampq4p1::su3_kp(double Q1, double Q2, int pa, int p1, int p2, int p3, 
		       int p4, int p5, su3_kp_i1 *res, double al) const
  {
    double cc, xi;
    int p[4] = {p1,p2,p3,p4};
    
    res->tree = res->pa = res->loop = 0.0;
    for(unsigned i = 0; i < 4; i++)
      for(unsigned j = i+1; j < 4; j++) {
	cc = su3_cc(Q1, Q2, p[i], p[j], p1, p2, p3, p4, p5);
	
	//----- loop contribution -----
	xi = 2.0*Xq(Sij(p[i],p[j]), 1.0);
	res->loop += xi*cc;
		
	//----- finite contributions -----
	if(p[i] == pa || p[j] == pa) {
	  int q = p[i] == pa ? p[j] : p[i];
	  xi = X(pa, q);
	  res->tree -= cc;
	  res->pa += xi*cc;
	}
      }
    
    res->tree /= Cf; res->pa /= Cf;
    res->ga = -Gq*res->tree;
    res->loop += (4.0*(Gq+Kq(al)) - 2.0*Cf + Ca/3.0)*(res->tree);
  }
  
  //   1-loop amplitude
#define F0(sij, skl) (Log(-(sij))-Log(-(skl)))
#define F1(sij, skl) (L0(-(sij), -(skl))/(skl))
#define F2(sij, skl) (L1(-(sij), -(skl))/((skl)*(skl)))
#define F3(sij, skl) (L2(-(sij), -(skl))/((skl)*(skl)*(skl)))
#define F(smn, sij, sjk) Ls_1(-(sij), -(smn), -(sjk), -(smn))
#define H(i,j,k,l) std::pow(a##i##k*a##j##l/(a##i##l*a##j##k), 2)

  std::complex<double> 
  ampq4p1::u0ppp(int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD a13 = A(1,3), a15 = A(1,5), a25 = A(2,5), a34 = A(3,4);    
    return -ICPLX*a13*a13/(a34*a15*a25);
  }    

  std::complex<double> 
  ampq4p1::u0pmp(int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD a14 = A(1,4), a15 = A(1,5), a25 = A(2,5), a34 = A(3,4);
    return ICPLX*a14*a14/(a34*a15*a25);
  }
  
  std::complex<double> 
  ampq4p1::u1lppp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s13 = S(1,3), s15 = S(1,5), s24 = S(2,4), s25 = S(2,5),
      s34 = S(3,4);
    
    _ComplexD a12 = A(1,2), a13 = A(1,3), a15 = A(1,5), a23 = A(2,3),
      a25 = A(2,5), a34 = A(3,4), b24 = B(2,4), b25 = B(2,5); 
    
    _ComplexD u0 = -ICPLX*a13*a13/(a34*a15*a25);
    
    return -ICPLX*u0*(2.0*Log(-s34)/3.0 - 29.0/18.0)
      + a15*a23*a23*b25*b25/(2.0*a25*a34)*F2(s15,s34)
      - a13*a23*b25/(a25*a34)*F1(s15, s34) + a12*a13*b24/(2.0*a15*a25*s34)
      + a13*a23*b25/(2.0*a25*a34*s34)
      - ICPLX*u0*(F(s25,s13,s34)+F(s15,s24,s34)+F(s24,s13,s15)+F(s13,s24,s25));
  }

  std::complex<double> 
  ampq4p1::u1lpmp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s13 = S(1,3), s15 = S(1,5), s24 = S(2,4), s25 = S(2,5),
      s34 = S(3,4);
    
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5),
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a34 = A(3,4), a35 = A(3,5),
      a45 = A(4,5), a53 = A(5,3), a54 = A(5,4), b23 = B(2,3), b25 = B(2,5),
      b35 = B(3,5); 
    
    _ComplexD h1234 = H(1,2,3,4), h1534 = H(1,5,3,4); 
    _ComplexD u0 = ICPLX*a14*a14/(a34*a15*a25);
    
    return -ICPLX*u0*(2.0*Log(-s34)/3.0 - 29.0/18.0)
      - a15*a24*a24*b25*b25/(2.0*a25*a34)*F2(s15,s34) 
      - a12*a13*a24*b23/(a15*a23*a25)*F1(s15,s34)
      - a14*a14/(a15*a25*a34)*F0(s15,s34)
      - a13*a34*b35/(a23*a35)*F1(s15,s24)
      - a12*a14/(a15*a23*a25)*F0(s15,s24)
      + a13*a45*b35/(a25*a35)*F1(s13,s24)
      + a12*a14*b23/(2.0*a15*a25*s34)
      - a14*a24*b25/(2.0*a25*a34*s34)
      - ICPLX*u0*(F(s25,s13,s34)+h1234*F(s15,s24,s34)
		  + h1534*F(s24,s13,s15) + F(s13,s24,s25));
  }
  
  std::complex<double> 
  ampq4p1::uBppp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s12 = S(1,2), s15 = S(1,5), s25 = S(2,5), s34 = S(3,4);
    
    _ComplexD a12 = A(1,2), a13 = A(1,3), a15 = A(1,5), a21 = A(2,1),
      a23 = A(2,3), a25 = A(2,5), a31 = A(3,1), a34 = A(3,4), a35 = A(3,5),
      b24 = B(2,4), b25 = B(2,5), b34 = B(3,4), b45 = B(4,5);
    
    _ComplexD h3251 = H(3,2,5,1);
    _ComplexD u0 = -ICPLX*a13*a13/(a34*a15*a25);
    
    return ICPLX*u0*(-3.0*Log(-s34) + 6.5)
      - a12*a12*a34*b24*b24/(2.0*a15*a25)*F2(s34,s15)
      - a12*a35*b25*b45/a25*F2(s12,s34) 
      + a12*a23*a35*b25/(a25*a25*a34)*(F1(s34,s12) + F1(s34,s15))
      + a13*a13/(2.0*a15*a25*a34)*F0(s34,s15) 
      + a12*a12*b24*b24/(2.0*a15*a25*b34*s15) - a13*b45/(a25*s34)
      + ICPLX*u0*(F(s34,s12,s25) + h3251*F(s34,s12,s15));
  }

  std::complex<double> 
  ampq4p1::uAppp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s13 = S(1,3), s14 = S(1,4), s15 = S(1,5), s23 = S(2,3), 
      s24 = S(2,4), s25 = S(2,5), s34 = S(3,4);
    
    _ComplexD a12 = A(1,2), a13 = A(1,3), a14 = A(1,4), a15 = A(1,5), 
      a23 = A(2,3), a24 = A(2,4), a25 = A(2,5), a31 = A(3,1), a32 = A(3,2),
      a34 = A(3,4), a35 = A(3,5), a41 = A(4,1), a42 = A(4,2), a45 = A(4,5),
      b25 = B(2,5), b45 = B(4,5);

    _ComplexD h3421 = H(3,4,2,1), h3451 = H(3,4,5,1);
    _ComplexD u0 = -ICPLX*a13*a13/(a34*a15*a25);
    
    return -a12*a23*b25/(a24*a25)*F1(s15,s34)
      - a14*a34*b45/(a24*a45)*F1(s15,s23)
      - a14*a35*b45/(a25*a45)*F1(s14,s23)
      + a12*a13/(a15*a24*a25)*F0(s34,s23)
      - ICPLX*u0*(F(s25,s13,s34) - h3421*F(s15,s23,s34) + F(s24,s13,s15)
		  - F(s14,s23,s25) + F(s15,s34,s24) - F(s25,s34,s14)
		  + F(s13,s24,s25) - h3451*F(s23,s14,s15));
  }
  
  void ampq4p1::uppp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _PrimAmp *res) const
  {
    double __nf = (double) nf;
    _ComplexD u1l = u1lppp(p1,p2,p3,p4,p5), uA = uAppp(p1,p2,p3,p4,p5), 
      uB = uBppp(p1,p2,p3,p4,p5), u0 = u0ppp(p1,p2,p3,p4,p5), 
      uNf = ICPLX*(10.0/9.0 - 2.0*Log(-S(3,4))/3.0)*u0;
 
    res->u0 = u0;
    res->u1 = ICPLX*Nc*(u1l + (uA+uB)/Nc2 - __nf*uNf/Nc);
    res->u2 = ICPLX*Nc*(uA-u1l - (2.0*uA+uB)/Nc2 + __nf*uNf/Nc);
  }
    
  void ampq4p1::upmp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _PrimAmp *res) const
  {
    double __nf = (double) nf;
    _ComplexD u1l = u1lpmp(p1,p2,p3,p4,p5), uA = uApmp(p1,p2,p3,p4,p5), 
      uB = uBpmp(p1,p2,p3,p4,p5), u0 = u0pmp(p1,p2,p3,p4,p5), 
      uNf = ICPLX*(10.0/9.0 - 2.0*Log(-S(3,4))/3.0)*u0;
 
    res->u0 = u0;
    res->u1 = ICPLX*Nc*(u1l + (uA+uB)/Nc2 - __nf*uNf/Nc);
    res->u2 = ICPLX*Nc*(uA-u1l - (2.0*uA+uB)/Nc2 + __nf*uNf/Nc);
  }

  void ampq4p1::umpp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _PrimAmp *res) const
  {
    upmp(nf, p2,p1,p4,p3,p5, res);
    res->u0 *= -1.0; res->u1 *= -1.0; res->u2 *= -1.0;
  }
  
  void ampq4p1::ummp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _PrimAmp *res) const
  {
    uppp(nf, p2,p1,p4,p3,p5, res);
    res->u0 *= -1.0; res->u1 *= -1.0; res->u2 *= -1.0;
  }
  
#define _PrimitiveAmplitudes(h1,h2)		\
  _PrimAmp u, d;				\
  u##h1##h2##p(nf, p1,p2,p3,p4,p5,&u);		\
  u##h2##h1##p(nf, p3,p4,p1,p2,p5,&d);		\
  res->A0 = Q1*u.u0 + Q2*d.u0;			\
  res->B1 = Q1*u.u1 + Q2*d.u1;			\
  res->A1 = Q1*u.u2 + Q2*d.u2

  
  void ampq4p1::
  matrix_1loop_pmpmp(double Q1, double Q2, unsigned int nf, int p1,
		     int p2, int p3, int p4, int p5, amp_1loop *res) const 
  {
    _PrimitiveAmplitudes(m,m);
  }

  void ampq4p1::
  matrix_1loop_pmmpp(double Q1, double Q2, unsigned int nf, int p1,
		    int p2, int p3, int p4, int p5, amp_1loop *res) const 
  {
    _PrimitiveAmplitudes(m,p);
  }

  void ampq4p1::
  matrix_1loop_mpmpp(double Q1, double Q2, unsigned int nf, int p1,
		     int p2, int p3, int p4, int p5, amp_1loop *res) const 
  {
    _PrimitiveAmplitudes(p,p);
  }

  void ampq4p1::
  matrix_1loop_mppmp(double Q1, double Q2, unsigned int nf, int p1,
		     int p2, int p3, int p4, int p5, amp_1loop *res) const 
  {
    _PrimitiveAmplitudes(p,m);
  }

  double ampq4p1::su3_1loop(double Q1, double Q2, unsigned int nf,
			    int p1, int p2, int p3, int p4, int p5) const
  {
    double out;
    amp_1loop M12;
    matrix_1loop_pmpmp(Q1, Q2, nf, p1,p2,p3,p4,p5, &M12); out  = amptree(M12);
    matrix_1loop_pmmpp(Q1, Q2, nf, p1,p2,p3,p4,p5, &M12); out += amptree(M12); 
    matrix_1loop_mppmp(Q1, Q2, nf, p1,p2,p3,p4,p5, &M12); out += amptree(M12);
    matrix_1loop_mpmpp(Q1, Q2, nf, p1,p2,p3,p4,p5, &M12); out += amptree(M12);
    
    swap();
    matrix_1loop_pmpmp(Q1, Q2, nf, p1,p2,p3,p4,p5, &M12); out += amptree(M12);
    matrix_1loop_pmmpp(Q1, Q2, nf, p1,p2,p3,p4,p5, &M12); out += amptree(M12); 
    matrix_1loop_mppmp(Q1, Q2, nf, p1,p2,p3,p4,p5, &M12); out += amptree(M12);
    matrix_1loop_mpmpp(Q1, Q2, nf, p1,p2,p3,p4,p5, &M12); out += amptree(M12);
    swap();
    
    return out;
  }
  

  void ampq4p1::su3_1loop(double Q, unsigned int nf, int p1, int p2, 
			  int p3, int p4, int p5, double *out) const
  {
    double tmp;
    amp_1loop M12, M14;
    
    matrix_1loop_pmpmp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
    matrix_1loop_pmpmp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
    out[0] = amptree(M12); out[1] = amptree(M12,M14);  
    
    matrix_1loop_pmmpp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
    matrix_1loop_pmmpp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
    out[0] += (tmp = amptree(M12)); out[1] += tmp + amptree(M14);  
    
    matrix_1loop_mppmp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
    matrix_1loop_mppmp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
    out[0] += (tmp = amptree(M12)); out[1] += tmp + amptree(M14);  
    
    matrix_1loop_mpmpp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
    matrix_1loop_mpmpp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
    out[0] += amptree(M12); out[1] += amptree(M12,M14);  

    swap();
    matrix_1loop_pmpmp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
    matrix_1loop_pmpmp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
    out[0] += amptree(M12); out[1] += amptree(M12,M14);  
    
    matrix_1loop_pmmpp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
    matrix_1loop_pmmpp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
    out[0] += (tmp = amptree(M12)); out[1] += tmp + amptree(M14);  
    
    matrix_1loop_mppmp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
    matrix_1loop_mppmp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
    out[0] += (tmp = amptree(M12)); out[1] += tmp + amptree(M14);  
    
    matrix_1loop_mpmpp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
    matrix_1loop_mpmpp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
    out[0] += amptree(M12); out[1] += amptree(M12,M14);  
    swap();
  }
  
  double ampq4p1::
  su3_1loop_mch(double Q1, double Q2, unsigned int nf,
		int p1, int p2, int p3, int p4, int p5) const
  {
    amp_1loop M12;
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (4*_M_rng());
    
    if(hpm == 1) swap();
    switch(hel) {
    case 0: matrix_1loop_pmpmp(Q1, Q2, nf, p1, p2, p3, p4, p5, &M12); break;
    case 1: matrix_1loop_mpmpp(Q1, Q2, nf, p1, p2, p3, p4, p5, &M12); break;
    case 2: matrix_1loop_pmmpp(Q1, Q2, nf, p1, p2, p3, p4, p5, &M12); break;
    case 3: matrix_1loop_mppmp(Q1, Q2, nf, p1, p2, p3, p4, p5, &M12); break;
    }
    if(hpm == 1) swap();
    
    return 8.0*amptree(M12);
  }
  
  void ampq4p1::su3_1loop_mch(double Q, unsigned int nf, int p1, int p2, 
			      int p3, int p4, int p5, double *out) const
  {
    amp_1loop M12, M14;
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (4*_M_rng());
    
    if(hpm == 1) swap();
    switch(hel) {
    case 0: 
      matrix_1loop_pmpmp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
      matrix_1loop_pmpmp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
      break;
    case 1:
      matrix_1loop_mpmpp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
      matrix_1loop_mpmpp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
      break;
    case 2:
      matrix_1loop_pmmpp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
      matrix_1loop_pmmpp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
      break;
    case 3:
      matrix_1loop_mppmp(Q, Q, nf, p1, p2, p3, p4, p5, &M12);
      matrix_1loop_mppmp(Q, Q, nf, p1, p4, p3, p2, p5, &M14);
      break;
    }
    if(hpm == 1) swap();

    out[0] = 8.0*amptree(M12);
    if(hel < 2) out[1] = 8.0*amptree(M12,M14);
    else out[1] = out[0] + 8.0*amptree(M14);
  }
}


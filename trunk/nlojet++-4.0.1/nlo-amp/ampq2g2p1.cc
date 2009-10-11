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
#include "ampq2g2p1.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {

  std::complex<double> 
  ampq2g2p1::Apm(int pi, int p1, int p2, int p3, int p4, int p5) const
  {
    std::complex<double> res(0.0, -1.0);
    
    res *= A(1,i)*pow(A(2,i),3);
    res /= A(2,3)*A(3,4)*A(4,1)*A(2,5)*A(5,1);
    return res;
  }

  std::complex<double> 
  ampq2g2p1::Amp(int pi, int p1, int p2, int p3, int p4, int p5) const
  {
    std::complex<double> res(0.0, 1.0);
    
    res *= A(2,i)*pow(A(1,i),3);
    res /= A(2,3)*A(3,4)*A(4,1)*A(2,5)*A(5,1);
    return res;
  }


#define AmpTree (Na*real(A1*conj(A1) + A2*conj(A2)) - 2.0*real(A1*conj(A2)))

  
  double ampq2g2p1::su3_tree(int p1, int p2, int p3, int p4, int p5) const
  {
    int p[3] = {p3,p4,p5};
    double ret_val = 0.0;
    _ComplexD A1, A2;
    
    for(unsigned i = 0; i < 3; i++) {
      A1 = Apm(p[i], p1,p2, p3,p4, p5);
      A2 = Apm(p[i], p1,p2, p4,p3, p5);
      ret_val += AmpTree;
      
      A1 = Amp(p[i], p1,p2, p3,p4, p5);
      A2 = Amp(p[i], p1,p2, p4,p3, p5);
      ret_val += AmpTree;
    }
    
    return 4.0*Na*ret_val/Nc;
  }


  double ampq2g2p1::su3_tree_mch(int p1, int p2, int p3, int p4, int p5) const
  {
    int p[6] = {p3,p4,p5, p3,p4,p5};
    unsigned int hel = (unsigned int) (6*_M_rng());
    _ComplexD A1, A2;
    
    if(hel < 4) {
      A1 = Apm(p[hel], p1,p2, p3,p4, p5);
      A2 = Apm(p[hel], p1,p2, p4,p3, p5);
    } else {
      A1 = Amp(p[hel], p1,p2, p3,p4, p5);
      A2 = Amp(p[hel], p1,p2, p4,p3, p5);
    }
    
    return 24.0*Na*AmpTree/Nc;
  }


  double ampq2g2p1::ampcc(double c1, double c2, double c3,
			  int p1, int p2, int p3, int p4, int p5) const
  {
    int p[3] = {p3,p4,p5};
    double S1 = 0.0, S2 = 0.0, S3 = 0.0;
    _ComplexD A1, A2;
    
    for(unsigned i = 0; i < 3; i++) {
      A1 = Apm(p[i], p1,p2, p3,p4, p5);
      A2 = Apm(p[i], p1,p2, p4,p3, p5);
      S1 += real(A1*conj(A1)); 
      S2 += real(A2*conj(A2));
      S3 += 2.0*real(A1*conj(A2)); 
            
      A1 = Amp(p[i], p1,p2, p3,p4, p5);
      A2 = Amp(p[i], p1,p2, p4,p3, p5);
      S1 += real(A1*conj(A1)); 
      S2 += real(A2*conj(A2));
      S3 += 2.0*real(A1*conj(A2)); 
    }
    
    return 2.0*Na*(c1*S1+c2*S2+c3*S3);
  }


  std::complex<double> 
  ampq2g2p1::amphcc(double c1, double c2, double c3, int p, 
		    int p1, int p2, int p3, int p4, int p5) const
  {
    int pi, pj;
    if(p == p3) pi = p4, pj = p5;
    else if(p == p4) pi = p3, pj = p5;
    else throw("Error in ampq2g2p1::amphtree");

    _ComplexD a23 = A(2,3), a34 = A(3,4), a41 = A(4,1), a24 = A(2,4), 
      a13 = A(1,3), a15 = A(1,5), a25 = A(2,5), a1i = A(1,i), a2i = A(2,i),
      a1j = A(1,j), a2j = A(2,j);
    
    _ComplexD num = 2.0*a1i*a2i*a1j*a2j*(a1i*a1i*a2j*a2j+a1j*a1j*a2i*a2i)/(a25*a25*a15*a15);
    _ComplexD m1 = 1.0/(a23*a34*a41), m2 = 1.0/(a24*a34*a13);
    
    return -Na*num*(c1*m1*m1 + c2*m2*m2 + 2.0*c3*m1*m2);
  }
   

#define __cc_switch12 (pi == p1 && pj == p2) || (pi == p2 && pj == p1)
#define __cc_switch34 (pi == p3 && pj == p4) || (pi == p4 && pj == p3)

#define __cc_switch13                                   \
  (pi == p1 && pj == p3) || (pi == p2 && pj == p4) ||   \
  (pi == p3 && pj == p1) || (pi == p4 && pj == p2)

#define __cc_switch14                                   \
  (pi == p1 && pj == p4) || (pi == p2 && pj == p3) ||   \
  (pi == p4 && pj == p1) || (pi == p3 && pj == p2)
  
  std::pair<double, std::complex<double> > ampq2g2p1::
  su3_cc(int pi,int pj, int p1, int p2, int p3, int p4, int p5) const
  {
    double c1, c2, c3;
    _ComplexD hcc(0.0);
    
    if(__cc_switch12) { c1 = c2 = Na/Nc2-1.0; c3 = Na/Nc2-2.0;}
    else if(__cc_switch13) { c1 = c3 = 1.0; c2 = -Na;} 
    else if(__cc_switch14) { c2 = c3 = 1.0; c1 = -Na;}
    else if(__cc_switch34) { c1 = c2 = -Nc2; c3 = 0.0;}
    else throw("Error in ampq2g2p1::su3_cc");     
    
    if(pj == p3 || pj == p4) hcc = amphcc(c1,c2,c3, pj, p1,p2,p3,p4,p5);
    return _Pair(ampcc(c1,c2,c3, p1,p2,p3,p4,p5), hcc);
  }

  double ampq2g2p1::
  ampcc(int pi,int pj, int p1, int p2, int p3, int p4, int p5) const
  {
    double c1, c2, c3;
    if(__cc_switch12) { c1 = c2 = Na/Nc2-1.0; c3 = Na/Nc2-2.0;}
    else if(__cc_switch13) { c1 = c3 = 1.0; c2 = -Na;} 
    else if(__cc_switch14) { c2 = c3 = 1.0; c1 = -Na;}
    else if(__cc_switch34) { c1 = c2 = -Nc2; c3 = 0.0;}
    else throw("Error in ampq2g2p1::ampcc");     
    return ampcc(c1,c2,c3, p1,p2,p3,p4,p5);
  }
  

#define X(i,j) (-std::log(std::abs(Sij(i,j))))
  
  void ampq2g2p1::su3_kp(unsigned int nf, int pa, int p1, int p2, int p3,
			 int p4, int p5, su3_kp_i1 *res, double al) const
  {
    double cc, xi, xj, xq = Gq/Cf, xg = Gg(nf)/Nc;
    int p[4] = {p1,p2,p3,p4};
    
    double Ta = (pa == p1 || pa == p2 ? Cf : Ca);
    double x[4] = {xq, xq, xg, xg};
    
    res->tree = res->pa = res->ga = res->loop = 0.0;
    for(unsigned i = 0; i < 4; i++)
      for(unsigned j = i+1; j < 4; j++) {
	cc = ampq2g2p1::ampcc(p[i], p[j], p1, p2, p3, p4, p5);
	
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
	  res->ga += x[iq]*cc;
	}
      }
    
    res->tree /= Ta; res->pa /= Ta; 
    res->loop += (2.0*(Gg(nf)+Kg(nf,al) + Gq+Kq(al)) - Cf)*(res->tree);
  }

#define __PrimitiveAmplitudes(h3,h4,h5)		\
  A1##h3##h4##h5(p1,p2, p3,p4,p5, A[0]);	\
  A1##h4##h5##h3(p1,p2, p4,p5,p3, A[3]);	\
  A1##h5##h3##h4(p1,p2, p5,p3,p4, A[2]);	\
  			       	                \
  A1##h3##h5##h4(p1,p2, p3,p5,p4, A[1]);	\
  A1##h5##h4##h3(p1,p2, p5,p4,p3, A[4]);	\
  A1##h4##h3##h5(p1,p2, p4,p3,p5, A[5]);	\
						\
						\
  A2##h3##h4##h5(p1, p3, p2, p4,p5, B[0]);	\
  A2##h4##h5##h3(p1, p4, p2, p5,p3, B[3]);	\
  A2##h5##h3##h4(p1, p5, p2, p3,p4, B[2]);	\
  		               	                \
  A2##h3##h5##h4(p1, p3, p2, p5,p4, B[1]);	\
  A2##h5##h4##h3(p1, p5, p2, p4,p3, B[4]);	\
  A2##h4##h3##h5(p1, p4, p2, p3,p5, B[5])

#define Qu  0.66666666666666666666 
#define Qd -0.33333333333333333333
 
  void ampq2g2p1::
  matrix_1loop(unsigned int nu, unsigned int nd, const _AmpPrim *A, const _AmpPrim *B, amp_1loop *res)
  { 
    unsigned int idx, nf = nu + nd;
    double TrQ = (nu*Qu + nd*Qd);
    _ComplexD Tsum[2] = {0.0, 0.0}, Csum[2] = {0.0, 0.0}, Asum[2] = {0.0, 0.0}, Bsum(0.0);
    
    for(unsigned i = 0; i < 6; i++) {
      idx = i < 3 ? 0 : 1;
      
      Tsum[idx] += A[i].Atree;
      Asum[idx] += A[i].Asusy - A[i].AL;
      Csum[idx] += A[i].Asf;
      
      if(i == 2 || i == 4) Bsum -= B[i].AL;
      else Bsum += B[i].Asusy - B[i].AL + B[i].Asf;
    }

    res->U0[0] = Qu*Tsum[0];
    res->U0[1] = Qu*Tsum[1];
    
    res->D0[0] = Qd*Tsum[0];
    res->D0[1] = Qd*Tsum[1];
    
    res->U1[0] = -Qu*B[2].AL + (nf*Qu-TrQ)*B[2].Asf/Nc - Qu*Asum[0]/Nc2 - (Qu+Nc*TrQ)*Csum[0]/Nc2;
    res->U1[1] = -Qu*B[4].AL + (nf*Qu-TrQ)*B[4].Asf/Nc - Qu*Asum[1]/Nc2 - (Qu+Nc*TrQ)*Csum[1]/Nc2;
    
    res->D1[0] = -Qd*B[2].AL + (nf*Qd-TrQ)*B[2].Asf/Nc - Qd*Asum[0]/Nc2 - (Qd+Nc*TrQ)*Csum[0]/Nc2;
    res->D1[1] = -Qd*B[4].AL + (nf*Qd-TrQ)*B[4].Asf/Nc - Qd*Asum[1]/Nc2 - (Qd+Nc*TrQ)*Csum[1]/Nc2;
    
    res->U3 = Qu*(Asum[0]+Asum[1] + Bsum) + (Qu+TrQ/Nc)*(Csum[0]+Csum[1]);
    res->D3 = Qd*(Asum[0]+Asum[1] + Bsum) + (Qd+TrQ/Nc)*(Csum[0]+Csum[1]);
  } 

  void ampq2g2p1::
  matrix_1loop_mpmpp(unsigned int nu, unsigned int nd, int p1, int p2,
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(m,p,p); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g2p1::
  matrix_1loop_mppmp(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(p,m,p); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g2p1::
  matrix_1loop_mpppm(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(p,p,m); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g2p1::
  matrix_1loop_mppmm(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(p,m,m); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g2p1::
  matrix_1loop_mpmpm(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(m,p,m); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g2p1::
  matrix_1loop_mpmmp(unsigned int nu, unsigned int nd, int p1, int p2, 
		     int p3, int p4, int p5, amp_1loop *res) const
  {
    static _AmpPrim A[6], B[6];
    __PrimitiveAmplitudes(m,m,p); 
    matrix_1loop(nu, nd, A, B, res);
  }

  void ampq2g2p1::amp1loop(const amp_1loop& A, double *out)
  {
    out[0] += 2.0*Na*(real(conj(A.U0[0])*(Na*A.U1[0] - A.U1[1] + A.U3))
		      + real(conj(A.U0[1])*(Na*A.U1[1] - A.U1[0] + A.U3)));
    
    out[1] += 2.0*Na*(real(conj(A.D0[0])*(Na*A.D1[0] - A.D1[1] + A.D3))
		      + real(conj(A.D0[1])*(Na*A.D1[1] - A.D1[0] + A.D3)));
  }

  void ampq2g2p1::
  su3_1loop(unsigned int nu, unsigned int nd, int p1, int p2,
	    int p3, int p4, int p5, double *out) const
  {
    static amp_1loop M;
    out[0] = out[1] = 0.0;
    
    matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    		  
    swap();	  
    matrix_1loop_mpppm(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mppmp(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mpmpp(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mppmm(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mpmpm(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    matrix_1loop_mpmmp(nu, nd, p1,p2,p3,p4,p5, &M); amp1loop(M, out);
    swap();
  }
    
  void ampq2g2p1::
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
    
    out[0] = out[1] = 0.0;
    amp1loop(M, out);
    
    out[0] *= 12.0;
    out[1] *= 12.0;
  }
}  //  namespace nlo

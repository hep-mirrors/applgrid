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
#include "ampq4.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {


  void ampq4::su3_tree(int p1, int p2, int p3, int p4, double *out) const
  {
    static double cfac = 2.0*(Nc2-1.0);
    double s = S(1,2), t = S(1,3), u = S(1,4);
    out[0] = cfac*(t*t+u*u)/(s*s);
    out[1] = out[0] + cfac*((t*t+s*s)/(u*u) - 2.0*t*t/(Nc*u*s));
  }

  void ampq4::ampcc12(int p1, int p2, int p3, int p4, double *out) const
  {
    static double cfac = -2.0*(Nc2-1.0)/Nc;
    double s = S(1,3), t = S(1,2), u = S(1,4);
    double At = (s*s+u*u)/(t*t), Au = (s*s+t*t)/(u*u), B = s*s/(t*u);
    
    out[0] = -0.5*cfac*At;
    out[1] = cfac*(B - Nc*(0.5*At + Au) + 0.5*Nc*Nc2*Au)/Nc;
  }

  void ampq4::ampcc13(int p1, int p2, int p3, int p4, double *out) const
  {
    static double cfac = -2.0*(Nc2-1.0)/Nc;
    double s = S(1,3), t = S(1,2), u = S(1,4);
    double At = (s*s+u*u)/(t*t), Au = (s*s+t*t)/(u*u), B = s*s/(t*u);
    
    out[0] = cfac*At;
    out[1] = cfac*(-B + Nc*(At + Au) - Nc2*B)/Nc;
  }
 
  void ampq4::ampcc14(int p1, int p2, int p3, int p4, double *out) const
  {
    static double cfac = -2.0*(Nc2-1.0)/Nc;
    double s = S(1,3), t = S(1,2), u = S(1,4);
    double At = (s*s+u*u)/(t*t), Au = (s*s+t*t)/(u*u), B = s*s/(t*u);
    
    out[0] = 0.5*(Nc2-2.0)*cfac*At;
    out[1] = cfac*(B - Nc*(0.5*Au + At) + 0.5*Nc*Nc2*At)/Nc;
  }
 
#define __cc_switch12 					\
  (pi == p1 && pj == p2) || (pi == p2 && pj == p1) ||	\
  (pi == p3 && pj == p4) || (pi == p4 && pj == p3)
  
#define __cc_switch13 					\
  (pi == p1 && pj == p3) || (pi == p2 && pj == p4) ||	\
  (pi == p3 && pj == p1) || (pi == p4 && pj == p2)

#define __cc_switch14 					\
  (pi == p1 && pj == p4) || (pi == p2 && pj == p3) ||	\
  (pi == p4 && pj == p1) || (pi == p3 && pj == p2)
  
  void ampq4::
  su3_cc(int pi,int pj, int p1, int p2, int p3, int p4, double *cc) const
  {
    if(__cc_switch12) ampcc12(p1, p2, p3, p4, cc);
    else if(__cc_switch13) ampcc13(p1, p2, p3, p4, cc);
    else if(__cc_switch14) ampcc14(p1, p2, p3, p4, cc);
    else throw("Error in ampq2g2::amp_cc");	
  }


#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampq4::su3_kp(int pa, int pb, int p1, int p2,
		     int p3, int p4, su3_kp_i2 *res, double al) const
  {
    double cc[2], xi;
    int p[4] = {p1,p2,p3,p4};

    res[0].tree=res[0].pa=res[0].pb=res[0].loop = 0.0;
    res[1].tree=res[1].pa=res[1].pb=res[1].loop = 0.0;
    for(unsigned i = 0; i < 4; i++)
      for(unsigned j = i+1; j < 4; j++) {
	su3_cc(p[i], p[j], p1, p2, p3, p4, cc);
	
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

	  if(q == pb) {
	    res[0].cca = res[0].ccb = -cc[0]/Cf;
	    res[1].cca = res[1].ccb = -cc[1]/Cf;
	  }
	}
	
	if(p[i] == pb || p[j] == pb) {
	  int q = p[i] == pb ? p[j] : p[i];

	  xi = X(pb, q);
	  res[0].pb += xi*cc[0];
	  res[1].pb += xi*cc[1];
	}
      }
    
    res[0].tree /= Cf; res[0].pa /= Cf; res[0].pb /= Cf;
    res[1].tree /= Cf; res[1].pa /= Cf; res[1].pb /= Cf;
    
    res[0].ga = res[0].gb = -Gq*(res[0].tree-res[0].cca);
    res[1].ga = res[1].gb = -Gq*(res[1].tree-res[1].cca);
    
    xi = 4.0*(Gq+Kq(al)); // - 2.0*Cf + Ca/3.0; I use the KST amps. No need conversion!
    res[0].loop += xi*res[0].tree;
    res[1].loop += xi*res[1].tree;
  }

#define NLO_PI2 11.31404884553380306325

  void ampq4::
  matrix_1loop_pmpm(unsigned int __nf, int p4, int p1, int p3, int p2, amp_1loop& amp) const
  {
    double nf = (double) __nf, s12 = S(1,2), s13 = S(1,3), s14 = S(1,4);
    _ComplexD l12 = Log(-s12), l13 = Log(-s13), l14 = Log(-s14);
    
    _ComplexD Vmm = -2.0*l14/3.0 + std::pow(l14-l13,2) + NLO_PI2;
    _ComplexD Vmp = -2.0*l14/3.0 + std::pow(l14-l12,2) + NLO_PI2;
    _ComplexD Vsf = -2.0*l14/3.0 + 1.11111111111111111111;
    _ComplexD Vsl = 3.0*l14 - 8.0;
    _ComplexD F = 0.5*s14/s12*(1.0-s13/s12)*(std::pow(l14-l13,2) + NLO_PI2)
      + s14/s12*(l14-l13);

    amp.A0 = ICPLX*A(1,2)*B(3,4)/s14;
    amp.A1 =  amp.A0*(Vmm+F - nf*Vsf/Nc - 2.0*(Vmm-Vmp+F)/Nc2 - Vsl/Nc2);
    amp.A2 = -amp.A0*(Vmp - nf*Vsf/Nc - (Vmm-Vmp+Vsl+F)/Nc2);
  }   
  
  void ampq4::
  matrix_1loop_pmmp(unsigned int __nf, int p4, int p1, int p3, int p2, amp_1loop& amp) const
  {
    double nf = (double) __nf, s12 = S(1,2), s13 = S(1,3), s14 = S(1,4);
    _ComplexD l12 = Log(-s12), l13 = Log(-s13), l14 = Log(-s14);

    _ComplexD Vmm = -2.0*l14/3.0 + std::pow(l14-l13,2) + NLO_PI2;
    _ComplexD Vmp = -2.0*l14/3.0 + std::pow(l14-l12,2) + NLO_PI2;
    _ComplexD Vsf = -2.0*l14/3.0 + 1.11111111111111111111;
    _ComplexD Vsl = 3.0*l14 - 8.0;
    _ComplexD F = 0.5*s14/s13*(1.0-s12/s13)*(std::pow(l14-l12,2) + NLO_PI2)
      + s14/s13*(l14-l12);

    amp.A0 = ICPLX*A(1,3)*B(2,4)/s14;
    amp.A1 =  amp.A0*(Vmm - nf*Vsf/Nc - 2.0*(Vmm-Vmp-F)/Nc2 - Vsl/Nc2);
    amp.A2 = -amp.A0*(Vmp+F - nf*Vsf/Nc - (Vmm-Vmp+Vsl-F)/Nc2);
  }   
  
  void ampq4::
  su3_1loop(unsigned int nf, int p1, int p2, int p3, int p4, double *amp) const
  {
    amp_1loop A, B;
    
    matrix_1loop_pmpm(nf, p1,p2,p3,p4, A); 
    matrix_1loop_pmpm(nf, p1,p4,p3,p2, B);
    amp[0] = real(A.A1*conj(A.A0));
    amp[1] = real((A.A1-B.A2/Nc)*conj(A.A0) + (B.A1-A.A2/Nc)*conj(B.A0));
    
    matrix_1loop_pmmp(nf, p1,p2,p3,p4, A); 
    matrix_1loop_pmmp(nf, p1,p4,p3,p2, B);
    amp[0] += real(A.A1*conj(A.A0));
    amp[1] += real(A.A1*conj(A.A0)) + real(B.A1*conj(B.A0));;
    
    swap();
    matrix_1loop_pmpm(nf, p1,p2,p3,p4, A); 
    matrix_1loop_pmpm(nf, p1,p4,p3,p2, B);
    amp[0] += real(A.A1*conj(A.A0));
    amp[1] += real((A.A1-B.A2/Nc)*conj(A.A0) + (B.A1-A.A2/Nc)*conj(B.A0));

    matrix_1loop_pmmp(nf, p1,p2,p3,p4, A); 
    matrix_1loop_pmmp(nf, p1,p4,p3,p2, B);
    amp[0] += real(A.A1*conj(A.A0));
    amp[1] += real(A.A1*conj(A.A0)) + real(B.A1*conj(B.A0));;
    swap();
    
    amp[0] *= Nc*Na;
    amp[1] *= Nc*Na;
  }
  
  void ampq4::
  su3_1loop_mch(unsigned int nf, int p1, int p2, int p3, int p4, double *amp) const
  {
    amp_1loop A, B;
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (2*_M_rng());

    if(hpm == 1) swap();
    switch(hel){
    case 0:
      matrix_1loop_pmpm(nf, p1,p2,p3,p4, A); 
      matrix_1loop_pmpm(nf, p1,p4,p3,p2, B);
      amp[0] = 4.0*Nc*Na*real(A.A1*conj(A.A0));
      amp[1] = 4.0*Nc*Na*real((A.A1-B.A2/Nc)*conj(A.A0) + (B.A1-A.A2/Nc)*conj(B.A0));
      break;
    case 1:
      matrix_1loop_pmmp(nf, p1,p2,p3,p4, A); 
      matrix_1loop_pmmp(nf, p1,p4,p3,p2, B);
      amp[0] = amp[1] = 4.0*Nc*Na*real(A.A1*conj(A.A0));
      amp[1] += 4.0*Nc*Na*real(B.A1*conj(B.A0));
      break;
    }
    if(hpm == 1) swap();
  }
}

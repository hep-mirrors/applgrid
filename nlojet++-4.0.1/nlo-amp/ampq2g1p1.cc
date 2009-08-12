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
#include "ampq2g1p1.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {

   
  double ampq2g1p1::su3_tree(int p1, int p2, int p3, int p4) const
  {
    double s13 = S(1,3), s14 = S(1,4), s23 = S(2,3), s24 = S(2,4); 
    return 2.0*Na*(s13*pow(s23,3)+s23*pow(s13,3)
		   +s14*pow(s24,3)+s24*pow(s14,3))/(s23*s13*s14*s24);
  }

  
#define __cc_switch12					\
(pi == p1 && pj == p2) || (pi == p2 && pj == p1)

#define __cc_switch13					\
(pi == p1 && pj == p3) || (pi == p3 && pj == p1) ||	\
(pi == p3 && pj == p2) || (pi == p2 && pj == p3)


  double ampq2g1p1::
  su3_cc(int pi,int pj, int p1, int p2, int p3, int p4) const
  {
    double s = 0.0;
    if(__cc_switch12)      s = 0.5/Nc;
    else if(__cc_switch13) s = -0.5*Nc;
    else throw("Error in ampq2g1p1::su3_cc(...)");   
    return s*su3_tree(p1, p2, p3, p4);
  }

#define X(i,j) (-std::log(std::abs(Sij(i,j))))

  void ampq2g1p1::su3_kp(unsigned int nf, int pa,  int p1, 
			 int p2, int p3, int p4, su3_kp_i1 *res, 
			 double al) const
  {
    double cc, xi, xj, xq = Gq/Cf, xg = Gg(nf)/Nc;
    int pi, pj, p[3] = {p1,p2,p3};
    
    double Ta = (pa == p1 || pa == p2 ? Cf : Ca);
    double x[3] = {xq, xq, xg};
    
    res->tree = su3_tree(p1,p2, p3, p4);
    res->pa = res->ga = res->loop = 0.0;
    
    for(unsigned i = 0; i < 3; i++)
      for(unsigned j = i+1; j < 3; j++) {
	pi = p[i]; pj = p[j];
	
	if(__cc_switch12)      cc = 0.5*(res->tree)/Nc;
	else if(__cc_switch13) cc = -0.5*Nc*(res->tree);
	else throw("Error in ampq2g1p1::su3_kp(...)");
	
	//----- loop contribution -----
	if(i < 2) {
	  xi = Xq(Sij(pi,pj), 1.0);
	  xj = j < 2 ? xi : Xg(Sij(pi,pj), 1.0, nf);
	} else {
	  xi = Xg(Sij(pi,pj), 1.0, nf);
	  xj = j < 2 ? Xq(Sij(pi,pj), 1.0) : xi;
	}
	
	res->loop += (xi+xj)*cc;
	
	//----- finite contributions -----
	if(pi == pa || pj == pa) {
	  int iq = pi == pa ? j : i;
	  res->pa += X(pa, p[iq])*cc;
	  res->ga += x[iq]*cc;
	}
      }
    
    res->pa /= Ta;
    res->loop += (Gg(nf)+Kg(nf,al) + 2.0*(Gq+Kq(al)) - Cf)*(res->tree);
  }


#define __PrimitiveAmplitudes(h3,h4)		\
  A1##h3##h4(p1,p2, p3,p4, A[0]);		\
  A1##h4##h3(p1,p2, p4,p3, A[1]);		\
  						\
  A2##h4##h3(p1, p4, p2, p3, B)

#define __amp_1loop				\
  (-2.0*Nc*Na*(B.AL + (A[0].Asusy+A[1].Asusy	\
		       -A[0].AL-A[1].AL)/Nc2)	\
   *(A[0].Atree + A[1].Atree))
  

  double ampq2g1p1::matrix_1loop_pmpm(int p1, int p2, int p3, int p4) const
  {
    static _AmpPrim A[2], B;
    __PrimitiveAmplitudes(m,p); 
    return real(__amp_1loop);
  }
  
  double ampq2g1p1::matrix_1loop_ppmm(int p1, int p2, int p3, int p4) const
  {
    static _AmpPrim A[2], B;
    __PrimitiveAmplitudes(p,m); 
    return real(__amp_1loop);
  }
  
  double ampq2g1p1::su3_1loop(int p1, int p2, int p3, int p4) const
  {
    double ret_val = 0.0;

    ret_val  = matrix_1loop_ppmm(p1,p2,p3,p4);
    ret_val += matrix_1loop_pmpm(p1,p2,p3,p4);
				      
    swap();			      
    ret_val += matrix_1loop_ppmm(p1,p2,p3,p4);
    ret_val += matrix_1loop_pmpm(p1,p2,p3,p4);
    swap();

    return ret_val;
  }

  double ampq2g1p1::su3_1loop_mch(int p1, int p2, int p3, int p4) const
  {
    double ret_val = 0.0;
    unsigned int hpm = (unsigned int) (2*_M_rng());
    unsigned int hel = (unsigned int) (2*_M_rng());
    
    if(hpm == 1) swap();
    switch(hel){
    case 0: ret_val = matrix_1loop_ppmm(p1,p2,p3,p4); break;
    case 1: ret_val = matrix_1loop_pmpm(p1,p2,p3,p4); break;
    }
    if(hpm == 1) swap();

    return 4.0*ret_val;
  }

  



}  //  namespace nlo

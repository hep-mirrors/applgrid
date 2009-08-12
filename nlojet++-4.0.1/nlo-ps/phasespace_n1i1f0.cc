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

//   Standard includes
#include <algorithm>

//   nlo includes
#include "bits/psg-phasespace_n1i1f0.h"



#define __DEB_PI    3.14159265358979323846
#define __DEB_2PI   6.28318530717958647692
#define __DEB_4PI2  157.91367041742973790108
#define __DEB_1_8PI 0.03978873577297383395



namespace nlo {
  
  typedef hadronic_event<lorentzvector<double>, hadronic_event_traits<1U,1U,0U> > _Event110;
  
  
  double basic_phasespace<_Event110>::_M_generate_event(double x, double y, event_type& p)
  {
    unsigned int up = p.upper();
    double s = 4.0*_M_el*_M_eh;
    
    if(up < 1) throw "unable to generate dis event";
        
    //----- weight of the leptonic phase space -----
    double weight = s*y/__DEB_4PI2;
    
    //----- incoming and outgoing lepton -----
    const random_generator& __rng = *_M_rng;

    double aa = _M_eh*x*y, bb = _M_el*(1.0 - y);
    double El = aa + bb, pz = aa - bb;
    double pt = std::sqrt(4.0*aa*bb), phi = __DEB_2PI*__rng();
    
    p[-2] = lorentzvector<double>(pt*std::cos(phi), pt*std::sin(phi), pz, El);
    p[-1] = lorentzvector<double>(0.0, 0.0, -_M_el, _M_el);
    
    //----- incoming hadron -----
    p[hadron(0)] = lorentzvector<double>(0.0, 0.0, _M_eh, _M_eh);
    
    //----- incoming and outgoing partons
    if(up == 1) {
      p[0] = x*p[hadron(0)];
      p[1] = p[0] + p[-1] - p[-2];
      
      weight *= __DEB_2PI/(s*y);
    } else {
      //----- generate the momentum fraction eta -----
      double lxb = -std::log(x);
      double eta = x*std::exp(lxb*__rng());
      weight *= lxb*eta;
      
      //----- incoming parton -----
      p[0] = eta*p[hadron(0)];
      
      //----- generate the outgoing partons ----
      lorentzvector<double> qq = p[0] + p[-1] - p[-2];
      threevector<double>   bb = qq.boostVector();
      
      event_type::iterator pbegin = p.begin() + 3;
      event_type::iterator pend   = p.end();
      
      weight *= _M_psgen -> operator()(s*y*(eta - x), pbegin, pend);
      
      //----- boost the outgoing partons to the labor frame -----
      for(event_type::iterator pi = pbegin; pi < pend; pi++) 
		pi -> boost(bb);
    }
    
    return weight;
  }
  
  double basic_phasespace<_Event110>::
  operator()(const event_type& q, event_type& p)
  {
    unsigned int up = p.upper(), uq = q.upper();
    if(uq < 1 || up < uq)  throw "unable to generate DIS event";
    
    //----- copy from q to p -----
    event_type::iterator pi;
    event_type::const_iterator qi;
    
    p[hadron(0)] = q[hadron(0)];
    for(pi = p.begin(), qi = q.begin(); qi < q.end(); pi++, qi++) 
      *pi = *qi; 
    
    if(up == uq) return 1.0;
    
    //----- dipole emmisions -----
    __fe_clear_exception();
    
    unsigned int n;
    int emit, omit, thd;
    double jac, weight = 1.0;
    const random_generator& __rng = *_M_rng;
    
    event_type::const_iterator p1 = p.begin()+3, pn;
    lorentzvector<double> pab = p[-1]/(p[-1]*p[hadron(0)]);
    
    for(n = uq+1, pn = p1 + uq+1; n <= up; n++, pn++) {
      //---- choose which one emit, omit, third ----    
      emit = (int) ((n-1)*__rng()) + 1;
      omit = (int) ((n-1)*__rng()) + 1;
      thd  = (int) (n*__rng()) + 1;
      if(emit == omit) omit = 0;
      
      //---- generate the dipole emmission ----
      if(omit != 0) dipole_emission::gendip_fff(__rng, _M_beta, _M_eps, p[emit], p[n], p[omit]);
      else dipole_emission::gendip_ffi(__rng, _M_beta, _M_eps, pab*p[0], p[emit], p[n], p[omit]);
      
      if(thd != (int) n) std::swap(p[thd], p[n]);
      
      //----- calculate the weight -----
      jac = dipole_emission::jacobi_fff(_M_beta, _M_eps, p1, pn)
	+   dipole_emission::jacobi_ffi(_M_beta, _M_eps, pab*p[0], p[0], p1, pn);
      weight *= n*(n-1)*(n-1)/jac;
      
      __fe_throw_exception();
    }
    
    return weight;
  }
}   //  namespace nlo

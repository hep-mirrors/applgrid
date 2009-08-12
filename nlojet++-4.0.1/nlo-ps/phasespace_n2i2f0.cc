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
#include "bits/psg-phasespace_n2i2f0.h"


#define __DEB_PI    3.14159265358979323846
#define __DEB_2PI   6.28318530717958647692
#define __DEB_4PI2  157.91367041742973790108
#define __DEB_1_8PI 0.03978873577297383395




namespace nlo {

  typedef hadronic_event<lorentzvector<double>, hadronic_event_traits<2U,2U,0U> > _Event220;
  
 
  void basic_phasespace<_Event220>::phasespace_cuts(double q2min, double q2max) 
  {
	if(q2min <= 0.0)  throw "basic_phasespace<...> : Q2min must be greater than zero";
	if(q2min > q2max) throw "basic_phasespace<...> : Q2min must be less than Q2max";
	if(q2min > _M_s)  throw "basic_phasespace<...> : Q2min must be less than s";
	if(q2max > _M_s) q2max = _M_s;

	_M_q2min = q2min; _M_q2max = q2max;
  }
    
  void basic_phasespace<_Event220>::_S_safety_cut(const event_type& p)
  {
    double s = p[-1]*p[0]; 
    int i, j, up = p.upper();
    for(i = -1; i < up; i++)  
      for(j = (i < 1 ? 1 : i+1); j <= up; j++)
        if(p[i]*p[j] < 1.0e-09*s) throw numeric_error();
  }
  
  
  double basic_phasespace<_Event220>::operator()(event_type& p)
  {
    unsigned int up = p.upper();
    if(up < 0)  throw "unable to generate hadron-hadron -> V+other event";
    
    //----- generate the variable Q2 -----
    const random_generator& __rng = *_M_rng;
    
    double xa, xb, weight = _M_q2max-_M_q2min;
    double Q2 = weight*__rng() + _M_q2min;
    double tau = Q2/_M_s;
    double ltau = -std::log(tau);
    
    //----- incoming hadrons -----
    double Eh = 0.5*std::sqrt(_M_s);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  Eh, Eh);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -Eh, Eh);
    
    //----- outgoing non-qcd particles (in their c.m. frame) -----
    event_type::iterator pbegin = p.begin();
    event_type::iterator pend   = p.end();
    weight *= _M_psgen -> operator()(Q2, pbegin, pbegin+2);
    
    if(up == 0) {
      xa = tau*std::exp(__rng()*ltau);
      xb = tau/xa;
      weight *= ltau/_M_s;
    } else {
      //----- incoming momentum fractions -----
      xa = tau*std::exp(std::sqrt(__rng())*ltau);
      xb = tau*std::exp(__rng()*std::log(xa/tau))/xa;
      weight *= 0.5*xa*xb*ltau*ltau;
      
      //----- generate n+1 massless momenta -----
      weight *= _M_psgen -> operator()(xa*xb*_M_s, pbegin+3, pend);
      double Ecm = std::sqrt(xa*xb*_M_s), E = p[0].T(), xi;
    
      if(up == 1) xi = 1.0-Q2/(Ecm*Ecm);
      else xi = (Ecm-E - std::sqrt(Q2*(1.0-2.0*E/Ecm)+E*E))/(Ecm-2.0*E);
      
      //----- rescale the momenta -----
      for(event_type::iterator pi = pbegin+3; pi < pend; pi++)
		pi -> operator*=(xi);
      
      double Eq = std::sqrt(Q2+xi*xi*E*E);
      weight *= Ecm*E/(Ecm*Eq-Q2)*std::pow(xi, (int) (2*up));
      
      //----- boost outgoing non-qcds to c.m. frame of incomings -----
      threevector<double> bvec(p[0]); bvec /= Eq; 
      p[-3].boost(bvec); p[-2].boost(bvec);
      
      weight /= __DEB_2PI;
    }
    
    //----- incoming partons -----
    p[-1] = xa*p[hadron(-1)];
    p[ 0] = xb*p[hadron( 0)];
    
    //---- boost to the original frame ----
    double bz = (xa - xb)/(xa + xb);
    if(bz != 0.0) {
      for(event_type::iterator pi = pbegin+4; pi < pend; pi++) 
		pi -> boost(0.0, 0.0, bz);
      
      p[-3].boost(0, 0, bz);
      p[-2].boost(0, 0, bz);
    }
    
    return weight;
  }


  double basic_phasespace<_Event220>::operator()(const event_type& q, event_type& p)
  {
    unsigned int up = p.upper(), uq = q.upper();
    if(uq < 0 || up < uq) throw "unable to generate the event";
    
    //----- copy from q to p -----
    event_type::iterator pi;
    event_type::const_iterator qi;
    
    p[hadron(-1)] = q[hadron(-1)];
    p[hadron( 0)] = q[hadron( 0)];
    for(pi = p.begin(), qi = q.begin(); qi < q.end(); pi++, qi++) 
      *pi = *qi; 
    
    if(up == uq) return 1.0;
    
    int emit, omit, thd;
    double weight = 1.0;
    const random_generator& __rng = *_M_rng;

    event_type::const_iterator p1 = p.begin()+4;
    event_type::const_iterator pn = p1+1;
    lorentzvector<double> pab(p[hadron(0)] + p[hadron(-1)]);
    pab /= p[hadron(0)]*p[hadron(-1)];
    
    if(uq == 0) {
      //---- generate the dipole emmissions ----
      __fe_clear_exception();
      
      dipole_emission::boost_ifi boost;
      int emit = -1, omit = (int) (2*__rng()) - 1;
      if(emit == omit) emit = 0;
      
      dipole_emission::gendip_ifi(__rng, _M_beta, _M_eps, pab*p[emit], p[emit], p[1], p[omit], boost);
      p[-2] = boost(p[-2]); p[-3] = boost(p[-3]);
      
      //----- calculate the weight -----
      weight /= dipole_emission::jacobi_ifi(_M_beta, _M_eps, pab*p[0], p[0], p[-1], p1, pn)
	+       dipole_emission::jacobi_ifi(_M_beta, _M_eps, pab*p[-1], p[-1], p[0], p1, pn);
      
      __fe_throw_exception();
      
      weight *= 2.0; uq++; 
      if(up == 1) return weight;
    }

    double jac;
    unsigned int n;
    
    for(n = uq+1, pn = p1+uq+1; n <= up; n++, pn++) {
      //---- choose which one emit, omit, third ----    
      emit = (int) ((n-1)*__rng()) + 1;
      omit = (int) (n*__rng()) - 1;
      thd  = (int) (n*__rng()) + 1;
      if(emit == omit) omit = n-1;
      
      //---- generate the dipole emmission ----
      if(omit > 0) dipole_emission::gendip_fff(__rng, _M_beta, _M_eps, p[emit], p[n], p[omit]);
      else dipole_emission::gendip_ffi(__rng, _M_beta, _M_eps, pab*p[omit], p[emit], p[n], p[omit]);
      
      if(thd != (int) n) std::swap(p[thd], p[n]);
      
      //----- calculate the weight -----
      jac = dipole_emission::jacobi_fff(_M_beta, _M_eps, p1, pn) 
	+   dipole_emission::jacobi_ffi(_M_beta, _M_eps, pab*p[ 0], p[ 0], p1, pn)
	+   dipole_emission::jacobi_ffi(_M_beta, _M_eps, pab*p[-1], p[-1], p1, pn);
      
      weight *= n*n*(n-1)/jac;
      __fe_throw_exception();
    }

    return weight;
  }
}  //   namespace nlo

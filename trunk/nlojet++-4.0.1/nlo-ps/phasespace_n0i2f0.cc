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
#include "bits/psg-phasespace_n0i2f0.h"



#define __eps__  0.01


namespace nlo {

  typedef hadronic_event<lorentzvector<double>, hadronic_event_traits<0U,2U,0U> > _Event020;



  double basic_phasespace<_Event020>::_S_fun_eta(double r1, double r2, double eps) 
  {
    double reps = std::sqrt(eps), leps = -std::log(eps);
    double N = leps + reps;
    
    return r2 < reps/N ? eps*r1*r1 : eps*std::exp(r1*leps);
  }
  
  double basic_phasespace<_Event020>::_S_jac_eta(double eta, double eps) {
    return (eta < eps ? 0.5/std::sqrt(eta) : 1.0/eta)/(std::sqrt(eps) - std::log(eps));
  }

  void basic_phasespace<_Event020>::_S_safety_cut(const event_type& p)
  {
    double s = p[-1]*p[0]; 
    int i, j, up = p.upper();
    for(i = -1; i < up; i++)  
      for(j = (i < 1 ? 1 : i+1); j <= up; j++)
        if(p[i]*p[j] < 1.0e-12*s) throw numeric_error();
  }


  double basic_phasespace<_Event020>::operator()(event_type& p)
  {
    unsigned int up = p.upper();
    if(up < 2)  throw "unable to generate hadron-hadron event";
  
    //----- incoming hadrons -----
    double E = 0.5*std::sqrt(_M_s);
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  E, E);
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -E, E);
    
    //----- incoming partons -----
    const random_generator& __rng = *_M_rng;

    double xa = _S_fun_eta(__rng(), __rng(), __eps__);
    double xb = _S_fun_eta(__rng(), __rng(), __eps__);
    double weight = 1.0/(_S_jac_eta(xa, __eps__)*_S_jac_eta(xb, __eps__));
    
    p[-1] = xa*p[hadron(-1)];
    p[ 0] = xb*p[hadron( 0)];
    
    //----- outgoing partons -----
    event_type::iterator pbegin = p.begin() + 2;
    event_type::iterator pend   = p.end();
    
    weight *= _M_psgen -> operator()(xa*xb*_M_s, pbegin, pend);
    
    double bz = (xa - xb)/(xa + xb);
    if(bz != 0.0)
      for(event_type::iterator pi = pbegin; pi < pend; pi++) 
        pi -> boost(0.0, 0.0, bz);
    
    return weight;
  }

  double basic_phasespace<_Event020>::
  operator()(const event_type& q, event_type& p)
  {
    unsigned int up = p.upper(), uq = q.upper();
    if(uq < 2 || up < uq)  throw "unable to generate hadron-hadron event";
    
    //----- copy from q to p -----
    event_type::iterator pi;
    event_type::const_iterator qi;
    
    p[hadron(-1)] = q[hadron(-1)];
    p[hadron( 0)] = q[hadron( 0)];
    for(pi = p.begin(), qi = q.begin(); qi < q.end(); pi++, qi++) 
      *pi = *qi; 

    if(up == uq) return 1.0;
    
    //----- dipole emmisions -----
    __fe_clear_exception();
    
    unsigned int n;
    int emit, omit, thd;
    double jac, weight = 1.0;
    const random_generator& __rng = *_M_rng;

    event_type::const_iterator p1 = p.begin()+2, pn;
    
    lorentzvector<double> pab(p[hadron(0)] + p[hadron(-1)]);
    pab /= p[hadron(0)]*p[hadron(-1)];
    
    for(n = uq+1, pn = p1 + uq+1; n <= up; n++, pn++) {
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
 
    //----- safety cuts -----
    _S_safety_cut(p);
    
    return weight;
  }
}   //  namespace nlo

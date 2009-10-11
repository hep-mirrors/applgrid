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
#include "bits/psg-phasespace_n0i0f0.h"


namespace nlo {
 
  typedef hadronic_event<lorentzvector<double>, hadronic_event_traits<0U,0U,0U> > _Event000;
  
  
  double basic_phasespace<_Event000>::operator()(_Event& p)
  {
    if(p.upper() < 2) throw("unable to generate e+e- event");
    
    //----- incoming leptons -----
    double E = 0.5*std::sqrt(_M_s);
    p[-1] = lorentzvector<double>(0.0, 0.0,  E, E);
    p[0]  = lorentzvector<double>(0.0, 0.0, -E, E);
    
    //----- outgoing partons -----
    return _M_psgen -> operator()(_M_s, p.begin() + 2, p.end());
  }
  
  double basic_phasespace<_Event000>::
  operator()(const event_type& q, event_type& p)
  {
    unsigned int up = p.upper(), uq = q.upper();
    if(uq < 2) throw("unable to generate e+e- event");
    
    //----- copy from q to p -----
    event_type::iterator pi;
    event_type::const_iterator qi;
    
    for(pi = p.begin(), qi = q.begin(); qi < q.end(); pi++, qi++) 
      *pi = *qi; 

    if(up == uq) return 1.0;
 
    //----- dipole emmisions -----
    __fe_clear_exception();    

    int emit, omit, thd;
    double weight = 1.0;
    const random_generator& __rng = *_M_rng;

    event_type::const_iterator p1 = p.begin()+2;
    
    for(unsigned int n = uq+1; n <= up; n++) {
      //---- choose which one emit, omit, third ----    
      emit = (int) ((n-1)*__rng()) + 1;
      omit = (int) ((n-2)*__rng()) + 1;
      thd  = (int) (n*__rng()) + 1;
      if(emit == omit) omit = n-1;
      
      //---- generate the dipole emmission ----
      dipole_emission::gendip_fff(__rng, _M_beta, _M_eps, p[emit], p[n], p[omit]);
      if(thd != (int) n) std::swap(p[thd], p[n]);
      
      weight *= n*(n-1)*(n-2)/dipole_emission::jacobi_fff(_M_beta, _M_eps, p1, p1+n);
      __fe_throw_exception();
    }
    
    return weight;
  }
 
}   //  namespace nlo

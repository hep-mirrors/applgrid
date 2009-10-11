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
#ifndef __NLO_PSG_PHASESPACE_N0I2F0_H__
#define __NLO_PSG_PHASESPACE_N0I2F0_H__ 1

#include <bits/nlo-event.h>
#include <bits/nlo-phasespace.h>

#include <bits/psg-phasespace.h>
#include <bits/psg-psgen.h>
#include <bits/psg-rambo.h>
#include <bits/psg-dipole.h>



namespace nlo {
  

  template<>
  class basic_phasespace<hadronic_event<lorentzvector<double>, hadronic_event_traits<0U,2U,0U> > >
    : public phasespace<hadronic_event<lorentzvector<double>, hadronic_event_traits<0U,2U,0U> > >
  {
    //  private types
    //   abstract base type
    typedef phasespace<hadronic_event<lorentzvector<double>, hadronic_event_traits<0U,2U,0U> > > _Base;
    
    //   event_type related type aliases 
    typedef _Base::event_type _Event;
    typedef event_type::iterator _Iterator;
    
    //   phase space generator related aliases
    typedef phasespace_generator<lorentzvector<double>, _Iterator> _PhaseSpaceGen;
    typedef rambo<_Iterator> _Rambo;
    
  public:
    typedef _Base::event_type event_type;

    //   constructors
    basic_phasespace(const random_generator *rng, double s, 
		     double beta = 2.0, double eps = 1.0)
      : _Base(rng), _M_s(s), _M_beta(beta), _M_eps(eps), 
	_M_intpsgen(true), _M_psgen(new _Rambo(rng)) {}
    
    basic_phasespace(const random_generator *rng, _PhaseSpaceGen *psg, 
		     double s, double beta = 2.0, double eps = 1.0)
      : _Base(rng), _M_s(s), _M_beta(beta), _M_eps(eps), 
	_M_intpsgen(false), _M_psgen(psg) {}
    
    ~basic_phasespace() {
      if(_M_intpsgen) delete _M_psgen;
    }
    
    double operator()(_Event& p);
    double operator()(const _Event& , _Event&);
        
  private:
    double _M_s, _M_beta, _M_eps;
    bool _M_intpsgen;
    _PhaseSpaceGen *_M_psgen;

    //   private static members
    static double _S_fun_eta(double, double, double);
    static double _S_jac_eta(double, double);
    static void   _S_safety_cut(const event_type&);
  };  
}

#endif

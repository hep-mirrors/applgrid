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
#include "epa5jet.h"
#include "ampq2g3l2.h"
#include "ampq4g1l2.h"


#define __PI_FAC5 (6185560.53048687015423266253*12.0*_M_ip.s(-1,0))


namespace nlo {


  epa5jet::epa5jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_epa(5, 3, al), _epa_jet_base(nu, nd), _M_mchel(mchel)
  {
    _M_q2g3 = new ampq2g3l2(_M_ip, rng);
    _M_q4g1 = new ampq4g1l2(_M_ip, rng);
  }
  
  epa5jet::~epa5jet() 
  {
    if(_M_q2g3) delete _M_q2g3;
    if(_M_q4g1) delete _M_q4g1;
  }

  void epa5jet::born_term(const event_type& p, weight_type& res) {
    _M_ip.calculate(p);
    res = __PI_FAC5*amp_tree(_M_q2g3, _M_q4g1);
  }
  
  void epa5jet::real_term(const event_type&, weight_type&){
    throw("Five jet NLO calculation is not implemented!");
  }
  

  void epa5jet::fini_term(const event_type&, weight_type&) {
    throw("Five jet NLO calculation is not implemented!");
  }
  
  void epa5jet::dipole_term(const event_type&, const event_type&,
			   int, int, int, weight_type&) {
    throw("Five jet NLO calculation is not implemented!");
  }
}   //  namespace nlo

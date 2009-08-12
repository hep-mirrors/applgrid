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
#include "dis4jet.h"
#include "ampq2g3l2.h"
#include "ampq4g1l2.h"


#define PI_FAC4 (4.0*77730046.08365014226587406347)


namespace nlo {
  

  dis4jet::dis4jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_dis(4, 3, nu, nd, al), _dis_jet_base(nu, nd), _M_mchel(mchel) 
  {
    _M_q2g3 = new ampq2g3l2(_M_ip, rng);
    _M_q4g1 = new ampq4g1l2(_M_ip, rng);
  }
  
  dis4jet::~dis4jet() 
  {
    if(_M_q2g3) delete _M_q2g3;
    if(_M_q4g1) delete _M_q4g1;
  }
    
  void dis4jet::born_term(const event_type& p, weight_type& res) 
  {
    double amp[3];
    _M_ip.calculate(p);
    
    if(_M_mchel) amp_tree_mch(_M_q2g3, _M_q4g1, amp);
    else amp_tree(_M_q2g3, _M_q4g1, amp);
    
    res[0] = PI_FAC4*amp[0];
    res[1] = PI_FAC4*amp[1];
    res[2] = PI_FAC4*amp[2];
  }
  
  void dis4jet::real_term(const event_type&, weight_type&) {
    throw("Four jet NLO calculation is not implemented");
  }

  void dis4jet::fini_term(double, double, const event_type&, weight_type *) 
  { throw("Four jet NLO calculation is not implemented");}
  
  void dis4jet::dipole_term(const event_type&, const event_type&,
			    int, int, int, weight_type&) 
  { throw("Four jet NLO calculation is not implemented");}
}      //   namespace nlo

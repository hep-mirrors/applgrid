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
#include "hhc2ph2jet.h"
#include "ampq2g2p2.h"
#include "ampq4p2.h"


//    PI_FACn = (2pi)^(2n+2) * 2^(n+1)
#define PI_FAC2 492231.26711055587177578917



namespace nlo {

 
 
  hhc2ph2jet::hhc2ph2jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_hhc2ph(2, 2, nu, nd, al), _hhc2ph_jet_base(nu, nd), _M_mchel(mchel) 
  {
    _M_q2g2p2 = new ampq2g2p2(_M_ip, rng);
    _M_q4p2   = new ampq4p2  (_M_ip, rng);
  }
  
  hhc2ph2jet::~hhc2ph2jet() 
  {
    if(_M_q2g2p2) delete _M_q2g2p2;
    if(_M_q4p2)   delete _M_q4p2;
  }
  
  void hhc2ph2jet::born_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    if(_M_mchel) amp_tree_mch(_M_q2g2p2, _M_q4p2, res.begin());
    else amp_tree(_M_q2g2p2, _M_q4p2, res.begin());
    res *= PI_FAC2;
  }
  
  void hhc2ph2jet::real_term(const event_type&, weight_type&) 
  {}
  
  void hhc2ph2jet::fini_term(double, double, double, double, const event_type&, weight_type *) 
  {}
  
  void hhc2ph2jet::dipole_term(const event_type&, const event_type&, int, int, int, weight_type&)
  {} 
}  //  namespace nlo

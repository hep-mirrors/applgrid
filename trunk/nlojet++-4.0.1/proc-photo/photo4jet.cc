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
#include "photo4jet.h"
#include "ampq2g3p1.h"
#include "ampq4g1p1.h"


//   PI_FACn = (2pi)^(2n-3) * 2^(n-1)
#define PI_FAC5 6185560.53048687015450831204



namespace nlo {


  photo4jet::photo4jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_photo(4, 3, nu, nd, al), _photo_jet_base(nu, nd), _M_mchel(mchel) 
  {
    _M_q2g3p1 = new ampq2g3p1(_M_ip, rng);
    _M_q4g1p1 = new ampq4g1p1(_M_ip, rng);
  }
  
  photo4jet::~photo4jet() 
  {
    if(_M_q2g3p1) delete _M_q2g3p1;
    if(_M_q4g1p1) delete _M_q4g1p1;
  }
  
  void photo4jet::born_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_q2g3p1, _M_q4g1p1, res.begin());
    res *= PI_FAC5;
  }
  
  void photo4jet::real_term(const event_type&, weight_type&) {}
  void photo4jet::fini_term(double, double, double, double, const event_type&, weight_type *){}
  void photo4jet::dipole_term(const event_type&, const event_type&, int, int, int, weight_type&) {}

  bool photo4jet::dipole_index(int i, int j, int k) {
    if(k > -1) return true;
    else return false;
  }
}

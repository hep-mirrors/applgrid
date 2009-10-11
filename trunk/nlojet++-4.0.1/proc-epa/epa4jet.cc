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
#include "epa4jet.h"
#include "ampq2g2l2.h"
#include "ampq2g3l2.h"
#include "ampq4l2.h"
#include "ampq4g1l2.h"


#define __PI_FAC4 (78341.03930503205203493218*12.0*_M_ip.s(-1,0))
#define __PI_FAC5 (6185560.53048687015423266253*12.0*_M_ip.s(-1,0))


namespace nlo {


  epa4jet::dipole_func epa4jet::_S_dipole[10] = 
  { &epa4jet::_M_d12, &epa4jet::_M_d13, &epa4jet::_M_d14, &epa4jet::_M_d15,
    &epa4jet::_M_d23, &epa4jet::_M_d13, &epa4jet::_M_d15, &epa4jet::_M_d34,
    &epa4jet::_M_d35, &epa4jet::_M_d35};


  epa4jet::epa4jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_epa(4, 2, al), _epa_jet_base(nu, nd), _M_mchel(mchel)
  {
    _M_q2g2 = new ampq2g2l2(_M_ip, rng);
    _M_q4   = new   ampq4l2(_M_ip, rng);
    _M_q2g3 = new ampq2g3l2(_M_ip, rng);
    _M_q4g1 = new ampq4g1l2(_M_ip, rng);
  }

  epa4jet::~epa4jet() 
  {
    if(_M_q4)   delete _M_q4;
    if(_M_q2g2) delete _M_q2g2;
    if(_M_q2g3) delete _M_q2g3;
    if(_M_q4g1) delete _M_q4g1;
  }

  void epa4jet::born_term(const event_type& p, weight_type& res) {
    _M_ip.calculate(p);
    if(_M_mchel) res = __PI_FAC4*amp_tree(_M_q2g2, _M_q4);
    else res = __PI_FAC4*amp_tree_mch(_M_q2g2, _M_q4);
  }
  
  void epa4jet::real_term(const event_type& p, weight_type& res) {
    _M_ip.calculate(p);
    res = __PI_FAC5*amp_tree(_M_q2g3, _M_q4g1);
  }
  
  void epa4jet::fini_term(const event_type& p, weight_type& res) {
    _M_ip.calculate(p);
    if(_M_mchel) res = __PI_FAC4*amp_1loop_mch(_M_q2g2, _M_q4, alpha());
    else res = __PI_FAC4*amp_1loop(_M_q2g2, _M_q4, alpha());
  }
 
  
  void epa4jet::dipole_term(const event_type& p, const event_type& dp, 
			   int i, int j, int k, weight_type& res) 
  {
    int kt  = (k == 5 ? j : k);
    int idx = 4*i-(i*i-i)/2 + j-6;
    
    _M_split.set(p[i], p[j], p[k]);
    _M_ip.calculate(dp);
    
    res = __PI_FAC5*((this->*_S_dipole[idx])(kt, i));
  }

#define __DEB_EE4_QA _M_split.Vqa()
#define __DEB_EE4_QG _M_split.Vqg()
#define __DEB_EE4_GG _M_split.Vgg()

#define __DEB_EE4_Q2G2(p1,p2,p3,p4) amp_cc(_M_q2g2, kt, i, p1,p2,p3,p4)
#define __DEB_EE4_Q4(p1,p2,p3,p4)   amp_cc(_M_q4,   kt, i, p1,p2,p3,p4)

  double epa4jet::_M_d12(int kt, int i) {
    return _M_nf*_M_Dijk(__DEB_EE4_QA, __DEB_EE4_Q2G2(3,4,1,2))/4.0;
  }
  
  double epa4jet::_M_d13(int kt, int i) {
    return _M_Dijk(__DEB_EE4_QG, __DEB_EE4_Q2G2(1,2,3,4))/6.0;
  }
  
  double epa4jet::_M_d14(int kt, int i) {
    return _M_Dijk(__DEB_EE4_QG, __DEB_EE4_Q2G2(1,2,3,4))/6.0 
      + _M_nf*_M_Dijk(__DEB_EE4_QA, __DEB_EE4_Q2G2(3,2,1,4))/4.0;
  }
  
  double epa4jet::_M_d15(int kt, int i) {
    std::pair<double, std::complex<double> > sp(__DEB_EE4_QG); 
    return _M_Dijk(sp, __DEB_EE4_Q2G2(1,2,3,4))/6.0
      + _M_Dijk(sp, __DEB_EE4_Q4(1,2,3,4))/4.0;
  }
  
  double epa4jet::_M_d23(int kt, int i) {
    return _M_Dijk(__DEB_EE4_QG, __DEB_EE4_Q2G2(1,2,3,4))/6.0 
      + _M_nf*_M_Dijk(__DEB_EE4_QA, __DEB_EE4_Q2G2(1,4,2,3))/4.0;
  }
  
  double epa4jet::_M_d34(int kt, int i) {
    std::pair<double, std::complex<double> > cc(__DEB_EE4_Q2G2(1,2,3,4));
    return _M_Dijk(__DEB_EE4_GG, cc)/6.0 + _M_nf*_M_Dijk(__DEB_EE4_QA, cc)/4.0;
  }
  
  double epa4jet::_M_d35(int kt, int i) {
    return _M_Dijk(__DEB_EE4_GG, __DEB_EE4_Q2G2(1,2,3,4))/6.0
      + _M_Dijk(__DEB_EE4_QG, __DEB_EE4_Q4(1,2,3,4))/4.0;
  }
}   // namespace nlo

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
#include "epa3jet.h"
#include "ampq2g1l2.h"
#include "ampq2g2l2.h"
#include "ampq4l2.h"

#define __PI_FAC3 (992.20085376959424561273*12.0*_M_ip.s(-1,0))
#define __PI_FAC4 (78341.03930503205203493218*12.0*_M_ip.s(-1,0))



namespace nlo {
  

  epa3jet::dipole_func epa3jet::_S_dipole[6] = 
  { &epa3jet::_M_d12, &epa3jet::_M_d13, &epa3jet::_M_d14, 
    &epa3jet::_M_d23, &epa3jet::_M_d13, &epa3jet::_M_d34};

  
  epa3jet::epa3jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_epa(3, 1, al), _epa_jet_base(nu, nd), _M_mchel(mchel)
  {
    _M_q2g1 = new ampq2g1l2(_M_ip, rng);
    _M_q2g2 = new ampq2g2l2(_M_ip, rng);
    _M_q4   = new   ampq4l2(_M_ip, rng);
  }
  
  epa3jet::~epa3jet() 
  {
    if(_M_q2g1) delete _M_q2g1;
    if(_M_q2g2) delete _M_q2g2;
    if(_M_q4)   delete _M_q4;
  }
  
  void epa3jet::born_term(const event_type& p, weight_type& res) {
    _M_ip.calculate(p);
    if(_M_mchel) res = __PI_FAC3*amp_tree_mch(_M_q2g1);
    else res = __PI_FAC3*amp_tree(_M_q2g1);
  }
  
  void epa3jet::real_term(const event_type& p, weight_type& res) {
    _M_ip.calculate(p);
    res = __PI_FAC4*amp_tree(_M_q2g2, _M_q4);
  }
  
  void epa3jet::fini_term(const event_type& p, weight_type& res) {
    _M_ip.calculate(p);
    if(_M_mchel) res = __PI_FAC3*amp_1loop_mch(_M_q2g1, alpha());
    else res = __PI_FAC3*amp_1loop(_M_q2g1, alpha());
  }
  
  void epa3jet::
  dipole_term(const event_type& p, const event_type& dp, 
	      int i,int j, int k, weight_type& res) 
  {
    int kt  = (k == 4 ? j : k);
    int idx = 3*i-(i*i-i)/2 + j-5;
        
    _M_split.set(p[i], p[j], p[k]);
    _M_ip.calculate(dp);
    
    res = __PI_FAC4*((this->*_S_dipole[idx])(kt, i));
  }
  
#define __DEB_EE3_QA _M_split.Vqa()
#define __DEB_EE3_QG _M_split.Vqg()
#define __DEB_EE3_GG _M_split.Vgg()
#define __DEB_EE3_Q2G1(p1,p2,p3) amp_cc(_M_q2g1, kt,i, p1,p2,p3)
  
  double epa3jet::_M_d12(int kt, int i) {
    return _M_nf*_M_Dijk(__DEB_EE3_QA, __DEB_EE3_Q2G1(3,2,1))/4.0;
  }
  
  double epa3jet::_M_d13(int kt, int i) {
    return _M_Dijk(__DEB_EE3_QG, __DEB_EE3_Q2G1(1,2,3))/2.0;
  }
  
  double epa3jet::_M_d14(int kt, int i) {
    return _M_Dijk(__DEB_EE3_QG, __DEB_EE3_Q2G1(1,2,3))/2.0 
      + _M_nf*_M_Dijk(__DEB_EE3_QA, __DEB_EE3_Q2G1(3,2,1))/4.0;
  }
  
  double epa3jet::_M_d23(int kt, int i) {
    return _M_Dijk(__DEB_EE3_QG, __DEB_EE3_Q2G1(1,2,3))/2.0
      + _M_nf*_M_Dijk(__DEB_EE3_QA, __DEB_EE3_Q2G1(1,3,2))/4.0;
  }
  
  double epa3jet::_M_d34(int kt, int i) {
    std::pair<double, std::complex<double> > cc(__DEB_EE3_Q2G1(1,2,3));
    return _M_Dijk(__DEB_EE3_GG, cc)/2.0 + _M_nf*_M_Dijk(__DEB_EE3_QA, cc)/4.0;
  }
}  //  namespace nlo

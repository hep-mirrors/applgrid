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
#include "bits/nlo-color.h"
#include "bits/proc-epa.h"

#include "ampq2g1l2.h"
#include "ampq2g2l2.h"
#include "ampq2g3l2.h"
#include "ampq4l2.h"
#include "ampq4g1l2.h"


namespace nlo {
  
  _epa_jet_base::_epa_jet_base(unsigned int __nu, unsigned int __nd) 
  {
    double __schr = 2.0*__nu - __nd;
    _M_nf = __nu + __nd;
    _M_xcharge = __schr*__schr/(4.0*__nu+__nd);
  }
  
  //
  //    Three parton contributions
  //
  double _epa_jet_base::amp_tree(ampq2g1l2 *amp) {
    return (amp -> su3_tree(1,3,2,0,-1))/Nc;
  }

  double _epa_jet_base::amp_tree_mch(ampq2g1l2 *amp) {
    return (amp -> su3_tree_mch(1,3,2,0,-1))/Nc;
  }
  
  double _epa_jet_base::amp_1loop(ampq2g1l2 *amp, double al) 
  {
    return (amp -> su3_1loop(_M_nf, 1,3,2,0,-1) 
	    + amp -> su3_ins(_M_nf, 1,3,2,0,-1, al))/Nc;
  }
  
  double _epa_jet_base::amp_1loop_mch(ampq2g1l2 *amp, double al) 
  {
    return (amp -> su3_1loop_mch(_M_nf, 1,3,2,0,-1)
	    + amp -> su3_ins_mch(_M_nf, 1,3,2,0,-1, al))/Nc;
  }
  
  std::pair<double, std::complex<double> > _epa_jet_base::
  amp_cc(ampq2g1l2 *amp, int i, int j, int p1, int p2, int p3) 
  {
    _Pair res(amp -> su3_cc(i,j, p1,p3,p2, 0,-1));
    res.first  /= Nc;
    res.second /= Nc;
    return res;
  }
  
  //
  //    Four parton contributions
  //
  double _epa_jet_base::amp_tree(ampq2g2l2 *amp1, ampq4l2 *amp2) 
  {
    double amp[10];
    double ret = 2.0*(amp1 -> su3_tree(1,3,4,2, 0,-1));
    amp2 -> su3_tree(1,4,3,2, 0,-1, amp);
    
    ret += _M_nf*(amp[0] + amp[1] + amp[3] + amp[4]);
    ret += amp[6] + amp[7] + amp[8] + amp[9];
    
    return 0.25*(ret + _M_xcharge*(amp[2] + amp[5]))/Nc;
  }

  double _epa_jet_base::amp_tree_mch(ampq2g2l2 *amp1, ampq4l2 *amp2) 
  {
    double amp[10];
    double ret = 2.0*(amp1 -> su3_tree_mch(1,3,4,2, 0,-1));
    amp2 -> su3_tree_mch(1,4,3,2, 0,-1, amp);
    
    ret += _M_nf*(amp[0] + amp[1] + amp[3] + amp[4]);
    ret += amp[6] + amp[7] + amp[8] + amp[9];
    
    return 0.25*(ret + _M_xcharge*(amp[2] + amp[5]))/Nc;
  }
  
  double _epa_jet_base::
  amp_1loop(ampq2g2l2 *amp1, ampq4l2 *amp2, double al) 
  {
    double ret, amp[10], ins[10];
    
    ret = 2.0*(amp1->su3_1loop(_M_nf, 1,3,4,2,0,-1)
	       + amp1->su3_ins(_M_nf, 1,3,4,2,0,-1, al));
    
    amp2 -> su3_1loop(_M_nf, 1,4,3,2, 0,-1, amp);
    amp2 -> su3_ins(1,4,3,2, 0,-1, ins, al);

    for(unsigned i = 0; i < 10; i++)
      amp[i] += ins[i];
    
    ret += _M_nf*(amp[0] + amp[1] + amp[3] + amp[4]);
    ret += amp[6] + amp[7] + amp[8] + amp[9];
    
    return 0.25*(ret + _M_xcharge*(amp[2] + amp[5]))/Nc;
  }
  
  double _epa_jet_base::
  amp_1loop_mch(ampq2g2l2 *amp1, ampq4l2 *amp2, double al) 
  {
    double ret, amp[10], ins[10];
    
    ret = 2.0*(amp1->su3_1loop_mch(_M_nf, 1,3,4,2,0,-1)
	       + amp1->su3_ins_mch(_M_nf, 1,3,4,2,0,-1, al));
    
    amp2 -> su3_1loop_mch(_M_nf, 1,4,3,2, 0,-1, amp);
    amp2 -> su3_ins_mch(1,4,3,2, 0,-1, ins, al);
    
    for(unsigned i = 0; i < 10; i++)
      amp[i] += ins[i];
    
    ret += _M_nf*(amp[0] + amp[1] + amp[3] + amp[4]);
    ret += amp[6] + amp[7] + amp[8] + amp[9];
    
    return 0.25*(ret + _M_xcharge*(amp[2] + amp[5]))/Nc;
  }
  
  std::pair<double, std::complex<double> > _epa_jet_base::
  amp_cc(ampq2g2l2 *amp, int i, int j, int p1, int p2, int p3, int p4) 
  {
    _Pair res(amp -> su3_cc(i,j, p1,p3,p4,p2, 0,-1));
    res.first  /= Nc;
    res.second /= Nc;
    return res;
  }
  
  std::pair<double, std::complex<double> > _epa_jet_base::
  amp_cc(ampq4l2 *amp2, int i, int j, int p1, int p2, int p3, int p4) 
  {
    double amp[10], ret;
    amp2 -> su3_cc(i,j, p1,p4,p3,p2, 0, -1, amp);
    
    ret = _M_nf*(amp[0] + amp[1] + amp[3] + amp[4]);
    ret += amp[6] + amp[7] + amp[8] + amp[9];
    
    double cc = (ret + _M_xcharge*(amp[2]+amp[5]))/Nc;
    
    return _Pair(cc, 0.0);
  }
  
  //
  //    Five parton contributions
  //
  double _epa_jet_base::amp_tree(ampq2g3l2 *amp1, ampq4g1l2 *amp2) 
  {
    double amp[10];
    double ret = 2.0*(amp1 -> su3_tree(1,3,4,5,2, 0,-1))/3.0;
    amp2 -> su3_tree(1,4,3,2,5,0,-1, amp);
    
    ret += _M_nf*(amp[0] + amp[1] + amp[3] + amp[4]);
    ret += amp[6] + amp[7] + amp[8] + amp[9];

    return 0.25*(ret + _M_xcharge*(amp[2] + amp[5]))/Nc;
  }
 
  double _epa_jet_base::amp_tree_mch(ampq2g3l2 *amp1, ampq4g1l2 *amp2) 
  {
    double amp[10];
    double ret = 2.0*(amp1 -> su3_tree_mch(1,3,4,5,2, 0,-1))/3.0;
    amp2 -> su3_tree_mch(1,4,3,2,5,0,-1, amp);
    
    ret += _M_nf*(amp[0] + amp[1] + amp[3] + amp[4]);
    ret += amp[6] + amp[7] + amp[8] + amp[9];

    return 0.25*(ret + _M_xcharge*(amp[2] + amp[5]))/Nc;
  }
}    //  namespace nlo

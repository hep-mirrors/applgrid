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
#ifndef __NLO_PROC_HHC2PH_H__
#define __NLO_PROC_HHC2PH_H__ 1


#include <complex>
#include <utility>
#include <bits/hhc2ph-weight.h>


namespace nlo {
  
  //    Foreward declarations
  class ampq2g1p2;
  class ampq2g2p2;
  class ampq4p2;
  struct su3_kp_i2;


  class _hhc2ph_jet_base
  {
  protected:
    //   constructor
    explicit _hhc2ph_jet_base(unsigned int nu=2U, unsigned int nd=3U) 
      : Nu(nu), Nd(nd), Nf(nu+nd) {}
    
    //  types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;
    
    //    three parton contributions    
    void amp_tree(ampq2g1p2 *, double *);
    void amp_1loop(ampq2g1p2 *, double *);
    void amp_1loop_mch(ampq2g1p2 *, double *);
    void amp_kp(double, ampq2g1p2 *, su3_kp_i2 *);

    _Pair amp_ccag(ampq2g1p2 *, int, int);
    _Pair amp_ccga(ampq2g1p2 *, int, int);
    _Pair amp_ccqa(ampq2g1p2 *, int, int);

    //    four parton contributions  
    void amp_tree(ampq2g2p2 *, ampq4p2 *, double *);
    void amp_tree_mch(ampq2g2p2 *, ampq4p2 *, double *);

    void convolutions(double e1, double x1, double xj1, double e2, double x2,
		      double xj2, double al, const su3_kp_i2 *kp, weight_hhc2ph *out)
    {
      __conv_x1(e1, x1, xj1, al, kp, out);
      __conv_x2(e2, x2, xj2, al, kp, out);
    }

    //   data members
    unsigned int Nu, Nd, Nf;
    
  private:
    //    Convolution
    void __conv_x1(double, double, double, double, const su3_kp_i2 *, weight_hhc2ph *);
    void __conv_x2(double, double, double, double, const su3_kp_i2 *, weight_hhc2ph *);  
  };
}

#endif


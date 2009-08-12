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
#ifndef __NLO_PROC_PHOTO_H__
#define __NLO_PROC_PHOTO_H__ 1


#include <complex>
#include <utility>

#include <bits/proc-hhc.h>
#include <bits/photo-weight.h>


namespace nlo {
  
  //    Foreward declarations
  class ampq2g1p1;
  class ampq2g2p1;
  class ampq2g3p1;
  class ampq4p1;
  class ampq4g1p1;

  struct su3_kp_i1;
  struct su3_kp_i2;

  class _photo_jet_base : protected _hhc_jet_base
  {
  protected:
    //   constructor
    explicit _photo_jet_base(unsigned int nu=2U, unsigned int nd=3U) 
      : _hhc_jet_base(nu+nd), Nu(nu), Nd(nd), Nf(nu+nd) {}
    
    //  types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

    //    three parton contributions    
    void amp_tree(ampq2g1p1 *, double *);
    void amp_1loop(ampq2g1p1 *, double *);
    void amp_1loop_mch(ampq2g1p1 *, double *);
    double amp_ccg(ampq2g1p1 *, int, int, int, int);
    double amp_ccq(ampq2g1p1 *, int, int, int, int);
    void amp_kp(double, ampq2g1p1 *, su3_kp_i1 *);
    void amp_kp(double, ampq2g2 *, ampq4 *, su3_kp_i2 *);

    //    four parton contributions    
    void amp_tree(ampq2g2p1 *, ampq4p1 *, double *);
    void amp_1loop(ampq2g2p1 *, ampq4p1 *, double *);
    void amp_1loop_mch(ampq2g2p1 *, ampq4p1 *, double *);
    void amp_kp(double, ampq2g2p1 *, ampq4p1 *, su3_kp_i1 *);
    void amp_kp(double, ampq2g3 *, ampq4g1 *, su3_kp_i2 *);

    _Pair amp_ccg(ampq2g2p1 *, int, int, int, int, int); 
    _Pair amp_ccq(ampq2g2p1 *, int, int, int, int, int); 
    void amp_ccq(ampq4p1 *, double, int, int, int, int, int, double *);
    double amp_ccq(ampq4p1 *, double, double, int, int, int, int, int);

    //    five parton contributions  
    void amp_tree(ampq2g3p1 *, ampq4g1p1 *, double *);
    void amp_tree_mch(ampq2g3p1 *, ampq4g1p1 *, double *);

    //    convolutions
    void conv_parton(double, double, double, double, const su3_kp_i1 *, weight_photo *);
    void conv_photon(double, double, double, double, const su3_kp_i2 *, weight_photo *);
    
    //   data members
    unsigned int Nu, Nd, Nf;
  };
}

#endif


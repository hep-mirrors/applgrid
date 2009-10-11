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
#ifndef __NLO_PROC_HHC_H__
#define __NLO_PROC_HHC_H__ 1


#include <complex>
#include <utility>
#include <bits/hhc-weight.h>


namespace nlo {
  
  //    Foreward declarations
  class ampg4;
  class ampg5;
  class ampg6;

  class ampq2g2;
  class ampq2g3;
  class ampq2g4;

  class ampq4;
  class ampq4g1;
  class ampq4g2;

  class ampq6;
  struct su3_kp_i2;

  class _hhc_jet_base
  {
  protected:
    //   constructor
    explicit _hhc_jet_base(unsigned int nf=5U) 
      : Nf(nf) {}
    
    //  types
    typedef std::pair<double, std::complex<double> > _Pair;

    //    Four parton contributions    
    void amp_tree(ampg4 *, ampq2g2 *, ampq4 *, double *, double *);
    void amp_1loop(ampg4 *, ampq2g2 *, ampq4 *, double *);
    void amp_1loop_mch(ampg4 *, ampq2g2 *, ampq4 *, double *);
    void amp_kp(double, ampg4 *, ampq2g2 *, ampq4 *, su3_kp_i2 *);

    void amp_ccgg(ampg4 *, int, int, int, int, double *);
    void amp_ccgg(ampq2g2 *, int, int, int, int, double *);
    void amp_ccqg(ampq2g2 *, int, int, int, int, double *);
    void amp_ccgq(ampq2g2 *, int, int, int, int, double *);
    void amp_ccag(ampq2g2 *, int, int, int, int, double *);
    void amp_ccga(ampq2g2 *, int, int, int, int, double *);
    void amp_ccqa(ampq2g2 *, int, int, int, int, double *);
    void amp_ccaq(ampq2g2 *, int, int, int, int, double *);
    void amp_cc(ampq4 *, int, int, int, int, int, int, double *);

    //    Five parton contributions    
    void amp_tree (ampg5 *, ampq2g3 *, ampq4g1 *, double *, double *);
    void amp_kp(double, ampg5 *, ampq2g3 *, ampq4g1 *, su3_kp_i2 *);

    void amp_ccgg(ampg5 *, int, int, int, int, int, _Pair *); 
    void amp_ccgg(ampq2g3 *, int, int, int, int, int, _Pair *); 
    void amp_ccqg(ampq2g3 *, int, int, int, int, int, _Pair *);
    void amp_ccgq(ampq2g3 *, int, int, int, int, int, _Pair *);
    void amp_ccag(ampq2g3 *, int, int, int, int, int, _Pair *);
    void amp_ccga(ampq2g3 *, int, int, int, int, int, _Pair *);
    void amp_ccqa(ampq2g3 *, int, int, int, int, int, _Pair *);
    void amp_ccaq(ampq2g3 *, int, int, int, int, int, _Pair *);
    void amp_cc(ampq4g1 *, int, int, int, int, int, int, int, _Pair *);

    void amp_1loop(ampg5 *, ampq2g3 *, ampq4g1 *, double *);
    void amp_1loop_mch(ampg5 *, ampq2g3 *, ampq4g1 *, double *);

    //    Six parton contributions
    void amp_tree(ampg6 *, ampq2g4 *, ampq4g2 *, ampq6 *, double *, double *);
    
    void convolutions(double e1, double x1, double xj1, double e2, double x2,
		      double xj2, double al, const su3_kp_i2 *kp, weight_hhc *out)
    {
      __conv_x1(e1, x1, xj1, al, kp, out);
      __conv_x2(e2, x2, xj2, al, kp, out);
    }
    
    //   data members
    unsigned int Nf;
    
  private:
    //    Convolution
    void __conv_x1(double, double, double, double, const su3_kp_i2 *, weight_hhc *); 
    void __conv_x2(double, double, double, double, const su3_kp_i2 *, weight_hhc *); 
  };
}

#endif


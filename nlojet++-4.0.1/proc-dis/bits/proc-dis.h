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
#ifndef __NLO_PROC_DIS_H__
#define __NLO_PROC_DIS_H__ 1

#include <complex>
#include <utility>

#include <bits/dis-weight.h>



namespace nlo {

  //    The amplitude classes
  class ampq2g1l2;
  class ampq2g2l2;
  class ampq2g3l2;
  class ampq4l2;
  class ampq4g1l2;
  struct su3_kp_i1;


  class _dis_jet_base 
  {
  protected:
    //   constructor
    explicit _dis_jet_base(unsigned int = 2U, unsigned int = 3U);
    
    //    typedefs
    typedef std::pair<double, std::complex<double> > _Pair;
    
    //   three parton contributions
    void amp_tree(ampq2g1l2 *, double *);
    void amp_tree_mch(ampq2g1l2 *, double *);

    void amp_1loop(ampq2g1l2 *, double *);
    void amp_1loop_mch(ampq2g1l2 *, double *);
    
    void amp_kp(double, ampq2g1l2 *, su3_kp_i1 *);
    void amp_kp_mch(double, ampq2g1l2 *, su3_kp_i1 *);
    
    void amp_ccg(ampq2g1l2 *, int, int, int, int, _Pair *);
    void amp_ccq(ampq2g1l2 *, int, int, int, int, _Pair *);
    void amp_cca(ampq2g1l2 *, int, int, int, int, _Pair *);
    
    //   four parton contributions
    void amp_tree(ampq2g2l2 *, ampq4l2 *, double *);
    void amp_tree_mch(ampq2g2l2 *, ampq4l2 *, double *);
    
    void amp_1loop(ampq2g2l2 *, ampq4l2 *, double *);
    void amp_1loop_mch(ampq2g2l2 *, ampq4l2 *, double *);
    
    void amp_kp(double, ampq2g2l2 *, ampq4l2 *, su3_kp_i1 *);
    void amp_kp_mch(double, ampq2g2l2 *, ampq4l2 *, su3_kp_i1 *);
    
    void amp_ccg(ampq2g2l2 *, int, int, int, int, int, _Pair *);
    void amp_ccq(ampq2g2l2 *, int, int, int, int, int, _Pair *);
    void amp_cca(ampq2g2l2 *, int, int, int, int, int, _Pair *);
    void amp_ccq(ampq4l2 *, int, int, int, int, int, _Pair *);
    void amp_cca(ampq4l2 *, int, int, int, int, int, _Pair *);
    
    //  five parton contributions
    void amp_tree(ampq2g3l2 *, ampq4g1l2 *, double *);
    void amp_tree_mch(ampq2g3l2 *, ampq4g1l2 *, double *);
    
    //  do the P & K convolution
    void convolution(double, double, double, double, const su3_kp_i1 *, weight_dis *);
    
    //  static members
    static double Dijk(const _Pair& v, const _Pair& cc){ 
      return v.first*cc.first + 2.0*real(v.second*cc.second);
    }
    
    //    data members
    unsigned int Nu, Nd, Nf;
    
  private:
    //    private data members
    double charge, charge_sqr;
  };
}    //   namespace nlo

#endif

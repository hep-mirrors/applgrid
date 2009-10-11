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
#ifndef __AMPQ4G1_H__
#define __AMPQ4G1_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq4g1 : private _Amp_base 
  {
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;
    
  public:
    //   constructor
    ampq4g1(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}

    struct color_subamp {
      _ComplexD A12, A34, A32, A14;
    };

    //   tree level helicity amplitudes
    typedef color_subamp amp_tree;

    void matrix_tree_pmpmp(int, int, int, int, int, amp_tree *) const;
    void matrix_tree_pmmpp(int, int, int, int, int, amp_tree *) const;
    void matrix_tree_mppmp(int, int, int, int, int, amp_tree *) const;
    void matrix_tree_mpmpp(int, int, int, int, int, amp_tree *) const;
    void matrix_tree_ppmmp(int, int, int, int, int, amp_tree *) const;
    void matrix_tree_mmppp(int, int, int, int, int, amp_tree *) const;

    //   1-loop level helicity amplitudes
    struct amp_1loop {
      color_subamp tree;
      color_subamp loop;
    };
      
    void matrix_1loop_pmpmp(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_pmmpp(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mppmp(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mpmpp(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_ppmmp(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mmppp(unsigned int, int, int, int, int, int, amp_1loop *) const;

    //   matrix element squared
    void su3_tree(int, int, int, int, int, double *) const;
    void su3_1loop(unsigned int, int, int, int, int, int, double *) const;
    void su3_1loop_mch(unsigned int, int, int, int, int, int, double *) const;

    //   color corraleted amplitude
    void su3_cc(int, int, int, int, int, int, int, _Pair *) const;
    
    //   amplitudes for the finite part
    void su3_kp(unsigned int, int, int, int, int, int, int, int, su3_kp_i2 *, double=1.0) const;

  private:
    //   private members
    void ampcc12(int, int, int, int, int, double *) const;
    void ampcc23(int, int, int, int, int, double *) const;
    void ampcc24(int, int, int, int, int, double *) const;
    void ampcc25(int, int, int, int, int, double *) const;
    void amphtree(int, int, int, int, int, _ComplexD *) const;
    void ampcc(int, int, int, int, int, int, int, double *) const;
    double amp1loop(const amp_1loop&) const;
    
    void A12ppp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A12mpp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A12pmp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A12mmp(unsigned int, int, int, int, int, int, _ComplexD *) const;

    void A14ppp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A14mpp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A14pmp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A14mmp(unsigned int, int, int, int, int, int, _ComplexD *) const;

    void A32ppp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A32mpp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A32pmp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A32mmp(unsigned int, int, int, int, int, int, _ComplexD *) const;

    void A34ppp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A34mpp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A34pmp(unsigned int, int, int, int, int, int, _ComplexD *) const;
    void A34mmp(unsigned int, int, int, int, int, int, _ComplexD *) const;

    _ComplexD D34ppp(unsigned int, int, int, int, int, int) const;
    _ComplexD D34mpp(unsigned int, int, int, int, int, int) const;
    _ComplexD D32ppp(unsigned int, int, int, int, int, int) const;
    _ComplexD D32mpp(unsigned int, int, int, int, int, int) const;
    _ComplexD D32pmp(unsigned int, int, int, int, int, int) const;
   
    _ComplexD E32ppp(int, int, int, int, int) const;
    _ComplexD E32mpp(int, int, int, int, int) const;
    _ComplexD E32pmp(int, int, int, int, int) const;
    _ComplexD E34ppp(int, int, int, int, int) const;
    _ComplexD E34pmp(int, int, int, int, int) const;
  };

  inline void ampq4g1::
  A34pmp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A34mpp(nf, p2,p1, p4,p3, p5, A);
    A[0] *= -1.0; A[1] *= -1.0;
  }

  inline void ampq4g1::
  A34mmp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A34ppp(nf, p2,p1, p4,p3, p5, A);
    A[0] *= -1.0; A[1] *= -1.0;
  }

  inline void ampq4g1::
  A32mmp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A32ppp(nf, p4,p3, p2,p1, p5, A);
    A[0] *= -1.0; A[1] *= -1.0;
  }


  inline void ampq4g1::
  A12ppp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A34ppp(nf, p3,p4, p1,p2, p5, A);
  }

  inline void ampq4g1::
  A12mpp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A34pmp(nf, p3,p4, p1,p2, p5, A);
  }

  inline void ampq4g1::
  A12pmp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A34mpp(nf, p3,p4, p1,p2, p5, A);
  }

  inline void ampq4g1::
  A12mmp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A34mmp(nf, p3,p4, p1,p2, p5, A);
  }


  inline void ampq4g1::
  A14ppp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A32ppp(nf, p3,p4, p1,p2, p5, A);
  }

  inline void ampq4g1::
  A14mpp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A32pmp(nf, p3,p4, p1,p2, p5, A);
  }

  inline void ampq4g1::
  A14pmp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A32mpp(nf, p3,p4, p1,p2, p5, A);
  }

  inline void ampq4g1::
  A14mmp(unsigned int nf, int p1, int p2, int p3, int p4, int p5, _ComplexD *A) const {
    A32mmp(nf, p3,p4, p1,p2, p5, A);
  }
}   //   namespace nlo

#endif

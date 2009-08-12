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
#ifndef __AMPQ2G1P2_H__
#define __AMPQ2G1P2_H__ 1

#include <ampq2g3.h>


namespace nlo {


  class ampq2g1p2 : private ampq2g3
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq2g1p2(const innerprod_type& __x, const random_generator& rng)
      : ampq2g3(__x, rng) {}
  
    //   tree level helicity amplitudes

    //   tree level matrix element squared
    double su3_tree(int, int, int, int, int) const;
   
    //  color corraleted amplitudes
    _Pair su3_cc(int, int, int, int, int, int, int) const;

    //   amplitudes for the finite part
    void su3_kp(unsigned int, int, int, int, int, int, int, int, su3_kp_i2 *, double=1.0) const;


    struct amp_1loop {
      _ComplexD U0, D0, U1, D1;
    };
    
    //   1-loop matrix elements (quark electric charges are included!)
    void matrix_1loop_mpmpp(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mppmp(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mpppm(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mppmm(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mpmpm(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_mpmmp(unsigned int, unsigned int, int, int, int, int, int, amp_1loop *) const;
    
    //   1-loop level matrix element squared (quark electric charges are included!)
    void su3_1loop(unsigned int, unsigned int, int, int, int, int, int, double *) const;
    void su3_1loop_mch(unsigned int, unsigned int, int, int, int, int, int, double *) const;
   
  private:
    //   private members
    _ComplexD amphtree(int, int, int, int, int) const;
    static void matrix_1loop(unsigned int, unsigned int, const _AmpPrim *, const _AmpPrim *, amp_1loop *);
  };
}   //   namespace nlo

#endif

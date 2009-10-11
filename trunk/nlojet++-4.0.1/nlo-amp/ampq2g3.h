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
#ifndef __AMPQ2G3_H__
#define __AMPQ2G3_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq2g3 : protected _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq2g3(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
  
    //   tree level helicity amplitudes
    void matrix_tree_pppmm(int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_ppmpm(int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmppm(int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_ppmmm(int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmpmm(int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmmpm(int, int, int, int, int, _ComplexD *) const;

    //   tree level matrix element squared
    double su3_tree(int, int, int, int, int) const;
   
    //  color corraleted amplitudes
    _Pair su3_cc(int, int, int, int, int, int, int) const;

    //   amplitudes for the finite part
    void su3_kp(unsigned int, int, int, int, int, int, int, int, su3_kp_i2 *, double=1.0) const;

    struct amp_1loop {
      _ComplexD A0, A1, A3, A4;
    };

    //   1-loop matrix elements
    void matrix_1loop_pmppm(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_ppmpm(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_pppmm(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_ppmmm(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_pmpmm(unsigned int, int, int, int, int, int, amp_1loop *) const;
    void matrix_1loop_pmmpm(unsigned int, int, int, int, int, int, amp_1loop *) const;

    //   1-loop level matrix element squared
    double su3_1loop(unsigned int, int, int, int, int, int) const;
    double su3_1loop_mch(unsigned int, int, int, int, int, int) const;
    
    //  protected:
    //   private members
    _ComplexD amp2p1m(int, int, int, int, int, int) const;
    _ComplexD amp1p2m(int, int, int, int, int, int) const;

    double ampcc12(int, int, int, int, int) const;
    double ampcc13(int, int, int, int, int) const;
    double ampcc34(int, int, int, int, int) const;
     
    double ampcc(int, int, int, int, int, int, int) const;
    _ComplexD amphtree(int, int, int, int, int, int) const;
    
    //   primitive amplitudes
    struct _AmpPrim {
      _ComplexD Atree, Asusy, Asf, AL;
    };

    //   permutations
    static const unsigned int s[5][6]; 
    
    _ComplexD Vg(int, int, int, int, int) const;
 
    void A1mpp(int, int, int, int, int, _AmpPrim&) const;
    void A1pmp(int, int, int, int, int, _AmpPrim&) const;
    void A1ppm(int, int, int, int, int, _AmpPrim&) const;

    void A1pmm(int, int, int, int, int, _AmpPrim&) const;
    void A1mpm(int, int, int, int, int, _AmpPrim&) const;
    void A1mmp(int, int, int, int, int, _AmpPrim&) const;

    void A2mpp(int, int, int, int, int, _AmpPrim&) const;
    void A2pmp(int, int, int, int, int, _AmpPrim&) const;
    void A2ppm(int, int, int, int, int, _AmpPrim&) const;

    void A2pmm(int, int, int, int, int, _AmpPrim&) const;
    void A2mpm(int, int, int, int, int, _AmpPrim&) const;
    void A2mmp(int, int, int, int, int, _AmpPrim&) const;

    double amp1loop(amp_1loop *) const;
    void matrix_1loop(unsigned int, const _AmpPrim *, 
		      const _AmpPrim *, amp_1loop *) const;
  };

  inline void ampq2g3::
  A1pmm(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const {
    swap(); A1ppm(p2, p1, p5, p4, p3, res); swap();
  }

  inline void ampq2g3::
  A1mpm(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const {
    swap(); A1pmp(p2, p1, p5, p4, p3, res); swap();
  }

  inline void ampq2g3::
  A1mmp(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const {
    swap(); A1mpp(p2, p1, p5, p4, p3, res); swap();
  }

  inline void ampq2g3::
  A2pmm(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const {
    swap(); A2mpp(p3, p2, p1, p5, p4, res); swap();
  }

  inline void ampq2g3::
  A2mpm(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const {
    swap(); A2ppm(p3, p2, p1, p5, p4, res); swap();
  }
 
  inline void ampq2g3::
  A2mmp(int p1, int p2, int p3, int p4, int p5, _AmpPrim& res) const {
    swap(); A2pmp(p3, p2, p1, p5, p4, res); swap();
  }
}   //   namespace nlo

#endif

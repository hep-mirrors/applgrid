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
#ifndef __AMPQ2G2_H__
#define __AMPQ2G2_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq2g2 : protected _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq2g2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
  
    //   tree level helicity amplitudes
    void matrix_tree_ppmm(int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmpm(int, int, int, int, _ComplexD *) const;

    struct amp_1loop {
      _ComplexD A0, A1, A3;
    };

    //   1-loop matrix elements
    void matrix_1loop_pmpm(int, int, int, int, amp_1loop *) const;
    void matrix_1loop_ppmm(int, int, int, int, amp_1loop *) const;

    //   matrix element squared
    double su3_tree(int, int, int, int) const;

    //  color corraleted amplitudes
    double su3_cc(int,int, int, int, int, int) const;

    //   amplitudes for the finite part
    void su3_kp(unsigned int, int, int, int, int, int, int, su3_kp_i2 *, double=1.0) const;

    //   1-loop level matrix element squared
    double su3_1loop(int, int, int, int) const;
    double su3_1loop_mch(int, int, int, int) const;

  private:
    //   private members
    double ampcc12(int, int, int, int) const;
    double ampcc13(int, int, int, int) const;
    double ampcc14(int, int, int, int) const;
    double ampcc34(int, int, int, int) const;

    _ComplexD amp2p1m(int, int, int, int, int, int) const;
    _ComplexD amp1p2m(int, int, int, int, int, int) const;
    
  protected:
    
    struct _AmpPrim {
      _ComplexD Atree, Asusy, AL;
    };
    
    void A1pm(int, int, int, int, _AmpPrim&) const;
    void A1mp(int, int, int, int, _AmpPrim&) const;
				
    void A2pm(int, int, int, int, _AmpPrim&) const;
    void A2mp(int, int, int, int, _AmpPrim&) const;

  private:
    double amp1loop(amp_1loop *) const;
    void matrix_1loop(const _AmpPrim *, const _AmpPrim *, amp_1loop *) const;
  };

  inline void 
  ampq2g2::A2pm(int p1, int p2, int p3, int p4, _AmpPrim& res) const {
    A2mp(p1,p4,p3,p2, res);
  }
  
  


}   //   namespace nlo

#endif

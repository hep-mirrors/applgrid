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
#ifndef __AMPQ4_H__
#define __AMPQ4_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq4 : private _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq4(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
  
    //   tree level helicity amplitudes
    
    //   1-loop helicity amplitudes
    struct amp_1loop {
      _ComplexD A0, A1, A2;
    };
    
    void matrix_1loop_pmpm(unsigned int, int, int, int, int, amp_1loop&) const;
    void matrix_1loop_pmmp(unsigned int, int, int, int, int, amp_1loop&) const;

    //   matrix element squared
    void su3_tree(int, int, int, int, double *) const;

    //  color corraleted amplitudes
    void su3_cc(int,int, int, int, int, int, double *) const;
    
    //   amplitudes for the finite part
    void su3_kp(int, int, int, int, int, int, su3_kp_i2 *, double=1.0) const;
    
    //   1-loop amplitudes
    void su3_1loop(unsigned int, int, int, int, int, double *) const;
    void su3_1loop_mch(unsigned int, int, int, int, int, double *) const;

  private:
    //   private members
    void ampcc12(int, int, int, int, double *) const;
    void ampcc13(int, int, int, int, double *) const;
    void ampcc14(int, int, int, int, double *) const;

  };
}   //   namespace nlo

#endif

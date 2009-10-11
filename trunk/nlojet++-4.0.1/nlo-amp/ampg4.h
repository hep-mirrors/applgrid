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
#ifndef __AMPG4_H__
#define __AMPG4_H__ 1

//   Standard includes
#include <utility>

//   nlo includes
#include <bits/amp-ampbase.h>


namespace nlo {


  class ampg4 : private _Amp_base
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampg4(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
    
    //  tree and 1-loop level helicity amplitudes  (SU(3) color)
    void matrix_tree(int, int, int, int, int, int , _ComplexD *) const;
    void matrix_1loop(unsigned, int, int, int, int, int, int, _ComplexD *) const;
     
    //  tree level matrix element squared
    double su3_tree(int, int, int, int) const;

    //  1-loop amplitudes
    double su3_1loop(unsigned, int, int, int, int) const;
    double su3_1loop_mch(unsigned, int, int, int, int) const;

    //  color corraleted amplitudes
    double su3_cc(int,int, int, int, int, int) const;

    //   amplitudes for the finite part
    void su3_kp(unsigned int, int, int, int, int, int, int, su3_kp_i2 *, double=1.0) const;

  private:
    //   private members 
    _ComplexD Atree(int, int, int, int, int, int) const;
    void A1mmpp(unsigned, int, int, int, int, _ComplexD *) const;
    void A1mpmp(unsigned, int, int, int, int, _ComplexD *) const;

    double ampcc12(int, int, int, int) const;
  };

  inline void ampg4::
  matrix_tree(int pi, int pj, int p1, int p2, int p3, int p4, _ComplexD *d) const {
    *d = Atree(p1,p2,p3,p4, pi,pj);
  }
}   //   namespace deb

#endif

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
#ifndef __AMPQ4G2_H__
#define __AMPQ4G2_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq4g2 : private _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;

  public:
    //   constructor
    ampq4g2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
  
    //   tree level helicity amplitudes

    //   matrix element squared
    void su3_tree(int, int, int, int, int, int, double *) const;
    
    //   other flavour decomposation
    void __su3_tree(int, int, int, int, int, int, double *) const;

  private:
    //   private members
    struct _Colmat;

    //  static members
    static const _Colmat _S_colmat;
    static const short _S_itab[4][8];

    //  calculate the color subamplitude
    void E20(int, int, int, int, int, int, const short *, _ComplexD *, double, bool) const;
    void E11(int, int, int, int, int, int, const short *, _ComplexD *, double, bool) const;
    void F20(int, int, int, int, int, int, const short *, _ComplexD *, double, bool) const;
    void F11(int, int, int, int, int, int, const short *, _ComplexD *, double, bool) const;
  };

  inline void ampq4g2::
  su3_tree(int p1, int p2, int p3, int p4, int p5, int p6, double *out) const
  {
    double tmp[3];
    __su3_tree(p1,p2,p3,p4,p5,p6, tmp);
    out[0] = tmp[0]; out[1] = tmp[0] + tmp[1] + tmp[2];
  }
    
}   //  namespace nlo

#endif

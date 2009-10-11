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
#ifndef __AMPQ4G1L2_H__
#define __AMPQ4G1L2_H__ 1

#include <bits/amp-ampbase.h>



namespace nlo {
  
  
  class ampq4g1l2 : private _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;

  public:
    //  constructors
    ampq4g1l2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}

    //   tree level helicity amplitudes (SU(3) color)
    void matrix_tree_pmpmp(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_pmpmm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_ppmmp(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_ppmmm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_pmmpp(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_pmmpm(int, int, int, int, int, int, int, _ComplexD *);

    //   tree level helicity amplitudes (color independent decomposation)
    void color_tree_pmpmp(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_pmpmm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_ppmmp(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_ppmmm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_pmmpp(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_pmmpm(int, int, int, int, int, int, int, _ComplexD *);

    //   matrix element squared
    void su3_tree(int, int, int, int, int, int, int, double *);
    void su3_tree_mch(int, int, int, int, int, int, int, double *);

  private:
    //   private members
    static double su3_tree_aa(_ComplexD *, _ComplexD *);
    static double su3_tree_ac(_ComplexD *, _ComplexD *);
    static void   su3_amptree(_ComplexD *, double *);

    _ComplexD A1ppp(int, int, int, int, int, int, int);
    _ComplexD A1ppm(int, int, int, int, int, int, int);
    _ComplexD A1pmp(int, int, int, int, int, int, int);
    _ComplexD A1pmm(int, int, int, int, int, int, int);
						       
    _ComplexD A2ppp(int, int, int, int, int, int, int);
    _ComplexD A2ppm(int, int, int, int, int, int, int);
    _ComplexD A2pmp(int, int, int, int, int, int, int);
    _ComplexD A2pmm(int, int, int, int, int, int, int);
						       
    _ComplexD A3ppp(int, int, int, int, int, int, int);
    _ComplexD A3pmm(int, int, int, int, int, int, int);
						       
    _ComplexD A4ppm(int, int, int, int, int, int, int);
    _ComplexD A4pmp(int, int, int, int, int, int, int);
  };

  inline void ampq4g1l2::
  color_tree_pmpmp(int p1, int p4, int p3, int p2, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pmpmp(p1, p4, p3, p2, p5, p6, p7, M);
  }
  
  inline void ampq4g1l2::
  color_tree_pmpmm(int p1, int p4, int p3, int p2, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pmpmm(p1, p4, p3, p2, p5, p6, p7, M);
  }
  
  inline void ampq4g1l2::
  color_tree_ppmmp(int p1, int p4, int p3, int p2, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_ppmmp(p1, p4, p3, p2, p5, p6, p7, M);
  }
  
  inline void ampq4g1l2::
  color_tree_ppmmm(int p1, int p4, int p3, int p2, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_ppmmm(p1, p4, p3, p2, p5, p6, p7, M);
  }
  
  inline void ampq4g1l2::
  color_tree_pmmpp(int p1, int p4, int p3, int p2, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pmmpp(p1, p4, p3, p2, p5, p6, p7, M);
  }
  
  inline void ampq4g1l2::
  color_tree_pmmpm(int p1, int p4, int p3, int p2, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pmmpm(p1, p4, p3, p2, p5, p6, p7, M);
  }
}   //   namespace nlo

#endif

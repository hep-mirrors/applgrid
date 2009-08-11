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
#ifndef __AMPQ2G3L2_H__
#define __AMPQ2G3L2_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {
  
  
  class ampq2g3l2 : private _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    
  public:
    //   constructor
    ampq2g3l2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
    
    //   tree level helicity amplitudes (SU(3) color)
    void matrix_tree_ppppm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_pmmmm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_pppmm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_ppmpm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_pmppm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_ppmmm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_pmpmm(int, int, int, int, int, int, int, _ComplexD *);
    void matrix_tree_pmmpm(int, int, int, int, int, int, int, _ComplexD *);

    //   tree level helicity amplitudes (color independent decomposation)
    void color_tree_ppppm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_pmmmm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_pppmm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_ppmpm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_pmppm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_ppmmm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_pmpmm(int, int, int, int, int, int, int, _ComplexD *);
    void color_tree_pmmpm(int, int, int, int, int, int, int, _ComplexD *);

    //   matrix element squared (SU(3) color)
    double su3_tree(int, int, int, int, int, int, int);
    double su3_tree_mch(int, int, int, int, int, int, int);

  private:
    //   private members
    static double su3_amptree(_ComplexD *);
    
    _ComplexD Appp(int, int, int, int, int, int, int); 
    _ComplexD Ammm(int, int, int, int, int, int, int); 
    _ComplexD Appm(int, int, int, int, int, int, int); 
    _ComplexD Apmp(int, int, int, int, int, int, int); 
    _ComplexD Ampp(int, int, int, int, int, int, int); 
    _ComplexD Apmm(int, int, int, int, int, int, int); 
    _ComplexD Ampm(int, int, int, int, int, int, int); 
    _ComplexD Ammp(int, int, int, int, int, int, int); 
  };
  
  inline void ampq2g3l2::
  color_tree_ppppm(int p1, int p2, int p3, int p4, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_ppppm(p1, p2, p3, p4, p5, p6, p7, M);
  }
  
  inline void ampq2g3l2::
  color_tree_pmmmm(int p1, int p2, int p3, int p4, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pmmmm(p1, p2, p3, p4, p5, p6, p7, M);
  }
  
  inline void ampq2g3l2::
  color_tree_pppmm(int p1, int p2, int p3, int p4, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pppmm(p1, p2, p3, p4, p5, p6, p7, M);
  }
  
  inline void ampq2g3l2::
  color_tree_ppmpm(int p1, int p2, int p3, int p4, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_ppmpm(p1, p2, p3, p4, p5, p6, p7, M);
  }
  
  inline void ampq2g3l2::
  color_tree_pmppm(int p1, int p2, int p3, int p4, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pmppm(p1, p2, p3, p4, p5, p6, p7, M);
  }
  
  inline void ampq2g3l2::
  color_tree_ppmmm(int p1, int p2, int p3, int p4, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_ppmmm(p1, p2, p3, p4, p5, p6, p7, M);
  }
  
  inline void ampq2g3l2::
  color_tree_pmpmm(int p1, int p2, int p3, int p4, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pmpmm(p1, p2, p3, p4, p5, p6, p7, M);
  }
    
  inline void ampq2g3l2::
  color_tree_pmmpm(int p1, int p2, int p3, int p4, int p5, int p6, int p7, _ComplexD *M) {
    matrix_tree_pmmpm(p1, p2, p3, p4, p5, p6, p7, M);
  }
}  //   namespace nlo

#endif

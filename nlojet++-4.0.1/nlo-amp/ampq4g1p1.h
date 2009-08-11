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
#ifndef __AMPQ4G1P1_H__
#define __AMPQ4G1P1_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq4g1p1 : private _Amp_base  
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq4g1p1(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
  
    //   tree level helicity amplitudes
    void matrix_tree_pmpmpp(double, double, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmmppp(double, double, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_mppmpp(double, double, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_mpmppp(double, double, int, int, int, int, int, int, _ComplexD *) const;
    
    void matrix_tree_pmpmmp(double, double, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmmpmp(double, double, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_mppmmp(double, double, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_mpmpmp(double, double, int, int, int, int, int, int, _ComplexD *) const;

    //   tree level matrix element squared
    void su3_tree(double, int, int, int, int, int, int, double *) const;
    double su3_tree(double, double, int, int, int, int, int, int) const;
    
    void su3_tree_mch(double, int, int, int, int, int, int, double *) const;
    double su3_tree_mch(double, double, int, int, int, int, int, int) const;

  private:
    //   private members
    _ComplexD g1(int, int, int, int, int, int) const;
    _ComplexD g2(int, int, int, int, int, int) const;
    _ComplexD f1(int, int, int, int, int, int) const;
    _ComplexD f2(int, int, int, int, int, int) const;
    _ComplexD f3(int, int, int, int, int, int) const;
    _ComplexD f4(int, int, int, int, int, int) const;
    void Amhv(double, double, const _ComplexD&, int, int, int, int, int, int, _ComplexD *) const;
    double amptree(_ComplexD *) const;
  };
  

  inline void ampq4g1p1::
  matrix_tree_pmpmpp(double Q1, double Q2, int p1, int p2, int p3, 
		     int p4, int p5, int p6, _ComplexD *E) const 
  {
    _ComplexD a24 = (_M_a->ptr)[p2][p4];
    Amhv(Q1,Q2, -a24*a24, p1,p2,p3,p4,p5,p6, E);
  }

  inline void ampq4g1p1::
  matrix_tree_mpmppp(double Q1, double Q2, int p1, int p2, int p3, 
		     int p4, int p5, int p6, _ComplexD *E) const 
  {
    _ComplexD a13 = (_M_a->ptr)[p1][p3];
    Amhv(Q1,Q2, -a13*a13, p1,p2,p3,p4,p5,p6, E);
  }

  inline void ampq4g1p1::
  matrix_tree_mppmpp(double Q1, double Q2, int p1, int p2, int p3, 
		     int p4, int p5, int p6, _ComplexD *E) const 
  {
    _ComplexD a14 = (_M_a->ptr)[p1][p4];
    Amhv(Q1,Q2, a14*a14, p1,p2,p3,p4,p5,p6, E);
  }

  inline void ampq4g1p1::
  matrix_tree_pmmppp(double Q1, double Q2, int p1, int p2, int p3, 
		     int p4, int p5, int p6, _ComplexD *E) const 
  {
    _ComplexD a23 = (_M_a->ptr)[p2][p3];
    Amhv(Q1,Q2, a23*a23, p1,p2,p3,p4,p5,p6, E);
  }
}   //   namespace nlo

#endif

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
#ifndef __AMPQ4L2_H__
#define __AMPQ4L2_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq4l2 : private _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;

  public:
    //  constructors
    ampq4l2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
    
    //   tree and 1-loop level helicity amplitudes (SU(3) color)
    void matrix_tree_ppmm(int, int, int, int, int, int, _ComplexD*) const;
    void matrix_tree_pmpm(int, int, int, int, int, int, _ComplexD*) const;
    void matrix_tree_pmmp(int, int, int, int, int, int, _ComplexD*) const;

    void matrix_1loop_ppmm(unsigned int, int, int, int, int, int, int, _ComplexD*) const;
    void matrix_1loop_pmpm(unsigned int, int, int, int, int, int, int, _ComplexD*) const;
    void matrix_1loop_pmmp(unsigned int, int, int, int, int, int, int, _ComplexD*) const;

    //  tree and 1-loop level helicity amplitudes 
    //  (color independent decomposation)
    void color_tree_ppmm(int, int, int, int, int, int, _ComplexD*) const;
    void color_tree_pmpm(int, int, int, int, int, int, _ComplexD*) const;
    void color_tree_pmmp(int, int, int, int, int, int, _ComplexD*) const;

    void color_1loop_ppmm(unsigned int, int, int, int, int, int, int, _ComplexD*) const;
    void color_1loop_pmpm(unsigned int, int, int, int, int, int, int, _ComplexD*) const;
    void color_1loop_pmmp(unsigned int, int, int, int, int, int, int, _ComplexD*) const;

    //  tree and 1-loop amplitudes (SU(3) color)
    void su3_tree(int, int, int, int, int, int, double *) const;
    void su3_tree_mch(int, int, int, int, int, int, double *) const;

    void su3_1loop(unsigned int, int, int, int, int, int, int, double *) const;
    void su3_1loop_mch(unsigned int, int, int, int, int, int, int, double *) const;

    //    color corraleted amplitudes
    void su3_cc(int, int, int, int, int, int, int, int, double *) const;
    
    //        give the remains of I term
    void su3_ins(int, int, int, int, int, int, double *, double = 1.0) const;
    void su3_ins_mch(int, int, int, int, int, int, double *, double = 1.0) const;

    //        give the K, P terms with one initial state parton
    void su3_kp(int, int, int, int, int, int, int, su3_kp_i1 *, double = 1.0) const;
    void su3_kp_mch(int, int, int, int, int, int, int, su3_kp_i1 *, double = 1.0) const;

  private:
    //  private members
    static void su3_amptree(_ComplexD *, double *);
    static void su3_ampcc  (_ComplexD *, double *); 
    static void su3_ampcc  (_ComplexD *, double, double, double, double *);
    static void su3_amploop(_ComplexD *, double *);

    void su3_m1_ppmm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    void su3_m1_pmpm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;

    void color_m1_ppmm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    void color_m1_pmpm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;

    _ComplexD App(int, int, int, int, int, int) const;
    _ComplexD Apm(int, int, int, int, int, int) const;
    _ComplexD Aax(int, int, int, int, int, int) const;
    _ComplexD Fpp(int, int, int, int, int, int) const;
    _ComplexD Fpm(int, int, int, int, int, int) const;
    _ComplexD Fsl(int, int, int, int, int, int) const;

    _ComplexD FApp(int, int, int, int, int, int) const;
    _ComplexD FApm(int, int, int, int, int, int) const;
    _ComplexD FCsl(int, int, int, int, int, int) const;
    _ComplexD F1sl(int, int, int, int, int, int) const;
  };

  inline void ampq4l2::
  matrix_tree_ppmm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const
  {
    M[0] = App(p1, p2, p3, p4, p5, p6);
    swap(); M[1] = App(p3, p4, p1, p2, p6, p5); swap();
    M[2] = M[3] = 0;
  } 
    
  inline void ampq4l2::
  matrix_tree_pmpm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const
  {
    M[0] = Apm(p1, p2, p3, p4, p5, p6);
    M[1] = Apm(p3, p4, p1, p2, p5, p6);
    M[2] = Apm(p1, p4, p3, p2, p5, p6);
    M[3] = Apm(p3, p2, p1, p4, p5, p6);
  }

  inline void ampq4l2::
  matrix_tree_pmmp(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const
  {
    M[0] = M[1] = 0;
    M[2] = App(p1, p4, p3, p2, p5, p6);
    swap(); M[3] = App(p3, p2, p1, p4, p6, p5); swap();
  } 
  
  inline void ampq4l2::
  color_tree_ppmm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const {
    matrix_tree_ppmm(p1, p2, p3, p4, p5, p6, M);
  } 
    
  inline void ampq4l2::
  color_tree_pmpm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const {
    matrix_tree_pmpm(p1, p2, p3, p4, p5, p6, M);
  }

  inline void ampq4l2::
  color_tree_pmmp(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const {
    matrix_tree_pmmp(p1, p2, p3, p4, p5, p6, M);
  } 
  
  inline std::complex<double> 
  ampq4l2::Fpp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD F = FApp(p1, p2, p3, p4, p5, p6);
    swap(); F -= FApp(p4, p3, p2, p1, p6, p5); swap();
    return F;
  }

  inline std::complex<double> 
  ampq4l2::Fpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD F = FApm(p1, p2, p3, p4, p5, p6);
    swap(); F -= FApm(p4, p3, p2, p1, p6, p5); swap();
    return F;
  }

  inline std::complex<double> 
  ampq4l2::Fsl(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return F1sl(p1, p2, p3, p4, p5, p6) + F1sl(p6, p5, p3, p4, p2, p1);
  }
}   //  namespace nlo

#endif

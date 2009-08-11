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
#ifndef __AMPQ2G1L2_H__
#define __AMPQ2G1L2_H__ 1

//   Standard includes
#include <utility>

//   nlo includes
#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq2g1l2 : private _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //  constructor
    ampq2g1l2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}

    //  tree and 1-loop level matrix elements (SU(3) color)
    void matrix_tree_ppm(int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmm(int, int, int, int, int, _ComplexD *) const;

    void matrix_1loop_ppm(int, int, int, int, int, _ComplexD *) const;
    void matrix_1loop_pmm(int, int, int, int, int, _ComplexD *) const;

    //  tree and 1-loop matrix elements (color independent decomposation)
    void color_tree_ppm(int, int, int, int, int, _ComplexD *) const;
    void color_tree_pmm(int, int, int, int, int, _ComplexD *) const;

    void color_1loop_ppm(int, int, int, int, int, _ComplexD *) const;
    void color_1loop_pmm(int, int, int, int, int, _ComplexD *) const;

    //  tree and 1-loop amplitudes (SU(3) color)
    double su3_tree(int, int, int, int, int) const;
    double su3_1loop(unsigned int, int, int, int, int, int) const;

    double su3_tree_mch(int, int, int, int, int) const;
    double su3_1loop_mch(unsigned int, int, int, int, int, int) const;

    //  color correlated amplitudes
    _Pair su3_cc(int, int, int, int, int, int, int) const;

    //        give the I terms
    double su3_ins(unsigned int, int, int, int, int, int, double = 1.0) const;
    double su3_ins_mch(unsigned int, int, int, int, int, int, double = 1.0) const;

    //        give the K, P terms
    void su3_kp(unsigned int, int, int, int, int, int, int, su3_kp_i1 *, double = 1.0) const;
    void su3_kp_mch(unsigned int, int, int, int, int, int, int, su3_kp_i1 *, double = 1.0) const;
    
  private:
    //   private members
    _ComplexD Atree1ppm(int, int, int, int, int) const;
    _ComplexD Fcc1ppm  (int, int, int, int, int) const;
    _ComplexD Fsc1ppm  (int, int, int, int, int) const;

    _ComplexD Atree2pmp(int, int, int, int, int) const;
    _ComplexD Fcc2pmp  (int, int, int, int, int) const;
    _ComplexD Fsc2pmp  (int, int, int, int, int) const;
    _ComplexD Fax2pmp  (int, int, int, int, int) const;

    _ComplexD Atree1pmm(int, int, int, int, int) const;  
    _ComplexD Fcc1pmm  (int, int, int, int, int) const;  
    _ComplexD Fsc1pmm  (int, int, int, int, int) const;    

    _ComplexD Atree2pmm(int, int, int, int, int) const;    
    _ComplexD Fcc2pmm  (int, int, int, int, int) const;    
    _ComplexD Fsc2pmm  (int, int, int, int, int) const;
    _ComplexD Fax2pmm  (int, int, int, int, int) const;
  };

  inline void ampq2g1l2::
  matrix_tree_ppm(int p1, int p2, int p3, int p5, int p6, _ComplexD *A) const { 
    *A = Atree1ppm(p1, p2, p3, p5, p6);
  }

  inline void ampq2g1l2::
  matrix_tree_pmm(int p1, int p2, int p3, int p5, int p6, _ComplexD *A) const {
    *A = Atree1pmm(p1, p2, p3, p5, p6);
  }

  inline void ampq2g1l2::
  color_tree_ppm(int p1, int p2, int p3, int p5, int p6, _ComplexD *A) const { 
    *A = Atree1ppm(p1, p2, p3, p5, p6);
  }

  inline void ampq2g1l2::
  color_tree_pmm(int p1, int p2, int p3, int p5, int p6, _ComplexD *A) const {
    *A = Atree1pmm(p1, p2, p3, p5, p6);
  }

  inline std::complex<double>
  ampq2g1l2::Atree1pmm(int p1, int p2, int p3, int p4, int p5) const {
    swap(); _ComplexD res = Atree1ppm(p3, p2, p1, p5, p4); swap();
    return res;
  }
  
  inline std::complex<double>
  ampq2g1l2::Fcc1pmm(int p1, int p2, int p3, int p4, int p5) const {
    swap(); _ComplexD res = Fcc1ppm(p3, p2, p1, p5, p4); swap();
    return res;
  }
  
  inline std::complex<double>
  ampq2g1l2::Fsc1pmm(int p1, int p2, int p3, int p4, int p5) const {
    swap(); _ComplexD res = Fsc1ppm(p3, p2, p1, p5, p4);swap();
    return res;
  }

  inline std::complex<double>
  ampq2g1l2::Atree2pmm(int p1, int p2, int p3, int p4, int p5) const {
    swap(); _ComplexD res = Atree2pmp(p2, p1, p3, p5, p4); swap();
    return res;
  }
  
  inline std::complex<double>
  ampq2g1l2::Fcc2pmm(int p1, int p2, int p3, int p4, int p5) const {
    swap(); _ComplexD res = Fcc2pmp(p2, p1, p3, p5, p4); swap();
    return res;
  }
  
  inline std::complex<double>
  ampq2g1l2::Fsc2pmm(int p1, int p2, int p3, int p4, int p5) const {
    swap(); _ComplexD res = Fsc2pmp(p2, p1, p3, p5, p4); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g1l2::Fax2pmm(int p1, int p2, int p3, int p4, int p5) const {
    swap(); _ComplexD res =  -Fax2pmp(p2, p1, p3, p5, p4); swap();
    return res;
  }
}   //   namespace nlo

#endif

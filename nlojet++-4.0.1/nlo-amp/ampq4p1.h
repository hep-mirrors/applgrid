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
#ifndef __AMPQ4P1_H__
#define __AMPQ4P1_H__ 1

#include <bits/amp-ampbase.h>
#include <bits/nlo-color.h>


namespace nlo {


  class ampq4p1 : private _Amp_base  
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq4p1(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
  
    //   tree level helicity amplitudes
    _ComplexD matrix_tree_pmpmp(double, double, int, int, int, int, int) const;
    _ComplexD matrix_tree_pmmpp(double, double, int, int, int, int, int) const;
    _ComplexD matrix_tree_mppmp(double, double, int, int, int, int, int) const;
    _ComplexD matrix_tree_mpmpp(double, double, int, int, int, int, int) const;

    //   tree level matrix element squared
    void su3_tree(double, int, int, int, int, int, double *) const;
    double su3_tree(double, double, int, int, int, int, int) const;
  
    //   color corraleted amplitude
    double su3_cc(double, double, int, int, int, int, int, int, int) const;
    void su3_cc(double, int, int, int, int, int, int, int, double *) const;
    
    //   amplitudes for the finite part
    void su3_kp(double, int, int, int, int, int, int, su3_kp_i1 *, double=1.0) const;
    void su3_kp(double, double, int, int, int, int, int, int, su3_kp_i1 *, double=1.0) const;
    
    //  1-loop level helicity amplitudes
    struct amp_1loop {
      _ComplexD A0, A1, B1;
    };
    
    void matrix_1loop_pmpmp(double, double, unsigned int, int, int, int, int, int, amp_1loop *) const; 
    void matrix_1loop_pmmpp(double, double, unsigned int, int, int, int, int, int, amp_1loop *) const; 
    void matrix_1loop_mppmp(double, double, unsigned int, int, int, int, int, int, amp_1loop *) const; 
    void matrix_1loop_mpmpp(double, double, unsigned int, int, int, int, int, int, amp_1loop *) const; 

    double su3_1loop(double, double, unsigned int, int, int, int, int, int) const;
    double su3_1loop_mch(double, double, unsigned int, int, int, int, int, int) const;

    void su3_1loop(double, unsigned int, int, int, int, int, int, double *) const;
    void su3_1loop_mch(double, unsigned int, int, int, int, int, int, double *) const;
 
  private:
    //   private members
    double ampcc(double, double, double, int, int, int, int, int) const;
    void ampcc(double, double, double, double, int, int, int, int, int, double *) const;
    _ComplexD Amhv(double, double, const _ComplexD&, int, int, int, int, int) const;

    //  1-loop primitive amplitudes
    _ComplexD u0ppp (int, int, int, int, int) const;
    _ComplexD u0pmp (int, int, int, int, int) const;
    
    _ComplexD u1lppp(int, int, int, int, int) const;
    _ComplexD u1lpmp(int, int, int, int, int) const;
    
    _ComplexD uAppp (int, int, int, int, int) const;
    _ComplexD uApmp (int, int, int, int, int) const;
    
    _ComplexD uBppp (int, int, int, int, int) const;
    _ComplexD uBpmp (int, int, int, int, int) const;

    struct _PrimAmp {
      _ComplexD u0, u1, u2;
    };
    
    void uppp(unsigned int, int, int, int, int, int, _PrimAmp *) const;
    void upmp(unsigned int, int, int, int, int, int, _PrimAmp *) const;
    void umpp(unsigned int, int, int, int, int, int, _PrimAmp *) const;
    void ummp(unsigned int, int, int, int, int, int, _PrimAmp *) const;

    static double amptree(const amp_1loop& M) {
      return 2.0*Na*real(M.A1*conj(M.A0));
    }
    
    static double amptree(const amp_1loop& M12, const amp_1loop& M14) {
      _ComplexD A1 = M12.A1 + M14.B1/Nc, B1 = M12.B1 + M14.A1*Nc;
      return 2.0*Na*real(A1*conj(M12.A0) + B1*conj(M14.A0)/Nc);
    }
  };
  
  inline std::complex<double> ampq4p1::
  matrix_tree_pmpmp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD a24 = (_M_a->ptr)[p2][p4];
    return Amhv(Q1,Q2, -a24*a24, p1,p2,p3,p4,p5);
  }
  
  inline std::complex<double> ampq4p1::
  matrix_tree_mpmpp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD a13 = (_M_a->ptr)[p1][p3];
    return Amhv(Q1,Q2, -a13*a13, p1,p2,p3,p4,p5);
  }
  
  inline std::complex<double> ampq4p1::
  matrix_tree_mppmp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD a14 = (_M_a->ptr)[p1][p4];
    return Amhv(Q1,Q2, a14*a14, p1,p2,p3,p4,p5);
  }
  
  inline std::complex<double> ampq4p1::
  matrix_tree_pmmpp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD a23 = (_M_a->ptr)[p2][p3];
    return Amhv(Q1,Q2, a23*a23, p1,p2,p3,p4,p5);
  }

  inline std::complex<double> 
  ampq4p1::uApmp(int p1, int p2, int p3, int p4, int p5) const {
    return -uAppp(p1,p2,p4,p3,p5);
  }
  
  inline std::complex<double> 
  ampq4p1::uBpmp(int p1, int p2, int p3, int p4, int p5) const {
    return uBppp(p1,p2,p4,p3,p5);
  }
 
}   //   namespace nlo

#endif

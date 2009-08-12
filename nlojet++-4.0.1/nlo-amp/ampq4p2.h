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
#ifndef __AMPQ4P2_H__
#define __AMPQ4P2_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq4p2 : private _Amp_base  
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq4p2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
  
    //   tree level helicity amplitudes

    //   tree level matrix element squared
    void su3_tree(double, int, int, int, int, int, int, double *) const;
    double su3_tree(double, double, int, int, int, int, int, int) const;
    
    void su3_tree_mch(double, int, int, int, int, int, int, double *) const;
    double su3_tree_mch(double, double, int, int, int, int, int, int) const;

  private:
    //   private members
    _ComplexD g1(int, int, int, int, int, int) const;
    _ComplexD g2(int, int, int, int, int, int) const;
    _ComplexD Amhv(double, double, int, int, int, int, int, int, int, int) const;
    
    _ComplexD Apmpmpp(double, double, int, int, int, int, int, int) const;
    _ComplexD Apmmppp(double, double, int, int, int, int, int, int) const;
    _ComplexD Amppmpp(double, double, int, int, int, int, int, int) const;
    _ComplexD Ampmppp(double, double, int, int, int, int, int, int) const;
    
    _ComplexD Apmpmmm(double, double, int, int, int, int, int, int) const;
    _ComplexD Apmmpmm(double, double, int, int, int, int, int, int) const;
    _ComplexD Amppmmm(double, double, int, int, int, int, int, int) const;
    _ComplexD Ampmpmm(double, double, int, int, int, int, int, int) const;

    _ComplexD Apmpmmp(double, double, int, int, int, int, int, int) const; 
    _ComplexD Apmmpmp(double, double, int, int, int, int, int, int) const; 
    _ComplexD Amppmmp(double, double, int, int, int, int, int, int) const; 
    _ComplexD Ampmpmp(double, double, int, int, int, int, int, int) const; 
    _ComplexD Apmpmpm(double, double, int, int, int, int, int, int) const; 
    _ComplexD Apmmppm(double, double, int, int, int, int, int, int) const; 
    _ComplexD Amppmpm(double, double, int, int, int, int, int, int) const; 
    _ComplexD Ampmppm(double, double, int, int, int, int, int, int) const; 
  };
  

  inline std::complex<double> ampq4p2::
  Apmpmpp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    return -Amhv(Q1,Q2,p2,p4,p1,p2,p3,p4,p5,p6);
  }

  inline std::complex<double> ampq4p2::
  Apmmppp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Amhv(Q1,Q2,p2,p3,p1,p2,p3,p4,p5,p6);
  }

  inline std::complex<double> ampq4p2::
  Amppmpp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Amhv(Q1,Q2,p1,p4,p1,p2,p3,p4,p5,p6);
  }
  
  inline std::complex<double> ampq4p2::
  Ampmppp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    return -Amhv(Q1,Q2,p1,p3,p1,p2,p3,p4,p5,p6);
  }



  inline std::complex<double> ampq4p2::
  Apmpmmm(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Ampmppp(Q1,Q2,p1,p2,p3,p4,p5,p6); swap();
    return Amp;
  }

  inline std::complex<double> ampq4p2::
  Apmmpmm(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Amppmpp(Q1,Q2,p1,p2,p3,p4,p5,p6); swap();
    return Amp;
  }

  inline std::complex<double> ampq4p2::
  Amppmmm(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Apmmppp(Q1,Q2,p1,p2,p3,p4,p5,p6); swap();
    return Amp;
  }

  inline std::complex<double> ampq4p2::
  Ampmpmm(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Apmpmpp(Q1,Q2,p1,p2,p3,p4,p5,p6); swap();
    return Amp;
  }
  


  inline std::complex<double> ampq4p2::
  Apmpmmp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Q1*Q1*g1(p1,p2,p3,p4,p5,p6) + Q2*Q2*g1(p3,p4,p1,p2,p5,p6)
      + Q1*Q2*(g2(p1,p2,p3,p4,p5,p6) + g2(p3,p4,p1,p2,p5,p6));
  }

  inline std::complex<double> ampq4p2::
  Apmmpmp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Q1*Q1*g1(p1,p2,p4,p3,p5,p6) + Q2*Q2*g1(p4,p3,p1,p2,p5,p6)
      - Q1*Q2*(g2(p1,p2,p4,p3,p5,p6) + g2(p4,p3,p1,p2,p5,p6));
  }
  
  inline std::complex<double> ampq4p2::
  Amppmmp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Apmmpmp(Q1,Q2,p2,p1,p4,p3,p5,p6); swap();
    return Amp;
  }
  
  inline std::complex<double> ampq4p2::
  Ampmpmp(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Apmpmmp(Q1,Q2,p2,p1,p4,p3,p5,p6); swap();
    return Amp;
  }

  inline std::complex<double> ampq4p2::
  Apmpmpm(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Ampmpmp(Q1,Q2,p1,p2,p3,p4,p5,p6); swap();
    return Amp;
  }

  inline std::complex<double> ampq4p2::
  Apmmppm(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Amppmmp(Q1,Q2,p1,p2,p3,p4,p5,p6); swap();
    return Amp;
  }

  inline std::complex<double> ampq4p2::
  Amppmpm(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Apmmpmp(Q1,Q2,p1,p2,p3,p4,p5,p6); swap();
    return Amp;
  }

  inline std::complex<double> ampq4p2::
  Ampmppm(double Q1, double Q2, int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD Amp = Apmpmmp(Q1,Q2,p1,p2,p3,p4,p5,p6); swap();
    return Amp;
  }
}   //   namespace nlo

#endif

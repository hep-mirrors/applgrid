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
#ifndef __AMPQ2G1P1_H__
#define __AMPQ2G1P1_H__ 1

#include <ampq2g2.h>



namespace nlo {


  class ampq2g1p1 : private ampq2g2
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq2g1p1(const innerprod_type& __x, const random_generator& rng)
      : ampq2g2(__x, rng) {}
  
    //   tree level helicity amplitudes

    //   tree level matrix element squared
    double su3_tree(int, int, int, int) const;
   
    //  color corraleted amplitudes
    double su3_cc(int, int, int, int, int, int) const;

    //   amplitudes for the finite part
    void su3_kp(unsigned int, int,  int, int, int, int, su3_kp_i1 *, double = 1.0) const;

    //   1-loop level matrix element squared
    double su3_1loop(int, int, int, int) const;
    double su3_1loop_mch(int, int, int, int) const;
     
  private:
    //   private members
    double matrix_1loop_pmpm(int, int, int, int) const;
    double matrix_1loop_ppmm(int, int, int, int) const;
   };
}   //   namespace nlo

#endif

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
#ifndef __AMPQ2G2P2_H__
#define __AMPQ2G2P2_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq2g2p2 : private _Amp_base  
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq2g2p2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
  
    //   tree level helicity amplitudes

    //   tree level matrix element squared
    double su3_tree(int, int, int, int, int, int) const;
    double su3_tree_mch(int, int, int, int, int, int) const;

  private:
    _ComplexD Amhv(int, int, int, int, int, int, int) const; 
    _ComplexD Apmmmpp(int, int, int, int, int, int) const;
    _ComplexD Apmpmmp(int, int, int, int, int, int) const;
    _ComplexD Apmmpmp(int, int, int, int, int, int) const; 
    
    _ComplexD Apmmppp(int, int, int, int, int, int) const;
    _ComplexD Apmpmpp(int, int, int, int, int, int) const;
    _ComplexD Apmppmp(int, int, int, int, int, int) const;
    _ComplexD Apmpppm(int, int, int, int, int, int) const;
    _ComplexD Apmmmmp(int, int, int, int, int, int) const;
    _ComplexD Apmmmpm(int, int, int, int, int, int) const;
    _ComplexD Apmmpmm(int, int, int, int, int, int) const;
    _ComplexD Apmpmmm(int, int, int, int, int, int) const;
    _ComplexD Apmppmm(int, int, int, int, int, int) const;
    _ComplexD Apmpmpm(int, int, int, int, int, int) const;
    _ComplexD Apmmppm(int, int, int, int, int, int) const;
  };

  
  inline std::complex<double> 
  ampq2g2p2::Apmmppp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Amhv(p3,p1,p2,p3,p4,p5,p6);
  }
  
  inline std::complex<double> 
  ampq2g2p2::Apmpmpp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Amhv(p4,p1,p2,p3,p4,p5,p6);
  }
  
  inline std::complex<double> 
  ampq2g2p2::Apmppmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Amhv(p5,p1,p2,p3,p4,p5,p6);
  }
  
  inline std::complex<double> 
  ampq2g2p2::Apmpppm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Amhv(p6,p1,p2,p3,p4,p5,p6);
  }
  
  inline std::complex<double> 
  ampq2g2p2::Apmmmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Apmpppm(p2,p1,p4,p3,p5,p6);
  }

  inline std::complex<double> 
  ampq2g2p2::Apmmmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Apmppmp(p2,p1,p4,p3,p5,p6);
  }

  inline std::complex<double> 
  ampq2g2p2::Apmmpmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Apmmppp(p2,p1,p4,p3,p5,p6);
  }

  inline std::complex<double> 
  ampq2g2p2::Apmpmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Apmpmpp(p2,p1,p4,p3,p5,p6);
  }

  inline std::complex<double> 
  ampq2g2p2::Apmppmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Apmmmpp(p2,p1,p4,p3,p5,p6);
  }
  
  inline std::complex<double> 
  ampq2g2p2::Apmpmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Apmpmmp(p2,p1,p4,p3,p5,p6);
  }
  
  inline std::complex<double> 
  ampq2g2p2::Apmmppm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return Apmmpmp(p2,p1,p4,p3,p5,p6);
  }
}   //   namespace nlo

#endif

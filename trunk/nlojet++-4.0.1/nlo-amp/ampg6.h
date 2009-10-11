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
#ifndef __AMPG6_H__
#define __AMPG6_H__ 1

//   nlo includes
#include <bits/amp-ampbase.h>


namespace nlo {


  class ampg6 : private _Amp_base
  {
    //   private types
    typedef std::complex<double> _ComplexD;

  public:
    //  constructors
    ampg6(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}

    //  tree level matrix element squared
    double su3_tree(int, int, int, int, int, int);

  private:
    //  private members 
    static const short perm[4][24];
    static const short hel[10][24];
  
    //  calculate the leading and subleading color contributions
    static double Leading   (_ComplexD *);
    static double SubLeading(_ComplexD *);
    
    //  color subamplitudes; 3+3-
    _ComplexD APPPMMM(int, int, int, int, int, int) const;
    _ComplexD APPMPMM(int, int, int, int, int, int) const;
    _ComplexD APMPMPM(int, int, int, int, int, int) const;
    
    _ComplexD APPMMPM(int p1, int p2, int p3, int p4, int p5, int p6) const { 
      return APPMPMM(p2,p1,p6,p5,p4,p3);
    }
    
    _ComplexD APPMMMP(int p1, int p2, int p3, int p4, int p5, int p6) const { 
      return APPPMMM(p6,p1,p2,p3,p4,p5);
    }
    
    _ComplexD APMPPMM(int p1, int p2, int p3, int p4, int p5, int p6) const { 
      return APPMPMM(p4,p3,p2,p1,p6,p5);
    }
    
    _ComplexD APMMMPP(int p1, int p2, int p3, int p4, int p5, int p6) const { 
      return APPPMMM(p1,p6,p5,p4,p3,p2);
    }
    
    _ComplexD APMPMMP(int p1, int p2, int p3, int p4, int p5, int p6) const { 
      return APPMPMM(p6,p1,p2,p3,p4,p5);
    }
    
    _ComplexD APMMPPM(int p1, int p2, int p3, int p4, int p5, int p6) const { 
      return APPMPMM(p4,p5,p6,p1,p2,p3);
    }
    
    _ComplexD APMMPMP(int p1, int p2, int p3, int p4, int p5, int p6) const { 
      return APPMPMM(p1,p6,p5,p4,p3,p2);
    }
  };
}   //  namespace nlo

#endif
  

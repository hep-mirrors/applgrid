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
#ifndef __AMPQ6_H__
#define __AMPQ6_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq6 : private _Amp_base
  {
    //   private types
    typedef std::complex<double> _ComplexD;

  public:
    //  constructors
    ampq6(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
    
    //  tree level matrix element squared
    void su3_tree(int, int, int, int, int, int, const char *, double *);

  private:
    //  private members 
    class _Colmat;
    static const _Colmat _S_colmat;

    static const short perm[36][6];
    static const short hel[36][7];
    static const short iq6[5][36];   
    static const short mapcol[4][36];
    
    //   color subamplitudes
    void diag(short, int, int, int, int, int, int, _ComplexD *) const;
  };
}   //  namespace nlo

#endif
  

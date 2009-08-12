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
#ifndef __AMPQ2G4_H__
#define __AMPQ2G4_H__ 1

#include <bits/amp-ampbase.h>


namespace nlo {


  class ampq2g4 : private _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;

  public:
    //   constructor
    ampq2g4(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}
    
    //   tree level helicity amplitudes

    //   matrix element squared
    double su3_tree(int, int, int, int, int, int);

  private:
    //   private members
    class _Colmat;
    static const _Colmat _S_colmat;

    static const short perm[4][24]; 
    static const short hel[6][24];
    
    //  calculate the color subamplitudes
    _ComplexD DPPMM(int, int, int, int, int, int) const;
    _ComplexD DMMPP(int, int, int, int, int, int) const;
    _ComplexD DMPMP(int, int, int, int, int, int) const;
    _ComplexD DPMPM(int, int, int, int, int, int) const;
    _ComplexD DPMMP(int, int, int, int, int, int) const;
    _ComplexD DMPPM(int, int, int, int, int, int) const;
  };
}    //  namespace nlo

#endif

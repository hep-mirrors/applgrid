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
#ifndef __NLO_RNG_CMRG_H__
#define __NLO_RNG_CMRG_H__ 1

#include <bits/hep-rng.h>


namespace nlo {

  class rng_cmrg : public random_generator
  {
  public:
    //   store the state of the generator
    struct state_type {
      unsigned long mt[624];
      int mti;
    };
    
    //   constructors
    explicit rng_cmrg(unsigned long int s = 4357) 
      : random_generator(0, 0xffffffffUL) { this -> set(s);}
    
    //   seeding the generator
    void srandom(unsigned long int);
    
    //   generate unsigned long and double random number
    unsigned long int get() const;
    double get_double() const;
    
  private:
    mutable state_type _M_state;
  };
  


}

#endif

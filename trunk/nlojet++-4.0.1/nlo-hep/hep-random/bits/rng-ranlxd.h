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
#ifndef __NLO_RNG_RANLXD_H__
#define __NLO_RNG_RANLXD_H__ 1

#include <bits/hep-rng.h>


namespace nlo {


  class rng_ranlxd : public random_generator
  {
  public:
    //   store the state of the generator
    struct state_type {
      double xdbl[12], carry;
      unsigned int ir, jr, ir_old, pr;
    };
    
    //   luxury levels
    enum luxury { ranlxd1 = 202, ranlxd2 = 397};
    
    //   constructors
    explicit rng_ranlxd(unsigned long int s = 1UL, luxury lux = ranlxd1) 
      : random_generator(0, 0xffffffffUL), _M_lux(lux) { this -> set(s);}
    
    //   seeding the generator
    void set(unsigned long int);
    
    //   generate unsigned long and double random number
    unsigned long int get() const;
    double get_double() const;
    
  private:
    //  state of the generator
    luxury _M_lux;
    mutable state_type _M_state;
    
    //  private members
    void _M_increment_state() const;
    void _M_set_lux(unsigned long int, unsigned int);

    //  static members
    static const int next[12];
    static const double one_bit;
  };
}

#endif

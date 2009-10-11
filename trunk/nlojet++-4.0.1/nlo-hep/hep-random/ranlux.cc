/* rng/ranlux.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* This is a lagged fibonacci generator with skipping developed by Luescher.
   The sequence is a series of 24-bit integers, x_n, 

   x_n = d_n + b_n

   where d_n = x_{n-10} - x_{n-24} - c_{n-1}, b_n = 0 if d_n >= 0 and
   b_n = 2^24 if d_n < 0, c_n = 0 if d_n >= 0 and c_n = 1 if d_n < 0,
   where after 24 samples a group of p integers are "skipped", to
   reduce correlations. By default p = 199, but can be increased to
   365.

   The period of the generator is around 10^171. 

   From: M. Luescher, "A portable high-quality random number generator
   for lattice field theory calculations", Computer Physics
   Communications, 79 (1994) 100-110.

   Available on the net as hep-lat/9309020 at http://xxx.lanl.gov/

   See also,

   F. James, "RANLUX: A Fortran implementation of the high-quality
   pseudo-random number generator of Luscher", Computer Physics
   Communications, 79 (1994) 111-114

   Kenneth G. Hamilton, F. James, "Acceleration of RANLUX", Computer
   Physics Communications, 101 (1997) 241-248

   Kenneth G. Hamilton, "Assembler RANLUX for PCs", Computer Physics
   Communications, 101 (1997) 249-253  
*/

#include "bits/rng-ranlux.h"


namespace nlo {
    

  
  unsigned long int rng_ranlux::_M_increment_state() const 
  {
    static const unsigned long int mask_lo = 0x00ffffffUL;
    static const unsigned long int mask_hi = ~0x00ffffffUL;
        
    unsigned int i = _M_state.i;
    unsigned int j = _M_state.j;
    long int delta = _M_state.u[j] - _M_state.u[i] - _M_state.carry;
    
    if(delta & mask_hi) {
      _M_state.carry = 1;
      delta &= mask_lo;
    } else _M_state.carry = 0;
    
    _M_state.u[i] = delta;
    if(i == 0) i = 23;
    else i--;
    
    _M_state.i = i;
    if(j == 0) j = 23;
    else j--;
    _M_state.j = j;
    
    return delta;
  }
    
  unsigned long int rng_ranlux::get() const 
  {
    const unsigned int skip = _M_state.skip;
    unsigned long int r = _M_increment_state();
      
    _M_state.n++;
    if(_M_state.n == 24) {
      unsigned int i;
      _M_state.n = 0;
      for(i = 0; i < skip; i++)
	_M_increment_state();
    }
      
    return r;
  }
  
  double rng_ranlux::get_double() const {
    return get()/16777216.0;
  }
    
  void rng_ranlux::_M_set_lux(unsigned long int s, unsigned int luxury)
  {
    static const unsigned long int mask_hi = ~0x00ffffffUL;
    static const unsigned long int two24 = 16777216;
    
    int i;
    long int seed;
    
    if(s == 0) s = 314159265;
    seed = s;
    
    // This is the initialization algorithm of F. James, 
    // widely in use for RANLUX.
    
    for(i = 0; i < 24; i++) {
      unsigned long int k = seed / 53668;
      seed = 40014 * (seed - k * 53668) - k * 12211;
      if(seed < 0) seed += 2147483563;
      _M_state.u[i] = seed % two24;
    }
    
    _M_state.i = 23;
    _M_state.j = 9;
    _M_state.n = 0;
    _M_state.skip = luxury - 24;
    
    if(_M_state.u[23] & mask_hi) _M_state.carry = 1;
    else _M_state.carry = 0;
  }
  
  void rng_ranlux::set(unsigned long int s) {
    _M_set_lux(s, _M_lux);
  }
}

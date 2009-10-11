/* rng/cmrg.c
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

/* This is a combined multiple recursive generator. The sequence is,

   z_n = (x_n - y_n) mod m1

   where the two underlying generators x and y are,

   x_n = (a_{1} x_{n-1} + a_{2} x_{n-2} + a_{3} x_{n-3}) mod m1
   y_n = (b_{1} y_{n-1} + b_{2} y_{n-2} + b_{3} y_{n-3}) mod m2

   with coefficients a11 ... a23,

   a_{1} = 0,     a_{2} = 63308, a_{3} = -183326
   b_{1} = 86098, b_{2} = 0,     b_{3} = -539608

   and moduli m1, m2,

   m1 = 2^31 - 1 = 2147483647
   m2 = 2^31 - 2000169 = 2145483479

   We initialize the generator with 

   x_1 = s_1 MOD m1, x_2 = s_2 MOD m1, x_3 = s_3 MOD m1
   y_1 = s_4 MOD m2, y_2 = s_5 MOD m2, y_3 = s_6 MOD m2

   where s_n = (69069 * s_{n-1}) mod 2^32 and s_0 = s is the
   user-supplied seed.

   NOTE: According to the paper the initial values for x_n must lie in
   the range 0 <= x_n <= (m1 - 1) and the initial values for y_n must
   lie in the range 0 <= y_n <= (m2 - 1), with at least one non-zero
   value -- our seeding procedure satisfies these constraints.

   We then use 7 iterations of the generator to "warm up" the internal
   state.

   The theoretical value of z_{10008} is 719452880. The subscript 10008
   means (1) seed the generator with s=1, (2) do the seven warm-up
   iterations that are part of the seeding process, (3) then do 10000
   actual iterations.

   The period of this generator is about 2^205.

   From: P. L'Ecuyer, "Combined Multiple Recursive Random Number
   Generators," Operations Research, 44, 5 (1996), 816--822.

   This is available on the net from L'Ecuyer's home page,

   http://www.iro.umontreal.ca/~lecuyer/myftp/papers/combmrg.ps
   ftp://ftp.iro.umontreal.ca/pub/simulation/lecuyer/papers/combmrg.ps */


#include "bits/rng-cmrg.h"

namespace nlo {
    
    
  unsigned long int rng_cmrg::get() const 
  {
    static const long int m1 = 2147483647, m2 = 2145483479;
    static const long int a2 = 63308, qa2 = 33921, ra2 = 12979;
    static const long int a3 = -183326, qa3 = 11714, ra3 = 2883;
    static const long int b1 = 86098, qb1 = 24919, rb1 = 7417;
    static const long int b3 = -539608, qb3 = 3976, rb3 = 2071;
      
    {
      long int h3 = _M_state.x3 / qa3;
      long int p3 = -a3 * (_M_state.x3 - h3 * qa3) - h3 * ra3;
	
      long int h2 = _M_state.x2 / qa2;
      long int p2 = a2 * (_M_state.x2 - h2 * qa2) - h2 * ra2;
	
      if(p3 < 0) p3 += m1;
      if(p2 < 0) p2 += m1;
	
      _M_state.x3 = _M_state.x2;
      _M_state.x2 = _M_state.x1;
      _M_state.x1 = p2 - p3;
      if (_M_state.x1 < 0)
	_M_state.x1 += m1;
    }
      
    {
      long int h3 = _M_state.y3 / qb3;
      long int p3 = -b3 * (_M_state.y3 - h3 * qb3) - h3 * rb3;
	
      long int h1 = _M_state.y1 / qb1;
      long int p1 = b1 * (_M_state.y1 - h1 * qb1) - h1 * rb1;
	
      if (p3 < 0) p3 += m2;
      if (p1 < 0) p1 += m2;
	
      _M_state.y3 = _M_state.y2;
      _M_state.y2 = _M_state.y1;
      _M_state.y1 = p1 - p3;
      if(_M_state.y1 < 0)
	_M_state.y1 += m2;
    }
      
    if(_M_state.x1 < _M_state.y1)
      return (_M_state.x1 - _M_state.y1 + m1);
    else return (_M_state.x1 - _M_state.y1);
  }
    
  double rng_cmrg::get_double() const {
    return get()/2147483647.0 ;
  }
  
  void rng_cmrg::set(unsigned long int s)
  {
    // An entirely adhoc way of seeding! 
    // This does **not** come from L'Ecuyer et al
    static const long int m1 = 2147483647, m2 = 2145483479;
    
    if(s == 0) s = 1;	// default seed is 1
      
#define LCG(n) ((69069 * n) & 0xffffffffUL)
    s = LCG(s); _M_state.x1 = s % m1;
    s = LCG(s); _M_state.x2 = s % m1;
    s = LCG(s); _M_state.x3 = s % m1;
      
    s = LCG (s); _M_state.y1 = s % m2;
    s = LCG (s); _M_state.y2 = s % m2;
    s = LCG (s); _M_state.y3 = s % m2;
      
    // "warm it up"
    for(unsigned int i = 0; i < 7; i++)
      this -> get();
  }
}

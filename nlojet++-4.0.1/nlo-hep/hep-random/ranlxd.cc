/* rng/ranlxd.c
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


/* This is an implementation of Martin Luescher's second generation
   double-precision (48-bit) version of the RANLUX generator. 

   Thanks to Martin Luescher for providing information on this
   generator.

*/

#include "bits/rng-ranlxd.h"


namespace nlo {

  
  const int rng_ranlxd::next[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0};
  const double rng_ranlxd::one_bit = 1.0/281474976710656.0;	/* 1/2^48 */
    
#define RANLUX_STEP(x1,x2,i1,i2,i3)		\
  x1=xdbl[i1] - xdbl[i2];			\
  if (x2 < 0)					\
  {						\
    x1-=one_bit;				\
    x2+=1;					\
  }						\
  xdbl[i3]=x2
  
  
  void rng_ranlxd::_M_increment_state() const
  {
    int k, kmax;
    double y1, y2, y3;
    
    double *xdbl = _M_state.xdbl;
    double carry = _M_state.carry;
    unsigned int ir = _M_state.ir;
    unsigned int jr = _M_state.jr;
    
    for(k = 0; ir > 0; ++k) {
      y1 = xdbl[jr] - xdbl[ir];
      y2 = y1 - carry;
    
      if(y2 < 0) {
	carry = one_bit;
	y2 += 1;
      } else carry = 0;
      
      xdbl[ir] = y2;
      ir = next[ir];
      jr = next[jr];
    }
    
    kmax = _M_state.pr - 12;  
    for(; k <= kmax; k += 12) {
      y1 = xdbl[7] - xdbl[0];
      y1 -= carry;
      
      RANLUX_STEP(y2, y1, 8, 1, 0);
      RANLUX_STEP(y3, y2, 9, 2, 1);
      RANLUX_STEP(y1, y3, 10, 3, 2);
      RANLUX_STEP(y2, y1, 11, 4, 3);
      RANLUX_STEP(y3, y2, 0, 5, 4);
      RANLUX_STEP(y1, y3, 1, 6, 5);
      RANLUX_STEP(y2, y1, 2, 7, 6);
      RANLUX_STEP(y3, y2, 3, 8, 7);
      RANLUX_STEP(y1, y3, 4, 9, 8);
      RANLUX_STEP(y2, y1, 5, 10, 9);
      RANLUX_STEP(y3, y2, 6, 11, 10);
      
      if(y3 < 0) {
	carry = one_bit;
	y3 += 1;
      } else carry = 0;
	    
      xdbl[11] = y3;
    }
    
    kmax = _M_state.pr;
    
    for(; k < kmax; ++k) {
      y1 = xdbl[jr] - xdbl[ir];
      y2 = y1 - carry;
      
      if(y2 < 0) {
	carry = one_bit;
	y2 += 1;
      } else carry = 0;
	    
      xdbl[ir] = y2;
      ir = next[ir];
      jr = next[jr];
    }
    
    _M_state.ir = ir;
    _M_state.ir_old = ir;
    _M_state.jr = jr;
    _M_state.carry = carry;
  }
    
  double rng_ranlxd::get_double() const
  {
    int ir = _M_state.ir;
    _M_state.ir = next[ir];
    
    if(_M_state.ir == _M_state.ir_old)
      _M_increment_state();
    
    return _M_state.xdbl[_M_state.ir];
  }
  
  unsigned long int rng_ranlxd::get() const {
    return (unsigned long int) (get_double()*4294967296.0); // 2^32
  }
    
  void rng_ranlxd::_M_set_lux (unsigned long int s, unsigned int luxury)
  {
    int ibit, jbit, i, k, l, xbit[31];
    double x, y;
    long int seed;
      
    if(s == 0) s = 1;
    seed = s;
    
    i = seed & 0xFFFFFFFFUL;
    for(k = 0; k < 31; ++k) {
      xbit[k] = i % 2;
      i /= 2;
    }
    
    ibit = 0; jbit = 18;
    for(k = 0; k < 12; ++k) {
      x = 0;
      
      for(l = 1; l <= 48; ++l) {
	y = (double) ((xbit[ibit] + 1) % 2);
	x += x + y;
	xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2;
	ibit = (ibit + 1) % 31;
	jbit = (jbit + 1) % 31;
      }
      _M_state.xdbl[k] = one_bit * x;
    }
      
    _M_state.carry = 0;
    _M_state.ir = 11;
    _M_state.jr = 7;
    _M_state.ir_old = 0;
    _M_state.pr = luxury;
  }
  
  void rng_ranlxd::set(unsigned long int s) {
    _M_set_lux(s, _M_lux);
  }
}


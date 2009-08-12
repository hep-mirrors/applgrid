//  Copyright (C) 2005 Zoltan Nagy
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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA


#include <math.h>
#include <cmath>
#include "bits/ran-gamma.h"

/* The Gamma distribution of order a>0 is defined by:

  p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

  for x>0.  If X and Y are independent gamma-distributed random
  variables of order a1 and a2 with the same scale parameter b, then
  X+Y has gamma distribution of order a1+a2.

  The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. \
*/



namespace nlo {

  
  double gamma_distribution<double>::operator()(double a, double b) const
  {
	/* assume a > 0 */
	unsigned int na = (unsigned int) std::floor(a);
	
	if(a == na) return b*_S_gamma_int(*_M_rng, na);
	else if(na == 0) return b*_S_gamma_frac(*_M_rng, a);
	else return b*(_S_gamma_int(*_M_rng, na) + _S_gamma_frac(*_M_rng, a - na));
  }

  
  double gamma_distribution<double>::
  _S_gamma_int(const random_generator& r, unsigned int a)
  {
	if(a < 12) {
	  double prod = 1.0;	
	  for(unsigned int i = 0; i < a; i++) prod *= r();
	
	/* Note: for 12 iterations we are safe against underflow, since
	  the smallest positive random number is O(2^-32). This means
	  the smallest possible product is 2^(-12*32) = 10^-116 which
	  is within the range of double precision. */
	
	  return -std::log(prod);
	} else return _S_gamma_large(r, (double) a);
  }


#define NLO_PI 3.14159265358979323846264338327950288
  
  double gamma_distribution<double>::
  _S_gamma_large(const random_generator& r, double a)
  {
	/* Works only if a > 1, and is most efficient if a is large
  
	This algorithm, reported in Knuth, is attributed to Ahrens.  A
	faster one, we are told, can be found in: J. H. Ahrens and
	U. Dieter, Computing 12 (1974) 223-246.  */
  
	double sqa = std::sqrt(2.0*a-1.0), x, y, v;

	do {
	  do {
		y = std::tan(NLO_PI*r());
		x = sqa*y + a-1.0;
	  } while (x <= 0.0);
	  v = r();
	} while (v > (1.0 + y*y)*std::exp((a-1.0)*std::log(x/(a-1.0)) - sqa*y));
  
	return x;
  }

  
#define NLO_E 2.71828182845904523536028747135266250
  
  double gamma_distribution<double>::
  _S_gamma_frac(const random_generator& r, double a)
  {
	/* This is exercise 16 from Knuth; see page 135, and the solution is
	on page 551.  */
	double p = NLO_E/(a + NLO_E), q, x, u, v;

	do {
	  u = r();
	  v = r();
	
	  if(u < p) {
		x = std::exp((1.0/a)*std::log(v));
		q = std::exp(-x);
	  } else {
		x = 1.0 - std::log(v);
		q = std::exp((a - 1.0)*std::log(x));
	  }
	} while(r() >= q);
  
	return x;
  }

  double gamma_distribution<double>::pdf(double x, double a, double b)
  {
	if(x < 0.0) return 0;
	else if(x == 0.0) return (a == 1.0 ? 1.0/b : 0.0);
	else if(a == 1.0) return std::exp(-x/b)/b;
	else return std::exp((a-1.0)*std::log(x/b) - x/b - ::lgamma(a))/b;
  }

}



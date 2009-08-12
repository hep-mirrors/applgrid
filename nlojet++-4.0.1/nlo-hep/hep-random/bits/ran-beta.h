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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#ifndef __NLO_RAN_BETA_H__
#define __NLO_RAN_BETA_H__ 1

//   C includes
#include <math.h>

//   C++ includes
#include <cmath>

//   nlo includes
#include <bits/ran-gamma.h>


namespace nlo {
  
  template<class _Float>
  class beta_distribution;
  
  
  template<>
  class beta_distribution<double>
  {
  public:
	//   public types
	typedef double value_type;
	
	//   constructors
	explicit beta_distribution(const random_generator& r)
	  : _M_gamma(r) {} 
	
	//   returns the random number
	value_type operator()(double a, double b) const {
	  double x1 = _M_gamma(a);
	  return x1/(x1+_M_gamma(b));
	}
	
	//   probability distribution function 
	static value_type pdf(value_type x, double a, double b) {
	  if(x < 0.0 || x > 1.0) return 0.0;
	  else return std::exp(::lgamma(a+b)-::lgamma(a)-::lgamma(b))*std::pow(x, a-1.0)*std::pow(1.0-x, b-1.0);
	}
	
  private:
	//   random engine
	gamma_distribution<double> _M_gamma;
  };
  
}
#endif


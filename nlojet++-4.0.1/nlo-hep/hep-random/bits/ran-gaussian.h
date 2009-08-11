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
#ifndef __NLO_RAN_GAUSSIAN_H__
#define __NLO_RAN_GAUSSIAN_H__ 1

//   C++ includes
#include <cmath>

//   nlo includes
#include <bits/hep-rng.h>


namespace nlo {
  
  template<class _Float>
  class gaussian_distribution;
  
  
  template<>
  class gaussian_distribution<double>
  {
  public:
	//   public types
	typedef double value_type;
	
	//   constructors
	explicit gaussian_distribution(const random_generator& r)
	  : _M_rng(&r), _M_new(true) {} 
	
	//   returns the random number
	value_type operator()(double = 1.0) const;
	
	//   probability distribution function 
	static value_type pdf(value_type x, double sigma = 1.0) 
	{
	  double u = x/sigma;
	  return std::exp(-0.5*u*u)/(2.50662827463100050241576528481104525*std::fabs(sigma));
	}

  private:
	//   random engine
	const random_generator *_M_rng;
	
	//   keep the second random number for the next call
	mutable bool _M_new;
	mutable double _M_x2, _M_sigma;
  };
  
  
}
#endif

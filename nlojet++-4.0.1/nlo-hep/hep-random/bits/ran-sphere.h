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
#ifndef __NLO_RAN_SPHERE_H__
#define __NLO_RAN_SPHERE_H__ 1

//   nlo includes
#include <bits/hep-rng.h>


namespace nlo {
  
  template<class _Tp>
  class sphere_distribution;
  
  
  template<>
  class sphere_distribution<double>
  {
  public:
	//   public types
	typedef double value_type;
	
	//   constructors
	explicit sphere_distribution(const random_generator& r)
	  : _M_rng(&r) {} 
	
	//   returns the random number
	void operator()(unsigned int, double *) const;
		
  private:
	//   random engine
	const random_generator *_M_rng;
  };
  
  
}
#endif



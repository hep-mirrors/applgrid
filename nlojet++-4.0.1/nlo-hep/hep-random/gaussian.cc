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


#include "bits/ran-gaussian.h"



namespace nlo {
  
#define NLO_2PI 6.28318530717958647692528676655900577
  
  
  double gaussian_distribution<double>::operator()(double sigma) const
  {
	if(_M_new || _M_sigma != sigma) {
	  double r = sigma*std::sqrt(-2.0*std::log(_M_rng -> get_double()));
	  double phi = NLO_2PI*(_M_rng -> get_double());
	  double x1 = r*std::cos(phi);
	  
	  _M_x2 = r*std::sin(phi);
	  _M_new = false;
	  _M_sigma = sigma;
	  
	  return x1;
	} else {
	  _M_new = true;
	  return _M_x2;
	}
  }

}


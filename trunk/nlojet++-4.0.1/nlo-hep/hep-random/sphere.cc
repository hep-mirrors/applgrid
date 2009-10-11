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


#include <cmath>
#include "bits/ran-sphere.h"
#include "bits/ran-gaussian.h"


namespace nlo {
  
#define NLO_2PI 6.28318530717958647692528676655900577
 
  void sphere_distribution<double>::operator()(unsigned int dim, double *x) const
  {
	switch(dim) {
	case 0: break;
	case 1: { x[0] = _M_rng->get_double();} break;	
	case 2: case 3: 
	  {
		double phi = NLO_2PI*(_M_rng->get_double());
		x[0] = std::cos(phi);
		x[1] = std::sin(phi);
		
		if(dim == 3U) {
		  x[2] = 1.0 - 2.0*(_M_rng->get_double());
		  double s = std::sqrt(1.0-x[2]*x[2]);
		  x[0] *= s; x[1] *= s; 
		}
	  } break;
	default:
	  {
		gaussian_distribution<double> g(*_M_rng);
		
		double s2 = 0.0;
		for(unsigned int i = 0; i < dim; i++) {
		  x[i] = g();
		  s2 += x[i]*x[i];
		}
		
		double s = 1.0/std::sqrt(s2);
		for(unsigned int i = 0; i < dim; i++)
		  x[i] *= s;
	  } break;
	}
	
	return;
  }
	
	
}


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

//   Standard includes
#include <cmath>

//   nlo++ includes
#include "bits/amp-ampbase.h"


#define __DEB_PI   3.14159265358979323846  
#define __DEB_PI_3 3.28986813369645287293
#define __DEB_PI_6 1.64493406684822643646

namespace nlo {


  double __specfunc_cl2(double x) 
  {
    static double a[] 
      = {0.0279528319735756613, 1.763088743898116e-4, 1.2662741461157e-6,
	 1.17171818134e-8, 1.230064129e-10,1.3952729e-12, 1.66908e-14,
	 2.076e-16,2.7e-18 };
    static double b[] 
      = {0.639097088857265341,-.054980569301851716, -9.61261945950606e-4,
	 -3.205468682255e-5, -1.329461695426e-6, -6.2093601824e-8,
	 -3.129600656e-9, -1.66351954e-10, -9.196527e-12, -5.24004e-13,
	 -3.058e-14, -1.82e-15, -1.1e-16, -7e-18 };
    double d, alfa, h, s, u, v, b0, b1, b2;
    short i;
    
    d = std::abs(x);
    v = d-2.0*__DEB_PI*std::floor(d/(2.0*__DEB_PI));
    s = ( x >= 0.0 ? 1.0 : -1.0);
    
    if(v > __DEB_PI) {
      v = 2.0*__DEB_PI-v;
      s = -s;
    }
    
    if (v == 0.0 || v == __DEB_PI) h = 0.0; 
    else if(v < __DEB_PI/2.0) {
      u = 2.0*v/__DEB_PI;
      d = u;
      h = d*d*2.0-1.0;
      alfa = 2.0*h;
      b1 = b2=0.0;
      
      for(i = 8; i >= 0; --i) {
	b0 = a[i]+alfa*b1-b2;
	b2 = b1;
	b1 = b0;
      }
      d = v;
      h = v*(1-std::log(v)+d*d*0.5*(b0-h*b2));
    } else {
      u = 2.0*v/__DEB_PI-2.0;
      d = u;
      h = d*d*2.0-1.0;
      alfa = 2.0*h ;
      b1 = b2 = 0.0;
      
      for(i = 13; i >= 0; --i) {
	b0 = b[i]+alfa*b1-b2;
	b2 = b1;
	b1 = b0;
      }
      h = (__DEB_PI-v)*(b0-h*b2);
    }
    return s*h;
  }
}

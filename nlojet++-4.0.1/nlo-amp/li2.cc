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
  
  double __specfunc_li2(double x) 
  {    
    static double c[] 
      = {0.429966935608137, 0.4097598753307711, -0.018588436650146,
	 0.0014575108406227, -1.430418444234e-4, 1.58841554188e-5,
	 -1.9078495939e-6, 2.419518085e-7, -3.19334127e-8, 4.3454506e-9,
	 -6.057848e-10, 8.6121e-11, -1.24433e-11, 1.8226e-12, -2.701e-13,
	 4.04e-14, -6.1e-15, 9e-16, -1e-16 };
    double d1, d2, alfa, a, h, s, t, y, b0, b1, b2;
    short i;

    if(x == 1.0)  return __DEB_PI_6;
    if(x == -1.0) return -0.5*__DEB_PI_6;
  
    t = -x;
    if(t <= -2.0) {
      y  = -1.0/(1.0 + t);
      s  = 1.0;
      d1 = std::log(-t);
      d2 = std::log(1.0 + 1./t);
      a = -__DEB_PI_3 + 0.5*(d1*d1 - d2*d2);
    } else if(t < -1.0) {
      y = -1.0 - t;
      s = -1.0;
      a = std::log(-t);
      a = -__DEB_PI_6 + a*(a + std::log(1.0 + 1.0/t));
    } else if(t <= -0.5) {
      y = (-1.0 - t)/t;
      s = 1.0;
      a = std::log(-t);
      a = -__DEB_PI_6 + a*(-0.5*a + std::log(1.0 + t));
    } else if(t < 0.0) {
      y = -t / (1.0 + t);
      s = -1.0;
      d1 = std::log(1.0 + t);
      a = 0.5*(d1 * d1);
    } else if(t <= 1.0) {
      y = t;
      s = 1.0;
      a = 0.0;
    } else {
      y = 1.0/t;
      s = -1.0;
      d1 = std::log(t);
      a = __DEB_PI_6 + 0.5*d1*d1;
    }
  
    h = 2.0*y - 1.0;
    alfa = 2.0*h;
    b1 = b2 = 0.0;
  
    for(i = 18; i >= 0; --i) {
      b0 = c[i] + alfa*b1 - b2;
      b2 = b1;
      b1 = b0;
    }
  
    return -(s*(b0 - h*b2) + a);
  } 

}

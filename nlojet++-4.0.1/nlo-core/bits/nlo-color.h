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
#ifndef __NLO_NLO_COLOR_H__
#define __NLO_NLO_COLOR_H__ 1

//   Standard includes
#include <cmath>


namespace nlo {


  //----- color factors (SU(N) and Tr = 0.5) -----
  extern const double Nc;              //   Nc  = 3;
  extern const double Na;              //   Na  = Nc^2-1;
  extern const double Nc2;             //   Nc2 = Nc^2;
  extern const double Na2;             //   Na2 = Na^2;
  extern const double Cf;              //   Cf  = Na/(2*Nc)
  extern const double Ca;              //   Ca  = Nc;
  extern const double Gq;              //   Gq  = 1.5*Cf 

  inline double Gg(unsigned int nf) {
    return (11.0*Ca - 2.0*nf)/6.0;
  }

  inline double B0(unsigned int nf) {
    return (11.0*Ca - 2.0*nf)/3.0;
  }

  inline double B1(unsigned int nf) {
    return (17.0*Ca*Ca - 3.0*Cf*nf - 5.0*Ca*nf)/3.0;
  }


  inline double Kq() {
    return 2.47342124420236475137;
  }

  inline double Kq(double al) 
  {
    double la = std::log(al);
    return 2.47342124420236475137 - Cf*la*la + Gq*(al - 1.0 - la);
  }

  inline double Kg(unsigned int nf, double al) 
  {
    double la = std::log(al);
    return 2.07728815537399578575*Ca - 0.55555555555555555555*nf
      -    Ca*la*la + Gg(nf)*(al - 1.0 - la); 
  }

  inline double Kg(unsigned int nf) {
    return 2.07728815537399578575*Ca - 0.55555555555555555555*nf;
  }
}   //   namespace nlo   

#endif

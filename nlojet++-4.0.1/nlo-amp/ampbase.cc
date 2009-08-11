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

//   nlo includes
#include "bits/nlo-color.h"
#include "bits/amp-ampbase.h"


#define DEB_PI  3.14159265358979323846  
#define DEB_PI3 3.28986813369645287293
#define DEB_PI6 1.64493406684822643646


namespace nlo {


  double _Amp_base::Xq(double s, double Q) {
    return 1.5*std::log(std::abs(s/Q)) + DEB_PI6*(s > 0.0 ? -1.0 : 2.0);
  }
  
  double _Amp_base::Xg(double s, double Q, unsigned int nf) { 
    return Gg(nf)*std::log(std::abs(s/Q))/Nc + DEB_PI6*(s > 0.0 ? -1.0 : 2.0);
  }

  std::complex<double> _Amp_base::Log(double x) {
    return _ComplexD(std::log(std::abs(x)), (x < 0.0 ? -DEB_PI : 0.0));
  }

  std::complex<double> _Amp_base::Li2(double x) 
  {
    if (x <= 1.0) return __specfunc_li2(x);
    else {
      double lx = std::log(x);
      return _ComplexD(DEB_PI6 - __specfunc_li2(1.0-x) - std::log(x-1)*lx, -DEB_PI*lx);
    }
  }
  
  std::complex<double> _Amp_base::Li2(double x, double y) 
  {
    double r = x/y;
    if (r >= 0.0) return __specfunc_li2(1.0-r);
    else return DEB_PI6 - __specfunc_li2(r) - std::log(1.0-r)*(Log(x) - Log(y));
  }

  std::complex<double> 
  _Amp_base::Li2(double x1, double x2, double y1, double y2) 
  {
    double r = (x1*x2)/(y1*y2);
    if (r >= 0.0) return __specfunc_li2(1.0-r);
    else return DEB_PI6 - __specfunc_li2(r) 
	   - std::log(1.0-r)*(Log(x1) + Log(x2) - Log(y1) - Log(y2));
  }

  std::complex<double> _Amp_base::L0(double x, double y) {
    return (Log(x) - Log(y))/(1.0 - x/y);
  }
  
  std::complex<double> _Amp_base::L1(double x, double y) {
    return (L0(x, y) + 1.0)/(1.0 - x/y);
  }

  std::complex<double> _Amp_base::L2(double x, double y) 
  {
    double r = x/y;
    double omr = 1-r;
    double omr3 = omr*omr*omr;

    return (Log(x)-Log(y)-(r-1.0/r)*0.5)/omr3;
  }

  std::complex<double> 
  _Amp_base::Ls_1(double x, double y, double z, double w) {
    return Li2(x, y) + Li2(z, w) + (Log(x)-Log(y))*(Log(z)-Log(w)) - DEB_PI6;
  }

  std::complex<double> 
  _Amp_base::Ls0(double x, double y, double z, double w) {
    return Ls_1(x,y,z,w)/(1.0-x/y-z/w);
  }

  std::complex<double> 
  _Amp_base::Ls1(double x, double y, double z, double w) {
    return (Ls0(x,y,z,w) + L0(x,y) +  L0(z,w))/(1.0-x/y-z/w);
  }

  std::complex<double> 
  _Amp_base::Ls2(double x, double y, double z, double w) {
    return (Ls1(x,y,z,w) + 0.5*(L1(x,y) +  L1(z,w)))/(1.0-x/y-z/w);
  }

  std::complex<double> 
  _Amp_base::Ls3(double x, double y, double z, double w) {
    return (Ls2(x,y,z,w) + (L2(x,y) +  L2(z,w))/3.0 - (y/x + w/z)/6.0)/(1.0-x/y-z/w);
  }

  std::complex<double> 
  _Amp_base::tLs_1_2mh(double s, double t, double m1, double m2) 
  {
    _ComplexD ls  = Log(-s);
    _ComplexD lst = ls - Log(-t);
    _ComplexD ls1 = ls - Log(-m1);
    _ComplexD ls2 = ls - Log(-m2);  

    return -Li2(-m1, -t) - Li2(-m2, -t) - 0.5*lst*lst + 0.5*ls1*ls2;
  }

  std::complex<double> 
  _Amp_base::Ls_1_2me(double s, double t, double m1, double m3)
  {
    _ComplexD lst = Log(-s)-Log(-t);  
    return -Li2(-m1, -s) - Li2(-m1, -t) - Li2(-m3, -s) 
      -     Li2(-m3, -t) + Li2(-m1, -m3, -s, -t) - 0.5*lst*lst;
  }

  std::complex<double>
  _Amp_base::Ls_1_2mh(double s, double t, double m1, double m2) {
    return tLs_1_2mh(s, t, m1, m2) + (0.5*(s-m1-m2)+m1*m2/t)*I3_3m(s,m1,m2);
  }

  std::complex<double> _Amp_base::I3m3(double s12, double s34, double s56) 
  {
    double rootD3, x, y, ro;
    double re, im = 0.0;

    double D3 = s12*s12+s34*s34+s56*s56-2.0*s12*s34-2.0*s34*s56-2.0*s56*s12;
    double d12 = s12-s34-s56, d34 = s34-s56-s12, d56 = s56-s12-s34;

    if(D3 < 0.0) {
      rootD3 = std::sqrt(-D3);
      re = 2.0*(  __specfunc_cl2(2.0*std::atan(rootD3/d12))
		+ __specfunc_cl2(2.0*std::atan(rootD3/d34))
	        + __specfunc_cl2(2.0*std::atan(rootD3/d56)))/rootD3;
    } else {
      rootD3 = std::sqrt(D3);
      x  = s12/s56;
      y  = s34/s56;
      ro = 2.0*s56/(d56+rootD3);
      
      re = -real(2.0*(Li2(-ro*x) + Li2(-ro*y)) + Log(ro*x)*Log(ro*y)
		 + Log(y/x)*Log((1.0+ro*y)/(1.0+ro*x)) + DEB_PI3)/rootD3;
      if(s34 > 0.0) {
	im  = std::log((d12+rootD3)*(d56+rootD3)/((d12-rootD3)*(d56-rootD3)));
	im *= -DEB_PI/rootD3;
      }
    }

    return _ComplexD(re, im);
  }

  std::complex<double> _Amp_base::I3_3m(double s1, double s2, double s3) 
  {
    bool do_conj = false;
    double tmp, min, max, thd, coef = 1.0;
    unsigned int npos = (s1 > 0.0 ? 1U : 0U) + (s2 > 0.0 ? 1U : 0U)
      + (s3 > 0.0 ? 1U : 0U);
    
    //---- when all the three arguments are negative ----
    if(npos == 0) 
      return I3m3(s1, s2, s3);
    
    //---- when all the three arguments are positive ----
    if(npos == 3) 
      return -I3m3(-s1, -s2, -s3);
    
    //---- when two arguments are positive
    if(npos == 2) {
      s1 *= -1.0; s2 *= -1.0; s3 *= -1.0; 
      coef = -1.0; do_conj = true;
    }
    
    tmp = (s1 > s2 ? s1 : s2); 
    max = (tmp > s3 ? tmp : s3);
    tmp = (s1 < s2 ? s1 : s2);
    min = (tmp < s3 ? tmp : s3);
    thd = s1 + s2 + s3 - max - min;
    _ComplexD retval = coef*I3m3(min, max, thd);
    
    return (do_conj ? conj(retval) : retval);
  }
  
  std::complex<double> 
  _Amp_base::f(double mt, double s12, double s34, double s56)
  {
    double D3 = s12*s12+s34*s34+s56*s56-2.0*s12*s34-2.0*s34*s56-2.0*s56*s12;
    double d12 = s12-s34-s56, d34 = s34-s56-s12, d56 = s56-s12-s34;
    
    if(mt != 0.0)
      return(_ComplexD(1.0/(24.0*mt*mt)+(2.0*s34+s12+s56)/(360.0*mt*mt*mt*mt), 0.0));
    else
      return((3.0*s34*d34/(D3*D3)-1.0/D3)*s12*s56*I3_3m(s12,s34,s56)
	     +(3.0*s56*d56/(D3*D3)-1.0/(2.0*D3))*s12*(Log(-s12)-Log(-s34))
	     +(3.0*s12*d12/(D3*D3)-1.0/(2.0*D3))*s56*(Log(-s56)-Log(-s34))
	     -d34/(2*D3));
  }
}    //  namespace nlo

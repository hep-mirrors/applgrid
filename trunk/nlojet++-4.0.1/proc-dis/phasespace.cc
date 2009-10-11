//  Copyright (C) 2003 Zoltan Nagy
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
#include "bits/dis-phasespace.h"


namespace nlo {
  
  void basic_phasespace_dis::
  phasespace_cuts(double q2min, double q2max, double xmin, double xmax, double ymin, double ymax)
  {
	double s = this -> _M_get_cm_energy_square();
	
	_M_q2min = q2min; _M_q2max = q2max; 
	_M_xmin  = xmin;  _M_xmax  = xmax; 
	_M_ymin  = ymin;  _M_ymax  = ymax;
	
	if(_M_q2min <= 0.0)  throw "Q2min must be greater than zero";
	if(_M_q2min > _M_q2max) throw "Q2min must be less than Q2max";
	if(_M_xmin > _M_xmax)   throw "xmin must be less than xmax";
	if(_M_ymin > _M_ymax)   throw "ymin must be less than ymax";
	
	if(_M_xmax == _M_xmin && _M_ymax == _M_ymin && _M_q2max == _M_q2min)
	  throw "at most two variables can be constrained";
	
	if(s*_M_xmax*_M_ymax < _M_q2min || s*_M_xmin*_M_ymin > _M_q2max)
	  throw "no phase space avaible";
	
	if(s*_M_xmin*_M_ymin > _M_q2min) _M_q2min = s*_M_xmin*_M_ymin;
	if(s*_M_xmax*_M_ymax < _M_q2max) _M_q2max = s*_M_xmax*_M_ymax;
  }


  double basic_phasespace_dis::operator()(event_dis& p)
  {
	double q2, y, x, yjac, qjac, s = this -> _M_get_cm_energy_square();
	
	if(_M_q2min == _M_q2max) {
	  q2 = _M_q2min;
	  qjac = 1.0;
	} else {
	  qjac  = std::log(_M_q2max/_M_q2min);
	  q2 = _M_q2min*std::exp((_M_rng -> operator()())*qjac);
	  qjac *= q2;
	}
	
	double yn = _M_ymin;
	double yx = _M_ymax;
	
	if(_M_ymin*_M_xmax*s < q2) yn = q2/(s*_M_xmax);
	if(_M_ymax*_M_xmin*s > q2) yx = q2/(s*_M_xmin);
	
	if(std::fabs(yn - yx) < 1e-10) {
	  y = yn;
	  yjac = 1.0;
	} else if(yn < yx) {
	  yjac = std::log(yx/yn);
	  y = yn*std::exp((_M_rng -> operator()())*yjac);
	} else throw "no phase space avaible";
	
	x = q2/(s*y);
	
	return qjac*yjac/s*(this -> _M_generate_event(x, y, p));
  }
  
}


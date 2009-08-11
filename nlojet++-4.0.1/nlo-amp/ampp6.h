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
#ifndef __AMPP6_H__
#define __AMPP6_H__ 1

//   nlo includes
#include <bits/amp-ampbase.h>


namespace nlo {
  
  
  class ampp6 : private _Amp_base
  {
	//   private types
	typedef std::complex<double> _ComplexD;
  
  public:
	//  constructors
	ampp6(const innerprod_type& __x, const random_generator& rng)
	  : _Amp_base(__x, rng) {}
  
	//  tree level matrix element
	_ComplexD matrix_element_1loop_mmpppp(int, int, int, int, int, int) const;
	
	_ComplexD matrix_element_1loop_ppmmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
	  swap();
	  _ComplexD res = matrix_element_1loop_mmpppp(p1,p2,p3,p4,p5,p6);
	  swap();
	  return res;
	}
  
  //private:
	//  private members 
	_ComplexD Lambda4(int, int, int, int, int, int) const; 
	_ComplexD Lambda5(int, int, int, int, int, int) const; 
	_ComplexD Lambda6(int, int, int, int, int, int) const; 
	_ComplexD Ammpppp(int, int, int, int, int, int) const; 
  };

}   //  namespace nlo

#endif

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
#include "ampp6.h"
#include "defmacros.h"


#define __NLO_64PI 201.0619298297467672616092


namespace nlo {
  
  std::complex<double> ampp6::Lambda4(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
	double q3 = S3(2,3,4), q4 = S(2,4), r3 = S(2,3), r4 = 0.0, s34 = S(3,4);
	_ComplexD l34 = Log(r3) - Log(q4);
	
	return (Li2(q3,r4, q4,r3) - Li2(q3,q4) - Li2(r4,r3) 
			- Li2(q3,r3) - Li2(r4,q4) - 0.5*l34*l34)/(s34*s34);
  }

  std::complex<double> ampp6::Lambda5(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
	double q3 = S(1,6), q4 = S3(2,4,5), r3 = S3(2,3,4), r4 = S(2,4), s35 = S(3,5);
	_ComplexD l34 = Log(r3) - Log(q4);
  
	return (Li2(q3,r4, q4,r3) - Li2(q3,q4) - Li2(r4,r3) 
			- Li2(q3,r3) - Li2(r4,q4) - 0.5*l34*l34)/(s35*s35);
  }

  std::complex<double> ampp6::Lambda6(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
	double q3 = 0.0, q4 = S(1,3), r3 = S(1,6), r4 = S3(2,4,5), s36 = S(3,6);
	_ComplexD l34 = Log(r3) - Log(q4);
  
	return (Li2(q3,r4, q4,r3) - Li2(q3,q4) - Li2(r4,r3) 
			- Li2(q3,r3) - Li2(r4,q4) - 0.5*l34*l34)/(s36*s36);
  }

  std::complex<double> ampp6::Ammpppp(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
	_ComplexD a1342 = A4(1,3,4,2), a1352 = A4(1,3,5,2), a1362 = A4(1,3,6,2);
	
	return (a1342*a1342*Lambda4(p1,p2,p3,p4,p5,p6) + a1352*a1352*Lambda5(p1,p2,p3,p4,p5,p6) 
			+ a1362*a1362*Lambda6(p1,p2,p3,p4,p5,p6))/(A(3,4)*A(4,5)*A(5,6)*A(6,3));
  }
  
  std::complex<double> 
  ampp6::matrix_element_1loop_mmpppp(int p1, int p2, int p3, int p4, int p5, int p6) const 
  {
	static const unsigned int perms[24][4] = {
	  {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, 
	  {0, 3, 1, 2}, {0, 3, 2, 1}, {1, 0, 2, 3}, {1, 0, 3, 2}, 
	  {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2}, {1, 3, 2, 0}, 
	  {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3}, {2, 1, 3, 0}, 
	  {2, 3, 0, 1}, {2, 3, 1, 0}, {3, 0, 1, 2}, {3, 0, 2, 1}, 
	  {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}
	};
	
	int p[4] = {p3,p4,p5,p6};
	_ComplexD res = 0.0;
	
	for(unsigned int i = 0; i < 24; i++)
	  res -= Ammpppp(p1,p2, p[perms[i][0]], p[perms[i][1]], p[perms[i][2]], p[perms[i][3]]);
	
	return __NLO_64PI*res;
  }
  
  
}


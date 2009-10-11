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
#include "ampq6.h"
#include "defmacros.h"


namespace nlo {


  const short ampq6::perm[36][6] = 
  {{0,1,2,3,4,5},{0,1,2,5,4,3},{0,3,2,1,4,5},{0,3,2,5,4,1},{0,5,2,3,4,1},
   {0,5,2,1,4,3},{0,1,4,3,2,5},{0,1,4,5,2,3},{0,3,4,1,2,5},{0,3,4,5,2,1},
   {0,5,4,3,2,1},{0,5,4,1,2,3},{2,1,0,3,4,5},{2,1,0,5,4,3},{2,3,0,1,4,5},
   {2,3,0,5,4,1},{2,5,0,3,4,1},{2,5,0,1,4,3},{2,1,4,3,0,5},{2,1,4,5,0,3},
   {2,3,4,1,0,5},{2,3,4,5,0,1},{2,5,4,3,0,1},{2,5,4,1,0,3},{4,1,2,3,0,5},
   {4,1,2,5,0,3},{4,3,2,1,0,5},{4,3,2,5,0,1},{4,5,2,3,0,1},{4,5,2,1,0,3},
   {4,1,0,3,2,5},{4,1,0,5,2,3},{4,3,0,1,2,5},{4,3,0,5,2,1},{4,5,0,3,2,1},
   {4,5,0,1,2,3}};

  const short ampq6::hel[36][7] = 
  {{ 6, 7, 8,11,12,13, 1},{ 4, 9, 8,11,10,15,-1},{ 6, 2, 3,16,17,13,-1},
   { 0, 9, 3,16,10,19, 1},{ 0, 7, 1,18,12,19,-1},{ 4, 2, 1,18,17,15, 1},
   { 9, 4, 8,11,15,10,-1},{ 7, 6, 8,11,13,12, 1},{ 9, 0, 3,16,19,10, 1},
   { 2, 6, 3,16,13,17,-1},{ 2, 4, 1,18,15,17, 1},{ 7, 0, 1,18,19,12,-1},
   { 6,16,17, 2, 3,13,-1},{ 4,18,17, 2, 1,15, 1},{ 6,11,12, 7, 8,13, 1},
   { 0,18,12, 7, 1,19,-1},{ 0,16,10, 9, 3,19, 1},{ 4,11,10, 9, 8,15,-1},
   {18, 4,17, 2,15, 1, 1},{16, 6,17, 2,13, 3,-1},{18, 0,12, 7,19, 1,-1},
   {11, 6,12, 7,13, 8, 1},{11, 4,10, 9,15, 8,-1},{16, 0,10, 9,19, 3, 1},
   {18, 7,19, 0,12, 1,-1},{16, 9,19, 0,10, 3, 1},{18, 2,15, 4,17, 1, 1},
   {11, 9,15, 4,10, 8,-1},{11, 7,13, 6,12, 8, 1},{16, 2,13, 6,17, 3,-1},
   { 9,16,19, 0, 3,10, 1},{ 7,18,19, 0, 1,12,-1},{ 9,11,15, 4, 8,10,-1},
   { 2,18,15, 4, 1,17, 1},{ 2,16,13, 6, 3,17,-1},{ 7,11,13, 6, 8,12, 1}};

  const short ampq6::iq6[5][36] =
  {{1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1},
   {1,0,1,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,1},
   {1,0,0,0,1,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,0,0,0,1},
   {1,1,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,0,0,1},
   {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
   
  const short ampq6::mapcol[4][36] =   
  {{5,1,2,0,3,4,2,4,5,3,0,1,3,0,4,1,5,2,4,2,3,5,1,0,1,5,0,2,4,3,0,3,1,4,2,5},
   {3,0,4,1,5,2,0,3,1,4,2,5,5,1,2,0,3,4,1,5,0,2,4,3,4,2,3,5,1,0,2,4,5,3,0,1},
   {1,5,0,2,4,3,4,2,3,5,1,0,0,3,1,4,2,5,2,4,5,3,0,1,5,1,2,0,3,4,3,0,4,1,5,2},
   {0,3,1,4,2,5,3,0,4,1,5,2,1,5,0,2,4,3,5,1,2,0,3,4,2,4,5,3,0,1,4,2,3,5,1,0}};

  class ampq6::_Colmat
  {
  public:
    //    constructor
    _Colmat();

    //    calculate the color sum
    double square(std::complex<double> *m) const;

  private:
    double _M_colmat[6][6];
  };

  ampq6::_Colmat::_Colmat()
  {
    const double c1 = Nc/16.0, c2 = Nc2/16.0, c3= Nc2*Nc/16.0; 
    double colfac[3] = {c1, c2, c3};
    short i, j, ord[6][6] = {{2,1,1,1,0,0},{1,2,0,0,1,1},{1,0,2,0,1,1},
			     {1,0,0,2,1,1},{0,1,1,1,2,0},{0,1,1,1,0,2}};

    for(i = 0; i < 6; i++)
      for(j = 0; j < 6; j++)
	_M_colmat[i][j] = colfac[ord[i][j]];
  }

  double ampq6::_Colmat::square(std::complex<double> *m) const
  {
    double ret_val = 0.0;
    std::complex<double> tmp;

    for(short i = 0; i < 6; i++) {
      tmp = 0.5*_M_colmat[i][i]*m[i];
      for(short j = i + 1; j < 6; j++)
        tmp += _M_colmat[i][j]*m[j];
      ret_val += 2.0*real(tmp*conj(m[i]));
    }
    
    return ret_val;
  }
  
  //    the color matrix
  const ampq6::_Colmat ampq6::_S_colmat;


#define PM1(i1,i2,i3,i4,i5,i6)						\
  (P1*A(i1,i3)*(A(i5,i1)*B(i2,i1)+A(i5,i3)*B(i2,i3))*B(i6,i4))

#define MP1(i1,i2,i3,i4,i5,i6)						\
    (P1*B(i2,i3)*(A(i1,i2)*B(i6,i2)+A(i1,i3)*B(i6,i3))*A(i5,i4))

#define DRIE(i1,i2,i3,i4,i5,i6)						\
    (P2*(A(i1,i5)*A(i1,i3)*B(i1,i2)*B(i6,i4)				\
	 + A(i1,i5)*A(i3,i5)*B(i2,i4)*B(i6,i5)				\
	 + A(i1,i3)*A(i3,i5)*B(i2,i6)*B(i4,i3)))


  void ampq6::
  diag(short ip, int q1, int q2, int q3, int q4, int q5, int q6, _ComplexD *H) const
  {
    int p[6] = {q1, q2, q3, q4, q5, q6};
    double sign = -4.0*hel[ip][6];
    
    int p1 = p[perm[ip][0]], p2 = p[perm[ip][1]], p3 = p[perm[ip][2]],
      p4 = p[perm[ip][3]], p5 = p[perm[ip][4]], p6 = p[perm[ip][5]];
    
    short h1 = hel[ip][0], h2 = hel[ip][1], h3 = hel[ip][2],
      h4 = hel[ip][3], h5 = hel[ip][4], h6 = hel[ip][5];
    
    double P1 = sign/(S3(1,2,3)*S(1,2)*S(5,6));
    
    H[5]  = PM1(1,2,3,4,5,6);
    H[h1] = (h1 < 10 ? PM1(1,2,3,4,6,5) : 0.0);
    H[h2] = (h2 < 10 ? MP1(1,2,3,4,5,6) : 0.0);
    H[h3] = (h3 < 10 ? MP1(1,2,3,4,6,5) : 0.0);
    H[h4] = (h4 < 10 ? PM1(2,1,3,4,5,6) : 0.0);
    H[h5] = (h5 < 10 ? PM1(2,1,3,4,6,5) : 0.0);
    H[h6] = (h6 < 10 ? MP1(2,1,3,4,5,6) : 0.0);
    
    if(p3 == 3) {
      double P2 = sign/(S(3,4)*S(1,2)*S(5,6));
      H[5] += DRIE(1,2,3,4,5,6);
      
      if(h1 < 10) H[h1] += DRIE(1,2,3,4,6,5);
      if(h2 < 10) H[h2] += DRIE(1,2,4,3,5,6);
      if(h3 < 10) H[h3] += DRIE(1,2,4,3,6,5);
      if(h4 < 10) H[h4] += DRIE(2,1,3,4,5,6);
      if(h5 < 10) H[h5] += DRIE(2,1,3,4,6,5);
      if(h6 < 10) H[h6] += DRIE(2,1,4,3,5,6);
    }
  }
  

  void ampq6::su3_tree(int p1, int p2, int p3, int p4, int p5, int p6, 
		       const char *which, double *out)
  {
    static short idiag[36], ido[5];
    static std::complex<double> m[36][10], q[6];
    short i, j, k, l;
    
    for(i = 0; i < 5; i++)
      ido[i] = (which[i] == '1' ? 1 : 0);
    
    if(ido[4] == 0) {
      for(i = 0; i < 36; i++) {
	idiag[i] = 0;
	for(j = 0; j < 5; j++)
	  idiag[i] += ido[j]*iq6[j][i];
      }
    }

    for(i = 0; i < 36; i++)
      if(idiag[i] != 0 || ido[4] == 1)
	diag(i, p1,p2,p3,p4,p5,p6, m[i]);

    for(i = 0; i < 5; i++) {
      out[i] = 0.0;
      if(ido[i] == 1)
	for(k = 0; k < 10; k++) {
	  for(l = 0; l < 6; l++)
	    q[l] = 0.0;
	  
	  for(j = 0; j < 36; j++)
	    if(iq6[i][j] == 1) {
	      q[mapcol[0][j]] += m[j][k];
	      q[mapcol[1][j]] -= m[j][k]/Nc;
	      q[mapcol[2][j]] -= m[j][k]/Nc;
	      q[mapcol[3][j]] += m[j][k]/Nc2;
	    }
	  
	  out[i] += 2.0*_S_colmat.square(q);
	}
    }
  }
}   //   namespace nlo

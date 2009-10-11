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
#include "ampq2g4.h"
#include "defmacros.h"

#define ICPLX std::complex<double>(0,1)

namespace nlo {

  
  const short ampq2g4::perm[4][24] = 
  {{0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3},
   {1,1,2,2,3,3,2,2,3,3,0,0,3,3,0,0,1,1,0,0,1,1,2,2},
   {2,3,3,1,1,2,3,0,0,2,2,3,0,1,1,3,3,0,1,2,2,0,0,1},
   {3,2,1,3,2,1,0,3,2,0,3,2,1,0,3,1,0,3,2,1,0,2,1,0}};

  const short ampq2g4::hel[6][24] = 
  {{1,1,2,2,3,3,4,4,5,5,1,1,6,6,2,2,4,4,3,3,5,5,6,6},
   {2,3,3,1,1,2,5,1,1,4,4,5,2,4,4,6,6,2,5,6,6,3,3,5},
   {3,2,1,3,2,1,1,5,4,1,5,4,4,2,6,4,2,6,6,5,3,6,5,3},
   {4,5,6,4,5,6,6,2,3,6,2,3,3,5,1,3,5,1,1,2,4,1,2,4},
   {5,4,4,6,6,5,2,6,6,3,3,2,5,3,3,1,1,5,2,1,1,4,4,2},
   {6,6,5,5,4,4,3,3,2,2,6,6,1,1,5,5,3,3,4,4,2,2,1,1}};
  
  
  class ampq2g4::_Colmat
  {
  public:
    //   constructor
    _Colmat();
    
    //   calculate the color sums
    static double Leading(std::complex<double> *);
    double SubLeading(std::complex<double> *) const;

  private:
    double _M_colmat[24][24];
  };

  ampq2g4::_Colmat::_Colmat() 
  {
    short i, j, ord1[6][24] = 
    {{-3,-1,0,-1,0,1,0,0,0,0,-1,0,1,1,0,0,0,1,0,0,1,0,1,0},
     {-1,-3,1,0,-1,0,0,0,0,0,0,-1,1,0,0,0,1,0,0,0,0,1,1,1},
     {0,1,-3,-1,0,-1,1,0,1,0,0,0,0,0,0,-1,0,0,0,0,1,1,1,0},
     {-1,0,-1,-3,1,0,0,1,1,1,0,0,0,0,-1,0,0,0,0,0,0,1,0,1},
     {0,-1,0,1,-3,-1,1,1,1,0,0,0,0,1,0,0,0,1,-1,0,0,0,0,0},
     {1,0,-1,0,-1,-3,0,1,0,1,0,0,1,0,0,0,1,1,0,-1,0,0,0,0}};
    
    short ord2[6][24] = 
    {{3,2,1,2,1,0,0,1,0,-1,2,1,-1,-2,1,0,-1,0,0,-1,-2,-1,-2,-3},
     {2,3,0,1,2,1,-1,0,1,0,1,2,-2,-3,0,-1,-2,-1,1,0,-1,0,-1,-2},
     {1,0,3,2,1,2,-2,-1,-2,-3,0,-1,1,0,1,2,-1,0,0,1,-2,-1,0,-1},
     {2,1,2,3,0,1,-1,0,-1,-2,1,0,0,-1,2,1,0,1,-1,0,-3,-2,-1,-2},
     {1,2,1,0,3,2,-2,-1,0,-1,0,1,-1,-2,-1,0,-3,-2,2,1,0,1,0,-1},
     {0,1,2,1,2,3,-3,-2,-1,-2,-1,0,0,-1,0,1,-2,-1,1,2,-1,0,1,0}};
    
    for(i = 0; i < 6; i++)
      for(j = 0; j < 24; j++)
	_M_colmat[i][j] = Nc2*Nc2*ord1[i][j] + Nc2*ord2[i][j] - 1;
    
    for(short irot = 1; irot <= 3; irot++)
      for(i = 0; i < 6; i++)
	for(j = 0; j < 24; j++)
	  _M_colmat[i+6*irot][j] = _M_colmat[i][(j+24-6*irot)%24];
  }

  double ampq2g4::_Colmat::Leading(std::complex<double> *df)
  {
    double ret_val = 0.0;
    for(short i = 0; i < 24; i++)
      ret_val += real(df[i]*conj(df[i]));
    
    return ret_val;
  }
 
  double ampq2g4::_Colmat::SubLeading(std::complex<double> *df) const
  {
    double ret_val = 0.0;
    std::complex<double> tmp;
    
    for(short i = 0; i < 24; i++) {
      tmp = 0.5*_M_colmat[i][i]*df[i];
      for(short j = i+1; j < 24; j++)
	tmp += _M_colmat[i][j]*df[j];
      ret_val += 2.0*real(tmp*conj(df[i]));
    }
    return ret_val;
  }

  //    the color matrix
  const ampq2g4::_Colmat ampq2g4::_S_colmat;


  std::complex<double> 
  ampq2g4::DPPMM(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    std::complex<double> Z1, Z3, Z4;
    std::complex<double> a43 = A(4,3), a31 = A(3,1), a23 = A(2,3), a26 = A(2,6), 
      b16 = B(1,6), b26 = B(2,6), b65 = B(6,5), b31 = B(3,1), 
      c4562 = C(4,5,6,2), c2561 = C(2,5,6,1), c1265 = C(1,2,6,5); 
     
    double s12 = S(1,2), s13 = S(1,3), s26 = S(2,6), s34 = S(3,4), 
      s45 = S(4,5), s56 = S(5,6), s126 = S3(1,2,6), s456 = S3(4,5,6);   
    
    Z1 = a43*a43*b16*b26*c1265*c1265;
    Z3 = a31*a23*b65*c4562*b65*c4562;
    Z4 = b26*b65*a43*a31
      *(c1265*c2561*c4562 - s34*a26*b65*c4562 - s56*a43*b31*c1265);
    
    return Z1/(s12*s126*s26*s45*s34) + Z3/(s45*s456*s56*s13*s12) 
      +    Z4/(s12*s26*s56*s45*s34*s13);
  }

  std::complex<double> 
  ampq2g4::DMMPP(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    std::complex<double> Z1, Z3, Z4;
    std::complex<double> a16 = A(1,6), a62 = A(6,2), a65 = A(6,5), b13 = B(1,3),
      b23 = B(2,3), b43 = B(4,3), c1564 = C(1,5,6,4), c1562 = C(1,5,6,2),  
      c5162 = C(5,1,6,2);  
    
    double s12 = S(1,2), s26 = S(2,6), s45 = S(4,5), s34 = S(3,4), 
      s56 = S(5,6), s13 = S(1,3), s126 = S3(1,2,6), s456 = S3(4,5,6);    
        
    Z1 = a16*a62*b43*c5162*b43*c5162;
    Z3 = a65*a65*b13*b23*c1564*c1564;
    Z4 = a62*a65*b43*b13*c5162*c1562*c1564;
    
    return Z1/(s12*s126*s26*s45*s34) + Z3/(s45*s456*s56*s13*s12) 
      +    Z4/(s12*s26*s56*s45*s34*s13);
  }

  std::complex<double> 
  ampq2g4::DPMPM(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    std::complex<double> Z1, Z2, Z3, Z4;
    std::complex<double> a53 = A(5,3), a13 = A(1,3), a32 = A(3,2), a43 = A(4,3),
      b16 = B(1,6), b26 = B(2,6), b64 = B(6,4), b65 = B(6,5), 
      c1264 = C(1,2,6,4), c3256 = C(3,2,5,6), c5264 = C(5,2,6,4),
      c5642 = C(5,6,4,2), c5462 = C(5,4,6,2);
    
    double s12 = S(1,2), s26 = S(2,6), s45 = S(4,5), s34 = S(3,4),
      s56 = S(5,6), s13 = S(1,3), s126 = S3(1,2,6), s256 = S3(2,5,6),
      s456 = S3(4,5,6);    
    
    Z1 = a53*a53*b16*b26*c1264*c1264;
    Z2 = a13*b26*c3256*c5264*c5264;
    Z3 = a13*a32*b64*c5642*b64*c5642;
    Z4 =-s126*a13*b64*c3256*c5264*c5462 - s256*a53*b64*c3256*c1264*c5462
      +  s456*a53*b26*c3256*c1264*c5264 - a53*a43*b65*b64*c1264*c5264*c5462;
    
    return Z1/(s12*s126*s26*s45*s34) + Z2/(s26*s256*s56*s34*s13)
      +    Z3/(s45*s456*s56*s13*s12) + Z4/(s12*s26*s56*s45*s34*s13);
  }

  std::complex<double> 
  ampq2g4::DMPPM(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    std::complex<double> Z1, Z2, Z3, Z4;
    std::complex<double> a14 = A(1,4), a43 = A(4,3), a54 = A(5,4), b16 = B(1,6),
      b26 = B(2,6), b13 = B(1,3), b23 = B(2,3), b12 = B(1,2), b63 = B(6,3),
      c1263 = C(1,2,6,3), c4256 = C(4,2,5,6), c5263 = C(5,2,6,3), 
      c1456 = C(1,4,5,6); 
    
    double s12 = S(1,2), s26 = S(2,6), s45 = S(4,5), s34 = S(3,4),
      s56 = S(5,6), s13 = S(1,3), s126 = S3(1,2,6), s256 = S3(2,5,6),
      s456 = S3(4,5,6);

    Z1 = a54*a54*b16*b26*c1263*c1263;
    Z2 = a14*b26*c4256*c5263*c5263;
    Z3 = a54*a54*b13*b23*c1456*c1456;
    Z4 =-s126*a14*a54*b16*b23*c5263*c1456 - s256*a54*a54*b16*b23*c1263*c1456
      +  s456*a14*a54*b16*b26*c1263*c5263 + a54*a43*b12*b63*c1263*c5263*c1456;
    
    return Z1/(s12*s126*s26*s45*s34) + Z2/(s26*s256*s56*s34*s13)
      +    Z3/(s45*s456*s56*s13*s12) + Z4/(s12*s26*s56*s45*s34*s13);
  }

  std::complex<double> 
  ampq2g4::DPMMP(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    std::complex<double> Z1, Z2, Z3, Z4;
    std::complex<double> a23 = A(2,3), a16 = A(1,6), a26 = A(2,6), a13 = A(1,3),
      a12 = A(1,2), a63 = A(6,3), b45 = B(4,5), b25 = B(2,5), b54 = B(5,4),
      b65 = B(6,5), c3162 = C(3,1,6,2), c3265 = C(3,2,6,5), c6452 = C(6,4,5,2),
      c6254 = C(6,2,5,4);

    double s12 = S(1,2), s26 = S(2,6), s45 = S(4,5), s34 = S(3,4), 
      s56 = S(5,6), s13 = S(1,3), s256 = S3(2,5,6), s126 = S3(1,2,6),
      s456 = S3(4,5,6);
    
    Z1 = -b45*b45*a16*a26*c3162*c3162;
    Z2 =  a13*b25*c3265*c6254*c6254;
    Z3 = -b45*b45*a13*a23*c6452*c6452;
    Z4 = s126*a13*a23*b25*b54*c6254*c6452 + s256*a16*a23*b45*b45*c3162*c6452
      -  s456*a16*a23*b25*b54*c3162*c6254 + a12*a63*b65*b54*c3162*c6254*c6452;
    
    return Z1/(s12*s126*s26*s45*s34) + Z2/(s26*s256*s56*s34*s13)
      +    Z3/(s45*s456*s56*s13*s12) + Z4/(s12*s26*s56*s45*s34*s13);
  }

  std::complex<double> 
  ampq2g4::DMPMP(int p1, int p2, int p3, int p4, int p5, int p6) const
  {
    std::complex<double> Z1, Z2, Z3, Z4;
    std::complex<double> a62 = A(6,2), a24 = A(2,4), a56 = A(5,6), a12 = A(1,2),
      a46 = A(4,6), a64 = A(6,4), a14 = A(1,4), a16 = A(1,6), a26 = A(2,6),
      b53 = B(5,3), b25 = B(2,5), b13 = B(1,3), b23 = B(2,3), b15 = B(1,5),
      b34 = B(3,4), b12 = B(1,2), b35 = B(3,5), c1352 = C(1,3,5,2), 
      c4352 = C(4,3,5,2), c4265 = C(4,2,6,5), c4162 = C(4,1,6,2), 
      c6253 = C(6,2,5,3), c1645 = C(1,6,4,5), c1465 = C(1,4,6,5); 
    
    double s13 = S(1,3), s26 = S(2,6), s12 = S(1,2), s45 = S(4,5), 
      s25 = S(2,5), s14 = S(1,4), s34 = S(3,4), s56 = S(5,6), 
      s256 = S3(2,5,6), s126 = S3(1,2,6), s456 = S3(4,5,6);
    
    Z1 = -a16*a26*b53*c4162*b53*c4162;
    Z2 =  a14*b25*c4265*c6253*c6253;
    Z3 =  a64*a64*b13*b23*c1645*c1645;
    Z4 = b15*b23*b25*b34*a12*a26*a14*a46*c4162
      -  b12*b13*b25*b35*a14*a16*a24*a56*c1465
      +  s13*s26*b35*a46*(b23*a46*c1465 - b35*a46*c1352 + b35*a16*c4352)
      +  s12*b35*b25*a46*a14
      *   ((s13-s45+s26)*c6253 - s25*b13*a16 + s14*b23*a26 - b13*a62*c1352);
    
    return Z1/(s12*s126*s26*s45*s34) + Z2/(s26*s256*s56*s34*s13)
      +    Z3/(s45*s456*s56*s13*s12) + Z4/(s12*s26*s56*s45*s34*s13);
  }


  double ampq2g4::su3_tree(int p1, int p2, int q3, int q4, int q5, int q6)
  {
    static const double coll = 2.0*Na*Nc2*Nc;
    static const double cols = 2.0*Na/(Nc2*Nc);
    static std::complex<double> df[7][24];

    double ret_val, num, s1i, s2i;
    int i, p3, p4, p5, p6, p[4]={q3,q4,q5,q6};

    for(i = 0; i < 24; i++) {
      p3 = p[perm[0][i]];
      p4 = p[perm[1][i]];
      p5 = p[perm[2][i]];
      p6 = p[perm[3][i]];
 
      df[0][i] = 1.0/(B(1,3)*B(3,4)*B(4,5)*B(5,6)*B(6,2));

      df[hel[0][i]][i] = DPPMM(p1,p2,p3,p4,p5,p6);
      df[hel[1][i]][i] = DPMPM(p1,p2,p3,p4,p5,p6);
      df[hel[2][i]][i] = DPMMP(p1,p2,p3,p4,p5,p6);
      df[hel[3][i]][i] = DMPPM(p1,p2,p3,p4,p5,p6);
      df[hel[4][i]][i] = DMPMP(p1,p2,p3,p4,p5,p6);
      df[hel[5][i]][i] = DMMPP(p1,p2,p3,p4,p5,p6);
   }

    s1i = Sij(p1,q3), s2i = Sij(p2,q3); num  = s1i*s2i*(s1i*s1i+s2i*s2i);
    s1i = Sij(p1,q4), s2i = Sij(p2,q4); num += s1i*s2i*(s1i*s1i+s2i*s2i);
    s1i = Sij(p1,q5), s2i = Sij(p2,q5); num += s1i*s2i*(s1i*s1i+s2i*s2i);
    s1i = Sij(p1,q6), s2i = Sij(p2,q6); num += s1i*s2i*(s1i*s1i+s2i*s2i);
    num /= S(1,2);
    
    ret_val = num*(coll*_Colmat::Leading(df[0]) + cols*_S_colmat.SubLeading(df[0]));
    for(i = 1; i <= 6; i++)
      ret_val += coll*_Colmat::Leading(df[i]) + cols*_S_colmat.SubLeading(df[i]);
    
    return ret_val;
  }
}  //  namespace nlo

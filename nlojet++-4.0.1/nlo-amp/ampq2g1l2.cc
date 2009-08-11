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
#include "ampq2g1l2.h"
#include "defmacros.h"



namespace nlo {
  
  std::complex<double> 
  ampq2g1l2::Atree1ppm(int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD a34 = A(3,4);
    return _ComplexD(0.0,-1.0)*a34*a34/(A(1,2)*A(2,3)*A(4,5));
  }

  std::complex<double> 
  ampq2g1l2::Fcc1ppm(int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD  a34 = A(3,4);
    double s23 = S(2,3), s45 = S(4,5);

    _ComplexD Fcc = a34*a34/(A(1,2)*A(2,3)*A(4,5))
      *(Ls_1(-S(1,2), -s45, -s23, -s45) 
	- 2.0*A4(3,1,5,4)/a34*L0(-s23, -s45)/s45);
    
    return Fcc;
  }

  std::complex<double> 
  ampq2g1l2::Fsc1ppm(int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD  a3154 = A4(3,1,5,4);
    double s23 = S(2,3), s45 = S(4,5);
    
    _ComplexD Fsc = a3154/(A(1,2)*A(2,3)*A(4,5))
      *(A(3,4)*L0(-s23, -s45)/s45 + 0.5*a3154*L1(-s23, -s45)/(s45*s45));

    return Fsc;
  }

  std::complex<double> 
  ampq2g1l2::Atree2pmp(int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD a24 = A(2,4);
    return _ComplexD(0.0,1.0)*a24*a24/(A(2,3)*A(3,1)*A(4,5));
  }

  std::complex<double> 
  ampq2g1l2::Fcc2pmp(int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD a24 = A(2,4), a23 = A(2,3), a31 = A(3,1), a45 = A(4,5),
      a12 = A(1,2), a34 = A(3,4), a14 = A(1,4), a13 = A(1,3), b13 = B(1,3); 

    double s12 = S(1,2), s45 = S(4,5), s13 = S(1,3), s23 = S(2,3); 

    _ComplexD Fcc = -a24*a24/(a23*a31*a45)*Ls_1(-s12, -s45, -s13, -s45)
      + a24*(a12*a34 - a14*a23)/(a23*a13*a13*a45)*Ls_1(-s12, -s45, -s23, -s45)
      + 2.0*b13*a14*a24/(a13*a45)*L0(-s23, -s45)/s45;

    return Fcc;
  }

  std::complex<double> 
  ampq2g1l2::Fsc2pmp(int p1, int p2, int p3, int p4, int p5) const 
  {
    double s12 = S(1,2), s45 = S(4,5), s23 = S(2,3);
    _ComplexD  a14 = A(1,4), a23 = A(2,3), a13 = A(1,3), a45 = A(4,5);

    _ComplexD a413 = A3(4,1,3), a231 = A3(2,3,1), a213 = A3(2,1,3),
      a435 = A3(4,3,5), a2134 = A4(2,1,3,4);

    _ComplexD b35 = B(3,5), b13 = B(1,3), b25 = B(2,5), b23 = B(2,3),
      b15 = B(1,5), b12 = B(1,2), b45 = B(4,5); 

    _ComplexD Fsc = a14*a14*a23/(a13*a13*a13*a45)*Ls_1(-s12, -s45, -s23, -s45)
      - 0.5*a413*a413*a23/(a13*a45*s23*s23)*L1(-s45, -s23)
      + a14*a14*a231/(a13*a13*a45*s23)*L0(-s45, -s23)
      - a213*a435/(a13*s12*s12)*L1(-s45, -s12)
      - a2134*a14/(a13*a13*a45*s12)*L0(-s45, -s12)
      - 0.5*b35*(b13*b25 + b23*b15)/(b12*b23*a13*b45);

    return Fsc;
  }

  std::complex<double> 
  ampq2g1l2::Fax2pmp(int p1, int p2, int p3, int p4, int p5) const
  {
    double s45 = S(4,5);
    
    return -B(5,3)*B(3,1)*A(2,4)
      * (L1(-S(1,2), -s45)/(s45*s45) - 1.0/(12*s45*mtsq));
  }


  void ampq2g1l2::
  matrix_1loop_ppm(int p1, int p2, int p3, int p5, int p6, _ComplexD *M) const 
  {
    _ComplexD  I(0.0, 1.0), a123, a132, atree;
    _ComplexD  V1 = -1.5*(Log(-S(5,6))-Log(-S(2,3)))-3.0;
    double V2 = -3.5;

    a132 = V2*Atree2pmp(p1, p3, p2, p5, p6)
      + I*(Fcc2pmp(p1, p3, p2, p5, p6) + Fsc2pmp(p1, p3, p2, p5, p6));

    atree = Atree1ppm(p1, p2, p3, p5, p6);
    a123  = V1*atree 
      + I*(Fcc1ppm(p1, p2, p3, p5, p6) + Fsc1ppm(p1, p2, p3, p5, p6));

    M[0] = atree;
    M[1] = a123+a132/Nc2;
  }

  void ampq2g1l2::
  matrix_1loop_pmm(int p1, int p2, int p3, int p5, int p6, _ComplexD *M) const 
  {
    _ComplexD  I(0.0, 1.0), a123, a132, atree;
    _ComplexD  V1 = -1.5*(Log(-S(5,6))-Log(-S(1,2)))-3.0;
    double V2 = -3.5;

    a132 = V2*Atree2pmm(p1, p3, p2, p5, p6)
      + I*(Fcc2pmm(p1, p3, p2, p5, p6) + Fsc2pmm(p1, p3, p2, p5, p6));

    atree = Atree1pmm(p1, p2, p3, p5, p6);
    a123 = V1*atree
      + I*(Fcc1pmm(p1, p2, p3, p5, p6) + Fsc1pmm(p1, p2, p3, p5, p6));

    M[0] = atree;
    M[1] = a123+a132/Nc2;
  }

  void ampq2g1l2::
  color_1loop_ppm(int p1, int p2, int p3, int p5, int p6, _ComplexD *M) const 
  {
    _ComplexD  I(0.0, 1.0), a123, a132, atree, ax132;
    _ComplexD  V1 = -1.5*(Log(-S(5,6))-Log(-S(2,3)))-3.0;
    double V2 = -3.5;

    a132 = V2*Atree2pmp(p1, p3, p2, p5, p6)
      + I*(Fcc2pmp(p1, p3, p2, p5, p6) + Fsc2pmp(p1, p3, p2, p5, p6));

    atree = Atree1ppm(p1, p2, p3, p5, p6);
    a123  = V1*atree 
      + I*(Fcc1ppm(p1, p2, p3, p5, p6) + Fsc1ppm(p1, p2, p3, p5, p6));

    ax132 = I*Fax2pmp(p1, p3, p2, p5, p6);
    
    M[0] = atree;
    M[1] = -a132;
    M[2] = a123 + a132;
    M[3] = ax132;
  }

  void ampq2g1l2::
  color_1loop_pmm(int p1, int p2, int p3, int p5, int p6, _ComplexD *M) const 
  {
    _ComplexD  I(0.0, 1.0), a123, a132, atree, ax132;
    _ComplexD  V1 = -1.5*(Log(-S(5,6))-Log(-S(1,2)))-3.0;
    double V2 = -3.5;
    
    a132 = V2*Atree2pmm(p1, p3, p2, p5, p6)
      + I*(Fcc2pmm(p1, p3, p2, p5, p6) + Fsc2pmm(p1, p3, p2, p5, p6));
    
    atree = Atree1pmm(p1, p2, p3, p5, p6);
    a123 = V1*atree
      + I*(Fcc1pmm(p1, p2, p3, p5, p6) + Fsc1pmm(p1, p2, p3, p5, p6));
    
    ax132 = I*Fax2pmm(p1, p3, p2, p5, p6);
    
    M[0] = atree;
    M[1] = - a132;
    M[2] = a123 + a132;
    M[3] = ax132;
  }
  


  //  TREE LEVEL AMPLITUDES
  double ampq2g1l2::su3_tree(int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD  P, M;
    double A;

    matrix_tree_ppm(p1, p2, p3, p4, p5, &P);
    matrix_tree_pmm(p1, p2, p3, p4, p5, &M);
    A = real(P*conj(P)+M*conj(M));

    matrix_tree_ppm(p1, p2, p3, p5, p4, &P);
    matrix_tree_pmm(p1, p2, p3, p5, p4, &M);
    A += real(P*conj(P)+M*conj(M));

    return 2.0*Na*A;
  }

  //  ONE-LOOP LEVEL AMPLITUDES
  double ampq2g1l2::
  su3_1loop(unsigned int nf, int p1, int p2, int p3, int p4, int p5) const
  {
    _ComplexD  P[2], M[2];

    matrix_1loop_ppm(p1, p2, p3, p4, p5, P);
    matrix_1loop_pmm(p1, p2, p3, p4, p5, M);
    double L = real(P[1]*conj(P[0]) + M[1]*conj(M[0]));

    matrix_1loop_ppm(p1, p2, p3, p5, p4, P);
    matrix_1loop_pmm(p1, p2, p3, p5, p4, M);
    L += real(P[1]*conj(P[0]) + M[1]*conj(M[0]));

    //-----the opposit helicity states (e+e- : means factor 2) -----
    swap();
    matrix_1loop_ppm(p1, p2, p3, p4, p5, P);
    matrix_1loop_pmm(p1, p2, p3, p4, p5, M);
    L += real(P[1]*conj(P[0]) + M[1]*conj(M[0]));
    
    matrix_1loop_ppm(p1, p2, p3, p5, p4, P);
    matrix_1loop_pmm(p1, p2, p3, p5, p4, M);
    L += real(P[1]*conj(P[0]) + M[1]*conj(M[0]));
    swap();

    return Na*Nc*L;
  }

  double ampq2g1l2::su3_tree_mch(int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD M;
    int hel = (int) (4*_M_rng());
    
    switch(hel) {
    case 0: matrix_tree_ppm(p1, p2, p3, p4, p5, &M); break;
    case 1: matrix_tree_pmm(p1, p2, p3, p4, p5, &M); break;
    case 2: matrix_tree_ppm(p1, p2, p3, p5, p4, &M); break;
    case 3: matrix_tree_pmm(p1, p2, p3, p5, p4, &M); break;
    }
    
    return 8.0*Na*real(M*conj(M));
  }
  
  double ampq2g1l2::
  su3_1loop_mch(unsigned int nf, int p1, int p2, int p3, int p4, int p5) const 
  {
    _ComplexD  M[2];
    int hpm = (int) (2*_M_rng());
    int hel = (int) (4*_M_rng());
    
    if(hpm == 1) swap();
    switch(hel) {
    case 0: matrix_1loop_ppm(p1, p2, p3, p4, p5, M); break;
    case 1: matrix_1loop_pmm(p1, p2, p3, p4, p5, M); break;
    case 2: matrix_1loop_ppm(p1, p2, p3, p5, p4, M); break;
    case 3: matrix_1loop_pmm(p1, p2, p3, p5, p4, M); break;
    }
    if(hpm == 1) swap();
       
    return 8.0*Na*Nc*real(M[1]*conj(M[0]));
  }


  //  COLOR CORRALETED AMPLITUDES
#define Cond13(i,j)			        	\
(i == p1 && j == p3) || (i == p3 && j == p1)

#define Cond12(i,j)			        	\
(i == p1 && j == p2) || (i == p2 && j == p1) ||		\
(i == p3 && j == p2) || (i == p2 && j == p3)


  std::pair<double, std::complex<double> >
  ampq2g1l2::su3_cc(int i, int j, int p1, int p2, int p3, int p4, int p5) const 
  {
    double amp, s = 0.0;
    _ComplexD hda(0.0), P, M;;
    
    if(Cond13(i,j))      s = 1.0/Nc;
    else if(Cond12(i,j)) s = -Nc;
    else throw("Error in ampq2g1l2::su3_cc(...)");
    
    matrix_tree_ppm(p1, p2, p3, p4, p5, &P);
    matrix_tree_pmm(p1, p2, p3, p4, p5, &M);
    amp = real(P*conj(P) + M*conj(M));
    if(j == p2) hda = P*conj(M);

    matrix_tree_ppm(p1, p2, p3, p5, p4, &P);
    matrix_tree_pmm(p1, p2, p3, p5, p4, &M);
    amp += real(P*conj(P) + M*conj(M));
    if(j == p2) hda += P*conj(M);

    return _Pair(s*Na*amp, s*Na*hda);
  }

#define X(i,j) std::log(std::abs(s/S(i,j)))
  
  void ampq2g1l2::
  su3_kp(unsigned int nf, int pa, int p1, int p2,
	 int p3, int p4, int p5, su3_kp_i1 *res, double al) const 
  {
    double amp, cc12, cc13, xq = Gq/Cf, xg = Gg(nf)/Ca, s = S(4,5);
    _ComplexD P, M;
    
    matrix_tree_ppm(p1, p2, p3, p4, p5, &P);
    matrix_tree_pmm(p1, p2, p3, p4, p5, &M);
    amp = real(P*conj(P) + M*conj(M));
    
    matrix_tree_ppm(p1, p2, p3, p5, p4, &P);
    matrix_tree_pmm(p1, p2, p3, p5, p4, &M);
    amp += real(P*conj(P) + M*conj(M));
    
    cc12 = -Nc*Na*amp, cc13 = Na*amp/Nc;
    res->tree = 2.0*Na*amp;
    
    if(pa==p1 || pa==p3) res->ga = xg*cc12+xq*cc13;
    else if(pa==p2)      res->ga = 2.0*xq*cc12;
    else throw "Error in ampq2g1l2::su3_kp(...)";
    
    if(pa==p1) res->pa = (X(1,2)*cc12+X(1,3)*cc13)/Cf;
    if(pa==p2) res->pa = (X(1,2)+X(2,3))*cc12/Nc;
    if(pa==p3) res->pa = (X(2,3)*cc12+X(1,3)*cc13)/Cf;

    //   1-loop contribution (log terms)
    double l12 = Xq(S(1,2),s)+Xq(S(2,3),s)+Xg(S(1,2),s,nf)+Xg(S(2,3),s,nf);
    double l13 = 2.0*Xq(S(1,3),s);

    res->loop = (Gg(nf)+Kg(nf,al) + 2.0*(Gq+Kq(al)) - Cf)*(res->tree) 
      +         l13*cc13 + l12*cc12;
  }

  void ampq2g1l2::
  su3_kp_mch(unsigned int nf, int pa, int p1, int p2, 
	     int p3, int p4, int p5, su3_kp_i1 *res, double al) const 
  {
    double xq = Gq/Cf, xg = Gg(nf)/Ca, s = S(4,5);
    int hel = (int) (4*_M_rng());
    _ComplexD M;
    
    switch(hel) {
    case 0: matrix_tree_ppm(p1, p2, p3, p4, p5, &M); break;
    case 1: matrix_tree_pmm(p1, p2, p3, p4, p5, &M); break;
    case 2: matrix_tree_ppm(p1, p2, p3, p5, p4, &M); break;
    case 3: matrix_tree_pmm(p1, p2, p3, p5, p4, &M); break;
    }
    
    double amp = 4.0*real(M*conj(M));
    double cc12 = -Nc*Na*amp, cc13 = Na*amp/Nc;
    
    res->tree = 2.0*Na*amp;
    
    if(pa==p1 || pa==p3) res->ga = xg*cc12+xq*cc13;
    else if(pa==p2)      res->ga = 2.0*xq*cc12;
    else throw "Error in ampq2g1l2::su3_kp_mch(...)";
    
    if(pa==p1) res->pa = (X(1,2)*cc12+X(1,3)*cc13)/Cf;
    if(pa==p2) res->pa = (X(1,2)+X(2,3))*cc12/Nc;
    if(pa==p3) res->pa = (X(2,3)*cc12+X(1,3)*cc13)/Cf;

    //   1-loop contribution (log terms)
    double l12 = Xq(S(1,2),s)+Xq(S(2,3),s)+Xg(S(1,2),s,nf)+Xg(S(2,3),s,nf);
    double l13 = 2.0*Xq(S(1,3),s);
    
    res->loop = (Gg(nf)+Kg(nf,al) + 2.0*(Gq+Kq(al)) - Cf)*(res->tree) 
      +         l13*cc13 + l12*cc12;
  }

  double ampq2g1l2::
  su3_ins(unsigned int nf, int p1, int p2,
	  int p3, int p4, int p5, double al) const 
  {
    double amp, tree, cc12, cc13, s = S(4,5);
    _ComplexD P, M;
    
    matrix_tree_ppm(p1, p2, p3, p4, p5, &P);
    matrix_tree_pmm(p1, p2, p3, p4, p5, &M);
    amp = real(P*conj(P) + M*conj(M));
    
    matrix_tree_ppm(p1, p2, p3, p5, p4, &P);
    matrix_tree_pmm(p1, p2, p3, p5, p4, &M);
    amp += real(P*conj(P) + M*conj(M));
    
    cc12 = -Nc*Na*amp, cc13 = Na*amp/Nc;
    tree = 2.0*Na*amp;
 
   //   1-loop contribution (log terms)
    double l12 = Xq(S(1,2),s)+Xq(S(2,3),s)+Xg(S(1,2),s,nf)+Xg(S(2,3),s,nf);
    double l13 = 2.0*Xq(S(1,3),s);

    return (Gg(nf)+Kg(nf,al)+2.0*(Gq+Kq(al)) - Cf)*tree + l13*cc13+l12*cc12;
  }

  double ampq2g1l2::
  su3_ins_mch(unsigned int nf, int p1, int p2, 
	      int p3, int p4, int p5, double al) const 
  {
    double tree, s = S(4,5);
    int hel = (int) (4*_M_rng());
    _ComplexD M;
    
    switch(hel) {
    case 0: matrix_tree_ppm(p1, p2, p3, p4, p5, &M); break;
    case 1: matrix_tree_pmm(p1, p2, p3, p4, p5, &M); break;
    case 2: matrix_tree_ppm(p1, p2, p3, p5, p4, &M); break;
    case 3: matrix_tree_pmm(p1, p2, p3, p5, p4, &M); break;
    }
    
    double amp = 4.0*real(M*conj(M));
    double cc12 = -Nc*Na*amp, cc13 = Na*amp/Nc;
    
    tree = 2.0*Na*amp;
    
    //   1-loop contribution (log terms)
    double l12 = Xq(S(1,2),s)+Xq(S(2,3),s)+Xg(S(1,2),s,nf)+Xg(S(2,3),s,nf);
    double l13 = 2.0*Xq(S(1,3),s);
    
    return (Gg(nf)+Kg(nf,al)+2.0*(Gq+Kq(al)) - Cf)*tree + l13*cc13+l12*cc12;
  }
}   //  namespace nlo

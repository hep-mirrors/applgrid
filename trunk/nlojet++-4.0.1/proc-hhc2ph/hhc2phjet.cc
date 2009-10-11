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
#include <cstring>

//   nlo includes
#include "bits/proc-hhc2ph.h"
#include "bits/nlo-color.h"
#include "bits/nlo-flkern.h"

#include "ampq2g1p2.h"
#include "ampq2g2p2.h"
#include "ampq4p2.h"

#define Qu  0.66666666666666666666 
#define Qd -0.33333333333333333333
#define Qu4 0.19753086419753086419
#define Qd4 0.01234567901234567901


namespace nlo {



  void _hhc2ph_jet_base::amp_tree(ampq2g1p2 *amp, double *out)
  {
    double tmp;
    std::memset(out, 0, 17*sizeof(double));

    //  gg -> 
    
    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    out[1] = Qu4*(tmp = -(amp -> su3_tree(-1,1, 0, -2,-3))/(Nc*Na)); 
    out[2] = Qd4*tmp;
    
    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    out[3] = Qu4*(tmp = -(amp -> su3_tree(0,1, -1, -2,-3))/(Nc*Na));
    out[4] = Qd4*tmp;
    
    //  uU ->
    //  uu ->
    //  dD ->
    //  dd ->
    //  ud ->
    //  du ->
    //  uUb ->
    
    //  uub -> g   (1/Nc * 1/Nc)
    out[12] = Qu4*(tmp = (amp -> su3_tree(-1,0, 1, -2,-3))/Nc2);

    //  dDb ->
    //  ddb -> g   (1/Nc * 1/Nc)
    out[14] = Qd4*tmp;

    //  udb ->
    //  dub ->
  }
 
  std::pair<double, std::complex<double> > _hhc2ph_jet_base::
  amp_ccag(ampq2g1p2 *amp, int i, int j) 
  {
    _Pair cc(amp -> su3_cc(i,j, 1, -1, 0, -2,-3));
    cc.first /= -Na*Nc;
    cc.second /= -Na*Nc;
    return cc;
  }

  std::pair<double, std::complex<double> > _hhc2ph_jet_base::
  amp_ccga(ampq2g1p2 *amp, int i, int j) 
  {
    _Pair cc(amp -> su3_cc(i,j, 1, 0, -1, -2,-3));
    cc.first /= -Na*Nc;
    cc.second /= -Na*Nc;
    return cc;
  }

  std::pair<double, std::complex<double> > _hhc2ph_jet_base::
  amp_ccqa(ampq2g1p2 *amp, int i, int j) 
  {
    _Pair cc(amp -> su3_cc(i,j, -1, 0, 1, -2,-3));
    cc.first /= Nc2;
    cc.second /= Nc2;
    return cc;
  }

  void _hhc2ph_jet_base::amp_kp(double al, ampq2g1p2 *amp, su3_kp_i2 *out)
  {
    static su3_kp_i2 tmp;

    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    amp -> su3_kp(Nf, -1,0, -1,1, 0, -2,-3, &tmp, al);
    out[0] = (-Qu4/(Nc*Na))*tmp; 
    out[1] = (-Qd4/(Nc*Na))*tmp;
    
    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    amp -> su3_kp(Nf, -1,0, 0,1, -1, -2,-3, &tmp, al);
    out[2] = (-Qu4/(Nc*Na))*tmp;
    out[3] = (-Qd4/(Nc*Na))*tmp;
    
    //  uub -> g   (1/Nc * 1/Nc)
    //  ddb -> g   (1/Nc * 1/Nc)
    amp -> su3_kp(Nf, -1,0, -1,0, 1, -2,-3, &tmp, al);
    out[4] = (Qu4/Nc2)*tmp;
    out[5] = (Qd4/Nc2)*tmp;
  }

  void _hhc2ph_jet_base::amp_1loop(ampq2g1p2 *amp, double *out)
  {
    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    amp -> su3_1loop(Nu, Nd, -1,1, 0, -2,-3, out);
    out[0] /= Nc*Na; 
    out[1] /= Nc*Na;
    
    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    amp -> su3_1loop(Nu, Nd, 0,1, -1, -2,-3, out+2);
    out[2] /= Nc*Na; 
    out[3] /= Nc*Na;
    
    //  uub -> g   (1/Nc * 1/Nc)
    //  ddb -> g   (1/Nc * 1/Nc)
    amp -> su3_1loop(Nu, Nd, -1,0, 1, -2,-3, out+4);
    out[4] /= Nc2;
    out[5] /= Nc2;
  }

  void _hhc2ph_jet_base::amp_1loop_mch(ampq2g1p2 *amp, double *out)
  {
    //  ug -> u    (1/Nc * 1/Na)
    //  dg -> d    (1/Nc * 1/Na)
    amp -> su3_1loop_mch(Nu, Nd, -1,1, 0, -2,-3, out);
    out[0] /= Nc*Na; 
    out[1] /= Nc*Na;
    
    //  gu -> u    (1/Na * 1/Nc)
    //  gd -> d    (1/Na * 1/Nc)
    amp -> su3_1loop_mch(Nu, Nd, 0,1, -1, -2,-3, out+2);
    out[2] /= Nc*Na; 
    out[3] /= Nc*Na;
    
    //  uub -> g   (1/Nc * 1/Nc)
    //  ddb -> g   (1/Nc * 1/Nc)
    amp -> su3_1loop_mch(Nu, Nd, -1,0, 1, -2,-3, out+4);
    out[4] /= Nc2;
    out[5] /= Nc2;
  }


  void _hhc2ph_jet_base::amp_tree(ampq2g2p2 *amp1, ampq4p2 *amp2, double *out)
  {
    double tmp[2];

    //          TWO QUARKS TWO GLUONS SUBPROCESSES
    //  gg -> qqb   (1/Na * 1/Na)
    out[0] = (Nu*Qu4 + Nd*Qd4)*(amp1 -> su3_tree(2,1, -1, 0, -2,-3))/Na2;
    
    //  ug -> ug    (1/Nc * 1/Na)
    //  dg -> dg    (1/Nc * 1/Na)
    out[1] = Qu4*(tmp[0] = (amp1 -> su3_tree(-1,1, 0, 2, -2,-3))/(Nc*Na)); 
    out[2] = Qd4*tmp[0];
    
    //  gu -> ug    (1/Na * 1/Nc)
    //  gd -> dg    (1/Na * 1/Nc)
    out[3] = Qu4*(tmp[0] = (amp1 -> su3_tree(0,1, -1, 2, -2,-3))/(Nc*Na));
    out[4] = Qd4*tmp[0];
    
    //  uub -> gg   (1/Nc * 1/Nc * 1/2)
    out[12] = Qu4*(tmp[0] = 0.5*(amp1 -> su3_tree(-1,0, 1, 2, -2,-3))/Nc2);

    //  ddb -> gg   (1/Nc * 1/Nc * 1/2)
    out[14] = Qd4*tmp[0];

 
    //          FOUR QUARKS SUBPROCESSES
    //  uU -> uU    (1/Nc * 1/Nc)
    //  uu -> uu    (1/Nc * 1/Nc * 1/2)
    //  dD -> dD    (1/Nc * 1/Nc)
    //  dd -> dd    (1/Nc * 1/Nc * 1/2)
    amp2 -> su3_tree(Qu, -1,1, 0, 2, -2,-3, tmp);
    out[5] = tmp[0]/Nc2;
    out[6] = tmp[1]/(2.0*Nc2);
    out[7] = out[5]/16.0;
    out[8] = out[6]/16.0;
    
    //  ud -> ud    (1/Nc * 1/Nc)
    out[9] = (amp2 -> su3_tree(Qu, Qd, -1,1, 0,2, -2,-3))/Nc2;

    //  du -> du    (1/Nc * 1/Nc)
    out[10] = (amp2 -> su3_tree(Qd, Qu, -1,1, 0, 2, -2,-3))/Nc2;

    //  uUb -> uUb    (1/Nc * 1/Nc)
    //  uub -> uub    (1/Nc * 1/Nc)
    //  dDb -> dDb    (1/Nc * 1/Nc)
    //  ddb -> ddb    (1/Nc * 1/Nc)
    amp2 -> su3_tree(Qu, -1,1, 2, 0, -2,-3, tmp);
    out[11]  = tmp[0]/Nc2;
    out[12] += tmp[1]/Nc2;
    out[13]  = out[11]/16.0;
    out[14] += tmp[1]/(16.0*Nc2);

    //  uub -> UUb    (1/Nc * 1/Nc * (nu-1)) 
    //  uub -> ddb    (1/Nc * 1/Nc * nd) 
    tmp[0] = (amp2 -> su3_tree(Qu, Qu, -1,0, 2, 1, -2,-3))/Nc2;
    out[12] += (Nu-1)*tmp[0]
      + Nd*(amp2 -> su3_tree(Qu, Qd, -1,0, 2, 1, -2,-3))/Nc2;

    //  ddb -> DDb    (1/Nc * 1/Nc * (nd-1)) 
    //  ddb -> uub    (1/Nc * 1/Nc * nu) 
    out[14] += (Nd-1)*tmp[0]/16.0
      + Nu*(amp2 -> su3_tree(Qd, Qu, -1,0, 2, 1, -2,-3))/Nc2;

    //  udb -> udb    (1/Nc * 1/Nc)
    out[15] = (amp2 -> su3_tree(Qu, Qd, -1,1, 2, 0, -2,-3))/Nc2;
    
    //  dub -> dub    (1/Nc * 1/Nc)
    out[16] = (amp2 -> su3_tree(Qd, Qu, -1,1, 2, 0, -2,-3))/Nc2;
  }
  
  
  void _hhc2ph_jet_base::
  amp_tree_mch(ampq2g2p2 *amp1, ampq4p2 *amp2, double *out)
  {
    double tmp[2];

    //          TWO QUARKS TWO GLUONS SUBPROCESSES
    //  gg -> qqb   (1/Na * 1/Na)
    out[0] = (Nu*Qu4 + Nd*Qd4)*(amp1 -> su3_tree_mch(2,1, -1, 0, -2,-3))/Na2;
    
    //  ug -> ug    (1/Nc * 1/Na)
    //  dg -> dg    (1/Nc * 1/Na)
    out[1] = Qu4*(tmp[0] = (amp1 -> su3_tree_mch(-1,1, 0, 2, -2,-3))/(Nc*Na)); 
    out[2] = Qd4*tmp[0];
    
    //  gu -> ug    (1/Na * 1/Nc)
    //  gd -> dg    (1/Na * 1/Nc)
    out[3] = Qu4*(tmp[0] = (amp1 -> su3_tree_mch(0,1, -1, 2, -2,-3))/(Nc*Na));
    out[4] = Qd4*tmp[0];
    
    //  uub -> gg   (1/Nc * 1/Nc * 1/2)
    out[12] = Qu4*(tmp[0] = 0.5*(amp1 -> su3_tree_mch(-1,0, 1, 2, -2,-3))/Nc2);

    //  ddb -> gg   (1/Nc * 1/Nc * 1/2)
    out[14] = Qd4*tmp[0];

 
    //          FOUR QUARKS SUBPROCESSES
    //  uU -> uU    (1/Nc * 1/Nc)
    //  uu -> uu    (1/Nc * 1/Nc * 1/2)
    //  dD -> dD    (1/Nc * 1/Nc)
    //  dd -> dd    (1/Nc * 1/Nc * 1/2)
    amp2 -> su3_tree_mch(Qu, -1,1, 0, 2, -2,-3, tmp);
    out[5] = tmp[0]/Nc2;
    out[6] = tmp[1]/(2.0*Nc2);
    out[7] = out[5]/16.0;
    out[8] = out[6]/16.0;
    
    //  ud -> ud    (1/Nc * 1/Nc)
    out[9] = (amp2 -> su3_tree_mch(Qu, Qd, -1,1, 0,2, -2,-3))/Nc2;

    //  du -> du    (1/Nc * 1/Nc)
    out[10] = (amp2 -> su3_tree_mch(Qd, Qu, -1,1, 0, 2, -2,-3))/Nc2;

    //  uUb -> uUb    (1/Nc * 1/Nc)
    //  uub -> uub    (1/Nc * 1/Nc)
    //  dDb -> dDb    (1/Nc * 1/Nc)
    //  ddb -> ddb    (1/Nc * 1/Nc)
    amp2 -> su3_tree_mch(Qu, -1,1, 2, 0, -2,-3, tmp);
    out[11]  = tmp[0]/Nc2;
    out[12] += tmp[1]/Nc2;
    out[13]  = out[11]/16.0;
    out[14] += tmp[1]/(16.0*Nc2);

    //  uub -> UUb    (1/Nc * 1/Nc * (nu-1)) 
    //  uub -> ddb    (1/Nc * 1/Nc * nd) 
    tmp[0] = (amp2 -> su3_tree_mch(Qu, Qu, -1,0, 2, 1, -2,-3))/Nc2;
    out[12] += (Nu-1)*tmp[0]
      + Nd*(amp2 -> su3_tree_mch(Qu, Qd, -1,0, 2, 1, -2,-3))/Nc2;

    //  ddb -> DDb    (1/Nc * 1/Nc * (nd-1)) 
    //  ddb -> uub    (1/Nc * 1/Nc * nu) 
    out[14] += (Nd-1)*tmp[0]/16.0
      + Nu*(amp2 -> su3_tree_mch(Qd, Qu, -1,0, 2, 1, -2,-3))/Nc2;

    //  udb -> udb    (1/Nc * 1/Nc)
    out[15] = (amp2 -> su3_tree_mch(Qu, Qd, -1,1, 2, 0, -2,-3))/Nc2;
    
    //  dub -> dub    (1/Nc * 1/Nc)
    out[16] = (amp2 -> su3_tree_mch(Qd, Qu, -1,1, 2, 0, -2,-3))/Nc2;
  }
  

#define A(i) kp[i].tree
#define C(i) kp[i].cca
#define P(i) kp[i].pa
#define G(i) kp[i].ga


  void _hhc2ph_jet_base::
  __conv_x1(double eta, double x, double xjac, double al, const su3_kp_i2 *kp, weight_hhc2ph *S) 
  {
    double lne = std::log(1.0-eta);
    double lie = lne*(lne - 2.0*std::log(eta)) - 2.0*__specfunc_li2(eta);
        
    //----- K term -----
    double k[4][2];
    Kgg(x, xjac, Nf, al, k[0]); Kgq(x, xjac, al, k[1]); 
    Kqg(x, xjac, al, k[2]); Kqq(x, xjac, al, k[3]); 
    k[0][1] += Ca*lie; k[3][1] += Cf*lie;
    
    S[0][0]  =                2.0*k[1][0]*(Nu*A(0) + Nd*A(1));
    S[0][1]  = k[3][0]*A(0);
    S[0][2]  = k[3][0]*A(1);
    S[0][3]  = k[0][0]*A(2) + k[1][0]*A(4);
    S[0][4]  = k[0][0]*A(3) + k[1][0]*A(5);
    S[0][12] = k[3][0]*A(4) + k[2][0]*A(2);
    S[0][14] = k[3][0]*A(5) + k[2][0]*A(3);
    S[0][5] = S[0][6] = S[0][10] = S[0][11] = S[0][16] = k[2][0]*A(2);
    S[0][7] = S[0][8] = S[0][ 9] = S[0][13] = S[0][15] = k[2][0]*A(3);
   
    //----- Ktilde term -----
    double tmp1, tmp2, t[4][2];
    Tgg(x, xjac, al, t[0]); Tgq(x, xjac, al, t[1]); 
    Tqg(x, xjac, al, t[2]); Tqq(x, xjac, al, t[3]); 
    t[0][1] += Ca*lne*lne; t[3][1] += Cf*lne*lne; 
    
    S[0][0]  +=                2.0*t[1][0]*(Nu*C(0) + Nd*C(1));
    S[0][1]  += t[3][0]*C(0);
    S[0][2]  += t[3][0]*C(1);
    S[0][3]  += t[0][0]*C(2) + t[1][0]*C(4);
    S[0][4]  += t[0][0]*C(3) + t[1][0]*C(5);
    S[0][12] += t[3][0]*C(4) + t[2][0]*C(2);
    S[0][14] += t[3][0]*C(5) + t[2][0]*C(3);

    tmp1 = t[2][0]*C(2);
    tmp2 = t[2][0]*C(3);
    
    S[0][5] += tmp1; S[0][6] += tmp1; S[0][10] += tmp1; S[0][11] += tmp1;
    S[0][16] += tmp1; S[0][7] += tmp2; S[0][8] += tmp2; S[0][9] += tmp2;
    S[0][13] += tmp2; S[0][15] += tmp2; 

    //----- P term -----
    double p[4][2];
    Pgg(x, xjac, Nf, p[0]); Pgq(x, xjac, p[1]); 
    Pqg(x, xjac, p[2]); Pqq(x, xjac, p[3]); 
    p[0][1] += 2.0*Ca*lne; p[3][1] += 2.0*Cf*lne; 

    S[0][0]  +=                2.0*p[1][0]*(Nu*P(0) + Nd*P(1));
    S[0][1]  += p[3][0]*P(0);
    S[0][2]  += p[3][0]*P(1);
    S[0][3]  += p[0][0]*P(2) + p[1][0]*P(4);
    S[0][4]  += p[0][0]*P(3) + p[1][0]*P(5);
    S[0][12] += p[3][0]*P(4) + p[2][0]*P(2);
    S[0][14] += p[3][0]*P(5) + p[2][0]*P(3);

    tmp1 = p[2][0]*P(2);
    tmp2 = p[2][0]*P(3);
    
    S[0][5] += tmp1; S[0][6] += tmp1; S[0][10] += tmp1; S[0][11] += tmp1;
    S[0][16] += tmp1; S[0][7] += tmp2; S[0][8] += tmp2; S[0][9] += tmp2;
    S[0][13] += tmp2; S[0][15] += tmp2; 

    //----- G term -----
    double g[2];
    g[0] = (x > 1.0-al ? xjac/(x-x*x) : 0.0);
    g[1] = -xjac/(1.0-x) + al-std::log(al) + lne;

    S[0][1]  += g[0]*G(0); 
    S[0][2]  += g[0]*G(1); 
    S[0][3]  += g[0]*G(2); 
    S[0][4]  += g[0]*G(3); 
    S[0][12] += g[0]*G(4); 
    S[0][14] += g[0]*G(5);

    //----- dirac delta and + description terms -----
    S[2][1]  = g[1]*G(0) + k[3][1]*A(0) + t[3][1]*C(0) + p[3][1]*P(0);
    S[2][2]  = g[1]*G(1) + k[3][1]*A(1) + t[3][1]*C(1) + p[3][1]*P(1);
    S[2][3]  = g[1]*G(2) + k[0][1]*A(2) + t[0][1]*C(2) + p[0][1]*P(2);
    S[2][4]  = g[1]*G(3) + k[0][1]*A(3) + t[0][1]*C(3) + p[0][1]*P(3);
    S[2][12] = g[1]*G(4) + k[3][1]*A(4) + t[3][1]*C(4) + p[3][1]*P(4);
    S[2][14] = g[1]*G(5) + k[3][1]*A(5) + t[3][1]*C(5) + p[3][1]*P(5);
    S[2][5] = S[2][6] = S[2][10] = S[2][11] = S[2][16] = 0.0;
    S[2][7] = S[2][8] = S[2][ 9] = S[2][13] = S[2][15] = 0.0;
    
    //---- factorization scale dependent term ----
    S[3][0]  =               -2.0*p[1][0]*(Nu*A(0) + Nd*A(1));
    S[3][1]  = -p[3][0]*A(0);
    S[3][2]  = -p[3][0]*A(1);
    S[3][3]  = -p[0][0]*A(2) - p[1][0]*A(4);
    S[3][4]  = -p[0][0]*A(3) - p[1][0]*A(5);
    S[3][12] = -p[3][0]*A(4) - p[2][0]*A(2);
    S[3][14] = -p[3][0]*A(5) - p[2][0]*A(3);
    S[3][5] = S[3][6] = S[3][10] = S[3][11] = S[3][16] = -p[2][0]*A(2);
    S[3][7] = S[3][8] = S[3][ 9] = S[3][13] = S[3][15] = -p[2][0]*A(3);
   
    S[5][1]  = -p[3][1]*A(0);
    S[5][2]  = -p[3][1]*A(1);
    S[5][3]  = -p[0][1]*A(2);
    S[5][4]  = -p[0][1]*A(3);
    S[5][12] = -p[3][1]*A(4);
    S[5][14] = -p[3][1]*A(5);
    S[5][5] = S[5][6] = S[5][10] = S[5][11] = S[5][16] = 0.0;
    S[5][7] = S[5][8] = S[5][ 9] = S[5][13] = S[5][15] = 0.0;
  }

#undef C
#undef P
#undef G

#define C(i) kp[i].ccb
#define P(i) kp[i].pb
#define G(i) kp[i].gb

 
  void _hhc2ph_jet_base::
  __conv_x2(double eta, double x, double xjac, double al, const su3_kp_i2 *kp, weight_hhc2ph *S) 
  {
    double lne = std::log(1.0-eta);
    double lie = lne*(lne - 2.0*std::log(eta)) - 2.0*__specfunc_li2(eta);
        
    //----- K term -----
    double k[4][2];
    Kgg(x, xjac, Nf, al, k[0]); Kgq(x, xjac, al, k[1]); 
    Kqg(x, xjac, al, k[2]); Kqq(x, xjac, al, k[3]); 
    k[0][1] += Ca*lie; k[3][1] += Cf*lie;
    
    S[1][0]  =                2.0*k[1][0]*(Nu*A(2) + Nd*A(3));
    S[1][1]  = k[0][0]*A(0) + k[1][0]*A(4);
    S[1][2]  = k[0][0]*A(1) + k[1][0]*A(5);
    S[1][3]  = k[3][0]*A(2);
    S[1][4]  = k[3][0]*A(3);
    S[1][12] = k[3][0]*A(4) + k[2][0]*A(0);
    S[1][14] = k[3][0]*A(5) + k[2][0]*A(1);
    S[1][5] = S[1][6] = S[1][ 9] = S[1][11] = S[1][15] = k[2][0]*A(0);
    S[1][7] = S[1][8] = S[1][10] = S[1][13] = S[1][16] = k[2][0]*A(1);
    
    //----- Ktilde term -----
    double tmp1, tmp2, t[4][2];
    Tgg(x, xjac, al, t[0]); Tgq(x, xjac, al, t[1]); 
    Tqg(x, xjac, al, t[2]); Tqq(x, xjac, al, t[3]); 
    t[0][1] += Ca*lne*lne; t[3][1] += Cf*lne*lne; 
    
    S[1][0]  +=                2.0*t[1][0]*(Nu*C(2) + Nd*C(3));
    S[1][1]  += t[0][0]*C(0) + t[1][0]*C(4);
    S[1][2]  += t[0][0]*C(1) + t[1][0]*C(5);
    S[1][3]  += t[3][0]*C(2);
    S[1][4]  += t[3][0]*C(3);
    S[1][12] += t[3][0]*C(4) + t[2][0]*C(0);
    S[1][14] += t[3][0]*C(5) + t[2][0]*C(1);
    
    tmp1 = t[2][0]*C(0);
    tmp2 = t[2][0]*C(1);

    S[1][5] += tmp1; S[1][6] += tmp1; S[1][9] += tmp1; S[1][11] += tmp1; 
    S[1][15] += tmp1; S[1][7] += tmp2; S[1][8] += tmp2; S[1][10] += tmp2;
    S[1][13] += tmp2; S[1][16] += tmp2; 
    
    //----- P term -----
    double p[4][2];
    Pgg(x, xjac, Nf, p[0]); Pgq(x, xjac, p[1]); 
    Pqg(x, xjac, p[2]); Pqq(x, xjac, p[3]); 
    p[0][1] += 2.0*Ca*lne; p[3][1] += 2.0*Cf*lne; 

    S[1][0]  +=                2.0*p[1][0]*(Nu*P(2) + Nd*P(3));
    S[1][1]  += p[0][0]*P(0) + p[1][0]*P(4);
    S[1][2]  += p[0][0]*P(1) + p[1][0]*P(5);
    S[1][3]  += p[3][0]*P(2);
    S[1][4]  += p[3][0]*P(3);
    S[1][12] += p[3][0]*P(4) + p[2][0]*P(0);
    S[1][14] += p[3][0]*P(5) + p[2][0]*P(1);
    
    tmp1 = p[2][0]*P(0);
    tmp2 = p[2][0]*P(1);

    S[1][5] += tmp1; S[1][6] += tmp1; S[1][9] += tmp1; S[1][11] += tmp1; 
    S[1][15] += tmp1; S[1][7] += tmp2; S[1][8] += tmp2; S[1][10] += tmp2;
    S[1][13] += tmp2; S[1][16] += tmp2; 

    //----- G term -----
    double g[2];
    g[0] = (x > 1.0-al ? xjac/(x-x*x) : 0.0);
    g[1] = -xjac/(1.0-x) + al-std::log(al) + lne;

    S[1][1]  += g[0]*G(0); 
    S[1][2]  += g[0]*G(1); 
    S[1][3]  += g[0]*G(2); 
    S[1][4]  += g[0]*G(3); 
    S[1][12] += g[0]*G(4); 
    S[1][14] += g[0]*G(5);

    //----- dirac delta and + description terms -----
    S[2][1]  += g[1]*G(0) + k[0][1]*A(0) + t[0][1]*C(0) + p[0][1]*P(0);
    S[2][2]  += g[1]*G(1) + k[0][1]*A(1) + t[0][1]*C(1) + p[0][1]*P(1);
    S[2][3]  += g[1]*G(2) + k[3][1]*A(2) + t[3][1]*C(2) + p[3][1]*P(2);
    S[2][4]  += g[1]*G(3) + k[3][1]*A(3) + t[3][1]*C(3) + p[3][1]*P(3);
    S[2][12] += g[1]*G(4) + k[3][1]*A(4) + t[3][1]*C(4) + p[3][1]*P(4);
    S[2][14] += g[1]*G(5) + k[3][1]*A(5) + t[3][1]*C(5) + p[3][1]*P(5);

    //---- factorization scale dependent term ----
    S[4][0]  =               -2.0*p[1][0]*(Nu*A(2) + Nd*A(3));
    S[4][1]  = -p[0][0]*A(0) - p[1][0]*A(4);
    S[4][2]  = -p[0][0]*A(1) - p[1][0]*A(5);
    S[4][3]  = -p[3][0]*A(2);
    S[4][4]  = -p[3][0]*A(3);
    S[4][12] = -p[3][0]*A(4) - p[2][0]*A(0);
    S[4][14] = -p[3][0]*A(5) - p[2][0]*A(1);
    S[4][5] = S[1][6] = S[1][ 9] = S[1][11] = S[1][15] = -p[2][0]*A(0);
    S[4][7] = S[1][8] = S[1][10] = S[1][13] = S[1][16] = -p[2][0]*A(1);

    S[5][1]  -= p[0][1]*A(0);
    S[5][2]  -= p[0][1]*A(1);
    S[5][3]  -= p[3][1]*A(2);
    S[5][4]  -= p[3][1]*A(3);
    S[5][12] -= p[3][1]*A(4);
    S[5][14] -= p[3][1]*A(5);
  }
}

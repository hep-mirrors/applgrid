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
#include "bits/proc-photo.h"
#include "bits/nlo-color.h"
#include "bits/nlo-flkern.h"

#include "ampq2g1p1.h"
#include "ampq2g2p1.h"
#include "ampq2g3p1.h"
#include "ampq4p1.h"
#include "ampq4g1p1.h"

#include "ampq2g2.h"
#include "ampq2g3.h"
#include "ampq4.h"
#include "ampq4g1.h"

#define Qu  0.66666666666666666666 
#define Qd -0.33333333333333333333
#define Qu2 0.44444444444444444444
#define Qd2 0.11111111111111111111


namespace nlo {


  //
  //      THREE PARTON AMPLITUDES
  //
  void _photo_jet_base::amp_tree(ampq2g1p1 *amp, double *out)
  {
    //   two quarks one gluon contributions
    //  g -> qqb  1/Na
    //  u -> ug   1/Nc
    //  d -> dg   1/Nc
    out[0] = (Nu*Qu2+Nd*Qd2)*(amp -> su3_tree(2,1, 0, -1))/Na;
    out[1] = -Qu2*(amp -> su3_tree(0,1, 2, -1))/Nc;
    out[2] = 0.25*out[1];
  }

  void _photo_jet_base::amp_1loop(ampq2g1p1 *amp, double *out)
  {
    //   two quarks one gluon contributions
    //  g -> qqb  1/Na
    //  u -> ug   1/Nc
    //  d -> dg   1/Nc
    out[0] = (Nu*Qu2+Nd*Qd2)*(amp -> su3_1loop(2,1, 0, -1))/Na;
    out[1] = Qu2*(amp -> su3_1loop(0,1, 2, -1))/Nc;
    out[2] = 0.25*out[1];
  }

  void _photo_jet_base::amp_1loop_mch(ampq2g1p1 *amp, double *out)
  {
    //   two quarks one gluon contributions
    //  g -> qqb  1/Na
    //  u -> ug   1/Nc
    //  d -> dg   1/Nc
    out[0] = (Nu*Qu2+Nd*Qd2)*(amp -> su3_1loop_mch(2,1, 0, -1))/Na;
    out[1] = Qu2*(amp -> su3_1loop_mch(0,1, 2, -1))/Nc;
    out[2] = 0.25*out[1];
  }

  void _photo_jet_base::amp_kp(double al, ampq2g1p1 *amp, su3_kp_i1 *out)
  {
    static su3_kp_i1 tmp[2];
    //   two quarks one gluon contributions
    //  g -> qqb  1/Na
    //  u -> ug   1/Nc
    //  d -> dg   1/Nc
    amp -> su3_kp(Nf, 0, 2,1, 0, -1, tmp, al);
    amp -> su3_kp(Nf, 0, 0,1, 2, -1, tmp+1, al);
    out[0] = (Nu*Qu2+Nd*Qd2)*tmp[0]/Na;
    out[1] = -Qu2*tmp[1]/Nc;
    out[2] = 0.25*out[1];
  }

  double _photo_jet_base::
  amp_ccg(ampq2g1p1 *amp, int i, int j, int p1, int p2) {
    return (amp -> su3_cc(i,j, p2, p1, 0, -1))/Na;
  }
  
  double _photo_jet_base::
  amp_ccq(ampq2g1p1 *amp, int i, int j, int p1, int p2) {
    return -(amp -> su3_cc(i,j, 0, p1, p2, -1))/Nc;
  }

  //
  //       FOUR PARTON AMPLITUDES
  //
  void _photo_jet_base::amp_tree(ampq2g2p1 *amp1, ampq4p1 *amp2, double *out)
  {
    //   two quarks two gluons contributions
    //  g -> qqbg  1/Na
    //  u -> ugg   1/Nc * 1/2  
    //  d -> dgg   1/Nc * 1/2  
    out[0] = (Nu*Qu2+Nd*Qd2)*(amp1 -> su3_tree(2,1, 3,0, -1))/Na;
    out[1] = Qu2*(amp1 -> su3_tree(0,1, 2,3, -1))/(2.0*Nc);
    out[2] = 0.25*out[1];
    
    //   four quarks contributions
    double tmp[2];

    //  u -> uUUb   1/Nc * (Nu-1)
    //  u -> uuub   1/Nc * 1/2
    //  u -> uddb   1/Nc * Nd
    amp2 -> su3_tree(Qu, 0,1, 3,2, -1, tmp);
    out[1] += (Nu-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc)
      + Nd*(amp2 -> su3_tree(Qu, Qd, 0,1, 3,2, -1))/Nc;
    
    //  d -> dDDb   1/Nc * (Nd-1)
    //  d -> dddb   1/Nc * 1/2
    //  d -> duub   1/Nc * Nu
    out[2] += 0.25*((Nd-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc))
      + Nu*(amp2 -> su3_tree(Qd, Qu, 0,1, 3,2, -1))/Nc;
  }
  
  void _photo_jet_base::
  amp_kp(double al, ampq2g2p1 *amp1, ampq4p1 *amp2, su3_kp_i1 *out)
  {
    static su3_kp_i1 tmp[3];

    //              TWO QUARKS TWO GLUONS CONTRIBUTIONS
    //  g -> qqbg  1/Na
    //  u -> ugg   1/Nc * 1/2  
    //  d -> dgg   1/Nc * 1/2
    amp1 -> su3_kp(Nf, 0, 2,1, 3,0, -1, tmp, al);  
    amp1 -> su3_kp(Nf, 0, 0,1, 2,3, -1, tmp+1, al);
    out[0] = (Nu*Qu2+Nd*Qd2)*tmp[0]/Na;    
    out[1] = Qu2*tmp[1]/(2.0*Nc);
    out[2] = 0.25*out[1];
    
    //                  FOUR QUARKS CONTRIBUTIONS
    //  u -> uUUb   1/Nc * (Nu-1)
    //  u -> uuub   1/Nc * 1/2
    //  u -> uddb   1/Nc * Nd
    amp2 -> su3_kp(Qu, 0, 0,1, 3,2, -1, tmp, al);
    amp2 -> su3_kp(Qu, Qd, 0, 0,1, 3,2, -1, tmp+2, al);
    out[1] += (Nu-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc) + Nd*tmp[2]/Nc;
    
    //  d -> dDDb   1/Nc * (Nd-1)
    //  d -> dddb   1/Nc * 1/2
    //  d -> duub   1/Nc * Nu
    amp2 -> su3_kp(Qd, Qu, 0, 0,1, 3,2, -1, tmp+2, al);
    out[2] += 0.25*((Nd-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc)) + Nu*tmp[2]/Nc;
  }
  
  void _photo_jet_base::amp_1loop(ampq2g2p1 *amp1, ampq4p1 *amp2, double *out)
  {
    double tmp[2];

    //   two quarks two gluons contributions
    //  g -> qqbg  1/Na
    //  u -> ugg   1/Nc * 1/2  
    //  d -> dgg   1/Nc * 1/2  
    amp1 -> su3_1loop(Nu,Nd, 2,1, 3,0, -1, tmp);
    out[0] = (Nu*tmp[0]+Nd*tmp[1])/Na;
    
    amp1 -> su3_1loop(Nu,Nd, 0,1, 2,3, -1, tmp);
    out[1] = tmp[0]/(2.0*Nc);
    out[2] = tmp[1]/(2.0*Nc);
    
    //   four quarks contributions
    //  u -> uUUb   1/Nc * (Nu-1)
    //  u -> uuub   1/Nc * 1/2
    //  u -> uddb   1/Nc * Nd
    amp2 -> su3_1loop(Qu, Nf, 0,1, 3,2, -1, tmp);
    out[1] += (Nu-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc)
      + Nd*(amp2 -> su3_1loop(Qu, Qd, Nf, 0,1, 3,2, -1))/Nc;
    
    //  d -> dDDb   1/Nc * (Nd-1)
    //  d -> dddb   1/Nc * 1/2
    //  d -> duub   1/Nc * Nu
    out[2] += 0.25*((Nd-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc))
      + Nu*(amp2 -> su3_1loop(Qd, Qu, Nf, 0,1, 3,2, -1))/Nc;
  }

  void _photo_jet_base::amp_1loop_mch(ampq2g2p1 *amp1, ampq4p1 *amp2, double *out)
  {
    double tmp[2];

    //   two quarks two gluons contributions
    //  g -> qqbg  1/Na
    //  u -> ugg   1/Nc * 1/2  
    //  d -> dgg   1/Nc * 1/2  
    amp1 -> su3_1loop_mch(Nu,Nd, 2,1, 3,0, -1, tmp);
    out[0] = (Nu*tmp[0]+Nd*tmp[1])/Na;
    
    amp1 -> su3_1loop_mch(Nu,Nd, 0,1, 2,3, -1, tmp);
    out[1] = tmp[0]/(2.0*Nc);
    out[2] = tmp[1]/(2.0*Nc);
    
    //   four quarks contributions
    //  u -> uUUb   1/Nc * (Nu-1)
    //  u -> uuub   1/Nc * 1/2
    //  u -> uddb   1/Nc * Nd
    amp2 -> su3_1loop_mch(Qu, Nf, 0,1, 3,2, -1, tmp);
    out[1] += (Nu-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc)
      + Nd*(amp2 -> su3_1loop_mch(Qu, Qd, Nf, 0,1, 3,2, -1))/Nc;
    
    //  d -> dDDb   1/Nc * (Nd-1)
    //  d -> dddb   1/Nc * 1/2
    //  d -> duub   1/Nc * Nu
    out[2] += 0.25*((Nd-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc))
      + Nu*(amp2 -> su3_1loop_mch(Qd, Qu, Nf, 0,1, 3,2, -1))/Nc;
  }

  std::pair<double, std::complex<double> > _photo_jet_base::
  amp_ccg(ampq2g2p1 *amp, int i, int j, int p1, int p2, int p3) 
  {
    _Pair cc(amp -> su3_cc(i,j, p2, p1, p3, 0, -1));
    cc.first /= Na;
    cc.second /= Na;
    return cc;
  }

  std::pair<double, std::complex<double> > _photo_jet_base::
  amp_ccq(ampq2g2p1 *amp, int i, int j, int p1, int p2, int p3) 
  {
    _Pair cc(amp -> su3_cc(i,j, 0, p1, p2, p3, -1));
    cc.first /= Nc;
    cc.second /= -Nc;
    return cc;
  }

  void _photo_jet_base::
  amp_ccq(ampq4p1 *amp, double Q1, int i, int j, int p1, int p2, int p3, double *cc) 
  {
    amp -> su3_cc(Q1, i,j, 0, p1, p3, p2, -1, cc);
    cc[0] /= Nc;
    cc[1] /= Nc;
  }

  double _photo_jet_base::
  amp_ccq(ampq4p1 *amp, double Q1, double Q2, int i, int j, int p1, int p2, int p3) 
  {
    return (amp -> su3_cc(Q1, Q2, i,j, 0, p1, p3, p2, -1))/Nc;
  }


  void _photo_jet_base::amp_tree(ampq2g3p1 *amp1, ampq4g1p1 *amp2, double *out)
  {
    //   TWO QUARKS THREE GLUONS CONTRIBUTIONS
    //  g -> qqbgg  1/Na * 1/2
    //  u -> uggg   1/Nc * 1/6  
    //  d -> dggg   1/Nc * 1/6  
    out[0] = (Nu*Qu2+Nd*Qd2)*(amp1 -> su3_tree(2,1, 3,4,0, -1))/(2.0*Na);
    out[1] = Qu2*(amp1 -> su3_tree(0,1, 2,3,4, -1))/(6.0*Nc);
    out[2] = 0.25*out[1];

    //   FOUR QUARKS ONE GLUONS CONTRIBUTIONS
    double tmp[2];

    //  g -> uubUUb   1/Na * 1/2*Nu*(Nu-1)
    //  g -> uubuub   1/Na * 1/2 * 1/2 * Nu
    //  g -> ddbDDb   1/Na * 1/2*Nd*(Nd-1)
    //  g -> ddbddb   1/Na * 1/2 * 1/2 * Nd
    //  g -> uubddb   1/Na * Nu*Nd
    amp2 -> su3_tree(Qu, 2,1, 4,3, 0, -1, tmp);
    out[0] += 0.5*Nu*(Nu-1)*tmp[0]/Na + Nu*tmp[1]/(4.0*Na)
      + 0.25*(0.5*Nd*(Nd-1)*tmp[0]/Na + Nd*tmp[1]/(4.0*Na))
      + Nu*Nd*(amp2 -> su3_tree(Qu, Qd, 2,1, 4,3, 0, -1))/Na;
    
    //  u -> uUUbg   1/Nc * (Nu-1)
    //  u -> uuubg   1/Nc * 1/2
    //  u -> uddbg   1/Nc * Nd
    amp2 -> su3_tree(Qu, 0,1, 3,2, 4, -1, tmp);
    out[1] += (Nu-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc)
      + Nd*(amp2 -> su3_tree(Qu, Qd, 0,1, 3,2, 4, -1))/Nc;
    
    //  d -> dDDbg   1/Nc * (Nd-1)
    //  d -> dddbg   1/Nc * 1/2
    //  d -> duubg   1/Nc * Nu
    out[2] += 0.25*((Nd-1)*tmp[0]/Nc + tmp[1]/(2.0*Nc))
      + Nu*(amp2 -> su3_tree(Qd, Qu, 0,1, 3,2, 4, -1))/Nc;
  }


#define A(i) kp[i].tree
#define C(i) kp[i].cca
#define P(i) kp[i].pa
#define G(i) kp[i].ga

  void _photo_jet_base::
  conv_parton(double eta, double x, double xjac, double al,
	      const su3_kp_i1 *kp, weight_photo *S)
  {
    double lne = std::log(1.0-eta);
    double lie = lne*(lne - 2.0*std::log(eta)) - 2.0*__specfunc_li2(eta);

    //---- K terms ----
    double kgq[2], kqg[2], kqq[2], kgg[2];
    Kgg(x, xjac, Nf, al, kgg); Kgq(x, xjac, al, kgq); 
    Kqg(x, xjac, al, kqg); Kqq(x, xjac, al, kqq); 
    kgg[1] += Ca*lie; kqq[1] += Cf*lie;

    S[1][0] = kgg[0]*A(0) + 2.0*Nu*kgq[0]*A(1) + 2.0*Nd*kgq[0]*A(2);
    S[1][1] = kqq[0]*A(1) + kqg[0]*A(0);
    S[1][2] = kqq[0]*A(2) + kqg[0]*A(0);
  
    //---- P terms ----
    double pgq[2],  pqg[2], pqq[2], pgg[2];
    Pgg(x, xjac, Nf, pgg); Pgq(x, xjac, pgq); 
    Pqg(x, xjac, pqg); Pqq(x, xjac, pqq); 
    pgg[1] += 2.0*Ca*lne; pqq[1] += 2.0*Cf*lne; 

    S[1][0] += pgg[0]*P(0) + 2.0*Nu*pgq[0]*P(1) + 2.0*Nd*pgq[0]*P(2);
    S[1][1] += pqq[0]*P(1) + pqg[0]*P(0);
    S[1][2] += pqq[0]*P(2) + pqg[0]*P(0);

    //---- G term ----
    double g[2];
    g[0] = (x > 1.0-al ? xjac/(x-x*x) : 0.0);
    g[1] = -xjac/(1.0-x) + al-std::log(al) + lne;

    S[1][0] += g[0]*G(0);    
    S[1][1] += g[0]*G(1);    
    S[1][2] += g[0]*G(2);

    //---- '+' description & dirac delta terms ----
    S[2][0] = g[1]*G(0) + kgg[1]*A(0) + pgg[1]*P(0);
    S[2][1] = g[1]*G(1) + kqq[1]*A(1) + pqq[1]*P(1);
    S[2][2] = g[1]*G(2) + kqq[1]*A(2) + pqq[1]*P(2);
    
    //---- factorization scale dependent term ----
    S[4][0] = -pgg[0]*A(0) - 2.0*Nu*pgq[0]*A(1) - 2.0*Nd*pgq[0]*A(2);    
    S[4][1] = -pqq[0]*A(1) - pqg[0]*A(0);
    S[4][2] = -pqq[0]*A(2) - pqg[0]*A(0);
    
    S[5][0] = -pgg[1]*A(0);
    S[5][1] = -pqq[1]*A(1);
    S[5][2] = -pqq[1]*A(2);
  }

  void _photo_jet_base::
  conv_photon(double eta, double x, double xjac, double al, const su3_kp_i2 *kp, weight_photo *S) 
  {
    //----- K, Ktilde & P term -----
    double kgq[2], tgq[2], pgq[2];
    Kgq(x, xjac, al, kgq); Tgq(x, xjac, al, tgq); Pgq(x, xjac, pgq); 
    
    S[0][0] = kgq[0]*A(0) + tgq[0]*C(0) + pgq[0]*P(0);
    S[0][1] = kgq[0]*A(1) + tgq[0]*C(1) + pgq[0]*P(1);
    S[0][2] = kgq[0]*A(2) + tgq[0]*C(2) + pgq[0]*P(2);
    
    //---- factorization scale dependent term ----
    S[3][0] = -pgq[0]*A(0);
    S[3][1] = -pgq[0]*A(1);
    S[3][2] = -pgq[0]*A(2);
    
    // conversion factor 2Nc
    S[0] *= 2.0*Nc; S[3] *= 2.0*Nc;
  }
 

  void _photo_jet_base::
  amp_kp(double al, ampq2g2 *amp2, ampq4 *amp3, su3_kp_i2 *out)
  {
    static su3_kp_i2 tmp[2];
    
    //----- TWO QUARKS TWO GLUONS PROCESSES -----
    //  qg -> qg    (1/Nc * 1/Na)
    amp2 -> su3_kp(Nf, -1,0, 1,-1,2,0, tmp, al);
    out[0] = 2.0*(Nu*Qu2+Nd*Qd2)*tmp[0]/(-Nc*Na);
    
    //  qbq -> gg   (1/2 * 1/Nc * 1/Nc)
    amp2 -> su3_kp(Nf, -1,0, -1,0,1,2, tmp, al);
    out[1] = Qu2*tmp[0]/(2.0*Nc2);
    out[2] = Qd2*tmp[0]/(2.0*Nc2);
    
    //----- FOUR QUARKS PROCESSES -----
    //  qr -> qr    (1/Nc * 1/Nc)
    //  qq -> qq    (1/2 * 1/Nc * 1/Nc)
    amp3 -> su3_kp(-1,0, 1,-1,2,0, tmp, al);
    out[1] += ((Nu-1)*Qu2 + Nd*Qd2)*tmp[0]/Nc2 + Qu2*tmp[1]/(2.0*Nc2);
    out[2] += ((Nd-1)*Qd2 + Nu*Qu2)*tmp[0]/Nc2 + Qd2*tmp[1]/(2.0*Nc2);
    
    //  qbq -> rrb  (1/Nc * 1/Nc * (nf-1))
    //  qbq -> qqb  (1/Nc * 1/Nc)
    amp3 -> su3_kp(-1,0, -1,0,1,2, tmp, al); 
    out[1] += Qu2*(tmp[0]*((Nf-1)/Nc2) + tmp[1]/Nc2);
    out[2] += Qd2*(tmp[0]*((Nf-1)/Nc2) + tmp[1]/Nc2);
    
    //  rbq -> qrb  (1/Nc * 1/Nc)
    amp3 -> su3_kp(-1,0, 1,0,-1,2, tmp, al);
    out[1] += ((Nu-1)*Qu2+Nd*Qd2)*tmp[0]/Nc2;
    out[2] += ((Nd-1)*Qd2+Nu*Qu2)*tmp[0]/Nc2;
  }
 
  void _photo_jet_base::
  amp_kp(double al, ampq2g3 *amp2, ampq4g1 *amp3, su3_kp_i2 *out)
  {
    static su3_kp_i2 tmp[2];
    
    //----- TWO QUARKS THREE GLUONS PROCESSES -----
    //  qg -> qgg    (1/2 * 1/Nc * 1/Na)
    amp2 -> su3_kp(Nf, -1,0, 1,-1,2,3,0, tmp, al);
    out[0] = 2.0*(Nu*Qu2+Nd*Qd2)*tmp[0]/(-2.0*Nc*Na);
    
    //  qbq -> ggg   (1/6 * 1/Nc * 1/Nc)
    amp2 -> su3_kp(Nf, -1,0, -1,0,1,2,3, tmp, al);
    out[1] = Qu2*tmp[0]/(6.0*Nc2);
    out[2] = Qd2*tmp[0]/(6.0*Nc2);
      
    //----- FOUR QUARKS ONE GLUON PROCESSES -----
    //  qg -> rrbq   (1/Nc * 1/Na * (nf-1))
    //  qg -> qqbq   (1/2 * 1/Nc * 1/Na)
    amp3 -> su3_kp(Nf, -1,0, 3,-1,1,2,0, tmp, al);
    out[0] -= 2.0*(Nu*Qu2+Nd*Qd2)*(tmp[0]*(Nf-1)/(Nc*Na) + tmp[1]/(2.0*Nc*Na));
    
    //  qr -> qrg    (1/Nc * 1/Nc)
    //  qq -> qqg    (1/2 * 1/Nc * 1/Nc)
    amp3 -> su3_kp(Nf, -1,0, 1,-1,2,0,3, tmp, al);
    out[1] += ((Nu-1)*Qu2+Nd*Qd2)*tmp[0]/Nc2 + Qu2*tmp[1]/(2.0*Nc2);
    out[2] += ((Nd-1)*Qd2+Nu*Qu2)*tmp[0]/Nc2 + Qd2*tmp[1]/(2.0*Nc2);
    
    //  qbq -> rrbg  (1/Nc * 1/Nc * (nf-1))
    //  qbq -> qqbg  (1/Nc * 1/Nc)
    amp3 -> su3_kp(Nf, -1,0, -1,0,1,2,3, tmp, al); 
    out[1] += Qu2*(tmp[0]*((Nf-1)/Nc2) + tmp[1]/Nc2);
    out[2] += Qd2*(tmp[0]*((Nf-1)/Nc2) + tmp[1]/Nc2);
    
    //  rbq -> qrbg  (1/Nc * 1/Nc)
    amp3 -> su3_kp(Nf, -1,0, 1,0,-1,2,3, tmp, al);
    out[1] += ((Nu-1)*Qu2+Nd*Qd2)*tmp[0]/Nc2;
    out[2] += ((Nd-1)*Qd2+Nu*Qu2)*tmp[0]/Nc2;
  }
}    //  namespace nlo


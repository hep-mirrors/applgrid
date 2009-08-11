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
#include "ampq2g1l2.h"
#include "ampq2g2l2.h"
#include "ampq2g3l2.h"
#include "ampq4l2.h"
#include "ampq4g1l2.h"

#include "bits/nlo-color.h"
#include "bits/nlo-flkern.h"
#include "bits/proc-dis.h"


#define Qu  0.666666666666666
#define Qd  (-0.333333333333333)
#define Qu2 0.444444444444444
#define Qd2 0.111111111111111


namespace nlo {
  

  _dis_jet_base::_dis_jet_base(unsigned int nu, unsigned int nd)
  {
    charge = (2.0*nu - nd)/3.0; 
    charge_sqr = (4.0*nu + nd)/9.0;
    Nu = nu; Nd = nd; Nf = nu + nd;
  }
  
  //
  //    Three parton contributions
  //
  //    amp[0] : -e(-2) + e(-1) + g(0) --> q(1) + qb(2)
  //    amp[1] : -e(-2) + e(-1) + u(0) --> q(1) + g(2)
  //    amp[2] : -e(-2) + e(-1) + d(0) --> q(1) + g(2)
  void _dis_jet_base::amp_tree(ampq2g1l2 *amp, double *res) 
  {
    res[0] = charge_sqr*(amp -> su3_tree(1,0,2,-1,-2))/Na;
    res[1] = Qu2*(amp -> su3_tree(1,2,0,-1,-2))/Nc;
    res[2] = 0.25*res[1];
  }

  void _dis_jet_base::amp_tree_mch(ampq2g1l2 *amp, double *res) 
  {
    res[0] = charge_sqr*(amp -> su3_tree_mch(1,0,2,-1,-2))/Na;
    res[1] = Qu2*(amp -> su3_tree_mch(1,2,0,-1,-2))/Nc;
    res[2] = 0.25*res[1];
  }
  
  void _dis_jet_base::
  amp_1loop(ampq2g1l2 *amp, double *res) 
  {
    res[0] = charge_sqr*(amp -> su3_1loop(Nf, 1,0,2,-1,-2))/Na;
    res[1] = Qu2*(amp -> su3_1loop(Nf, 1,2,0,-1,-2))/Nc;
    res[2] = 0.25*res[1];
  }

  void _dis_jet_base::
  amp_1loop_mch(ampq2g1l2 *amp, double *res) 
  {
    res[0] = charge_sqr*(amp -> su3_1loop_mch(Nf, 1,0,2,-1,-2))/Na;
    res[1] = Qu2*(amp -> su3_1loop_mch(Nf, 1,2,0,-1,-2))/Nc;
    res[2] = 0.25*res[1];
  }
  
  void _dis_jet_base::amp_kp(double al, ampq2g1l2 *amp, su3_kp_i1 *res) 
  {
    amp -> su3_kp(Nf, 0, 1,0,2,-1,-2, res, al);
    res[0] *= charge_sqr/Na;
    
    amp -> su3_kp(Nf, 0, 1,2,0,-1,-2, res+1, al);
    res[1] *= Qu2/Nc;
    res[2] = 0.25*res[1];
  }
  
  void _dis_jet_base::amp_kp_mch(double al, ampq2g1l2 *amp, su3_kp_i1 *res) 
  {
    amp -> su3_kp_mch(Nf, 0, 1,0,2,-1,-2, res, al);
    res[0] *= charge_sqr/Na;
    
    amp -> su3_kp_mch(Nf, 0, 1,2,0,-1,-2, res+1, al);
    res[1] *= Qu2/Nc;
    res[2] = 0.25*res[1];
  }
  
  void _dis_jet_base::
  amp_ccg(ampq2g1l2 *amp, int i, int j, int p1, int p2, _Pair *res) 
  {
    _Pair cc(amp -> su3_cc(i,j, p1,0,p2,-1,-2));
    
    res -> first  = cc.first*charge_sqr/Na;
    res -> second = cc.second*charge_sqr/Na;
  }
  
  void _dis_jet_base::
  amp_ccq(ampq2g1l2 *amp, int i, int j, int p1, int p2, _Pair *res)
  {
    _Pair cc(amp -> su3_cc(i,j, p1,p2,0,-1,-2));
    
    res[1].first  = Qu2*cc.first/Nc;
    res[1].second = Qu2*cc.second/Nc;
    
    res[2].first  = Qd2*cc.first/Nc;
    res[2].second = Qd2*cc.second/Nc;
  }
  
  void _dis_jet_base::
  amp_cca(ampq2g1l2 *amp, int i, int j, int p1, int p2, _Pair *res) 
  {
    _Pair cc(amp -> su3_cc(i,j, 0,p2,p1,-1,-2));
    
    res[1].first  = Qu2*cc.first/Nc;
    res[1].second = Qu2*cc.second/Nc;
    
    res[2].first  = Qd2*cc.first/Nc;
    res[2].second = Qd2*cc.second/Nc;
  }
  
  //
  //   Four parton
  //
  //    amp[0] : -e(-2) + e(-1) + g(0) --> q(1) + qb(2) + g(3)
  //    amp[1] : -e(-2) + e(-1) + u(0) --> q(1) + g(2)  + g(3)
  //    amp[2] : -e(-2) + e(-1) + d(0) --> q(1) + g(2)  + g(3)
  //
  //    Amp[0] : 0
  //    Amp[1] : -e(-2) + e(-1) + u(0) --> q(1) + Q(2)  + Qb(3)
  //    Amp[2] : -e(-2) + e(-1) + d(0) --> q(1) + Q(2)  + Qb(3)
  void _dis_jet_base::amp_tree(ampq2g2l2 *amp1, ampq4l2 *amp2, double *res)
  {
    res[0] = charge_sqr*(amp1 -> su3_tree(1,0,3,2,-1,-2))/Na;
    res[1] = 0.5*Qu2*(amp1 -> su3_tree(1,2,3,0,-1,-2))/Nc;
    res[2] = 0.25*res[1];
    
    double A1, A2, A3, amp[10];
    amp2 -> su3_tree(1,3,2,0,-1,-2, amp);
    
    A1  = charge_sqr*(amp[1] + amp[3]);
    A2  = charge*(amp[2] + amp[5]);
    A3  = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[1] += 0.5*(A1 + Qu*A2 + Qu2*A3)/Nc;
    res[2] += 0.5*(A1 + Qd*A2 + Qd2*A3)/Nc;
  }
  
  void _dis_jet_base::amp_tree_mch(ampq2g2l2 *amp1, ampq4l2 *amp2, double *res)
  {
    res[0] = charge_sqr*(amp1 -> su3_tree_mch(1,0,3,2,-1,-2))/Na;
    res[1] = 0.5*Qu2*(amp1 -> su3_tree_mch(1,2,3,0,-1,-2))/Nc;
    res[2] = 0.25*res[1];
    
    double A1, A2, A3, amp[10];
    amp2 -> su3_tree_mch(1,3,2,0,-1,-2, amp);
    
    A1  = charge_sqr*(amp[1] + amp[3]);
    A2  = charge*(amp[2] + amp[5]);
    A3  = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[1] += 0.5*(A1 + Qu*A2 + Qu2*A3)/Nc;
    res[2] += 0.5*(A1 + Qd*A2 + Qd2*A3)/Nc;
  }
  
  void _dis_jet_base::
  amp_1loop(ampq2g2l2 *amp1, ampq4l2 *amp2, double *res)
  {   
    res[0] = charge_sqr*(amp1 -> su3_1loop(Nf, 1,0,3,2,-1,-2))/Na;
    res[1] = 0.5*Qu2*(amp1 -> su3_1loop(Nf, 1,2,3,0,-1,-2))/Nc;
    res[2] = 0.25*res[1];
    
    double A1, A2, A3, amp[10];
    amp2 -> su3_1loop(Nf, 1,3,2,0,-1,-2, amp);
    
    A1 = charge_sqr*(amp[1] + amp[3]);
    A2 = charge*(amp[2] + amp[5]);
    A3 = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[1] += 0.5*(A1 + Qu*A2 + Qu2*A3)/Nc;
    res[2] += 0.5*(A1 + Qd*A2 + Qd2*A3)/Nc;
  }
  
  void _dis_jet_base::
  amp_1loop_mch(ampq2g2l2 *amp1, ampq4l2 *amp2, double *res)
  {
    res[0] = charge_sqr*(amp1 -> su3_1loop_mch(Nf, 1,0,3,2,-1,-2))/Na;
    res[1] = 0.5*Qu2*(amp1 -> su3_1loop_mch(Nf, 1,2,3,0,-1,-2))/Nc;
    res[2] = 0.25*res[1];
    
    double A1, A2, A3, amp[10];
    amp2 -> su3_1loop_mch(Nf, 1,3,2,0,-1,-2, amp);
    
    A1 = charge_sqr*(amp[1] + amp[3]);
    A2 = charge*(amp[2] + amp[5]);
    A3 = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[1] += 0.5*(A1 + Qu*A2 + Qu2*A3)/Nc;
    res[2] += 0.5*(A1 + Qd*A2 + Qd2*A3)/Nc;
  }
  
  void _dis_jet_base::amp_kp(double al, ampq2g2l2 *amp1, ampq4l2 *amp2, su3_kp_i1 *res)
  {
    static su3_kp_i1 kp[10], A1, A2, A3;
    
    amp1 -> su3_kp(Nf, 0, 1,0,3,2,-1,-2, res, al);
    res[0] *= charge_sqr/Na;
    
    amp1 -> su3_kp(Nf, 0, 1,2,3,0,-1,-2, res+1, al);
    res[1] *= 0.5*Qu2/Nc;
    res[2]  = 0.25*res[1];
    
    amp2 -> su3_kp(0, 1,3,2,0,-1,-2, kp, al);
    A1 = charge_sqr*(kp[1] + kp[3]); 
    A2 = charge*(kp[2] + kp[5]);	   
    A3 = Nf*(kp[0] + kp[4]) + kp[6] + kp[7] + kp[8] + kp[9];
    res[1] += 0.5*(A1 + Qu*A2 + Qu2*A3)/Nc;	    
    res[2] += 0.5*(A1 + Qd*A2 + Qd2*A3)/Nc;
  }

  void _dis_jet_base::amp_kp_mch(double al, ampq2g2l2 *amp1, ampq4l2 *amp2, su3_kp_i1 *res)
  {
    static su3_kp_i1 kp[10], A1, A2, A3;
    
    amp1 -> su3_kp_mch(Nf, 0, 1,0,3,2,-1,-2, res, al);
    res[0] *= charge_sqr/Na;
    
    amp1 -> su3_kp_mch(Nf, 0, 1,2,3,0,-1,-2, res+1, al);
    res[1] *= 0.5*Qu2/Nc;
    res[2]  = 0.25*res[1];
    
    amp2 -> su3_kp_mch(0, 1,3,2,0,-1,-2, kp, al);
    A1 = charge_sqr*(kp[1] + kp[3]); 
    A2 = charge*(kp[2] + kp[5]);	   
    A3 = Nf*(kp[0] + kp[4]) + kp[6] + kp[7] + kp[8] + kp[9];
    res[1] += 0.5*(A1 + Qu*A2 + Qu2*A3)/Nc;	    
    res[2] += 0.5*(A1 + Qd*A2 + Qd2*A3)/Nc;
  }

  void _dis_jet_base::
  amp_ccg(ampq2g2l2 *amp, int i, int j, int p1, int p2, int p3, _Pair *res)
  {
    _Pair cc(amp -> su3_cc(i,j, p1,0,p3,p2,-1,-2));
    
    res -> first  = cc.first*charge_sqr/Na;
    res -> second = cc.second*charge_sqr/Na;
  }
  
  void _dis_jet_base::
  amp_ccq(ampq2g2l2 *amp, int i, int j, int p1, int p2, int p3, _Pair *res)
  {
    _Pair cc(amp -> su3_cc(i,j, p1,p2,p3,0,-1,-2));
    
    res[1].first  = Qu2*cc.first/Nc;
    res[1].second = Qu2*cc.second/Nc;
    			
    res[2].first  = Qd2*cc.first/Nc;
    res[2].second = Qd2*cc.second/Nc;
  }
  
  void _dis_jet_base::
  amp_cca(ampq2g2l2 *amp, int i, int j, int p1, int p2, int p3, _Pair *res)
  {
    _Pair cc(amp -> su3_cc(i,j, 0,p2,p3,p1,-1,-2));
    
    res[1].first  = Qu2*cc.first/Nc;
    res[1].second = Qu2*cc.second/Nc;
    			
    res[2].first  = Qd2*cc.first/Nc;
    res[2].second = Qd2*cc.second/Nc;
  }
  
  void _dis_jet_base::
  amp_ccq(ampq4l2 *amp1, int i, int j, int p1, int p2, int p3, _Pair *res)
  {
    double A1, A2, A3, amp[10];
    amp1 -> su3_cc(i,j, p1,p3,p2,0, -1,-2, amp);
    
    A1 = charge_sqr*(amp[1] + amp[3]);
    A2 = charge*(amp[2] + amp[5]);
    A3 = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[1].first = (A1 + Qu*A2 + Qu2*A3)/Nc;
    res[2].first = (A1 + Qd*A2 + Qd2*A3)/Nc;
    
    res[1].second = res[2].second = 0.0;
  }
  
  void _dis_jet_base::
  amp_cca(ampq4l2 *amp1, int i, int j, int p1, int p2, int p3, _Pair *res)
  {
    double A1, A2, A3, amp[10];
    amp1 -> su3_cc(i,j, 0,p3,p2,p1, -1,-2, amp);
    
    A1 = charge_sqr*(amp[1] + amp[3]);
    A2 = charge*(amp[2] + amp[5]);
    A3 = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[1].first = (A1 + Qu*A2 + Qu2*A3)/Nc;
    res[2].first = (A1 + Qd*A2 + Qd2*A3)/Nc;
    
    res[1].second = res[2].second = 0.0;
  }
  
  
  //
  //   Five parton
  //
  //    amp[0] : -e(-2) + e(-1) + g(0) --> q(1) + qb(2) + g(3) + g(4)
  //    amp[1] : -e(-2) + e(-1) + u(0) --> q(1) + g(2)  + g(3) + g(4)
  //    amp[2] : -e(-2) + e(-1) + d(0) --> q(1) + g(2)  + g(3) + g(4)
  //
  //    Amp[0] : -e(-2) + e(-1) + g(0) --> q(1) + q(2)  + Q(3)  + Qb(4)
  //    Amp[1] : -e(-2) + e(-1) + u(0) --> q(1) + Q(2)  + Qb(3) + g(4)
  //    Amp[2] : -e(-2) + e(-1) + d(0) --> q(1) + Q(2)  + Qb(3) + g(4)
  void _dis_jet_base::
  amp_tree(ampq2g3l2 *amp1, ampq4g1l2 *amp2, double *res)
  {
    res[0] = 0.5*charge_sqr*(amp1 -> su3_tree(1,0,3,4,2,-1,-2))/Na;
    res[1] = Qu2*(amp1 -> su3_tree(1,2,3,4,0,-1,-2))/(6.0*Nc);
    res[2] = 0.25*res[1];
    
    double A1, A2, A3, amp[10];
    amp2 -> su3_tree(1,3,2,0,4,-1,-2, amp);
    
    A1 = charge_sqr*(amp[1] + amp[3]);
    A2 = charge*(amp[2] + amp[5]);
    A3 = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[1] += 0.5*(A1 + Qu*A2 + Qu2*A3)/Nc;
    res[2] += 0.5*(A1 + Qd*A2 + Qd2*A3)/Nc;
    
    amp2 -> su3_tree(1,4,3,2,0,-1,-2, amp);
    A1 = charge_sqr*(amp[1] + amp[3]);
    A2 = charge*(amp[2] + amp[5]);
    A3 = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[0] += 0.25*(Nf*A1 + charge*A2 + charge_sqr*A3)/Na;
  }

  void _dis_jet_base::
  amp_tree_mch(ampq2g3l2 *amp1, ampq4g1l2 *amp2, double *res)
  {
    res[0] = 0.5*charge_sqr*(amp1 -> su3_tree_mch(1,0,3,4,2,-1,-2))/Na;
    res[1] = Qu2*(amp1 -> su3_tree_mch(1,2,3,4,0,-1,-2))/(6.0*Nc);
    res[2] = 0.25*res[1];
    
    double A1, A2, A3, amp[10];
    amp2 -> su3_tree_mch(1,3,2,0,4,-1,-2, amp);
    A1 = charge_sqr*(amp[1] + amp[3]);
    A2 = charge*(amp[2] + amp[5]);
    A3 = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[1] += 0.5*(A1 + Qu*A2 + Qu2*A3)/Nc;
    res[2] += 0.5*(A1 + Qd*A2 + Qd2*A3)/Nc;
    
    amp2 -> su3_tree_mch(1,4,3,2,0,-1,-2, amp);
    A1 = charge_sqr*(amp[1] + amp[3]);
    A2 = charge*(amp[2] + amp[5]);
    A3 = Nf*(amp[0] + amp[4]) + amp[6] + amp[7] + amp[8] + amp[9];
    
    res[0] += 0.25*(Nf*A1 + charge*A2 + charge_sqr*A3)/Na;
  }



#define A(i) kp[i].tree
#define P(i) kp[i].pa
#define G(i) kp[i].ga

  void _dis_jet_base::
  convolution(double eta, double x, double xjac, double al,
	      const su3_kp_i1 *kp, weight_dis *S)
  {
    double lne = std::log(1.0-eta);
    double lie = lne*(lne - 2.0*std::log(eta)) - 2.0*__specfunc_li2(eta);

    //---- K terms ----
    double kgq[2], kqg[2], kqq[2], kgg[2];
    Kgg(x, xjac, Nf, al, kgg); Kgq(x, xjac, al, kgq); 
    Kqg(x, xjac, al, kqg); Kqq(x, xjac, al, kqq); 
    kgg[1] += Ca*lie; kqq[1] += Cf*lie;

    S[0][0] = kgg[0]*A(0) + 2.0*Nu*kgq[0]*A(1) + 2.0*Nd*kgq[0]*A(2);
    S[0][1] = kqq[0]*A(1) + kqg[0]*A(0);
    S[0][2] = kqq[0]*A(2) + kqg[0]*A(0);
  
    //---- P terms ----
    double pgq[2],  pqg[2], pqq[2], pgg[2];
    Pgg(x, xjac, Nf, pgg); Pgq(x, xjac, pgq); 
    Pqg(x, xjac, pqg); Pqq(x, xjac, pqq); 
    pgg[1] += 2.0*Ca*lne; pqq[1] += 2.0*Cf*lne; 

    S[0][0] += pgg[0]*P(0) + 2.0*Nu*pgq[0]*P(1) + 2.0*Nd*pgq[0]*P(2);
    S[0][1] += pqq[0]*P(1) + pqg[0]*P(0);
    S[0][2] += pqq[0]*P(2) + pqg[0]*P(0);

    //---- G term ----
    double g[2];
    g[0] = (x > 1.0-al ? xjac/(x-x*x) : 0.0);
    g[1] = -xjac/(1.0-x) + al-std::log(al) + lne;

    S[0][0] += g[0]*G(0);    
    S[0][1] += g[0]*G(1);    
    S[0][2] += g[0]*G(2);

    //---- '+' description & dirac delta terms ----
    S[1][0] = g[1]*G(0) + kgg[1]*A(0) + pgg[1]*P(0);
    S[1][1] = g[1]*G(1) + kqq[1]*A(1) + pqq[1]*P(1);
    S[1][2] = g[1]*G(2) + kqq[1]*A(2) + pqq[1]*P(2);
    
    //---- factorization scale dependent term ----
    S[2][0] = -pgg[0]*A(0) - 2.0*Nu*pgq[0]*A(1) - 2.0*Nd*pgq[0]*A(2);    
    S[2][1] = -pqq[0]*A(1) - pqg[0]*A(0);
    S[2][2] = -pqq[0]*A(2) - pqg[0]*A(0);
    
    S[3][0] = -pgg[1]*A(0);
    S[3][1] = -pqq[1]*A(1);
    S[3][2] = -pqq[1]*A(2);
  }
}   //   namespace nlo

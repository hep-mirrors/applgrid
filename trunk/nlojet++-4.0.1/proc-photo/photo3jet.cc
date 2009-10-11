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
#include "photo3jet.h"
#include "ampq2g2p1.h"
#include "ampq2g3p1.h"
#include "ampq4p1.h"
#include "ampq4g1p1.h"
#include "ampq4g1.h"
#include "ampq2g3.h"
#include "bits/nlo-color.h"


//   PI_FACn = (2pi)^(2n-3) * 2^(n-1)
#define PI_FAC4 78341.03930503205203742586
#define PI_FAC5 6185560.53048687015450831204



namespace nlo {

  
  photo3jet::dipole_func photo3jet::_S_dipole[14] = 
  { &photo3jet::_M_di1, &photo3jet::_M_di2, &photo3jet::_M_di3, 
    &photo3jet::_M_di4, &photo3jet::_M_d01, &photo3jet::_M_d02,
    &photo3jet::_M_d03, &photo3jet::_M_d04, &photo3jet::_M_d12, 
    &photo3jet::_M_d13, &photo3jet::_M_d14, &photo3jet::_M_d23,
    &photo3jet::_M_d24, &photo3jet::_M_d34};


  photo3jet::photo3jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_photo(3, 2, nu, nd, al), _photo_jet_base(nu, nd), _M_mchel(mchel) 
  {
    _M_q2g2p1 = new ampq2g2p1(_M_ip, rng);
    _M_q4p1   = new ampq4p1  (_M_ip, rng);
    _M_q2g3p1 = new ampq2g3p1(_M_ip, rng);
    _M_q4g1p1 = new ampq4g1p1(_M_ip, rng);

    _M_q2g3 = new ampq2g3(_M_ip, rng);
    _M_q4g1 = new ampq4g1(_M_ip, rng);
  }
  
  photo3jet::~photo3jet() 
  {
    if(_M_q2g2p1) delete _M_q2g2p1;
    if(_M_q4p1)   delete _M_q4p1;
    if(_M_q2g3p1) delete _M_q2g3p1;
    if(_M_q4g1p1) delete _M_q4g1p1;

    if(_M_q2g3) delete _M_q2g3;
    if(_M_q4g1) delete _M_q4g1;    
  }
  
  void photo3jet::born_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_q2g2p1, _M_q4p1, res.begin());
    res *= PI_FAC4;
  }
  
  void photo3jet::real_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_q2g3p1, _M_q4g1p1, res.begin());
    res *= PI_FAC5;
  }

  void photo3jet::fini_term(double x1, double xjac1, double x2, double xjac2, 
			    const event_type& p, weight_type *res)
  { 
    double al = alpha();
    static double loop[3] = {0.0, 0.0, 0.0};
    static su3_kp_i1 kp_parton[3];
    static su3_kp_i2 kp_photon[3];
    
    _M_ip.calculate(p);
    amp_kp(al, _M_q2g2p1, _M_q4p1, kp_parton);
    amp_kp(al, _M_q2g3,   _M_q4g1, kp_photon);

    if(_M_mchel) amp_1loop_mch(_M_q2g2p1, _M_q4p1, loop);
    else amp_1loop(_M_q2g2p1, _M_q4p1, loop);
    
    double s = p[hadron(-1)]*p[hadron(0)];
    double e1 = (p[-1]*p[hadron( 0)])/s;
    double e2 = (p[ 0]*p[hadron(-1)])/s;
    
    //---- convolutions ----
    conv_photon(e1, x1, xjac1, al, kp_photon, res);
    conv_parton(e2, x2, xjac2, al, kp_parton, res);
    
    for(unsigned int i = 0; i < 3; i++) {
      //---- 1-loop contributions ----
      res[2][i] += kp_parton[i].loop + loop[i];
      
      //---- renormalization scale dependent term ----
      res[6][i] = kp_parton[i].tree*2.0*Gg(Nf);
    }
    
    //---- overall Pi factors ----
    for(unsigned int i = 0; i < 7; i++)
      res[i] *= PI_FAC4;    
  }
    
  void photo3jet::dipole_term(const event_type& p, const event_type& dp,
			    int i, int j, int k, weight_type& res) 
  {
    typedef split_fin<lorentzvector<double> > _SplitF;
    typedef split_ini<lorentzvector<double> > _SplitI;
    
    if(i <= 0) {      
      _M_sini = (k <= 0 ? (_SplitI *) &_M_sifi : (_SplitI *) &_M_siff);
      _M_sini -> set(p[i], p[j], p[k]);
    } else {
      _M_sfin = (k <= 0 ? (_SplitF *)  &_M_sffi : (_SplitF *)  &_M_sfff);
      _M_sfin -> set(p[i], p[j], p[k]);
    }
    
    int kt = (k == 4 ? j : k);
    int idx = (i==-1 ? j-1 : 3*i-(i*i-i)/2 + j+3);

    _M_ip.calculate(dp);    
    (this ->* _S_dipole[idx])(kt, i, res);
    res *= PI_FAC5;
  }

  bool photo3jet::dipole_index(int i, int j, int k) {
    if(k > -1) return true;
    else return false;
  }


  //
  //             dipole functions
  //
  namespace _NLO_Pair_Tools {
    typedef std::pair<double, std::complex<double> > _Pair;
    
    inline double operator*(const _Pair& v, const _Pair& cc) { 
      return v.first*cc.first + 2.0*real(v.second*cc.second);
    }
    
    inline _Pair operator*(const _Pair& v, const double& cc) { 
      return _Pair(v.first*cc, v.second*cc);
    }
    
    inline _Pair operator*(const double& cc, const _Pair& v) { 
      return _Pair(v.first*cc, v.second*cc);
    }

    inline _Pair operator/(const _Pair& v, const double& cc) { 
      return _Pair(v.first/cc, v.second/cc);
    }
  }
  
  using namespace _NLO_Pair_Tools;

#define HHC_FQG Vqg(_M_sfin -> Vqg())
#define HHC_FQA Vqa(_M_sfin -> Vqa())
#define HHC_FGG Vgg(_M_sfin -> Vgg())
                   		     
#define HHC_IQG Vqg(_M_sini -> Vqg())
#define HHC_IQQ Vqq(_M_sini -> Vqq())
#define HHC_IGA Vga(_M_sini -> Vga())
#define HHC_IGG Vgg(_M_sini -> Vgg())

#define HHC_CCQG(p1,p2,p3) amp_ccqg(_M_q2g3, kt, i, p1, p2, p3, cc) 
#define HHC_CCAG(p1,p2,p3) amp_ccag(_M_q2g3, kt, i, p1, p2, p3, cc) 
#define HHC_CCAQ(p1,p2,p3) amp_ccaq(_M_q2g3, kt, i, p1, p2, p3, cc) 

#define HHC_CCQ4G1(p1,p2,p3,p4,p5)                      \
  amp_cc(_M_q4g1, kt, i, p1, p2, p3, p4, p5, q4g1) 

#define PHOTO_CCG(p1,p2,p3) amp_ccg(_M_q2g2p1, kt, i, p1, p2, p3) 
#define PHOTO_CCQ(p1,p2,p3) amp_ccq(_M_q2g2p1, kt, i, p1, p2, p3) 

#define SAME_CCQ(Q1, p1,p2,p3) amp_ccq(_M_q4p1, Q1, kt, i, p1, p2, p3, q4) 
#define DIFF_CCQ(Q1,Q2, p1,p2,p3) amp_ccq(_M_q4p1, Q1, Q2, kt, i, p1, p2, p3) 

#define Qu  0.66666666666666666666 
#define Qd -0.33333333333333333333
#define Qu2 0.44444444444444444444
#define Qd2 0.11111111111111111111



  void photo3jet::_M_di1(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp;
    _Pair cc[7], q4g1[2], HHC_IGA;

    HHC_CCAG(2,1,3);
    HHC_CCAQ(1,2,3);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*2.0*Nc*(Vga*cc[1])/2.0;
    d[1] = Qu2*(tmp = 2.0*Nc*(Vga*cc[5])/6.0);
    d[2] = Qd2*tmp;

    //   four quark and one gloun contributions 
    HHC_CCQ4G1(-1,2,3,1,0);
    ueq = 2.0*Nc*(Vga*q4g1[0])/2.0;
    eq  = 2.0*Nc*(Vga*q4g1[1])/4.0;
    
    d[0] += Nu*Qu2*((Nu+2*Nd-1)*ueq + eq) + Nd*Qd2*((Nd-1)*ueq + eq);

    HHC_CCQ4G1(-1,0,2,3,1);
    ueq = 2.0*Nc*(Vga*q4g1[0]);
    eq  = 2.0*Nc*(Vga*q4g1[1])/2.0;
    
    d[1] += Qu2*((Nf-1)*ueq + eq);
    d[2] += Qd2*((Nf-1)*ueq + eq);
  }

  void photo3jet::_M_di2(int kt, int i, weight_type& d) 
  { 
    double ueq, eq;
    _Pair cc[7], q4g1[2], HHC_IGA;

    HHC_CCQG(1,2,3);
   
    d[0] = (Nu*Qu2+Nd*Qd2)*2.0*Nc*(Vga*cc[1])/2.0;
    
    //   four quark and one gloun contributions 
    HHC_CCQ4G1(1,-1,3,2,0);
    ueq = 2.0*Nc*(Vga*q4g1[0])/2.0;
    eq  = 2.0*Nc*(Vga*q4g1[1])/4.0;
    
    d[0] += Nu*Qu2*((Nu+2*Nd-1)*ueq + eq) + Nd*Qd2*((Nd-1)*ueq + eq);

    HHC_CCQ4G1(1,0,-1,3,2);
    ueq = 2.0*Nc*(Vga*q4g1[0]);
    eq  = 2.0*Nc*(Vga*q4g1[1])/2.0;
    
    d[1] = Qu2*((Nu-1)*ueq + eq) + Nd*Qd2*ueq;
    d[2] = Qd2*((Nd-1)*ueq + eq) + Nu*Qu2*ueq;
  }

  void photo3jet::_M_di3(int kt, int i, weight_type& d) 
  { 
    double ueq, eq;
    _Pair q4g1[2], HHC_IGA;
    
    //   four quark and one gloun contributions 
    HHC_CCQ4G1(1,2,-1,3,0);
    ueq = 2.0*Nc*(Vga*q4g1[0])/2.0;
    eq  = 2.0*Nc*(Vga*q4g1[1])/4.0;
    
    d[0] = Nu*Qu2*((Nu-1)*ueq + eq) + Nd*Qd2*((Nd+2*Nu-1)*ueq + eq);

    HHC_CCQ4G1(1,0,2,-1,3);
    ueq = 2.0*Nc*(Vga*q4g1[0]);
    eq  = 2.0*Nc*(Vga*q4g1[1])/2.0;
    
    d[1] = Qu2*((Nu-1)*ueq + eq) + Nd*Qd2*ueq;
    d[2] = Qd2*((Nd-1)*ueq + eq) + Nu*Qu2*ueq;
  }

  void photo3jet::_M_di4(int kt, int i, weight_type& d) 
  {
    double ueq, eq;
    _Pair q4g1[2], HHC_IGA;
    
    //   four quark and one gloun contributions 
    HHC_CCQ4G1(1,2,3,-1,0);
    ueq = 2.0*Nc*(Vga*q4g1[0])/2.0;
    eq  = 2.0*Nc*(Vga*q4g1[1])/4.0;
    
    d[0] = Nu*Qu2*((Nu-1)*ueq + eq) + Nd*Qd2*((Nd+2*Nu-1)*ueq + eq);
    d[1] = d[2] = 0.0;
  }
  
  void photo3jet::_M_d01(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp, q4[2];
    _Pair cc[2], HHC_IGA, HHC_IQQ;
    
    cc[0] = PHOTO_CCQ(2,1,3);
    cc[1] = PHOTO_CCG(2,3,1);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vga*cc[0])/2.0;
    
    //   four quark and one gloun contributions 
    SAME_CCQ(Qu, 2,1,3);
    ueq = (Vga.first*q4[0])/2.0;
    eq  = (Vga.first*q4[1])/4.0;
    
    d[0] += Nu*(Nu-1)*ueq + Nu*eq + (Nd*(Nd-1)*ueq + Nd*eq)/4.0
      +     Nu*Nd*(Vga.first*DIFF_CCQ(Qu, Qd, 2,1,3));    
    
    d[1] = ((Nu-0.5)*Qu2 + Nd*Qd2)*(tmp = Vqq*cc[1]);
    d[2] = ((Nd-0.5)*Qd2 + Nu*Qu2)*tmp;
  }
  
  void photo3jet::_M_d02(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp, q4[2];
    _Pair cc[2], HHC_IGA, HHC_IQG, HHC_IQQ;
    
    cc[0] = PHOTO_CCQ(1,2,3);
    cc[1] = PHOTO_CCG(1,3,2);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vga*cc[0])/2.0;
    d[1] = Qu2*(tmp = Vqg*cc[0])/6.0;
    d[2] = Qd2*tmp/6.0;

    //   four quark and one gloun contributions 
    SAME_CCQ(Qu, 1,3,2);
    ueq = (Vga.first*q4[0])/2.0;
    eq  = (Vga.first*q4[1])/4.0;
    
    d[0] += Nu*(Nu-1)*ueq + Nu*eq + (Nd*(Nd-1)*ueq + Nd*eq)/4.0
      +     Nu*Nd*(Vga.first*DIFF_CCQ(Qu, Qd, 1,3,2));    
    
    d[1] += 0.5*Qu2*(Vqq*cc[1]);
    d[2] += 0.5*Qd2*(Vqq*cc[1]);
  }

  void photo3jet::_M_d03(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp, q4[2];
    _Pair cc[2], HHC_IGG, HHC_IGA, HHC_IQG;
    
    cc[0] = PHOTO_CCG(1,2,3);
    cc[1] = PHOTO_CCQ(1,2,3);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vgg*cc[0])/2.0;
    d[1] = Qu2*(tmp = Vqg*cc[1])/6.0;
    d[2] = Qd2*tmp/6.0;

    //   four quark and one gloun contributions 
    SAME_CCQ(Qu, 3,2,1);
    ueq = (Vga.first*q4[0])/2.0;
    eq  = (Vga.first*q4[1])/4.0;
    
    d[0] += Nu*(Nu-1)*ueq + Nu*eq + (Nd*(Nd-1)*ueq + Nd*eq)/4.0
      +     Nu*Nd*(Vga.first*DIFF_CCQ(Qd, Qu, 3,2,1));    
  }
  
  void photo3jet::_M_d04(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp, q4[2];
    _Pair cc[2], HHC_IGG, HHC_IGA, HHC_IQG;
    
    cc[0] = PHOTO_CCG(1,2,3);
    cc[1] = PHOTO_CCQ(1,2,3);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vgg*cc[0])/2.0;
    d[1] = Qu2*(tmp = Vqg*cc[1])/6.0;
    d[2] = Qd2*tmp/6.0;

    //   four quark and one gloun contributions 
    SAME_CCQ(Qu, 3,1,2);
    ueq = (Vga.first*q4[0])/2.0;
    eq  = (Vga.first*q4[1])/4.0;
    
    d[0] += Nu*(Nu-1)*ueq + Nu*eq + (Nd*(Nd-1)*ueq + Nd*eq)/4.0
      +     Nu*Nd*(Vga.first*DIFF_CCQ(Qd, Qu, 3,1,2));
    
    SAME_CCQ(Qu, 1,2,3);
    ueq = Vqg.first*q4[0];
    eq  = (Vqg.first*q4[1])/2.0;
    
    d[1] += (Nu-1)*ueq + eq + Nd*(Vqg.first*DIFF_CCQ(Qu, Qd, 1,2,3));
    d[2] += ((Nd-1)*ueq + eq)/4.0 + Nu*(Vqg.first*DIFF_CCQ(Qd, Qu, 1,2,3));
  }
  
  void photo3jet::_M_d12(int kt, int i, weight_type& d) 
  {
    double tmp;
    _Pair cc[2], HHC_FQG, HHC_FQA;
    
    cc[0] = PHOTO_CCQ(1,2,3);
    cc[1] = PHOTO_CCG(3,2,1);
    
    d[1] = Qu2*(tmp = Vqg*cc[0])/6.0;
    d[2] = Qd2*tmp/6.0;

    //   four quark and one gloun contributions
    d[0] = 0.5*(Nu*(Nu-0.5)*Qu2 + Nd*(Nd+2.0*Nu-0.5)*Qd2)*(Vqa*cc[1]);
  }

  void photo3jet::_M_d13(int kt, int i, weight_type& d) 
  {
    double tmp;
    _Pair cc[3], HHC_FQG, HHC_FQA;
    
    cc[0] = PHOTO_CCG(1,2,3);
    cc[1] = PHOTO_CCQ(1,2,3);
    cc[2] = PHOTO_CCQ(2,1,3);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vqg*cc[0])/2.0;
    d[1] = Qu2*(tmp = Vqg*cc[1])/6.0;
    d[2] = Qd2*tmp/6.0;    
    
    //   four quark and one gloun contributions 
    d[1] += 0.5*Qu2*Vqa*cc[2];
    d[2] += 0.5*Qd2*Vqa*cc[2];
  }

  void photo3jet::_M_d14(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp, q4[2];
    _Pair cc[3], HHC_FQG, HHC_FQA;
    
    cc[0] = PHOTO_CCG(1,2,3);
    cc[1] = PHOTO_CCQ(1,2,3);
    cc[2] = PHOTO_CCG(3,2,1);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vqg*cc[0])/2.0;
    d[1] = Qu2*(tmp = Vqg*cc[1])/6.0;
    d[2] = Qd2*tmp/6.0;        

    //   four quark and one gloun contributions 
    SAME_CCQ(Qu, 1,2,3);
    ueq = Vqg.first*q4[0];
    eq  = (Vqg.first*q4[1])/2.0; 

    d[0] += (Nu*Qu2 + Nd*Qd2)*(Vqa*cc[2])/4.0;    
    d[1] += (Nu-1)*ueq + eq + Nd*(Vqg.first*DIFF_CCQ(Qu, Qd, 1,2,3));
    d[2] += ((Nd-1)*ueq + eq)/4.0 + Nu*(Vqg.first*DIFF_CCQ(Qd, Qu, 1,2,3));
  }

  void photo3jet::_M_d23(int kt, int i, weight_type& d) 
  {
    double tmp;
    _Pair cc[3], HHC_FQG, HHC_FGG, HHC_FQA;
    
    cc[0] = PHOTO_CCG(1,2,3);
    cc[1] = PHOTO_CCQ(1,2,3);
    cc[2] = PHOTO_CCG(1,3,2);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vqg*cc[0])/2.0;
    d[1] = Qu2*(tmp = Vgg*cc[1])/6.0;
    d[2] = Qd2*tmp/6.0;        
    
    //   four quark and one gloun contributions 
    d[0] += (Nu*Qu2 + Nd*Qd2)*(Vqa*cc[2])/4.0;
    d[1] += (Nu-0.5 + Nd)*Qu2*Vqa*cc[1];
    d[2] += (Nd-0.5 + Nu)*Qd2*Vqa*cc[1];
  }

  void photo3jet::_M_d24(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp, q4[2];
    _Pair cc[2], HHC_FQG, HHC_FGG;
    
    cc[0] = PHOTO_CCG(1,2,3);
    cc[1] = PHOTO_CCQ(1,2,3);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vqg*cc[0])/2.0;
    d[1] = Qu2*(tmp = Vgg*cc[1])/6.0;
    d[2] = Qd2*tmp/6.0;        
    
    //   four quark and one gloun contributions 
    SAME_CCQ(Qu, 1,2,3);
    ueq = Vqg.first*q4[0];
    eq  = (Vqg.first*q4[1])/2.0; 
    
    d[1] += (Nu-1)*ueq + eq + Nd*(Vqg.first*DIFF_CCQ(Qu, Qd, 1,2,3));
    d[2] += ((Nd-1)*ueq + eq)/4.0 + Nu*(Vqg.first*DIFF_CCQ(Qd, Qu, 1,2,3));
  }

  void photo3jet::_M_d34(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp, q4[2];
    _Pair cc[2], HHC_FGG, HHC_FQG, HHC_FQA;
    
    cc[0] = PHOTO_CCG(1,2,3);
    cc[1] = PHOTO_CCQ(1,2,3);
	
    d[0] = (Nu*Qu2+Nd*Qd2)*(Vgg*cc[0])/2.0;
    d[1] = Qu2*(tmp = Vgg*cc[1])/6.0;
    d[2] = Qd2*tmp/6.0;        
    
    //   four quark and one gloun contributions 
    SAME_CCQ(Qu, 1,2,3); 
    ueq = Vqg.first*q4[0];
    eq  = (Vqg.first*q4[1])/2.0; 
    
    d[0] += 0.5*(Nu*(Nu+2.0*Nd-0.5)*Qu2 + Nd*(Nd-0.5)*Qd2)*(Vqa*cc[0]);
    d[1] += (Nu-1)*ueq + eq + Nd*(Vqg.first*DIFF_CCQ(Qu, Qd, 1,2,3));
    d[2] += ((Nd-1)*ueq + eq)/4.0 + Nu*(Vqg.first*DIFF_CCQ(Qd, Qu, 1,2,3));
  }
}  //  namespace nlo

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
#include "photo2jet.h"
#include "ampq2g1p1.h"
#include "ampq2g2p1.h"
#include "ampq4p1.h"
#include "ampq2g2.h"
#include "ampq4.h"
#include "bits/nlo-color.h"


//   PI_FACn = (2pi)^(2n-3) * 2^(n-1)
#define PI_FAC3 992.20085376959424563168
#define PI_FAC4 78341.03930503205203742586




namespace nlo {

  photo2jet::dipole_func photo2jet::_S_dipole[9] = 
  { &photo2jet::_M_di1, &photo2jet::_M_di2, &photo2jet::_M_di3, 
    &photo2jet::_M_d01, &photo2jet::_M_d02, &photo2jet::_M_d03,
    &photo2jet::_M_d12, &photo2jet::_M_d13, &photo2jet::_M_d23};

 
  photo2jet::photo2jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_photo(2, 1, nu, nd, al), _photo_jet_base(nu, nd), _M_mchel(mchel) 
  {
    _M_q2g1p1 = new ampq2g1p1(_M_ip, rng);
    _M_q2g2p1 = new ampq2g2p1(_M_ip, rng);
    _M_q4p1   = new ampq4p1  (_M_ip, rng);

    _M_q2g2 = new ampq2g2(_M_ip, rng);
    _M_q4   = new ampq4  (_M_ip, rng);
  }
  
  photo2jet::~photo2jet() 
  {
    if(_M_q2g1p1) delete _M_q2g1p1;
    if(_M_q2g2p1) delete _M_q2g2p1;
    if(_M_q4p1)   delete _M_q4p1;

    if(_M_q2g2) delete _M_q2g2;
    if(_M_q4)   delete _M_q4;
  }
  
  void photo2jet::born_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_q2g1p1, res.begin());
    res *= PI_FAC3;
  }
  
  void photo2jet::real_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_q2g2p1, _M_q4p1, res.begin());
    res *= PI_FAC4;
  }
  

  void photo2jet::fini_term(double x1, double xjac1, double x2, double xjac2, 
			    const event_type& p, weight_type *res)
  {
    double al = alpha();
    static double loop[3] = {0.0, 0.0, 0.0};
    static su3_kp_i1 kp_parton[3];
    static su3_kp_i2 kp_photon[3];
    
    _M_ip.calculate(p);
    amp_kp(al, _M_q2g1p1, kp_parton);
    amp_kp(al, _M_q2g2, _M_q4, kp_photon);
    
    if(_M_mchel) amp_1loop_mch(_M_q2g1p1, loop);
    else amp_1loop(_M_q2g1p1, loop);
    
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
      res[6][i] = kp_parton[i].tree*Gg(Nf);
    }
    
    //---- overall Pi factors ----
    for(unsigned int i = 0; i < 7; i++)
      res[i] *= PI_FAC3;    
  }
  
  void photo2jet::dipole_term(const event_type& p, const event_type& dp,
			    int i, int j, int k, weight_type& res) 
  {
    typedef split_fin<lorentzvector<double> > _SplitF;
    typedef split_ini<lorentzvector<double> > _SplitI;
    
    if(i <= 0) {      
      _M_sini = (k <= 0 ? (_SplitI *) &_M_sifi : (_SplitI *) &_M_siff);
      _M_sini -> set(p[i], p[j], p[k]);
    } else {
      _M_sfin = (k <= 0 ? (_SplitF *) &_M_sffi : (_SplitF *) &_M_sfff);
      _M_sfin -> set(p[i], p[j], p[k]);
    }
    
    int kt = (k == 3 ? j : k);
    int idx = (i==-1 ? j-1 : 2*i-(i*i-i)/2 + j+2);
    
    _M_ip.calculate(dp); 
    (this ->* _S_dipole[idx])(kt, i, res);
    res *= PI_FAC4;
  }
  
  bool photo2jet::dipole_index(int i, int j, int k) {
    if(k > -1) return true;
    else return false;
  }

	 
  //
  //             dipole functions
  //
#define HHC_FQG Vqg = (_M_sfin -> Vqg()).first
#define HHC_FQA Vqa = (_M_sfin -> Vqa()).first
#define HHC_FGG Vgg = (_M_sfin -> Vgg()).first
                   
#define HHC_IQG Vqg = (_M_sini -> Vqg()).first
#define HHC_IQQ Vqq = (_M_sini -> Vqq()).first
#define HHC_IGA Vga = (_M_sini -> Vga()).first
#define HHC_IGG Vgg = (_M_sini -> Vgg()).first

#define HHC_CCQG(p1,p2) amp_ccqg(_M_q2g2, kt, i, p1, p2, cc) 
#define HHC_CCAG(p1,p2) amp_ccag(_M_q2g2, kt, i, p1, p2, cc) 
#define HHC_CCAQ(p1,p2) amp_ccaq(_M_q2g2, kt, i, p1, p2, cc) 
#define HHC_CCQ4(p1,p2,p3,p4) amp_cc(_M_q4, kt, i, p1, p2, p3, p4, q4) 

#define PHOTO_CCG(p1,p2) amp_ccg(_M_q2g1p1, kt, i, p1, p2) 
#define PHOTO_CCQ(p1,p2) amp_ccq(_M_q2g1p1, kt, i, p1, p2) 

#define Qu2 0.44444444444444444444
#define Qd2 0.11111111111111111111

  void photo2jet::_M_di1(int kt, int i, weight_type& d) 
  {
    double ueq, eq, tmp;
    double cc[7], q4[2], HHC_IGA;

    HHC_CCAG(2,1);
    HHC_CCAQ(1,2);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*2.0*Nc*(Vga*cc[1]);
    d[1] = Qu2*(tmp = 2.0*Nc*(Vga*cc[5])/2.0);
    d[2] = Qd2*tmp;

    //   four quark and one gloun contributions 
    HHC_CCQ4(-1,0,2,1);
    ueq = 2.0*Nc*(Vga*q4[0]);
    eq  = 2.0*Nc*(Vga*q4[1])/2.0;
    
    d[1] += Qu2*((Nf-1)*ueq + eq);
    d[2] += Qd2*((Nf-1)*ueq + eq);
  }

  void photo2jet::_M_di2(int kt, int i, weight_type& d) 
  { 
    double ueq, eq;
    double cc[7], q4[2], HHC_IGA;

    HHC_CCQG(1,2);
   
    d[0] = (Nu*Qu2+Nd*Qd2)*2.0*Nc*(Vga*cc[1]);
    
    //   four quark and one gloun contributions 
    HHC_CCQ4(1,0,-1,2);
    ueq = 2.0*Nc*(Vga*q4[0]);
    eq  = 2.0*Nc*(Vga*q4[1])/2.0;
    
    d[1] = Qu2*((Nu-1)*ueq + eq) + Nd*Qd2*ueq;
    d[2] = Qd2*((Nd-1)*ueq + eq) + Nu*Qu2*ueq;
  }

  void photo2jet::_M_di3(int kt, int i, weight_type& d) 
  { 
    double ueq, eq;
    double q4[2], HHC_IGA;
    
    //   four quark and one gloun contributions 
    HHC_CCQ4(1,0,2,-1);
    ueq = 2.0*Nc*(Vga*q4[0]);
    eq  = 2.0*Nc*(Vga*q4[1])/2.0;
    
    d[0] = 0.0;
    d[1] = Qu2*((Nu-1)*ueq + eq) + Nd*Qd2*ueq;
    d[2] = Qd2*((Nd-1)*ueq + eq) + Nu*Qu2*ueq;
  }

  void photo2jet::_M_d01(int kt, int i, weight_type& d) 
  {
    double cc[2], HHC_IGA, HHC_IQQ;
    
    cc[0] = PHOTO_CCQ(2,1);
    cc[1] = PHOTO_CCG(2,1);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*Vga*cc[0];

    //  four quark dipoles
    d[1] = ((Nu-0.5)*Qu2 + Nd*Qd2)*Vqq*cc[1];
    d[2] = ((Nd-0.5)*Qd2 + Nu*Qu2)*Vqq*cc[1];
  }
  
  void photo2jet::_M_d02(int kt, int i, weight_type& d) 
  {
    double tmp, cc[2], HHC_IGA, HHC_IQG, HHC_IQQ;
    
    cc[0] = PHOTO_CCQ(1,2);
    cc[1] = PHOTO_CCG(1,2);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*Vga*cc[0];
    d[1] = Qu2*(tmp = Vqg*cc[0]/2.0);
    d[2] = Qd2*tmp;

    //  four quark dipoles
    d[1] += 0.5*Qu2*Vqq*cc[1];
    d[2] += 0.5*Qd2*Vqq*cc[1];
  }
  
  void photo2jet::_M_d03(int kt, int i, weight_type& d) 
  {
    double tmp, HHC_IGG, HHC_IQG;
        
    d[0] = (Nu*Qu2+Nd*Qd2)*Vgg*PHOTO_CCG(1,2);
    d[1] = Qu2*(tmp = Vqg*PHOTO_CCQ(1,2)/2.0);
    d[2] = Qd2*tmp;
  }
  
  void photo2jet::_M_d12(int kt, int i, weight_type& d) 
  { 
    double tmp, HHC_FQG;

    d[0] = 0.0;
    d[1] = Qu2*(tmp = Vqg*PHOTO_CCQ(1,2)/2.0);
    d[2] = Qd2*tmp;
  }
  
  void photo2jet::_M_d13(int kt, int i, weight_type& d) 
  {
    double tmp, cc[3], HHC_FQG, HHC_FQA;
    
    cc[0] = PHOTO_CCG(1,2);
    cc[1] = PHOTO_CCQ(1,2);
    cc[2] = PHOTO_CCQ(2,1);
   
    d[0] = (Nu*Qu2+Nd*Qd2)*Vqg*cc[0];    
    d[1] = Qu2*(tmp = Vqg*cc[1]/2.0);
    d[2] = Qd2*tmp;

    //  four quark dipoles
    d[1] += 0.5*Qu2*Vqa*cc[2];
    d[2] += 0.5*Qd2*Vqa*cc[2];
  }
  
  void photo2jet::_M_d23(int kt, int i, weight_type& d) 
  {
     double tmp, cc[2], HHC_FQG, HHC_FGG, HHC_FQA;
    
    cc[0] = PHOTO_CCG(1,2);
    cc[1] = PHOTO_CCQ(1,2);
    
    d[0] = (Nu*Qu2+Nd*Qd2)*Vqg*cc[0];    
    d[1] = Qu2*(tmp = Vgg*cc[1]/2.0);
    d[2] = Qd2*tmp;
   
    //  four quark dipoles
    d[1] += (Nu-0.5 + Nd)*Qu2*Vqa*cc[1];
    d[2] += (Nd-0.5 + Nu)*Qd2*Vqa*cc[1];
  }
}  //  namespace nlo

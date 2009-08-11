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
#include "dis3jet.h"
#include "ampq2g2l2.h"
#include "ampq2g3l2.h"
#include "ampq4l2.h"
#include "ampq4g1l2.h"
#include "bits/nlo-color.h"


#define PI_FAC3 (4.0*984462.53422111174351397439)
#define PI_FAC4 (4.0*77730046.08365014226587406347)


namespace nlo {
  
  
  dis3jet::dipole_func dis3jet::_S_dipole[10] = 
  { &dis3jet::_M_d01, &dis3jet::_M_d02, &dis3jet::_M_d03, &dis3jet::_M_d04,
    &dis3jet::_M_d12, &dis3jet::_M_d13, &dis3jet::_M_d14, &dis3jet::_M_d23,
    &dis3jet::_M_d24, &dis3jet::_M_d34};
    
  dis3jet::dis3jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_dis(3, 2, nu, nd, al), _dis_jet_base(nu, nd), _M_mchel(mchel) 
  {
    _M_q2g2 = new ampq2g2l2(_M_ip, rng);  
    _M_q4   = new   ampq4l2(_M_ip, rng);
    _M_q2g3 = new ampq2g3l2(_M_ip, rng);
    _M_q4g1 = new ampq4g1l2(_M_ip, rng);
  }

  
  dis3jet::~dis3jet() {
    if(_M_q2g2) delete _M_q2g2;
    if(_M_q4)   delete _M_q4;
    if(_M_q2g3) delete _M_q2g3;
    if(_M_q4g1) delete _M_q4g1;
  }
  
  void dis3jet::born_term(const event_type& p, weight_type& res) 
  {
    double amp[3];
    _M_ip.calculate(p);
    
    if(_M_mchel) amp_tree_mch(_M_q2g2, _M_q4, amp);
    else amp_tree(_M_q2g2, _M_q4, amp);
    
    res[0] = PI_FAC3*amp[0];
    res[1] = PI_FAC3*amp[1];
    res[2] = PI_FAC3*amp[2];
  }

  void dis3jet::real_term(const event_type& p, weight_type& res) 
  {
    double amp[3];
    _M_ip.calculate(p);
    amp_tree(_M_q2g3, _M_q4g1, amp);
    
    res[0] = PI_FAC4*amp[0];
    res[1] = PI_FAC4*amp[1];
    res[2] = PI_FAC4*amp[2];
  }
 
  void dis3jet::
  fini_term(double x, double xjac, const event_type& p, weight_type *res)
  {
    double loop[3];
    su3_kp_i1 kp[3];
    _M_ip.calculate(p);
    
    if(_M_mchel) {
      amp_kp_mch(alpha(), _M_q2g2, _M_q4, kp);
      amp_1loop_mch(_M_q2g2, _M_q4, loop);
    } else {
      amp_kp(alpha(), _M_q2g2, _M_q4, kp);
      amp_1loop(_M_q2g2, _M_q4, loop);
    }
   
    //---- convolutions ----
    double eta = (p[-1]*p[0])/(p[-1]*p[hadron(0)]);
    convolution(eta, x, xjac, alpha(), kp, res);

    for(unsigned int i = 0; i < 3; i++) {
      //---- 1-loop contributions ----
      res[1][i] += kp[i].loop + loop[i];
      
      //---- renormalization scale dependent term ----
      res[4][i] = kp[i].tree*2.0*Gg(Nf);
    }

    //---- overall Pi factors ----
    for(unsigned int i = 0; i < 5; i++)
      res[i] *= PI_FAC3;
  }   
  
  void dis3jet::dipole_term(const event_type& p, const event_type& dp,
			    int i, int j, int k, weight_type& res) 
  {
    _M_ip.calculate(dp);
 
    if(i == 0) _M_siff.set(p[i], p[j], p[k]);
    else {
      typedef split_fin<lorentzvector<double> > _SplitF;
      _M_sfin = (k == 0 ? (_SplitF *) &_M_sffi : (_SplitF *) &_M_sfff);
      _M_sfin -> set(p[i], p[j], p[k]);
    }
    
    int kt = (k == 4 ? j : k);
    int idx = (i == 0 ? j-1 : 3*i-(i*i-i)/2 + j-1);

    (this ->* _S_dipole[idx])(kt, i, res); 
    res *= PI_FAC4;
  }
  
  //
  //             dipole functions
  //
  typedef std::pair<double, std::complex<double> > _Pair;
  
  
#define DIS3_FQG Vqg(_M_sfin -> Vqg())
#define DIS3_FQA Vqa(_M_sfin -> Vqa())
#define DIS3_FGG Vgg(_M_sfin -> Vgg())
		       
#define DIS3_IQG Vqg(_M_siff.Vqg())
#define DIS3_IQQ Vqq(_M_siff.Vqq())
#define DIS3_IGA Vga(_M_siff.Vga())
#define DIS3_IGG Vgg(_M_siff.Vgg())

#define DIS3_CCG(amp,p1,p2,p3) amp_ccg(amp, kt, i, p1, p2, p3, cc) 
#define DIS3_CCQ(amp,p1,p2,p3) amp_ccq(amp, kt, i, p1, p2, p3, cc) 
#define DIS3_CCA(amp,p1,p2,p3) amp_cca(amp, kt, i, p1, p2, p3, cc) 



  void dis3jet::_M_d12(int kt, int i, weight_type& d) 
  {
    _Pair cc[3], DIS3_FQG, DIS3_FQA;
    
    DIS3_CCG(_M_q2g2, 3,2,1);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = Nf*Dijk(Vqa, cc[0])/4.0;
    d[1] = Dijk(Vqg, cc[1])/6.0;
    d[2] = Dijk(Vqg, cc[2])/6.0;
  }


  void dis3jet::_M_d13(int kt, int i, weight_type& d) 
  {
    _Pair cc[3], DIS3_FQG, DIS3_FQA;
    
    DIS3_CCG(_M_q2g2, 1,2,3);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = Dijk(Vqg, cc[0])/2.0;
    d[1] = Dijk(Vqg, cc[1])/6.0;
    d[2] = Dijk(Vqg, cc[2])/6.0;
    
    DIS3_CCQ(_M_q2g2, 2,1,3);
    
    d[1] += Nf*Dijk(Vqa, cc[1])/2.0;
    d[2] += Nf*Dijk(Vqa, cc[2])/2.0;
  }
  

  void dis3jet::_M_d14(int kt, int i, weight_type& d) 
  {
    _Pair cc[3], DIS3_FQG, DIS3_FQA;
    
    DIS3_CCG(_M_q2g2, 1,2,3);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = Dijk(Vqg, cc[0])/2.0;
    d[1] = Dijk(Vqg, cc[1])/6.0;
    d[2] = Dijk(Vqg, cc[2])/6.0;
    
    DIS3_CCG(_M_q2g2, 3,2,1);
    DIS3_CCQ(_M_q4,   1,2,3);
    
    d[0] += Nf*Dijk(Vqa, cc[0])/4.0;
    d[1] += Dijk(Vqg, cc[1])/2.0;
    d[2] += Dijk(Vqg, cc[2])/2.0;
  }

  
  void dis3jet::_M_d23(int kt, int i, weight_type& d) 
  {
    _Pair cc[3], DIS3_FQG, DIS3_FGG, DIS3_FQA;
    
    DIS3_CCG(_M_q2g2, 1,2,3);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = Dijk(Vqg, cc[0])/2.0;
    d[1] = Dijk(Vgg, cc[1])/6.0;
    d[2] = Dijk(Vgg, cc[2])/6.0;
    
    DIS3_CCG(_M_q2g2, 1,3,2);
    
    d[0] += Nf*Dijk(Vqa, cc[0])/4.0;
    d[1] += Nf*Dijk(Vqa, cc[1])/2.0;
    d[2] += Nf*Dijk(Vqa, cc[2])/2.0;
  }


  void dis3jet::_M_d24(int kt, int i, weight_type& d)
  {
    _Pair cc[3], DIS3_FQG, DIS3_FGG;
    
    DIS3_CCG(_M_q2g2, 1,2,3);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = Dijk(Vqg, cc[0])/2.0;
    d[1] = Dijk(Vgg, cc[1])/6.0;
    d[2] = Dijk(Vgg, cc[2])/6.0;
    
    DIS3_CCQ(_M_q4, 1,2,3);
    
    d[1] += Dijk(Vqg, cc[1])/2.0;
    d[2] += Dijk(Vqg, cc[2])/2.0;
  }


  void dis3jet::_M_d34(int kt, int i, weight_type& d) 
  {
    _Pair cc[3], DIS3_FQG, DIS3_FGG, DIS3_FQA;
    
    DIS3_CCG(_M_q2g2, 1,2,3);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = Dijk(Vgg, cc[0])/2.0;
    d[1] = Dijk(Vgg, cc[1])/6.0;
    d[2] = Dijk(Vgg, cc[2])/6.0;
    
    DIS3_CCQ(_M_q4, 1,2,3);
    
    d[0] += Nf*Dijk(Vqa, cc[0])/4.0;
    d[1] += Dijk(Vqg, cc[1])/2.0;
    d[2] += Dijk(Vqg, cc[2])/2.0;
  }
  
  
  void dis3jet::_M_d01(int kt, int i, weight_type& d) 
  {
    _Pair cc[3], DIS3_IGA, DIS3_IQQ;
    
    DIS3_CCG(_M_q2g2, 2,3,1);
    DIS3_CCA(_M_q2g2, 2,3,1);
    
    d[0] = (Nu*Dijk(Vga, cc[1]) + Nd*Dijk(Vga, cc[2]))/2.0;
    d[1] = d[2] = Dijk(Vqq, cc[0])/2.0;
    
    DIS3_CCA(_M_q4, 2,3,1);
    d[0] += (Nu*Dijk(Vga, cc[1]) + Nd*Dijk(Vga,cc[2]))/4.0;
  }
  
  
  void dis3jet::_M_d02(int kt, int i, weight_type& d)
  {
    _Pair cc[3], DIS3_IGA, DIS3_IQG, DIS3_IQQ;
    
    DIS3_CCG(_M_q2g2, 1,3,2);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = (Nu*Dijk(Vga, cc[1]) + Nd*Dijk(Vga, cc[2]))/2.0;
    d[1] = d[2] = Dijk(Vqq, cc[0])/2.0;
    d[1] += Dijk(Vqg, cc[1])/6.0;
    d[2] += Dijk(Vqg, cc[2])/6.0;
    
    DIS3_CCQ(_M_q4, 1,3,2);
    d[0] += (Nu*Dijk(Vga, cc[1]) + Nd*Dijk(Vga, cc[2]))/4.0;
  }
  
  
  void dis3jet::_M_d03(int kt, int i, weight_type& d) 
  {
    _Pair cc[3], DIS3_IQG, DIS3_IGG, DIS3_IGA;
    
    DIS3_CCG(_M_q2g2, 1,2,3);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = Dijk(Vgg, cc[0])/2.0;
    d[1] = Dijk(Vqg, cc[1])/6.0;
    d[2] = Dijk(Vqg, cc[2])/6.0;
    
    DIS3_CCA(_M_q4, 2,1,3);
    d[0] += (Nu*Dijk(Vga, cc[1]) + Nd*Dijk(Vga, cc[2]))/4.0;
  }
  
  
  void dis3jet::_M_d04(int kt, int i, weight_type& d) 
  {
    _Pair cc[3], DIS3_IQG, DIS3_IGG, DIS3_IGA;
    
    DIS3_CCG(_M_q2g2, 1,2,3);
    DIS3_CCQ(_M_q2g2, 1,2,3);
    
    d[0] = Dijk(Vgg, cc[0])/2.0;
    d[1] = Dijk(Vqg, cc[1])/6.0;
    d[2] = Dijk(Vqg, cc[2])/6.0;
    
    DIS3_CCQ(_M_q4, 1,3,2);
    d[0] += (Nu*Dijk(Vga, cc[1]) + Nd*Dijk(Vga, cc[2]))/4.0;
    
    DIS3_CCQ(_M_q4, 1,2,3);
    d[1] += Dijk(Vqg, cc[1])/2.0;
    d[2] += Dijk(Vqg, cc[2])/2.0;
  }
}  //  namespace nlo









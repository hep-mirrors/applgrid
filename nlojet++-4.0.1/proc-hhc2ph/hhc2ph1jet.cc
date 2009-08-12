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
#include "hhc2ph1jet.h"
#include "ampq2g1p2.h"
#include "ampq2g2p2.h"
#include "ampq4p2.h"
#include "bits/nlo-color.h"


//    PI_FACn = (2pi)^(2n+2) * 2^(n+1)
#define PI_FAC1 6234.18182617615598326995
#define PI_FAC2 492231.26711055587177578917



namespace nlo {

  hhc2ph1jet::dipole_func hhc2ph1jet::_S_dipole[5] = 
  { &hhc2ph1jet::_M_di1, &hhc2ph1jet::_M_di2, &hhc2ph1jet::_M_d01, 
    &hhc2ph1jet::_M_d02, &hhc2ph1jet::_M_d12};

 
  hhc2ph1jet::hhc2ph1jet(const random_generator& rng, bool mchel, unsigned int nu, unsigned int nd, double al)
    : process_hhc2ph(1, 1, nu, nd, al), _hhc2ph_jet_base(nu, nd), _M_mchel(mchel) 
  {
    _M_q2g1p2 = new ampq2g1p2(_M_ip, rng);
    _M_q2g2p2 = new ampq2g2p2(_M_ip, rng);
    _M_q4p2   = new ampq4p2  (_M_ip, rng);
  }
  
  hhc2ph1jet::~hhc2ph1jet() 
  {
    if(_M_q2g1p2) delete _M_q2g1p2;
    if(_M_q2g2p2) delete _M_q2g2p2;
    if(_M_q4p2)   delete _M_q4p2;
  }
  
  void hhc2ph1jet::born_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_q2g1p2, res.begin());
    res *= PI_FAC1;
  }
  
  void hhc2ph1jet::real_term(const event_type& p, weight_type& res) 
  {
    _M_ip.calculate(p);
    amp_tree(_M_q2g2p2, _M_q4p2, res.begin());
    res *= PI_FAC2;
  }

  void hhc2ph1jet::fini_term(double x1, double xjac1, double x2, double xjac2, 
			    const event_type& p, weight_type *res)
  {
    double al = alpha();
    static double loop[6];
    static su3_kp_i2 kp[6];

    _M_ip.calculate(p);
    amp_kp(al, _M_q2g1p2, kp);

    if(_M_mchel) amp_1loop_mch(_M_q2g1p2, loop);
    else amp_1loop(_M_q2g1p2, loop);

    double s = p[hadron(-1)]*p[hadron(0)];
    double e1 = (p[-1]*p[hadron( 0)])/s;
    double e2 = (p[ 0]*p[hadron(-1)])/s;

    //---- convolutions ----
    convolutions(e1, x1, xjac1, e2, x2, xjac2, al, kp, res);
    
    //---- 1-loop contributions ----
    res[2][1]  += kp[0].loop + loop[0];
    res[2][2]  += kp[1].loop + loop[1];
    res[2][3]  += kp[2].loop + loop[2];
    res[2][4]  += kp[3].loop + loop[3];
    res[2][12] += kp[4].loop + loop[4];
    res[2][14] += kp[5].loop + loop[5];
      
    //---- renormalization scale dependent term ----
    double gfac = Gg(Nf);
    res[6][1]  = gfac*kp[0].tree;
    res[6][2]  = gfac*kp[1].tree;
    res[6][3]  = gfac*kp[2].tree;
    res[6][4]  = gfac*kp[3].tree;
    res[6][12] = gfac*kp[4].tree;
    res[6][14] = gfac*kp[5].tree;
    res[6][0] = res[6][5] = res[6][6] = res[6][7] = 0.0;
    res[6][8] = res[6][9] = res[6][10] = res[6][11] = 0.0; 
    res[6][13] = res[6][15] = res[6][16] = 0.0;

    //---- overall Pi factors ----
    for(unsigned int i=0; i < 7; i++)
      res[i] *= PI_FAC1;
  }

  void hhc2ph1jet::dipole_term(const event_type& p, const event_type& dp,
			       int i, int j, int k, weight_type& res) 
  {
    typedef split_ini<lorentzvector<double> > _SplitI;
    
    if(i <= 0) {      
      _M_sini = (k <= 0 ? (_SplitI *) &_M_sifi : (_SplitI *) &_M_siff);
      _M_sini -> set(p[i], p[j], p[k]);
    } else _M_sffi.set(p[i], p[j], p[k]);
    
    int kt = (k == 2 ? j : k);
    int idx = (i==-1 ? j-1 : i+j+1);
    
    _M_ip.calculate(dp); 
    (this ->* _S_dipole[idx])(kt, i, res);
    res *= PI_FAC2;
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
  
#define HHC_FQG Vqg(_M_sffi.Vqg())
#define HHC_FQA Vqa(_M_sffi.Vqa())
#define HHC_FGG Vgg(_M_sffi.Vgg())
                                     
#define HHC_IQG Vqg(_M_sini -> Vqg())
#define HHC_IQQ Vqq(_M_sini -> Vqq())
#define HHC_IGA Vga(_M_sini -> Vga())
#define HHC_IGG Vgg(_M_sini -> Vgg())


#define HHC_CCAG amp_ccag(_M_q2g1p2, kt, i) 
#define HHC_CCGA amp_ccga(_M_q2g1p2, kt, i) 
#define HHC_CCQA amp_ccqa(_M_q2g1p2, kt, i) 

#define Qu4 0.19753086419753086419
#define Qd4 0.01234567901234567901

  void hhc2ph1jet::_M_di1(int kt, int i, weight_type& d) 
  {
    double tmp;
    _Pair cc[3], HHC_IGA, HHC_IQG, HHC_IQQ;
    
    cc[0] = HHC_CCAG;
    cc[1] = HHC_CCQA;
    cc[2] = HHC_CCGA;

    d[0] = (Nu*Qu4 + Nd*Qd4)*(Vga*cc[0]);
    d[1] = d[2] = 0.0;
    d[3] = Qu4*(tmp = Vga*cc[1]);
    d[4] = Qd4*tmp;
    d[12] = Qu4*(tmp = 0.5*(Vqg*cc[1]));
    d[14] = Qd4*tmp;

    //----- Four quark subprocesses -----
    d[5] = d[10] = d[11] = d[16] = Qu4*(tmp = Vqq*cc[2]);
    d[6] = Qu4*tmp/2.0;
    d[7] = d[9] = d[13] = d[15] = Qd4*tmp;
    d[8] = Qd4*tmp/2.0;
    d[12] += Qu4*tmp;
    d[14] += Qd4*tmp;
  }
  
  void hhc2ph1jet::_M_di2(int kt, int i, weight_type& d) 
  { 
    double tmp;
    _Pair cc[3], HHC_IGA, HHC_IQG, HHC_IGG, HHC_IQQ;
    
    cc[0] = HHC_CCAG;
    cc[1] = HHC_CCGA;
    cc[2] = HHC_CCQA;
    
    d[0] = (Nu*Qu4 + Nd*Qd4)*(Vga*cc[0]);
    d[1] = Qu4*(tmp = Vqg*cc[0]);
    d[2] = Qd4*tmp;
    d[3] = Qu4*(tmp = Vgg*cc[1]);
    d[4] = Qd4*tmp;
    d[12] = Qu4*(tmp = 0.5*(Vqg*cc[2]));
    d[14] = Qd4*tmp;

    //----- Four quark subprocesses -----
    d[5] = d[7] = d[9] = d[10] = d[11] = d[13] = d[15] = d[16] = 0.0;
    d[6] = Qu4*(tmp = Vqq*cc[1])/2.0;
    d[8] = Qd4*tmp/2.0;
  }
  
  void hhc2ph1jet::_M_d01(int kt, int i, weight_type& d) 
  {
    double tmp;
    _Pair cc[3], HHC_IGA, HHC_IQG, HHC_IQQ;
    
    cc[0] = HHC_CCGA;
    cc[1] = HHC_CCQA;
    cc[2] = HHC_CCAG;

    d[0] = (Nu*Qu4 + Nd*Qd4)*(Vga*cc[0]);
    d[1] = Qu4*(tmp = Vga*cc[1]);
    d[2] = Qd4*tmp;
    d[3] = d[4] = 0.0;
    d[12] = Qu4*(tmp = 0.5*(Vqg*cc[1]));
    d[14] = Qd4*tmp;

    //----- Four quark subprocesses -----
    d[5] = d[7] = d[9] = d[10] = d[11] = d[13] = d[15] = d[16] = 0.0;
    d[6] = Qu4*(tmp = Vqq*cc[2])/2.0;
    d[8] = Qd4*tmp/2.0;
  }
  
  void hhc2ph1jet::_M_d02(int kt, int i, weight_type& d) 
  { 
    double tmp;
    _Pair cc[3], HHC_IGA, HHC_IGG, HHC_IQG, HHC_IQQ;
    
    cc[0] = HHC_CCGA;
    cc[1] = HHC_CCAG;
    cc[2] = HHC_CCQA;

    d[0] = (Nu*Qu4 + Nd*Qd4)*(Vga*cc[0]);
    d[1] = Qu4*(tmp = Vgg*cc[1]);
    d[2] = Qd4*tmp;
    d[3] = Qu4*(tmp = Vqg*cc[0]);
    d[4] = Qd4*tmp;
    d[12] = Qu4*(tmp = 0.5*(Vqg*cc[2]));
    d[14] = Qd4*tmp;

    //----- Four quark subprocesses -----
    d[5] = d[9] = d[11] = d[15] = Qu4*(tmp = Vqq*cc[1]);
    d[6] = Qu4*tmp/2.0;
    d[7] = d[10] = d[13] = d[16] = Qd4*tmp;
    d[8] = Qd4*tmp/2.0;
    d[12] += Qu4*tmp;
    d[14] += Qd4*tmp;
  }
    
  void hhc2ph1jet::_M_d12(int kt, int i, weight_type& d) 
  {
    double tmp;
    _Pair cc[3], HHC_FQG, HHC_FGG, HHC_FQA;
    
    cc[0] = HHC_CCAG;
    cc[1] = HHC_CCGA;
    cc[2] = HHC_CCQA;

    d[0] = 0.0;
    d[1] = Qu4*(tmp = Vqg*cc[0]);
    d[2] = Qd4*tmp;
    d[3] = Qu4*(tmp = Vqg*cc[1]);
    d[4] = Qd4*tmp;
    d[12] = Qu4*(tmp = 0.5*(Vgg*cc[2]));
    d[14] = Qd4*tmp;

    //----- Four quark subprocesses -----
    d[5] = d[6] = d[7] = d[8] = d[9] = d[10] = 0.0;
    d[11] = d[13] = d[15] = d[16] = 0.0;
    d[12] += Nf*Qu4*(tmp = Vqa*cc[2]);
    d[14] += Nf*Qd4*tmp;
  }
}  //  namespace nlo

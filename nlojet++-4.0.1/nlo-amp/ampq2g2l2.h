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
#ifndef __AMPQ2G2L2_H__
#define __AMPQ2G2L2_H__ 1

//   Standard includes
#include <utility>

//   nlo includes
#include <bits/amp-ampbase.h>
#include <bits/nlo-color.h>



namespace nlo {


  class ampq2g2l2 : private _Amp_base 
  {
    //   private types
    typedef std::complex<double> _ComplexD;
    typedef std::pair<double, _ComplexD> _Pair;

  public:
    //   constructor
    ampq2g2l2(const innerprod_type& __x, const random_generator& rng)
      : _Amp_base(__x, rng) {}

    //  tree and 1-loop level helicity amplitudes (SU(3) color)
    void matrix_tree_pppm(int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_ppmm(int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmpm(int, int, int, int, int, int, _ComplexD *) const;
    void matrix_tree_pmmm(int, int, int, int, int, int, _ComplexD *) const;
   
    void matrix_1loop_pppm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_1loop_ppmm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_1loop_pmpm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    void matrix_1loop_pmmm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;

    //  tree and 1-loop level helicity amplitudes 
    //  (color independent decomposation)
    void color_tree_pppm(int, int, int, int, int, int, _ComplexD *) const;
    void color_tree_ppmm(int, int, int, int, int, int, _ComplexD *) const;
    void color_tree_pmpm(int, int, int, int, int, int, _ComplexD *) const;
    void color_tree_pmmm(int, int, int, int, int, int, _ComplexD *) const;
    
    void color_1loop_pppm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    void color_1loop_ppmm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    void color_1loop_pmpm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    void color_1loop_pmmm(unsigned int, int, int, int, int, int, int, _ComplexD *) const;
    
    //  tree and 1-loop amplitudes (SU(3) color)
    double su3_tree(int, int, int, int, int, int) const;
    double su3_tree_mch(int, int, int, int, int, int) const;

    double su3_1loop(unsigned int, int, int, int, int, int, int) const;
    double su3_1loop_mch(unsigned int, int, int, int, int, int, int) const;

    //  color corraleted amplitudes
    _Pair su3_cc(int, int, int, int, int, int, int, int) const;

    //   give the I terms
    double su3_ins(unsigned int, int, int, int, int, int, int, double = 1.0) const;
    double su3_ins_mch(unsigned int, int, int, int, int, int, int, double = 1.0) const;

    //   give the K, P terms
    void su3_kp(unsigned int, int, int, int, int, int, int, int, su3_kp_i1 *, double = 1.0) const;
    void su3_kp_mch(unsigned int, int, int, int, int, int, int, int, su3_kp_i1 *, double = 1.0) const;

  private:
    //        private members 
    static double su3_amptree(_ComplexD *);
    static void   su3_ampcc(_ComplexD *, double *);
    static void   su3_ampcc(_ComplexD *, _ComplexD *, double *);
    static void   su3_ampcc(_ComplexD *, _ComplexD *, unsigned int, _Pair *);
    static void   su3_ampcc(_ComplexD *, _ComplexD *, unsigned int, double *);
    static double su3_amploop(_ComplexD *);

    _ComplexD As1pppm   (int, int, int, int, int, int) const;
    _ComplexD As1pmmm   (int, int, int, int, int, int) const;

    _ComplexD Atree1pppm(int, int, int, int, int, int) const;
    _ComplexD Fcc1pppm  (int, int, int, int, int, int) const;
    _ComplexD Fsc1pppm  (int, int, int, int, int, int) const;

    _ComplexD Atree1ppmm(int, int, int, int, int, int) const;
    _ComplexD Fcc1ppmm  (int, int, int, int, int, int) const;
    _ComplexD Fsc1ppmm  (int, int, int, int, int, int) const;

    _ComplexD Atree1pmpm(int, int, int, int, int, int) const;
    _ComplexD Fcc1pmpm  (int, int, int, int, int, int) const;
    _ComplexD Fsc1pmpm  (int, int, int, int, int, int) const;

    _ComplexD Atree2ppmp(int, int, int, int, int, int) const;
    _ComplexD Fcc2ppmp  (int, int, int, int, int, int) const;
    _ComplexD Fsc2ppmp  (int, int, int, int, int, int) const;

    _ComplexD Atree2ppmm(int, int, int, int, int, int) const;
    _ComplexD Fcc2ppmm  (int, int, int, int, int, int) const;
    _ComplexD Fsc2ppmm  (int, int, int, int, int, int) const;

    _ComplexD Atree3pmpp(int, int, int, int, int, int) const;
    _ComplexD Fcc3pmpp  (int, int, int, int, int, int) const;
    _ComplexD Fsc3pmpp  (int, int, int, int, int, int) const;

    _ComplexD Atree3pmmp(int, int, int, int, int, int) const;
    _ComplexD Fcc3pmmp  (int, int, int, int, int, int) const;
    _ComplexD Fsc3pmmp  (int, int, int, int, int, int) const;

    _ComplexD Atree3pmpm(int, int, int, int, int, int) const;
    _ComplexD Fcc3pmpm  (int, int, int, int, int, int) const;
    _ComplexD Fsc3pmpm  (int, int, int, int, int, int) const;

    _ComplexD Atree1pmmm(int, int, int, int, int, int) const;
    _ComplexD Fcc1pmmm  (int, int, int, int, int, int) const;
    _ComplexD Fsc1pmmm  (int, int, int, int, int, int) const;

    _ComplexD Atree2pmmp(int, int, int, int, int, int) const;
    _ComplexD Fcc2pmmp  (int, int, int, int, int, int) const;
    _ComplexD Fsc2pmmp  (int, int, int, int, int, int) const;

    _ComplexD Atree2pmmm(int, int, int, int, int, int) const;
    _ComplexD Fcc2pmmm  (int, int, int, int, int, int) const;
    _ComplexD Fsc2pmmm  (int, int, int, int, int, int) const;

    _ComplexD Atree3pmmm(int, int, int, int, int, int) const;
    _ComplexD Fcc3pmmm  (int, int, int, int, int, int) const;
    _ComplexD Fsc3pmmm  (int, int, int, int, int, int) const;

    _ComplexD Fvs3pmpp  (int, int, int, int, int, int) const;
    _ComplexD Fvf3pmpp  (int, int, int, int, int, int) const;
    _ComplexD Fax3pmpp  (int, int, int, int, int, int) const;
    _ComplexD Faxsl3pmpp(int, int, int, int, int, int) const;
    	     	       
    _ComplexD Fvs3pmpm  (int, int, int, int, int, int) const;
    _ComplexD Fvf3pmpm  (int, int, int, int, int, int) const;
    _ComplexD Fax3pmpm  (int, int, int, int, int, int) const;
    _ComplexD Fax3pmmp  (int, int, int, int, int, int) const;
    _ComplexD Faxsl3pmpm(int, int, int, int, int, int) const;
    
    _ComplexD Fax3pmmm  (int, int, int, int, int, int) const;
    _ComplexD Faxsl3pmmm(int, int, int, int, int, int) const;
    _ComplexD Fvf3pmmm  (int, int, int, int, int, int) const;
    _ComplexD Fvs3pmmm  (int, int, int, int, int, int) const;
    		       
    _ComplexD Faxsl3pmmp(int, int, int, int, int, int) const;
    _ComplexD Fvf3pmmp  (int, int, int, int, int, int) const;
    _ComplexD Fvs3pmmp  (int, int, int, int, int, int) const;
    
    _ComplexD M1        (int, int, int, int, int, int) const;
    _ComplexD M2a       (int, int, int, int, int, int) const;
    _ComplexD M2        (int, int, int, int, int, int) const;
    _ComplexD M3a       (int, int, int, int, int, int) const;
    _ComplexD M3        (int, int, int, int, int, int) const;
    _ComplexD FAcc1ppmm (int, int, int, int, int, int) const;
    _ComplexD FAsc1ppmm (int, int, int, int, int, int) const;
    _ComplexD FAcc1pmpm (int, int, int, int, int, int) const;
    _ComplexD FAsc1pmpm (int, int, int, int, int, int) const;
    _ComplexD T2ppmm    (int, int, int, int, int, int) const;
    _ComplexD TT2ppmm   (int, int, int, int, int, int) const;
    _ComplexD FBcc3pmmp (int, int, int, int, int, int) const;
    _ComplexD FBsc3pmmp (int, int, int, int, int, int) const;
    _ComplexD FBcc3pmpm (int, int, int, int, int, int) const;
    _ComplexD FBsc3pmpm (int, int, int, int, int, int) const;
    _ComplexD Tc3pmpm   (int, int, int, int, int, int) const;
    _ComplexD TTc3pmpm  (int, int, int, int, int, int) const;
    _ComplexD Ts3pmpm   (int, int, int, int, int, int) const;
    _ComplexD TTs3pmpm  (int, int, int, int, int, int) const;

    _ComplexD F1vs3pmpp  (int, int, int, int, int, int) const;
    _ComplexD F2vs3pmpp  (int, int, int, int, int, int) const;
    _ComplexD F2ax3pmpp  (int, int, int, int, int, int) const;
    _ComplexD FCvs3pmpm  (int, int, int, int, int, int) const;
    _ComplexD FBvs3pmpm  (int, int, int, int, int, int) const;
    _ComplexD FCvf3pmpm  (int, int, int, int, int, int) const;
    _ComplexD FBvf3pmpm  (int, int, int, int, int, int) const;
    _ComplexD C1ax       (int, int, int, int, int, int) const;
    _ComplexD Cax        (int, int, int, int, int, int) const;
    _ComplexD FBax3pmpm  (int, int, int, int, int, int) const;
    _ComplexD FBax3pmmp  (int, int, int, int, int, int) const;
    _ComplexD FBaxsl3pmpm(int, int, int, int, int, int) const;
  };

  inline void ampq2g2l2::
  matrix_tree_pppm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const  {
    M[0] = Atree1pppm(p1, p2, p3, p4, p5, p6);
    M[1] = Atree1pppm(p1, p3, p2, p4, p5, p6);
  }
  
  inline void ampq2g2l2::
  matrix_tree_ppmm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const  {
    M[0] = Atree1ppmm(p1, p2, p3, p4, p5, p6);
    M[1] = Atree1pmpm(p1, p3, p2, p4, p5, p6);
  }

  inline void ampq2g2l2::
  matrix_tree_pmpm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const  {
    M[0] = Atree1pmpm(p1, p2, p3, p4, p5, p6);
    M[1] = Atree1ppmm(p1, p3, p2, p4, p5, p6);
  }

  inline void ampq2g2l2::
  matrix_tree_pmmm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const  {
    M[0] = Atree1pmmm(p1, p2, p3, p4, p5, p6);
    M[1] = Atree1pmmm(p1, p3, p2, p4, p5, p6);
  }

  inline void ampq2g2l2::
  color_tree_pppm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const  {
    M[0] = Atree1pppm(p1, p2, p3, p4, p5, p6);
    M[1] = Atree1pppm(p1, p3, p2, p4, p5, p6);
  }
  
  inline void ampq2g2l2::
  color_tree_ppmm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const  {
    M[0] = Atree1ppmm(p1, p2, p3, p4, p5, p6);
    M[1] = Atree1pmpm(p1, p3, p2, p4, p5, p6);
  }

  inline void ampq2g2l2::
  color_tree_pmpm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const  {
    M[0] = Atree1pmpm(p1, p2, p3, p4, p5, p6);
    M[1] = Atree1ppmm(p1, p3, p2, p4, p5, p6);
  }

  inline void ampq2g2l2::
  color_tree_pmmm(int p1, int p2, int p3, int p4, int p5, int p6, _ComplexD *M) const  {
    M[0] = Atree1pmmm(p1, p2, p3, p4, p5, p6);
    M[1] = Atree1pmmm(p1, p3, p2, p4, p5, p6);
  }

  inline std::complex<double> 
  ampq2g2l2::As1pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = As1pppm(p4, p3, p2, p1, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Atree1pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Atree1pppm(p4, p3, p2, p1, p6, p5); swap();
    return res;
  }
  
  inline std::complex<double> 
  ampq2g2l2::Fcc1ppmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fcc = FAcc1ppmm(p1, p2, p3, p4, p5, p6);
    swap(); Fcc += FAcc1ppmm(p4, p3, p2, p1, p6, p5); swap();
    return Fcc;
  }
  
  inline std::complex<double> 
  ampq2g2l2::Fsc1ppmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fsc = FAsc1ppmm(p1, p2, p3, p4, p5, p6);
    swap(); Fsc += FAsc1ppmm(p4, p3, p2, p1, p6, p5); swap();
    return Fsc;
  }
  
  inline std::complex<double> 
  ampq2g2l2::Fcc1pmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fcc = FAcc1pmpm(p1, p2, p3, p4, p5, p6);
    swap(); Fcc += FAcc1pmpm(p4, p3, p2, p1, p6, p5); swap();
    return Fcc;
  }

  inline std::complex<double> 
  ampq2g2l2::Fsc1pmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fsc = FAsc1pmpm(p1, p2, p3, p4, p5, p6);
    swap(); Fsc += FAsc1pmpm(p4, p3, p2, p1, p6, p5); swap();
    return Fsc;
  }
  
  inline std::complex<double> 
  ampq2g2l2::Fcc3pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fcc = FBcc3pmmp(p1, p2, p3, p4, p5, p6);    
    swap(); Fcc += FBcc3pmmp(p2, p1, p4, p3, p6, p5); swap();
    return Fcc;
  }
  
  inline std::complex<double> 
  ampq2g2l2::Fsc3pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fsc = FBsc3pmmp(p1, p2, p3, p4, p5, p6);
    swap(); Fsc += FBsc3pmmp(p2, p1, p4, p3, p6, p5); swap();
    return Fsc;
  }

  inline std::complex<double> 
  ampq2g2l2::Fcc3pmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fcc = FBcc3pmpm(p1, p2, p3, p4, p5, p6);
    swap(); Fcc += FBcc3pmpm(p2, p1, p4, p3, p6, p5); swap();
    return Fcc;
  }
  
  inline std::complex<double> 
  ampq2g2l2::Fsc3pmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fsc = FBsc3pmpm(p1, p2, p3, p4, p5, p6);
    swap(); Fsc += FBsc3pmpm(p2, p1, p4, p3, p6, p5); swap();
    return Fsc;
  }

  inline std::complex<double> 
  ampq2g2l2::Fcc1pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fcc1pppm(p4, p3, p2, p1, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fsc1pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fsc1pppm(p4, p3, p2, p1, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Atree2pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Atree2ppmm(p3, p2, p1, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fcc2pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fcc2ppmm(p3, p2, p1, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fsc2pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fsc2ppmm(p3, p2, p1, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Atree2pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const{
    swap(); _ComplexD res = Atree2ppmp(p3, p2, p1, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fcc2pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fcc2ppmp(p3, p2, p1, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fsc2pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fsc2ppmp(p3, p2, p1, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Atree3pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Atree3pmpp(p2, p1, p4, p3, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fcc3pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fcc3pmpp(p2, p1, p4, p3, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fsc3pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fsc3pmpp(p2, p1, p4, p3, p6, p5); swap();
    return res;
  }


  inline std::complex<double> 
  ampq2g2l2::Fax3pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = -Fax3pmpp(p2, p1, p4, p3, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Faxsl3pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = -Faxsl3pmpp(p2, p1, p3, p4, p6, p5); swap();
    return res;
  }
  
  inline std::complex<double> 
  ampq2g2l2::Faxsl3pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = -Faxsl3pmpm(p2, p1, p3, p4, p6, p5); swap();
    return res;
  }
  
  inline std::complex<double> 
  ampq2g2l2::Fvf3pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fvf3pmpp(p2, p1, p3, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fvf3pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fvf3pmpm(p2, p1, p3, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fvs3pmmm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fvs3pmpp(p2, p1, p3, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double> 
  ampq2g2l2::Fvs3pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    swap(); _ComplexD res = Fvs3pmpm(p2, p1, p3, p4, p6, p5); swap();
    return res;
  }

  inline std::complex<double>
  ampq2g2l2::Fvs3pmpp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    return F2vs3pmpp(p1, p2, p3, p4, p5, p6)
      +    F2vs3pmpp(p1, p2, p4, p3, p5, p6);
  }

  inline std::complex<double>
  ampq2g2l2::Fvs3pmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fvs = FBvs3pmpm(p1, p2, p3, p4, p5, p6);
    swap(); Fvs += FBvs3pmpm(p2, p1, p4, p3, p6, p5); swap();
    return Fvs;
  }

  inline std::complex<double>
  ampq2g2l2::Fvf3pmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fvf = FBvf3pmpm(p1, p2, p3, p4, p5, p6);
    swap(); Fvf += FBvf3pmpm(p2, p1, p4, p3, p6, p5); swap();
    return Fvf;
  }

  inline std::complex<double>
  ampq2g2l2::Fax3pmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fax = FBax3pmpm(p1, p2, p3, p4, p5, p6);
    swap(); Fax += FBax3pmpm(p2, p1, p4, p3, p6, p5); swap();
    return Fax;
  }

  inline std::complex<double>
  ampq2g2l2::Fax3pmmp(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fax = FBax3pmmp(p1, p2, p3, p4, p5, p6);
    swap(); Fax += FBax3pmmp(p2, p1, p4, p3, p6, p5); swap();
    return Fax;
  }

  inline std::complex<double> 
  ampq2g2l2::Faxsl3pmpm(int p1, int p2, int p3, int p4, int p5, int p6) const {
    _ComplexD Fax = FBaxsl3pmpm(p1, p2, p3, p4, p5, p6);
    swap(); Fax += FBaxsl3pmpm(p2, p1, p4, p3, p6, p5); swap();
    return Fax;
  }

  inline double ampq2g2l2::su3_amptree(_ComplexD *__m) {
    return Na*(Na*real(__m[0]*conj(__m[0]) + __m[1]*conj(__m[1])) 
	       - 2.0*real(__m[0]*conj(__m[1])))/Nc;
  }
}   //  namespace nlo

#endif

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
#ifndef __NLO_AMP_AMPBASE_H__
#define __NLO_AMP_AMPBASE_H__ 1


//   standard headers
#include <complex>

//   nlo headers
#include <bits/hep-lorentzvector.h>
#include <bits/nlo-innerprod.h>
#include <bits/hep-rng.h>



namespace nlo {


  class _Amp_base 
  {
    //    private types
    typedef std::complex<double> _ComplexD;

  protected:
    //   protected types
    typedef innerprod<lorentzvector<double> >   innerprod_type;
    typedef innerprod_type::lorentz_prod_type   lorentz_prod_type;
    typedef innerprod_type::spinor_prod_type    spinor_prod_type;

    
    //  constructors
    _Amp_base(const innerprod_type& __x, const random_generator& rng)
      : _M_s(&__x.s()), _M_a(&__x.a()), _M_b(&__x.b()), _M_rng(rng) {}
    
    //  operator :  <ij>  <--> [ij]
    void swap() const { 
      const spinor_prod_type *__tmp = _M_a; _M_a = _M_b; _M_b = __tmp;
    }
    
    //  help to calculate the loop contributions
    static double Xq(double, double);
    static double Xg(double, double, unsigned int);
  
    static _ComplexD Log(double);
    static _ComplexD Li2(double);
    static _ComplexD Li2(double, double);
    static _ComplexD Li2(double, double, double, double);
    
    //    BDK functions
    static _ComplexD L0  (double, double);
    static _ComplexD L1  (double, double);
    static _ComplexD L2  (double, double);
    static _ComplexD Ls_1(double, double, double, double);
    static _ComplexD Ls0 (double, double, double, double);
    static _ComplexD Ls1 (double, double, double, double);
    static _ComplexD Ls2 (double, double, double, double);
    static _ComplexD Ls3 (double, double, double, double);

    static _ComplexD I3_3m    (double, double, double);
    static _ComplexD tLs_1_2mh(double, double, double, double);
    static _ComplexD Ls_1_2mh (double, double, double, double);
    static _ComplexD Ls_1_2me (double, double, double, double);
    static _ComplexD f        (double, double, double, double);
    
    //   data members
    mutable const lorentz_prod_type *_M_s;          // Lorentz products
    mutable const spinor_prod_type  *_M_a, *_M_b;   // Spinor products
    
    //   random number generator to do Monte Carlo helicity sum
    const random_generator& _M_rng;

  private:
    //  private static functions
    static _ComplexD I3m3(double, double, double);  
  };

  //   Some global special function
  //  (I don't want to compute a special finction library. 
  //   That is another issue. The best way is to use the gsl library.)
  double __specfunc_li2(double);
  double __specfunc_cl2(double);
  
  //   struct for storing the result of the K+P term (one initial parton)
  struct su3_kp_i1 
  { 
    //  default constructor  
    su3_kp_i1() 
      : tree(0.0), loop(0.0), ga(0.0), pa(0.0) {}
    
    //  computed assignments
    su3_kp_i1& operator+=(const su3_kp_i1&);
    su3_kp_i1& operator-=(const su3_kp_i1&);
    su3_kp_i1& operator*=(const double&);
    su3_kp_i1& operator/=(const double&);

    double tree, loop, ga, pa;
  };

  inline  su3_kp_i1& su3_kp_i1::operator+=(const su3_kp_i1& x)
  {
    tree += x.tree; ga += x.ga;
    loop += x.loop; pa += x.pa;
    return *this;
  }

  inline  su3_kp_i1& su3_kp_i1::operator-=(const su3_kp_i1& x)
  {
    tree -= x.tree; ga -= x.ga;
    loop -= x.loop; pa -= x.pa;
    return *this;
  }

  inline  su3_kp_i1& su3_kp_i1::operator*=(const double& x)
  {
    tree *= x; ga *= x;
    loop *= x; pa *= x;
    return *this;
  }

  inline  su3_kp_i1& su3_kp_i1::operator/=(const double& x)
  {
    tree /= x; ga /= x; 
    loop /= x; pa /= x;
    return *this;
  }

  inline su3_kp_i1 operator-(const su3_kp_i1& x)
  { return su3_kp_i1(x) *= -1.0;}

  inline su3_kp_i1 operator+(const su3_kp_i1& x, const su3_kp_i1& y)
  { return su3_kp_i1(x) += y;}

  inline su3_kp_i1 operator-(const su3_kp_i1& x, const su3_kp_i1& y)
  { return su3_kp_i1(x) -= y;}

  inline su3_kp_i1 operator*(const su3_kp_i1& x, const double& y)
  { return su3_kp_i1(x) *= y;}

  inline su3_kp_i1 operator*(const double& y, const su3_kp_i1& x)
  { return su3_kp_i1(x) *= y;}
  
  inline su3_kp_i1 operator/(const su3_kp_i1& x, const double& y)
  { return su3_kp_i1(x) /= y;}
  
  
  //   struct for storing the result of the K+P term (two initial parton)
  struct su3_kp_i2 
  { 
    //  default constructor  
    su3_kp_i2() 
      : tree(0.0), loop(0.0), ga(0.0), gb(0.0), 
	cca(0.0), ccb(0.0), pa(0.0), pb(0.0) {}
    
    //  computed assignments
    su3_kp_i2& operator+=(const su3_kp_i2&);
    su3_kp_i2& operator-=(const su3_kp_i2&);
    su3_kp_i2& operator*=(const double&);
    su3_kp_i2& operator/=(const double&);

    double tree, loop, ga, gb, cca, ccb, pa, pb;
  };

  inline  su3_kp_i2& su3_kp_i2::operator+=(const su3_kp_i2& x)
  {
    tree += x.tree; ga += x.ga; cca += x.cca; pa += x.pa;
    loop += x.loop; gb += x.gb; ccb += x.ccb; pb += x.pb;
    return *this;
  }

  inline  su3_kp_i2& su3_kp_i2::operator-=(const su3_kp_i2& x)
  {
    tree -= x.tree; ga -= x.ga; cca -= x.cca; pa -= x.pa;
    loop -= x.loop; gb -= x.gb; ccb -= x.ccb; pb -= x.pb;
    return *this;
  }

  inline  su3_kp_i2& su3_kp_i2::operator*=(const double& x)
  {
    tree *= x; ga *= x; cca *= x; pa *= x;
    loop *= x; gb *= x; ccb *= x; pb *= x;
    return *this;
  }

  inline  su3_kp_i2& su3_kp_i2::operator/=(const double& x)
  {
    tree /= x; ga /= x; cca /= x; pa /= x;
    loop /= x; gb /= x; ccb /= x; pb /= x;
    return *this;
  }

  inline su3_kp_i2 operator-(const su3_kp_i2& x)
  { return su3_kp_i2(x) *= -1.0;}

  inline su3_kp_i2 operator+(const su3_kp_i2& x, const su3_kp_i2& y)
  { return su3_kp_i2(x) += y;}

  inline su3_kp_i2 operator-(const su3_kp_i2& x, const su3_kp_i2& y)
  { return su3_kp_i2(x) -= y;}

  inline su3_kp_i2 operator*(const su3_kp_i2& x, const double& y)
  { return su3_kp_i2(x) *= y;}

  inline su3_kp_i2 operator*(const double& y, const su3_kp_i2& x)
  { return su3_kp_i2(x) *= y;}
  
  inline su3_kp_i2 operator/(const su3_kp_i2& x, const double& y)
  { return su3_kp_i2(x) /= y;} 
}   // namespace nlo

#endif





  
  

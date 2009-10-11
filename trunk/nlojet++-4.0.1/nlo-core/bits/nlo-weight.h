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
#ifndef __NLO_NLO_WEIGHT_H__
#define __NLO_NLO_WEIGHT_H__ 1

//   standard headers
#include <cmath>
#include <cstring>
#include <iostream>

//   nlojet++ includes
#include <bits/hep-sample_traits.h>



namespace nlo {

  extern const char *weight_no_label[0];

  template<unsigned int _Size, const char **_Label = weight_no_label>
  class weight
  {
  public:
    //   constructors
    weight() {
      std::memset(_M_sub, 0, _Size*sizeof(double));
    }
    
    explicit weight(const double *arr) {
      std::memcpy(_M_sub, arr, _Size*sizeof(double));
    }

    //   copy constructor
    weight(const weight<_Size, _Label>& w) {
      std::memcpy(_M_sub, w._M_sub, _Size*sizeof(double));
    }
    
    //   assignment
    weight<_Size, _Label>& operator=(const weight<_Size, _Label>& w) {
      if(this != &w) std::memcpy(_M_sub, w._M_sub, _Size*sizeof(double));
      return *this;
    }
    
    //   computed assignments
    weight& operator*=(const double& x) { 
      for(double *i = _M_sub; i < _M_sub+_Size; *(i++) *= x);
       return *this;
    }
    
    weight& operator/=(const double& x) { 
      for(double *i = _M_sub; i < _M_sub+_Size; *(i++) /= x);
       return *this;
    }
    
    weight& operator*=(const weight& x) { 
      for(unsigned i = 0; i < _Size; i++) 
	_M_sub[i] *= x[i];
      return *this;
    }
 
   weight& operator/=(const weight& x) { 
      for(unsigned i = 0; i < _Size; i++) 
	_M_sub[i] /= x[i];
      return *this;
    }
 
    weight& operator+=(const weight& x) { 
      for(unsigned i = 0; i < _Size; i++) 
	_M_sub[i] += x[i];
      return *this;
    }

    weight& operator-=(const weight& x) { 
      for(unsigned i = 0; i < _Size; i++) 
	_M_sub[i] -= x[i];
      return *this;
    }

    //   iterators
    double * begin() { return _M_sub;} 
    double * end() { return _M_sub + _Size;} 

    const double * begin() const { return _M_sub;} 
    const double * end() const { return _M_sub + _Size;} 

    //   element access
    const double& operator[](unsigned int n) const { return _M_sub[n];}
    double& operator[](unsigned int n) { return _M_sub[n];}
    
    //  size of the weight
    unsigned size() const { return _Size;}
    
    //  get the labels
    const char *label(unsigned int i) const { 
      if(_Label == weight_no_label) return "no label defined";
      else return _Label[i];
    }
    
  private:
    double _M_sub[_Size];
  };

  
  template<unsigned int _Size, const char **_Label>
  inline weight<_Size, _Label> operator-(const weight<_Size, _Label>& a) {
    return weight<_Size, _Label>(a) *= -1.0;
  }
  
  template<unsigned int _Size, const char **_Label> inline 
  weight<_Size, _Label> operator+(const weight<_Size, _Label>& a, const weight<_Size, _Label>& b) {
    return weight<_Size, _Label>(a) += b;
  }

  template<unsigned int _Size, const char **_Label> inline 
  weight<_Size, _Label> operator-(const weight<_Size, _Label>& a, const weight<_Size, _Label>& b) {
    return weight<_Size, _Label>(a) -= b;
  }

  template<unsigned int _Size, const char **_Label> inline 
  weight<_Size, _Label> operator/(const weight<_Size, _Label>& a, const weight<_Size, _Label>& b) {
    return weight<_Size, _Label>(a) /= b;
  }

  template<unsigned int _Size, const char **_Label> inline 
  weight<_Size, _Label> operator*(const weight<_Size, _Label>& a, const weight<_Size, _Label>& b) {
    return weight<_Size, _Label>(a) *= b;
  }
  
  template<unsigned int _Size, const char **_Label> inline 
  weight<_Size, _Label> operator*(const weight<_Size, _Label>& a, double b) {
    return weight<_Size, _Label>(a) *= b;
  }
  
  template<unsigned int _Size, const char **_Label> inline 
  weight<_Size, _Label> operator*(double b, const weight<_Size, _Label>& a) {
    return weight<_Size, _Label>(a) *= b;
  }

  template<unsigned int _Size, const char **_Label> inline 
  weight<_Size, _Label> operator/(const weight<_Size, _Label>& a, double b) {
    return weight<_Size, _Label>(a) /= b;
  }
  
  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> operator/(double a, const weight<_Size, _Label>& b) 
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = a/b[i];
    return res;
  }

  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> apply(const weight<_Size, _Label>& x, double f(double)) 
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = f(x[i]);
    return res;
  }
  
  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> apply(const weight<_Size, _Label>& x, double f(const double&)) 
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = f(x[i]);
    return res;
  }
  
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> abs(const weight<_Size, _Label>& x) {
    return apply(x, std::abs);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> acos(const weight<_Size, _Label>& x) {
    return apply(x, std::acos);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> asin(const weight<_Size, _Label>& x) {
    return apply(x, std::asin);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> atan(const weight<_Size, _Label>& x) {
    return apply(x, std::atan);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> cos(const weight<_Size, _Label>& x) {
    return apply(x, std::cos);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> cosh(const weight<_Size, _Label>& x) {
    return apply(x, std::cosh);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> exp(const weight<_Size, _Label>& x) {
    return apply(x, std::exp);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> log(const weight<_Size, _Label>& x) {
    return apply(x, std::log);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> log10(const weight<_Size, _Label>& x) {
    return apply(x, std::log10);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> sin(const weight<_Size, _Label>& x) {
    return apply(x, std::sin);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> sinh(const weight<_Size, _Label>& x) {
    return apply(x, std::sinh);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> sqrt(const weight<_Size, _Label>& x) {
    return apply(x, std::sqrt);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> tan(const weight<_Size, _Label>& x) {
    return apply(x, std::tan);
  }
  
  template<unsigned int _Size, const char **_Label> 
  inline weight<_Size, _Label> tanh(const weight<_Size, _Label>& x) {
    return apply(x, std::tanh);
  }
  
  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> atan2(const weight<_Size, _Label>& x, const weight<_Size, _Label>& y) 
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::atan2(x[i], y[i]);
    return res;
  }

  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> atan2(const weight<_Size, _Label>& x, const double& y)
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::atan2(x[i], y);
    return res;
  }
  
  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> atan2(const double& x, const weight<_Size, _Label>& y) 
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::atan2(x, y[i]);
    return res;
  }
  
  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> pow(const weight<_Size, _Label>& x, const weight<_Size, _Label>& y)
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::pow(x[i], y[i]);
    return res;
  }
   
  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> pow(const weight<_Size, _Label>& x, const double& y) 
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::pow(x[i], y);
    return res;
  }
  
  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> pow(const double& x, const weight<_Size, _Label>& y)
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::pow(x, y[i]);
    return res;
  }
  
  template<unsigned int _Size, const char **_Label> 
  weight<_Size, _Label> pow(const weight<_Size, _Label>& x, int y) 
  { 
    weight<_Size, _Label> res;
    for(unsigned i = 0; i < _Size; i++) 
      res[i] = std::pow(x[i], y);
    return res;
  }

  template<unsigned int _Size, const char **_Label>
  inline bool IsFinite(const weight<_Size, _Label> & x) 
  {
    for(unsigned int i = 0; i < _Size; i++)
      if(!::finite(x[i])) return false;
    return true;
  }
  
  template<unsigned int _Size, const char **_Label>
  std::ostream& operator<<(std::ostream& os, const weight<_Size, _Label>& x) 
  {  
    os<<"[";
    for(unsigned int i = 0; i < _Size-1; i++) os<<x[i]<<", ";
    os<<x[_Size-1]<<"]";
    return os;
  }
  
  template<unsigned int _Size, const char **_Label>
  std::istream& operator>>(std::istream& is, weight<_Size, _Label>& x) 
  { 
    char dummy;
    is>>dummy;
    for(unsigned int i = 0; i < _Size; i++) is>>x[i]>>dummy;
    return is;
  }

  template<unsigned int _Size, const char **_Label>
  std::ostream& write(std::ostream& os, const weight<_Size, _Label>& x) { 
    os.write((const char *) (x.begin()), _Size*sizeof(double));
    return os;
  }
  
  template<unsigned int _Size, const char **_Label>
  std::istream& read(std::istream& is, weight<_Size, _Label>& x) { 
    is.read((char *) (x.begin()), _Size*sizeof(double));
    return is;
  }
}

namespace nlo {

  //   empty declaration of the conversion function of class weight
  template<class _From, class _To>
  struct weight_conversion;
  
  //   conversion : weight --> double (effect: sum of the subprocesses)
  template<unsigned int _Size, const char **_Label>
  struct weight_conversion<weight<_Size, _Label>, double> 
    : public std::unary_function<weight<_Size, _Label>, double>
  {
    double operator()(const weight<_Size, _Label>& x) { 
      double ret_val = 0.0;
      for(const double *i = x.begin(); i < x.end(); ret_val += *(i++));
      return ret_val;
    }
  };  
}


//
//      class weight<_Size> as sample type
//
namespace nlo {
  
  template<unsigned int _Size, const char **_Label>
  struct sample_traits<weight<_Size, _Label> >
  {
    //   Sample type
    typedef weight<_Size, _Label> sample_type;
    
    //   Neutral element (a + zero = a)
    static const sample_type zero;

    //   Assignable  (a = b)
    static void assign(sample_type& a, const sample_type& b) { a = b;}

    //   Computed assignments (a+=b, a-=b, a*=b)
    static void assadd(sample_type& a, const sample_type& b) { a += b;}
    static void asssub(sample_type& a, const sample_type& b) { a -= b;}
    static void assmul(sample_type& a, const sample_type& b) { a *= b;}
    static void assdiv(sample_type& a, const sample_type& b) { a /= b;}

    //   Computed assignment (a*=d, where d is double)
    static void assmul(sample_type& a, const double& d) { a *= d;};
    
    //   Other multiplicative assignments
    template<typename _Xp>
    static void assmul(sample_type& a, const _Xp& x) { a *= x;}
    
    //   Squareroot function
    static void fsqrt(sample_type& x) { x = sqrt(x);}
    
    //   True if the sample is finite (not Inf, not NaN)
    static bool finite(const sample_type& x) { return IsFinite(x);}
        
    //   I/O functions
    static void write_txt(std::ostream& os, const sample_type& x) { os<<x;}
    static void write_bin(std::ostream& os, const sample_type& x) { write(os, x);}
    
    static void read_txt(std::istream& is, sample_type& x) { is>>x;} 
    static void read_bin(std::istream& is, sample_type& x) { read(is, x);} 

    static void print(std::ostream&, const sample_type&, const sample_type&);  
  };

  template<unsigned int _Size, const char **_Label>
  const weight<_Size, _Label> sample_traits<weight<_Size, _Label> >::
  zero = weight<_Size, _Label>();



  template<unsigned int _Size, const char **_Label>
  void sample_traits<weight<_Size, _Label> >::
  print(std::ostream& os, const sample_type& mean, const sample_type& err) 
  {    
    double m, v;

    for(unsigned int l = 0; l < _Size; l++) {
      m = mean[l]; v = err[l];
      os<<"   "; if(m >= 0.0) os<<" "; os<<m;
      os<<"  ";  if(v >= 0.0) os<<" "; os<<v;
    }
  }





}   //  namespace nlo


#endif


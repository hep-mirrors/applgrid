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
#ifndef __NLO_HEP_SAMPLE_TRAITS_H__
#define __NLO_HEP_SAMPLE_TRAITS_H__ 1


//   Standard includes
#include <cmath>
#include <iostream>


namespace nlo {


  //   sample_traits template class
  template<typename _Tp>
  struct sample_traits 
  {
    //   Sample type
    typedef _Tp sample_type;
    
    //   Neutral element (a + zero = a)
    static const sample_type zero;

    //   Assignable  (a = b)
    static void assign(sample_type&, const sample_type&);

    //   Computed assignments (a+=b, a-=b, a*=b, a /= b)
    static void assadd(sample_type&, const sample_type&);
    static void asssub(sample_type&, const sample_type&);
    static void assmul(sample_type&, const sample_type&);
    static void assdiv(sample_type&, const sample_type&);

    //   Computed assignment (a*=d, where d is double)
    static void assmul(sample_type&, const double&);
    
    //   Other multiplicative assignments
    template<typename _Xp>
    static void assmul(sample_type&, const _Xp&);
    
    //   Squareroot function
    static void fsqrt(sample_type&);

    //   True if the sample is finite (not Inf, not NaN)
    bool finite(const sample_type&);
        
    //   I/O functions
    static void write_txt(std::ostream&, const sample_type&);
    static void write_bin(std::ostream&, const sample_type&);

    static void read_txt(std::istream&, sample_type&); 
    static void read_bin(std::istream&, sample_type&); 
    
    static void print(std::ostream&, const sample_type&, const sample_type&);
  };

  //
  //   Specialization : sample_traits<double>
  //
  template<>
  struct sample_traits<double>
  {
    //   Sample type
    typedef double sample_type;
    
    //   Neutral element (a + zero = a)
    static const sample_type zero;

    //   Assignable  (a = b)
    static void assign(sample_type&, const sample_type&);

    //   Computed assignments (a+=b, a-=b, a*=b, a/=b)
    static void assadd(sample_type&, const sample_type&);
    static void asssub(sample_type&, const sample_type&);
    static void assmul(sample_type&, const sample_type&);
    static void assdiv(sample_type&, const sample_type&);

    //   Other multiplicative assignments
    template<typename _Xp>
    static void assmul(sample_type&, const _Xp&);
    
    //   Squareroot function
    static void fsqrt(sample_type&);
        
    //   True if the sample is finite (not Inf, not NaN)
    static bool finite(const sample_type&);

    //   I/O functions
    static void write_txt(std::ostream&, const sample_type&);
    static void write_bin(std::ostream&, const sample_type&);

    static void read_txt(std::istream&, sample_type&); 
    static void read_bin(std::istream&, sample_type&); 

    static void print(std::ostream&, const sample_type&, const sample_type&);
  };
  
  inline void sample_traits<double>::
  assign(double& ths, const double& x) {
    ths = x;
  }
  
  inline void sample_traits<double>::
  assadd(double& ths, const double& x) { 
    ths += x;
  }

  inline void sample_traits<double>::
  asssub(double& ths, const double& x) { 
    ths -= x;
  }
 
  inline void sample_traits<double>::
  assmul(double& ths, const double& x) { 
    ths *= x;
  }
 
  inline void sample_traits<double>::
  assdiv(double& ths, const double& x) { 
    ths /= x;
  }
  
  inline void sample_traits<double>::
  fsqrt(double& x) { 
    x = std::sqrt(x);
  }
   
  inline bool sample_traits<double>::finite(const double& x) {
    return (::finite(x) ? true : false);
  }

  template<typename _Xp>
  inline void sample_traits<double>::
  assmul(double& ths, const _Xp& x) {
    ths *= x;
  }
  
  inline void sample_traits<double>::
  write_txt(std::ostream& os, const double& x) { 
    os<<x;
  }
  
  inline void sample_traits<double>::
  write_bin(std::ostream& os, const double& x) {
    os.write((const char *) &x, sizeof(double));
  }

  inline void sample_traits<double>::
  read_txt(std::istream& is, double& x) { 
    is>>x;
  }
  
  inline void sample_traits<double>::
  read_bin(std::istream& is, double& x) {
    is.read((char *) &x, sizeof(double));
  }

  inline void sample_traits<double>::
  print(std::ostream& os, const double& mean, const double& err) {
    if(mean >= 0.0) os<<" "; os<<mean<<"  ";    
    if(err >= 0.0)  os<<" "; os<<err;
  }

  //
  //   Specialization : sample_traits<float>
  //
  template<>
  struct sample_traits<float>
  {
    //   Sample type
    typedef float sample_type;
    
    //   Neutral element (a + zero = a)
    static const sample_type zero;

    //   Assignable  (a = b)
    static void assign(sample_type&, const sample_type&);

    //   Computed assignments (a+=b, a-=b, a*=b, a/=b)
    static void assadd(sample_type&, const sample_type&);
    static void asssub(sample_type&, const sample_type&);
    static void assmul(sample_type&, const sample_type&);
    static void assdiv(sample_type&, const sample_type&);

    //   Computed assignment (a*=d, where d is double)
    static void assmul(sample_type&, const double&);

    //   Other multiplicative assignments
    template<typename _Xp>
    static void assmul(sample_type&, const _Xp&);
    
    //   Squareroot function
    static void fsqrt(sample_type&);
        
    //   True if the sample is finite (not Inf, not NaN)
    static bool finite(const sample_type&);

    //   I/O functions
    static void write_txt(std::ostream&, const sample_type&);
    static void write_bin(std::ostream&, const sample_type&);

    static void read_txt(std::istream&, sample_type&); 
    static void read_bin(std::istream&, sample_type&); 

    static void print(std::ostream&, const sample_type&, const sample_type&);
  };

  inline void sample_traits<float>::
  assign(float& ths, const float& x) {
    ths = x;
  }
  
  inline void sample_traits<float>::
  assadd(float& ths, const float& x) { 
    ths += x;
  }

  inline void sample_traits<float>::
  asssub(float& ths, const float& x) { 
    ths -= x;
  }
 
  inline void sample_traits<float>::
  assmul(float& ths, const float& x) { 
    ths *= x;
  }
  
  inline void sample_traits<float>::
  assdiv(float& ths, const float& x) { 
    ths /= x;
  }
  
  inline void sample_traits<float>::
  assmul(float& ths, const double& x) { 
    ths *= x;
  }
  
  inline void sample_traits<float>::
  fsqrt(float& x) { 
    x = std::sqrt(x);
  }
   
  inline bool sample_traits<float>::finite(const float& x) {
    return (::finite(x) ? true : false);
  }

  template<typename _Xp>
  inline void sample_traits<float>::
  assmul(float& ths, const _Xp& x) {
    ths *= x;
  }
  
  inline void sample_traits<float>::
  write_txt(std::ostream& os, const float& x) { 
    os<<x;
  }
  
  inline void sample_traits<float>::
  write_bin(std::ostream& os, const float& x) {
    os.write((const char *) &x, sizeof(float));
  }

  inline void sample_traits<float>::
  read_txt(std::istream& is, float& x) { 
    is>>x;
  }
  
  inline void sample_traits<float>::
  read_bin(std::istream& is, float& x) {
    is.read((char *) &x, sizeof(float));
  }

  inline void sample_traits<float>::
  print(std::ostream& os, const float& mean, const float& err) {
    if(mean >= 0.0f) os<<" "; os<<mean<<"  ";    
    if(err >= 0.0f)  os<<" "; os<<err;
  }

  //
  //   Specialization : sample_traits<long double>
  //
  template<>
  struct sample_traits<long double>
  {
    //   Sample type
    typedef long double sample_type;
    
    //   Neutral element (a + zero = a)
    static const sample_type zero;

    //   Assignable  (a = b)
    static void assign(sample_type&, const sample_type&);

    //   Computed assignments (a+=b, a-=b, a*=b, a/=b)
    static void assadd(sample_type&, const sample_type&);
    static void asssub(sample_type&, const sample_type&);
    static void assmul(sample_type&, const sample_type&);
    static void assdiv(sample_type&, const sample_type&);

    //   Computed assignment (a*=d, where d is double)
    static void assmul(sample_type&, const double&);

    //   Other multiplicative assignments
    template<typename _Xp>
    static void assmul(sample_type&, const _Xp&);
    
    //   Squareroot function
    static void fsqrt(sample_type&);
        
    //   True if the sample is finite (not Inf, not NaN)
    static bool finite(const sample_type&);

    //   I/O functions
    static void write_txt(std::ostream&, const sample_type&);
    static void write_bin(std::ostream&, const sample_type&);

    static void read_txt(std::istream&, sample_type&); 
    static void read_bin(std::istream&, sample_type&); 

    static void print(std::ostream&, const sample_type&, const sample_type&);
  };

  
  inline void sample_traits<long double>::
  assign(long double& ths, const long double& x) {
    ths = x;
  }
  
  inline void sample_traits<long double>::
  assadd(long double& ths, const long double& x) { 
    ths += x;
  }

  inline void sample_traits<long double>::
  asssub(long double& ths, const long double& x) { 
    ths -= x;
  }
 
  inline void sample_traits<long double>::
  assmul(long double& ths, const long double& x) { 
    ths *= x;
  }
  
  inline void sample_traits<long double>::
  assdiv(long double& ths, const long double& x) { 
    ths /= x;
  }
  
  inline void sample_traits<long double>::
  assmul(long double& ths, const double& x) { 
    ths *= x;
  }
  
  inline void sample_traits<long double>::
  fsqrt(long double& x) { 
    x = std::sqrt(x);
  }
   
  inline bool sample_traits<long double>::finite(const long double& x) {
    return (::finite(x) ? true : false);
  }

  template<typename _Xp>
  inline void sample_traits<long double>::
  assmul(long double& ths, const _Xp& x) {
    ths *= x;
  }
  
  inline void sample_traits<long double>::
  write_txt(std::ostream& os, const long double& x) { 
    os<<x;
  }
  
  inline void sample_traits<long double>::
  write_bin(std::ostream& os, const long double& x) {
    os.write((const char *) &x, sizeof(long double));
  }

  inline void sample_traits<long double>::
  read_txt(std::istream& is, long double& x) { 
    is>>x;
  }
  
  inline void sample_traits<long double>::
  read_bin(std::istream& is, long double& x) {
    is.read((char *) &x, sizeof(long double));
  }

  inline void sample_traits<long double>::
  print(std::ostream& os, const long double& mean, const long double& err) {
    if(mean >= 0.0l) os<<" "; os<<mean<<"  ";    
    if(err >= 0.0l)  os<<" "; os<<err;
  }

}   //  namespace nlo

#endif

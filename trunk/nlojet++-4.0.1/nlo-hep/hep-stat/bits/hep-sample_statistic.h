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
#ifndef __NLO_HEP_SAMPLE_STATISTIC_H__
#define __NLO_HEP_SAMPLE_STATISTIC_H__ 1


//   nlo includes
#include <bits/hep-sample_traits.h>
#include <bits/hep-exception.h>


namespace nlo {
 
  //     forward declarations
  template<typename _Tp, class _Traits = sample_traits<_Tp> >
  class sample_statistic;
  
  template<typename _Tp, class _Traits> std::ostream& 
  operator<<(std::ostream&, const sample_statistic<_Tp, _Traits>&); 
  
  template<typename _Tp, class _Traits> 
  std::ostream& write(std::ostream&, const sample_statistic<_Tp, _Traits>&); 
  
  template<typename _Tp, class _Traits> 
  std::istream& operator>>(std::istream&, sample_statistic<_Tp, _Traits>&);
  
  template<typename _Tp, class _Traits> 
  std::istream& read(std::istream&, sample_statistic<_Tp, _Traits>&); 
  
  template<typename _Tp, class _Traits> 
  std::ostream& print(std::ostream&, const sample_statistic<_Tp, _Traits>&);
  
  
  //    template class sample_statistic
  template<typename _Tp, class _Traits>
  class sample_statistic 
  {
  public:
    //    types
    typedef _Traits  traits_type;
    typedef typename _Traits::sample_type sample_type;
    
    //   default constructor
    sample_statistic() 
      : _M_n(0UL), _M_dat(_Traits::zero), _M_dat2(_Traits::zero)
    {}
    
    //   constructor
    sample_statistic(unsigned long int n, const sample_type& dat, const sample_type& dat2)
      : _M_n(n), _M_dat(dat), _M_dat2(dat2)
    {} 
   
    //   add the samples
    sample_statistic& operator+=(const sample_type&); 
    sample_statistic& operator+=(const sample_statistic& x) 
    {
      _Traits::assadd(_M_dat,  x._M_dat);
      _Traits::assadd(_M_dat2, x._M_dat2);
      _M_n += x._M_n;
      return *this;
    } 
    
    sample_statistic& assadd_weighted(const sample_statistic&); 
    sample_statistic& add_samples_no_check(const sample_type&); 
      
    //    add zero sample
    unsigned long int operator++() { return ++_M_n;}
    unsigned long int operator++(int) { return _M_n++;}
    
    //   re-initialize
    void reset() 
    { 
      _M_n = 0UL; 
      _Traits::assign(_M_dat,  _Traits::zero);
      _Traits::assign(_M_dat2, _Traits::zero);
    }
    
    //   number of samples
    unsigned long int samples() const { return _M_n;}
    
    //   sum of the samples
    const sample_type& sum_samples() const { return _M_dat;}
    
    //   sum of the samples square
    const sample_type& sum_samples_sqr() const { return _M_dat2;}
    
    //   the mean of the samples
    sample_type mean() const 
    { 
      double on = 1.0/((double) _M_n);
      sample_type dat(_M_dat);
      _Traits::assmul(dat, on);
      return dat;
    }

    //   the mean of the squared samples
    sample_type mean2() const 
    { 
      double on = 1.0/((double) _M_n);
      sample_type dat2(_M_dat2);
      _Traits::assmul(dat2, on);
      return dat2;
    }
    
    //   the variance of the samples
    sample_type variance() const;
    
    //   estimated error
    sample_type error() const;
    
  protected:
    //     data members
    unsigned long int _M_n;
    sample_type _M_dat, _M_dat2;

    //    friend declarations
    friend std::ostream& operator<<<>(std::ostream&, const sample_statistic&); 
    friend std::ostream& write<>(std::ostream&, const sample_statistic&); 
    friend std::istream& operator>><>(std::istream&, sample_statistic&);
    friend std::istream& read<>(std::istream&, sample_statistic&); 
  };
  
  template<typename _Tp, class _Traits>
  sample_statistic<_Tp, _Traits>& 
  sample_statistic<_Tp, _Traits>::operator+=(const sample_type& x) 
  {
    sample_type xx(x), dat(_M_dat), dat2(_M_dat2);
    _Traits::assmul(xx, xx);
    
    _Traits::assadd(dat,  x); 
    _Traits::assadd(dat2, xx);
    
    if(_Traits::finite(dat)  == false) return *this;
    if(_Traits::finite(dat2) == false) return *this;
    
    _Traits::assign(_M_dat,  dat); 
    _Traits::assign(_M_dat2, dat2);
    
    ++_M_n;
    return *this;
  }

  template<typename _Tp, class _Traits>
  sample_statistic<_Tp, _Traits>& 
  sample_statistic<_Tp, _Traits>::
  add_samples_no_check(const sample_type& x) 
  {
    sample_type xx(x);
    _Traits::assmul(xx, xx);   
    _Traits::assadd(_M_dat,  x); 
    _Traits::assadd(_M_dat2, xx);
    
    ++_M_n;
    return *this;
  }

  template<typename _Tp, class _Traits>
  typename sample_statistic<_Tp, _Traits>::sample_type
  sample_statistic<_Tp, _Traits>::variance() const 
  { 
    sample_type m(this -> mean()), m2(this -> mean2());
    
    traits_type::assmul(m,  m);
    traits_type::asssub(m2, m);
    
    return m2;
  }
  
  template<typename _Tp, class _Traits>
  typename sample_statistic<_Tp, _Traits>::sample_type
  sample_statistic<_Tp, _Traits>::error() const 
  { 
    double on = 1.0/((double) _M_n);
    sample_type v(this -> variance());
    
    traits_type::assmul(v, on);
    traits_type::fsqrt(v);
    
    return v;
  }

  template<typename _Tp, class _Traits>
  inline sample_statistic<_Tp, _Traits> 
  operator+(const sample_statistic<_Tp, _Traits> a,
	    const sample_statistic<_Tp, _Traits> b) 
  { return sample_statistic<_Tp, _Traits>(a) += b;}
 

  template<typename _Tp, class _Traits> sample_statistic<_Tp, _Traits>& 
  sample_statistic<_Tp, _Traits>::assadd_weighted(const sample_statistic& x) 
  {
    sample_type m1(this->mean()), v1(this->variance());
    sample_type m2(x.mean()), v2(x.variance());
    sample_type v;
    
    try {
      __fe_clear_exception();
      //  m1 = m1/v1 + m2/v2
      _Traits::assdiv(m1, v1);
      _Traits::assdiv(m2, v2);
      _Traits::assadd(m1, m2);
      
      //  Var(f) = v = v1*v2/(v1+v2) 
      _Traits::assign(v,  v1);
      _Traits::assmul(v,  v2);
      _Traits::assadd(v1, v2);
      _Traits::assdiv(v,  v1);
      __fe_throw_exception();
    } catch(fp_exception) {
      this -> operator+=(x);
      return *this;
    }
    
    //  <f> = m1 = m1*v
    _Traits::assmul(m1, v);

    //  <f^2> = v = v + m1*m1
    _Traits::assign(m2, m1);
    _Traits::assmul(m2, m1);
    _Traits::assadd(v,  m2);
    
    
    _M_n = (unsigned long int) ((2.0*x._M_n*_M_n)/((double) (_M_n + x._M_n)));
    _Traits::assmul(m1, (double) _M_n);    
    _Traits::assmul(v,  (double) _M_n);    
    
    _Traits::assign(_M_dat, m1);
    _Traits::assign(_M_dat2, v);
    
    return *this;
  } 
 

  template<typename _Tp, class _Traits>
  inline sample_statistic<_Tp, _Traits> 
  add_weighted(const sample_statistic<_Tp, _Traits> a,
	       const sample_statistic<_Tp, _Traits> b) 
  { return sample_statistic<_Tp, _Traits>(a).assadd_weighted(b);}
  

  //        I/O operations
  template<typename _Tp, class _Traits> std::ostream& 
  operator<<(std::ostream& os, const sample_statistic<_Tp, _Traits>& x) 
  {
    os<<x._M_n<<"  ";
    _Traits::write_txt(os, x._M_dat);  os<<"  ";
    _Traits::write_txt(os, x._M_dat2); os<<"\n";
    return os;
  }
  
  template<typename _Tp, class _Traits> std::ostream& 
  write(std::ostream& os, const sample_statistic<_Tp, _Traits>& x) 
  {
    os.write((const char *) &(x._M_n), sizeof(unsigned long int));
    _Traits::write_bin(os, x._M_dat);
    _Traits::write_bin(os, x._M_dat2);
    return os;
  }

  template<typename _Tp, class _Traits> std::istream& 
  operator>>(std::istream& is, sample_statistic<_Tp, _Traits>& x) 
  {
    is>>x._M_n;
    _Traits::read_txt(is, x._M_dat);
    _Traits::read_txt(is, x._M_dat2);
    return is;
  }
  
  template<typename _Tp, class _Traits> std::istream& 
  read(std::istream& is, sample_statistic<_Tp, _Traits>& x) 
  {
    is.read((char *) &(x._M_n), sizeof(unsigned long int));
    _Traits::read_bin(is, x._M_dat);
    _Traits::read_bin(is, x._M_dat2);
    return is;
  }

  template<typename _Tp, class _Traits> std::ostream& 
  print(std::ostream& os, const sample_statistic<_Tp, _Traits>& x) 
  {
    _Traits::print(os, x.mean(), x.error());
    return os;
  }
}   //  namespace nlo

#endif

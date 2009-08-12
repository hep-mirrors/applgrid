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
#ifndef __NLO_HEP_DIST_H__
#define __NLO_HEP_DIST_H__ 1


//   Standard includes
#include <vector>
#include <iostream>

//   nlo includes
#include <bits/hep-func.h>
#include <bits/hep-sample_statistic.h>
#include <bits/hep-distpoint_traits.h>


namespace nlo {
  
  //   Forward declaration of the class distriburion
  template<typename _Tp, class _Point, 
	   class _Traits = sample_traits<_Tp>, 
	   class _PointTraits = distpoint_traits<_Point> >
  class distribution;

  //   I/O operations
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::ostream& operator<<(std::ostream&, const distribution<_Tp, _Point, _Traits, _PointTraits>&);   
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::ostream& print(std::ostream&, const distribution<_Tp, _Point, _Traits, _PointTraits>&);
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::ostream& write(std::ostream&, const distribution<_Tp, _Point, _Traits, _PointTraits>&);
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::istream& operator>>(std::istream&, distribution<_Tp, _Point, _Traits, _PointTraits>&);
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::istream& read(std::istream&, distribution<_Tp, _Point, _Traits, _PointTraits>&);

  //
  //   Specialization : class _Point = void
  //
  template<typename _Tp, class _Traits>
  class distribution<_Tp, void, _Traits, distpoint_traits<void> > 
    : public sample_statistic<_Tp, _Traits>
  {
    //    private types 
    typedef sample_statistic<_Tp, _Traits> _Base;
    
  public:
    //       types
    typedef _Base statistic_type;
    typedef typename _Base::traits_type sample_traits;
    typedef typename _Base::sample_type sample_type;
    
    //   constructors
    distribution()
      : _M_tmp(_Traits::zero) {}
    
    //   assignments
    distribution& operator=(const distribution& x) {
      if(this != &x) {
				_Base::operator=(x);
				_Traits::assign(_M_tmp, x._M_tmp);
      }
      return *this;
    }
    
    distribution& operator+=(const distribution& x) {
      _Base::operator+=(x);
      return *this;
    }

    distribution& assadd_weighted(const distribution& x) {
      _Base::assadd_weighted(x);
      return *this;
    }
    
    //  add the samples 
    unsigned long int operator++() {
      _Base::operator+=(_M_tmp);
      _Traits::assign(_M_tmp, _Traits::zero);
      return _Base::samples();
    }
    
    unsigned long int operator++(int) 
    { return operator++() - 1UL;}
    
    unsigned long int add_samples_no_check() {
      _Base::add_samples_no_check(_M_tmp);
      _Traits::assign(_M_tmp, _Traits::zero);
      return _Base::samples();
    }
    
    //        accumulate the weights
    void accumulate(const sample_type& weight) {
      _Traits::assadd(_M_tmp, weight);
    }

    //     reset 
    void reset_tmp() { _Traits::assign(_M_tmp, _Traits::zero);}
    void reset() { _Base::reset(); reset_tmp();}

    //   element access
    const sample_type& get_tmp() const { return _M_tmp;}

  private:
    sample_type _M_tmp;
  };
  
  //   I/O functions for _Point = void
  template<typename _Tp, class _Traits> inline std::ostream& 
  operator<<(std::ostream& os, 
	           const distribution<_Tp, void, _Traits, distpoint_traits<void> >& x) 
  {
    return os<<dynamic_cast<const sample_statistic<_Tp, _Traits>&>(x);
  }

  template<typename _Tp, class _Traits> inline std::istream& 
  operator>>(std::istream& is, 
						 distribution<_Tp, void, _Traits, distpoint_traits<void> >& x) 
  {
    is>>dynamic_cast<sample_statistic<_Tp, _Traits>&>(x);
    x.reset_tmp();
    return is;
  }
  
  template<typename _Tp, class _Traits> inline std::ostream& 
  write(std::ostream& os, 
	      const distribution<_Tp, void, _Traits, distpoint_traits<void> >& x) 
  {
    return write(os, dynamic_cast<const sample_statistic<_Tp, _Traits>&>(x));
  }

  template<typename _Tp, class _Traits> inline std::istream& 
	read(std::istream& is, 
	     distribution<_Tp, void, _Traits, distpoint_traits<void> >& x) 
  {
    read(is, dynamic_cast<sample_statistic<_Tp, _Traits>&>(x));
    x.reset_tmp();
    return is;
  }
  
  template<typename _Tp, class _Traits> inline std::ostream& 
  print(std::ostream& os, 
	      const distribution<_Tp, void, _Traits, distpoint_traits<void> >& d) {
    return print(os, dynamic_cast<const sample_statistic<_Tp, _Traits>&>(d));
  }


  //
  //    General template class 
  //
  template<typename _Tp, class _Point, class _Traits, class _PointTraits>
  class distribution
  {
  private:
    //    private types
    typedef distribution<_Tp, void, _Traits, distpoint_traits<void> > _Dist;
    typedef std::pair<_Point, _Dist> _Value;
    typedef std::vector<_Value> _Base;

  public:
    //   public types
    typedef _Tp sample_type;
    typedef _Traits sample_traits;
    typedef _PointTraits point_traits;

    typedef typename _Base::size_type size_type;
    typedef typename _Base::value_type value_type;
    typedef typename _Base::iterator iterator;
    typedef typename _Base::reverse_iterator reverse_iterator;
    typedef typename _Base::reference reference;
    typedef typename _Base::const_iterator const_iterator;
    typedef typename _Base::const_reference const_reference;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;

    //   constructors
    distribution() {}
    distribution(size_type n, const _Point *base) : _M_base(n) { 
      for(iterator i = _M_base.begin(); i < _M_base.end(); i++)
	i -> first = *(base++);
    }
    
    template<class _InputIterator>
    distribution(_InputIterator first, _InputIterator last)
      : _M_base(first, last) {}
    
    //   computed assignment
    distribution& operator+=(const distribution& x) {
      const_iterator ix = x._M_base.begin();
      for(iterator i = _M_base.begin(); i < _M_base.end(); i++, ix++)
	i -> second += ix -> second;
      return *this;
    }
    
    distribution& assadd_weighted(const distribution& x) {
      const_iterator ix = x._M_base.begin();
      for(iterator i = _M_base.begin(); i < _M_base.end(); i++, ix++)
	i -> second.assadd_weighted(ix -> second);
      return *this;
    }
    
    //   add the samples 
    unsigned long int operator++();
    unsigned long int operator++(int) { 
      return operator++()-1UL;
    }
    
    unsigned long int add_samples_no_check() {
      unsigned long stmp = 0UL;
      for(iterator i = _M_base.begin(); i < _M_base.end(); i++)
	      stmp += i -> second.add_samples_no_check();
      return stmp/_M_base.size();
    }
    
    //   reinitialize the distribution
    void reset() {
      for(iterator i = _M_base.begin(); i < _M_base.end(); i++)
	i -> second.reset();
    }
    
    void reset_tmp() {
      for(iterator i = _M_base.begin(); i < _M_base.end(); i++)
	      i -> second.reset_tmp();
    }
    
    //  capacity
    size_type size() const { return _M_base.size();}
    void resize(size_type n, const value_type& x = value_type()) {
      _M_base.resize(n, x);
    }

    //  element access
    iterator begin() { return _M_base.begin();}
    iterator end() { return _M_base.end();}
    
    const_iterator begin() const { return _M_base.begin();}
    const_iterator end() const { return _M_base.end();}
    
    reverse_iterator rbegin() { return _M_base.rbegin();}
    reverse_iterator rend() { return _M_base.rend();}
    
    const_reverse_iterator rbegin() const { return _M_base.rbegin();}
    const_reverse_iterator rend() const { return _M_base.rend();}
    
    reference operator[](size_type i) { return _M_base[i];}
    const_reference operator[](size_type i) const { return _M_base[i];}
    
    //   accumulate the weights
    template<class _Arg, class _Res> 
    void accumulate(const std::pointer_to_unary_function<_Arg, _Res>& func) {
      this -> _M_acc_unary(func);
    }
    
    template<class _Xp, class _Arg, class _Res> 
    void accumulate(const pointer_to_unary_member_function<_Xp, _Arg, _Res>& func) {
      this -> _M_acc_unary(func);
    }      
    
    template<class _Xp, class _Arg, class _Res> 
    void accumulate(const pointer_to_unary_const_member_function<_Xp, _Arg, _Res>& func) {
      this -> _M_acc_unary(func);
    }     
     
    template<class _Arg1, class _Arg2, class _Res> 
    void accumulate(const std::pointer_to_binary_function<_Arg1, _Arg2, _Res>& func, _Arg1 x, const _Tp& w) {
      this -> _M_acc_binary(func, x, w);
    }

    template<class _Arg1, class _Arg2> 
    void accumulate(const std::pointer_to_binary_function<_Arg1, _Arg2, double>& func, _Arg1 x, const _Tp& w) {
      this -> _M_acc_double_bin(func, x, w);
    }
   
    template<class _Xp, class _Arg1, class _Arg2, class _Res> 
    void accumulate(const pointer_to_binary_member_function<_Xp, _Arg1, _Arg2, _Res>& func, _Arg1 x, const _Tp& w) {
      this -> _M_acc_binary(func, x, w);
    }    
    
    template<class _Xp, class _Arg1, class _Arg2, class _Res> 
    void accumulate(const pointer_to_binary_const_member_function<_Xp, _Arg1, _Arg2, _Res>& func, _Arg1 x, const _Tp& w) {
      this -> _M_acc_binary(func, x, w);
    }    
  
    template<class _Xp, class _Arg1, class _Arg2> 
    void accumulate(const pointer_to_binary_member_function<_Xp, _Arg1, _Arg2, double>& func, _Arg1 x, const _Tp& w) {
      this -> _M_acc_double_bin(func, x, w);
    }
    
    template<class _Xp, class _Arg1, class _Arg2> 
    void accumulate(const pointer_to_binary_const_member_function<_Xp, _Arg1, _Arg2, double>& func, _Arg1 x, const _Tp& w) {
      this -> _M_acc_double_bin(func, x, w);
    }
    
  private:
    //   data members
    _Base _M_base;

    //   private members
    template<class _Func> 
    void _M_acc_unary(const _Func& func) {
      for(iterator i = _M_base.begin(); i < _M_base.end(); i++)
	i -> second.accumulate(func(i -> first));
    }
    
    template<class _Func, class _Arg1> 
    void _M_acc_binary(const _Func&, _Arg1, const sample_type&);
    
    template<class _Func, class _Arg1> 
    void _M_acc_double_bin(const _Func&, _Arg1, const sample_type&);
  };


  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  unsigned long int distribution<_Tp, _Point, _Traits, _PointTraits>::operator++()
  {
    sample_type tmp;
    for(iterator i = _M_base.begin(); i < _M_base.end(); i++) {
      if(!_Traits::finite(i->second.get_tmp())) {
	      this -> reset_tmp();
	      return 0;
      }
      
      _Traits::assign(tmp, i -> second.get_tmp());
      _Traits::assmul(tmp, tmp);
      
      if(!_Traits::finite(tmp)) {
	      this -> reset_tmp();
	      return 0;
      }
    }
    
    return this -> add_samples_no_check();
  }
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits>
  inline distribution<_Tp, _Point, _Traits, _PointTraits>
  operator+(const distribution<_Tp, _Point, _Traits, _PointTraits>& x, const distribution<_Tp, _Point, _Traits, _PointTraits>& y) {
    return distribution<_Tp, _Point, _Traits, _PointTraits>(x) += y;
  }

  template<typename _Tp, class _Point, class _Traits, class _PointTraits>
  inline distribution<_Tp, _Point, _Traits, _PointTraits>
  add_weighted(const distribution<_Tp, _Point, _Traits, _PointTraits>& x, const distribution<_Tp, _Point, _Traits, _PointTraits>& y) {
    return distribution<_Tp, _Point, _Traits, _PointTraits>(x).assadd_weighted(y);
  }
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  template<class _Func, class _Arg1> 
  void distribution<_Tp, _Point, _Traits, _PointTraits>::
  _M_acc_binary(const _Func& func, _Arg1 x, const sample_type& w)
  {
    sample_type tmp;
    for(iterator i = _M_base.begin(); i < _M_base.end(); i++) {
      _Traits::assign(tmp, w);
      _Traits::assmul(tmp, func(x, i -> first));
      i -> second.accumulate(tmp);
    }
  }
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  template<class _Func, class _Arg1> void distribution<_Tp, _Point, _Traits, _PointTraits>::
  _M_acc_double_bin(const _Func& func, _Arg1 x, const sample_type& w)
  {
    double func_tmp;
    sample_type tmp;
   
    for(iterator i = _M_base.begin(); i < _M_base.end(); i++) {
      func_tmp = func(x, i -> first);
      if(func_tmp != 0.0) {
	      _Traits::assign(tmp, w);
	      _Traits::assmul(tmp, func_tmp);
	      i -> second.accumulate(tmp);
      }
    }

  }
    
 
  //      I/O operations
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> std::ostream& 
  operator<<(std::ostream& os, const distribution<_Tp, _Point, _Traits, _PointTraits>& d)
  {
    typedef typename distribution<_Tp, _Point, _Traits, _PointTraits>::const_iterator const_iterator;
    
    os<<d.size()<<"\n";
    for(const_iterator i = d.begin(); i < d.end(); i++) {
      _PointTraits::write_txt(os, i->first);
      os<<"\n"<<i->second<<"\n";
    }
    
    return os;
  }
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::ostream& write(std::ostream& os, const distribution<_Tp, _Point, _Traits, _PointTraits>& d)
  {
    typedef typename distribution<_Tp, _Point, _Traits, _PointTraits>::const_iterator const_iterator;
    typedef typename distribution<_Tp, _Point, _Traits, _PointTraits>::size_type size_type;
    
    size_type n = d.size();
    os.write((const char *) &n, sizeof(size_type));
    
    for(const_iterator i = d.begin(); i < d.end(); i++) {
      _PointTraits::write_bin(os, i->first); write(os, i->second);
    }
    return os;
  }
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::istream& operator>>(std::istream& is, distribution<_Tp, _Point, _Traits, _PointTraits>& d)
  {
    typedef typename distribution<_Tp, _Point, _Traits, _PointTraits>::size_type size_type;
    typedef typename distribution<_Tp, _Point, _Traits, _PointTraits>::iterator iterator;
    
    size_type size;

    is>>size;
    d.resize(size);

    for(iterator i = d.begin(); i < d.end(); i++) {
      _PointTraits::read_txt(is, i->first);
      is>>i->second;
      i -> second.reset_tmp();
    }
    
    return is;
  }

  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::istream& read(std::istream& is, distribution<_Tp, _Point, _Traits, _PointTraits>& d)
  {
    typedef typename distribution<_Tp, _Point, _Traits, _PointTraits>::size_type size_type;
    typedef typename distribution<_Tp, _Point, _Traits, _PointTraits>::iterator iterator;

    size_type size;
    is.read((char *) &size, sizeof(size_type));
    d.resize(size);

    for(iterator i = d.begin(); i < d.end(); i++) {
      _PointTraits::read_bin(is, i->first);
      read(is, i->second);
      i -> second.reset_tmp();
    }
    
    return is;
  }

  template<class _Tp, typename _Point, class _Traits, class _PointTraits> 
  std::ostream& print(std::ostream& os, const distribution<_Tp, _Point, _Traits, _PointTraits>& d)
  {
    typedef typename distribution<_Tp, _Point, _Traits, _PointTraits>::const_iterator const_iterator;

    for(const_iterator i = d.begin(); i < d.end(); i++) {
      _PointTraits::print(os, i->first); os<<"    ";
      _Traits::print(os, i->second.mean(), i->second.error());
      os<<"\n";
    }
    
    return os;
  }  	
	
}  //  namespace nlo

#endif

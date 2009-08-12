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
#ifndef __NLO_HEP_DISTBOOK_H__
#define __NLO_HEP_DISTBOOK_H__ 1


//   Standard includes
#include <map>
#include <string>

//   nlo includes
#include <bits/hep-dist.h>


namespace nlo {

  
  template<typename _Tp, class _Point, 
	   class _Traits = sample_traits<_Tp>,
	   class _PointTraits = distpoint_traits<_Point> >
  class distbook
    : public std::map<int, std::pair<std::basic_string<char>, distribution<_Tp, _Point, _Traits, _PointTraits> > >
  {
  public:
    //   types
    typedef std::basic_string<char> string_type;
    typedef distribution<_Tp, _Point, _Traits, _PointTraits> dist_type;
    typedef _Point point_type;
    typedef _PointTraits point_traits_type;

  private:
    //   private types
    typedef std::map<int, std::pair<string_type, dist_type> > _Base;
    typedef typename _Base::value_type  _Value;

  public:
    //  inherited types
    typedef typename _Base::key_type               key_type;
    typedef typename _Base::mapped_type            mapped_type;
    typedef typename _Base::value_type             value_type;
    typedef typename _Base::key_compare            key_compare;
    typedef typename _Base::allocator_type         allocator_type;
    typedef typename _Base::reference              reference;
    typedef typename _Base::const_reference        const_reference;
    typedef typename _Base::iterator               iterator;
    typedef typename _Base::const_iterator         const_iterator; 
    typedef typename _Base::size_type              size_type;
    typedef typename _Base::difference_type        difference_type;
    typedef typename _Base::pointer                pointer;
    typedef typename _Base::const_pointer          const_pointer;
    typedef typename _Base::reverse_iterator       reverse_iterator;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;

    //  assignments
    distbook& operator+=(const distbook& __x) {
      iterator __iter = _Base::begin();
      const_iterator __xiter = __x.begin();
      while(__iter != _Base::end()) 
	(__iter++ -> second).second += (__xiter++ -> second).second;
      return *this;
    }

    distbook& assadd_weighted(const distbook& __x) {
      iterator __iter = _Base::begin();
      const_iterator __xiter = __x.begin();
      while(__iter != _Base::end()) 
	(__iter++ -> second).second.assadd_weighted((__xiter++ -> second).second);
      return *this;
    }

    //  add the samples 
    unsigned long operator++();
    unsigned long operator++(int) { return operator++() - 1UL;}

    //    create and fill the distributions
    void create(int id, const char *title, unsigned n, const _Point *base) {
      _Base::insert(_Value(id, mapped_type(title, dist_type(n, base))));
    }
    
    //   acumulate the weights 
    template<class _Arg, class _Res> 
    void accumulate(int id, const std::pointer_to_unary_function<_Arg, _Res>& f) {
      (_Base::find(id) -> second).second.accumulate(f);
    }
    
    template<class _Xp, class _Arg, class _Res> 
    void accumulate(int id, const pointer_to_unary_member_function<_Xp, _Arg, _Res>& f) {
      (_Base::find(id) -> second).second.accumulate(f);
    }
    
    template<class _Xp, class _Arg, class _Res> 
    void accumulate(int id, const pointer_to_unary_const_member_function<_Xp, _Arg, _Res>& f) {
      (_Base::find(id) -> second).second.accumulate(f);
    }
        
    template<typename _Arg1, typename _Arg2, typename _Res> 
    void accumulate(int id, const std::pointer_to_binary_function<_Arg1, _Arg2, _Res>& f, _Arg1 x, const _Tp& w) {
      (_Base::find(id) -> second).second.accumulate(f, x, w);
    }
    
    template<class _Xp, class _Arg1, class _Arg2, class _Res> 
    void accumulate(int id, const pointer_to_binary_member_function<_Xp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Tp& w) {
      (_Base::find(id) -> second).second.accumulate(f, x, w);
    }
    
    template<class _Xp, class _Arg1, class _Arg2, class _Res> 
    void accumulate(int id, const pointer_to_binary_const_member_function<_Xp, _Arg1, _Arg2, _Res>& f, _Arg1 x, const _Tp& w) {
      (_Base::find(id) -> second).second.accumulate(f, x, w);
    }
  };
  
 
  template<typename _Tp, class _Point, class _Traits, class _PointTraits>
  unsigned long distbook<_Tp, _Point, _Traits, _PointTraits>::operator++() 
  {
    unsigned long __ret_val = 0UL;
    iterator __iter = _Base::begin();
    while(__iter != _Base::end())
      __ret_val = ++((__iter++ -> second).second);
    return __ret_val;
  }

  template<typename _Tp, class _Point, class _Traits, class _PointTraits>
  inline distbook<_Tp, _Point, _Traits, _PointTraits>
  operator+(const distbook<_Tp, _Point, _Traits, _PointTraits>& __x, 
	    const distbook<_Tp, _Point, _Traits, _PointTraits>& __y) {
    return distbook<_Tp, _Point, _Traits, _PointTraits>(__x) += __y;
  }
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits>
  inline distbook<_Tp, _Point, _Traits, _PointTraits>
  add_weighted(const distbook<_Tp, _Point, _Traits, _PointTraits>& __x, 
	       const distbook<_Tp, _Point, _Traits, _PointTraits>& __y) {
    return distbook<_Tp, _Point, _Traits, _PointTraits>(__x).assadd_weighted(__y);
  }
  
  //   I/O operations
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> std::ostream& 
  operator<<(std::ostream& __os, const distbook<_Tp, _Point, _Traits, _PointTraits>& __x) 
  {
    __os<<__x.size()<<"\n";
    typename distbook<_Tp, _Point, _Traits, _PointTraits>::const_iterator __iter = __x.begin();
    while(__iter != __x.end()) {
      __os<<(__iter -> first)<<" "
	  <<(__iter -> second).first<<"\n"
	  <<(__iter -> second).second<<"\n";
      ++__iter;
    }
    return __os;
  }

  template<typename _Tp, class _Point, class _Traits, class _PointTraits> std::ostream& 
  write(std::ostream& __os, const distbook<_Tp, _Point, _Traits, _PointTraits>& __x) 
  {
    __os<<__x.size()<<"\n";
    typename distbook<_Tp, _Point, _Traits, _PointTraits>::const_iterator __iter = __x.begin();
    while(__iter != __x.end()) {
      __os<<(__iter -> first)<<" "
	  <<(__iter -> second).first<<"\n";
      write(__os, (__iter++ -> second).second);
    }
    return __os;
  }

  template<typename _Tp, class _Point, class _Traits, class _PointTraits> std::ostream& 
  print(std::ostream& __os, const distbook<_Tp, _Point, _Traits, _PointTraits>& __x) 
  {
    typename distbook<_Tp, _Point, _Traits, _PointTraits>::const_iterator __iter = __x.begin();
    while(__iter != __x.end()) {
      __os<<"#  "<<(__iter -> second).first<<"\n";
      print(__os, (__iter -> second).second)<<"\n";
      ++__iter;
    }
    return __os;
  }
  
  template<typename _Tp, class _Point, class _Traits, class _PointTraits> std::istream& 
  operator>>(std::istream& __is, distbook<_Tp, _Point, _Traits, _PointTraits>& __x)
  {
    typedef typename distbook<_Tp, _Point, _Traits, _PointTraits>::mapped_type mapped_type;
    typedef typename distbook<_Tp, _Point, _Traits, _PointTraits>::value_type value_type;

    int __key;
    unsigned int __size;
    typename distbook<_Tp, _Point, _Traits, _PointTraits>::string_type __title;
    typename distbook<_Tp, _Point, _Traits, _PointTraits>::dist_type __dist;

    __x.clear();
    __is>>__size;
    for(unsigned int __i = 0; __i < __size; __i++) { 
      __is>>__key; std::getline(__is, __title); __is>>__dist;
      __x.insert(value_type(__key, mapped_type(__title, __dist)));
    }
    
    return __is;
  }

  template<typename _Tp, class _Point, class _Traits, class _PointTraits> 
  std::istream& read(std::istream& __is, distbook<_Tp, _Point, _Traits, _PointTraits>& __x)
  {
    typedef typename distbook<_Tp, _Point, _Traits, _PointTraits>::mapped_type mapped_type;
    typedef typename distbook<_Tp, _Point, _Traits, _PointTraits>::value_type value_type;

    int __key;
    unsigned int __size;
    typename distbook<_Tp, _Point, _Traits, _PointTraits>::string_type __title;
    typename distbook<_Tp, _Point, _Traits, _PointTraits>::dist_type __dist;

    __x.clear();
    __is>>__size;
    for(unsigned int __i = 0; __i < __size; __i++) { 
      __is>>__key; std::getline(__is, __title); read(__is, __dist);
      __x.insert(value_type(__key, mapped_type(__title, __dist)));
    }
    
    return __is;
  }


  //   Specialization
  template<typename _Tp, class _Traits>
  class distbook<_Tp, void, _Traits, distpoint_traits<void> >
    : public std::map<int, std::pair<std::basic_string<char>, distribution<_Tp, void, _Traits, distpoint_traits<void> > > >
  {
  public:
    //   types
    typedef std::basic_string<char> string_type;
    typedef distribution<_Tp, void, _Traits, distpoint_traits<void> > dist_type;

  private:
    //   private types
    typedef std::map<int, std::pair<string_type, dist_type> > _Base;
    typedef typename _Base::value_type  _Value;
    
  public:
    //  inherited types
    typedef typename _Base::key_type               key_type;
    typedef typename _Base::mapped_type            mapped_type;
    typedef typename _Base::value_type             value_type;
    typedef typename _Base::key_compare            key_compare;
    typedef typename _Base::allocator_type         allocator_type;
    typedef typename _Base::reference              reference;
    typedef typename _Base::const_reference        const_reference;
    typedef typename _Base::iterator               iterator;
    typedef typename _Base::const_iterator         const_iterator; 
    typedef typename _Base::size_type              size_type;
    typedef typename _Base::difference_type        difference_type;
    typedef typename _Base::pointer                pointer;
    typedef typename _Base::const_pointer          const_pointer;
    typedef typename _Base::reverse_iterator       reverse_iterator;
    typedef typename _Base::const_reverse_iterator const_reverse_iterator;

    //  assignments
    distbook& operator+=(const distbook& __x) {
      iterator __iter = _Base::begin();
      const_iterator __xiter = __x.begin();
      while(__iter != _Base::end()) 
				(__iter++ -> second).second += (__xiter++ -> second).second;
      return *this;
    }
    
    distbook& assadd_weighted(const distbook& __x) {
      iterator __iter = _Base::begin();
      const_iterator __xiter = __x.begin();
      while(__iter != _Base::end()) 
				(__iter++ -> second).second.assadd_weighted((__xiter++ -> second).second);
      return *this;
    }
    
    //  add the samples 
    unsigned long operator++();
    unsigned long operator++(int) { return operator++() - 1UL;}
    
    //    create and fill the distributions
    void create(int id, const char *title) {
      _Base::insert(_Value(id, mapped_type(title, dist_type())));
    }
    
    void accumulate(int id, const _Tp& w) { 
      (_Base::find(id) -> second).second.accumulate(w);
    }
  };
 
  template<typename _Tp, class _Traits> unsigned long 
	distbook<_Tp, void, _Traits, distpoint_traits<void> >::operator++() 
  {
    unsigned long __ret_val = 0UL;
    iterator __iter = _Base::begin();
    while(__iter != _Base::end()) 
      __ret_val = ++((__iter++ -> second).second);
    return __ret_val;
  }
}   //   namespace nlo

#endif

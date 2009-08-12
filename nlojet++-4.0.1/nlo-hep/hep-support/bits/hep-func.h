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
#ifndef __HEP_FUNC_H__
#define __HEP_FUNC_H__

#include <functional>


namespace nlo {


  template <class _Tp, class _Arg, class _Res>
  class pointer_to_unary_member_function 
    : public std::unary_function<_Arg, _Res> 
  {
  public:
    pointer_to_unary_member_function(_Tp& ref, _Res (_Tp::*func)(_Arg))
      : _M_ptr(&ref), _M_func(func) {}
    
    pointer_to_unary_member_function(_Tp *ptr, _Res (_Tp::*func)(_Arg))
      : _M_ptr(ptr), _M_func(func) {}
    
    _Res operator()(_Arg x) const {
      return (_M_ptr ->* _M_func)(x);
    }
    
  private:
    _Tp *_M_ptr;
    _Res (_Tp::*_M_func)(_Arg);
  };

  template <class _Tp, class _Arg, class _Res>
  class pointer_to_unary_const_member_function 
    : public std::unary_function<_Arg, _Res> 
  {
  public:
    pointer_to_unary_const_member_function(const _Tp& ref,
					   _Res (_Tp::*func)(_Arg) const)
      : _M_ptr(&ref), _M_func(func) {}
    
    pointer_to_unary_const_member_function(const _Tp *ptr,
					   _Res (_Tp::*func)(_Arg) const)
      : _M_ptr(ptr), _M_func(func) {}
    
    _Res operator()(_Arg x) const {
      return (_M_ptr ->* _M_func)(x);
    }
    
  private:
    const _Tp *_M_ptr;
    _Res (_Tp::*_M_func)(_Arg) const;
  };


  template <class _Tp, class _Arg1, class _Arg2, class _Res>
  class pointer_to_binary_member_function 
    : public std::binary_function<_Arg1, _Arg2, _Res> 
  {
  public:
    pointer_to_binary_member_function(_Tp& ref, _Res (_Tp::*func)(_Arg1, _Arg2))
      : _M_ptr(&ref), _M_func(func) {}
    
    pointer_to_binary_member_function(_Tp *ptr, _Res (_Tp::*func)(_Arg1, _Arg2))
      : _M_ptr(ptr), _M_func(func) {}
    
    _Res operator()(_Arg1 x, _Arg2 y) const {
      return (_M_ptr ->* _M_func)(x, y);
    }
    
  private:
    _Tp *_M_ptr;
    _Res (_Tp::*_M_func)(_Arg1, _Arg2);
  };

  template <class _Tp, class _Arg1, class _Arg2, class _Res>
  class pointer_to_binary_const_member_function 
    : public std::binary_function<_Arg1, _Arg2, _Res> 
  {
  public:
    pointer_to_binary_const_member_function(const _Tp& ref,
					    _Res (_Tp::*func)(_Arg1, _Arg2) const)
      : _M_ptr(&ref), _M_func(func) {}
    
    pointer_to_binary_const_member_function(const _Tp *ptr,
					    _Res (_Tp::*func)(_Arg1, _Arg2) const)
      : _M_ptr(ptr), _M_func(func) {}
    
    _Res operator()(_Arg1 x, _Arg2 y) const {
      return (_M_ptr ->* _M_func)(x, y);
    }
    
  private:
    const _Tp *_M_ptr;
    _Res (_Tp::*_M_func)(_Arg1, _Arg2) const;
  };


  template <class _Tp, class _Arg>
  class pointer_to_unary_member_function<_Tp, _Arg, void>
    : public std::unary_function<_Arg, void> 
  {
  public:
    pointer_to_unary_member_function(_Tp& ref, void (_Tp::*func)(_Arg))
      : _M_ptr(&ref), _M_func(func) {}
    
    pointer_to_unary_member_function(_Tp *ptr, void (_Tp::*func)(_Arg))
      : _M_ptr(ptr), _M_func(func) {}
    
    void operator()(_Arg x) const {
      (_M_ptr ->* _M_func)(x);
    }
    
  private:
    _Tp *_M_ptr;
    void (_Tp::*_M_func)(_Arg);
  };
  
  template <class _Tp, class _Arg>
  class pointer_to_unary_const_member_function<_Tp, _Arg, void>
    : public std::unary_function<_Arg, void> 
  {
  public:
    pointer_to_unary_const_member_function(const _Tp& ref,
					   void (_Tp::*func)(_Arg) const)
      : _M_ptr(&ref), _M_func(func) {}
    
    pointer_to_unary_const_member_function(const _Tp *ptr,
					   void (_Tp::*func)(_Arg) const)
      : _M_ptr(ptr), _M_func(func) {}
    
    void operator()(_Arg x) const {
      (_M_ptr ->* _M_func)(x);
    }
    
  private:
    const _Tp *_M_ptr;
    void (_Tp::*_M_func)(_Arg) const;
  };
  
  
  template <class _Tp, class _Arg1, class _Arg2>
  class pointer_to_binary_member_function<_Tp, _Arg1, _Arg2, void> 
    : public std::binary_function<_Arg1, _Arg2, void> 
  {
  public:
    pointer_to_binary_member_function(_Tp& ref, void (_Tp::*func)(_Arg1, _Arg2))
      : _M_ptr(&ref), _M_func(func) {}
    
    pointer_to_binary_member_function(_Tp *ptr, void (_Tp::*func)(_Arg1, _Arg2))
      : _M_ptr(ptr), _M_func(func) {}
    
    void operator()(_Arg1 x, _Arg2 y) const {
      (_M_ptr ->* _M_func)(x, y);
    }
    
  private:
    _Tp *_M_ptr;
    void (_Tp::*_M_func)(_Arg1, _Arg2);
  };
  
  template <class _Tp, class _Arg1, class _Arg2>
  class pointer_to_binary_const_member_function<_Tp, _Arg1, _Arg2, void> 
    : public std::binary_function<_Arg1, _Arg2, void> 
  {
  public:
    pointer_to_binary_const_member_function(const _Tp& ref,
					    void (_Tp::*func)(_Arg1, _Arg2) const)
      : _M_ptr(&ref), _M_func(func) {}
    
    pointer_to_binary_const_member_function(const _Tp *ptr,
					    void (_Tp::*func)(_Arg1, _Arg2) const)
      : _M_ptr(ptr), _M_func(func) {}
    
    void operator()(_Arg1 x, _Arg2 y) const {
      (_M_ptr ->* _M_func)(x, y);
    }
    
  private:
    const _Tp *_M_ptr;
    void (_Tp::*_M_func)(_Arg1, _Arg2) const;
  };



  template<class Tp, class Arg, class Result>
  inline pointer_to_unary_member_function<Tp, Arg, Result>
  ptr_fun(Tp *ptr, Result (Tp::*f)(Arg)) {
    return pointer_to_unary_member_function<Tp, Arg, Result>(ptr, f);
  }

  template<class Tp, class Arg, class Result>
  inline pointer_to_unary_member_function<Tp, Arg, Result>
  ptr_fun(Tp& ref, Result (Tp::*f)(Arg)) {
    return pointer_to_unary_member_function<Tp, Arg, Result>(ref, f);
  }
  
  template<class Tp, class Arg, class Result>
  inline pointer_to_unary_const_member_function<Tp, Arg, Result>
  ptr_fun(const Tp *ptr, Result (Tp::*f)(Arg) const) {
    return pointer_to_unary_const_member_function<Tp, Arg, Result>(ptr, f);
  }
  
  template<class Tp, class Arg, class Result>
  inline pointer_to_unary_const_member_function<Tp, Arg, Result>
  ptr_fun(const Tp& ref, Result (Tp::*f)(Arg) const) {
    return pointer_to_unary_const_member_function<Tp, Arg, Result>(ref, f);
  }


  template<class Tp, class Arg1, class Arg2, class Result>
  inline pointer_to_binary_member_function<Tp, Arg1, Arg2, Result>
  ptr_fun(Tp *ptr, Result (Tp::*f)(Arg1, Arg2)) {
    return pointer_to_binary_member_function<Tp, Arg1, Arg2, Result>(ptr, f);
  }

  template<class Tp, class Arg1, class Arg2, class Result>
  inline pointer_to_binary_member_function<Tp, Arg1, Arg2, Result>
  ptr_fun(Tp& ref, Result (Tp::*f)(Arg1, Arg2)) {
    return pointer_to_binary_member_function<Tp, Arg1, Arg2, Result>(ref, f);
  }
  
  template<class Tp, class Arg1, class Arg2, class Result>
  inline pointer_to_binary_const_member_function<Tp, Arg1, Arg2, Result>
  ptr_fun(const Tp *ptr, Result (Tp::*f)(Arg1, Arg2) const) {
    return pointer_to_binary_const_member_function<Tp, Arg1, Arg2, Result>(ptr, f);
  }
  
  template<class Tp, class Arg1, class Arg2, class Result>
  inline pointer_to_binary_const_member_function<Tp, Arg1, Arg2, Result>
  ptr_fun(const Tp& ref, Result (Tp::*f)(Arg1, Arg2) const) {
    return pointer_to_binary_const_member_function<Tp, Arg1, Arg2, Result>(ref, f);
  }
}   //  namespace nlo

#endif


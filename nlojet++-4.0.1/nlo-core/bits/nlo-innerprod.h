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
#ifndef __NLO_NLO_INNERPROD_H__
#define __NLO_NLO_INNERPROD_H__ 1

//   Standard includes
#include <complex>

//   nlo includes
#include <bits/nlo-event.h>


namespace nlo {

  
  template<class _Lv>
  class innerprod
  {
  public:
    //   types
    typedef partonic_event<_Lv>             event_type;
    typedef typename _Lv::value_type        value_type;
    typedef typename event_type::index_type index_type;
    typedef typename event_type::size_type  size_type;
    
    //   the matrix type which help us to store the inner products
    template<typename _Tp>
    class matrix {
    public:
      //   constructors
      matrix() 
	: ptr(0), low(0), high(-1) {}
      
      matrix(index_type __low, index_type __high)
	: ptr(_S_allocate(__low, __high)), low(__low), high(__high) {}
 
      matrix(const matrix<_Tp>& __x)
	: ptr(_S_allocate(__x.low, __x.high)), low(__x.low), high(__x.high)
	{ this -> _M_copy(__x);} 
      
      ~matrix() { 
	_S_deallocate(ptr, low, high);
      }
      
      //   assignment
      matrix& operator=(const matrix& __x) {
	if(&__x != this) {
	  this -> resize(__x.low, __x.high);
	  this -> _M_copy(__x);
	}
	return *this;
      }
      
      //   element acces
      _Tp* operator[](index_type __r) { return ptr[__r];}
      const _Tp* operator[](index_type __r) const { return ptr[__r];}

      //   resize object
      void resize(index_type __low, index_type __high) {
	if(low != __low || __high > high) {
	  _S_deallocate(ptr, low, high);
	  ptr = _S_allocate(low = __low, high = __high);
	}
      }

      //   data members
      _Tp **ptr;
      index_type low, high;
      
    private:
      void _M_copy(const matrix&);
      static _Tp** _S_allocate(index_type, index_type);
      static void  _S_deallocate(_Tp **, index_type, index_type);
    };
    
    typedef matrix<value_type>                lorentz_prod_type;
    typedef matrix<std::complex<value_type> > spinor_prod_type;
    
    //   calculate the inner products form the partonic event
    void calculate(const event_type&);
    void calculate_s(const event_type&);
    void calculate_ab(const event_type&);

    //  elements access
    const value_type& s(index_type __i, index_type __j) const { 
      return _M_s.ptr[__i][__j];
    }
    
    const std::complex<value_type>& a(index_type __i, index_type __j) const { 
      return _M_a.ptr[__i][__j];
    }
    
    const std::complex<value_type>& b(index_type __i, index_type __j) const { 
      return _M_b.ptr[__i][__j];
    }
    
    const lorentz_prod_type& s() const { return _M_s;}
    const spinor_prod_type&  a() const { return _M_a;}
    const spinor_prod_type&  b() const { return _M_b;}
    
  protected:
    //    innner products
    lorentz_prod_type _M_s;     //  2*pi*pj
    spinor_prod_type  _M_a;     //  <pi|pj> 
    spinor_prod_type  _M_b;     //  [pi|pj]
  };
  

  template<class _Lv> template<typename _Tp>
  _Tp** innerprod<_Lv>::matrix<_Tp>::
  _S_allocate(index_type __low, index_type __up) 
  {
    size_type __n = __up - __low + 1;
    _Tp **__ptr = new _Tp*[__n] - __low;
    
    for(index_type __i = __low; __i <= __up; __i++)
      __ptr[__i] = new _Tp[__n] - __low;
    
    return __ptr;
  }

  template<class _Lv> template<typename _Tp> 
  void innerprod<_Lv>::matrix<_Tp>::
  _S_deallocate(_Tp **__ptr, index_type __low, index_type __up)
  {
    if(__ptr + __low) {
      for(index_type __i = __low; __i <= __up; __i++)
        if(__ptr[__i] + __low) delete [] (__ptr[__i] + __low);
      delete [] (__ptr + __low);
    }
  }

  template<class _Lv> template<typename _Tp> 
  void innerprod<_Lv>::matrix<_Tp>::
  _M_copy(const innerprod<_Lv>::matrix<_Tp>& __x)
  {
    for(index_type __i = low; __i <= high; __i++)
      for(index_type __j = low; __j <= high; __j++)
        ptr[__i][__j] = __x.ptr[__i][__j];
  }
  
  template<class _Lv> 
  void innerprod<_Lv>::calculate(const event_type& __p) 
  {
    //----- resize the objects -----
    index_type __new_low = __p.lower();
    index_type __new_up  = __p.upper();
    
    _M_s.resize(__new_low, __new_up);
    _M_a.resize(__new_low, __new_up);
    _M_b.resize(__new_low, __new_up);
    
    //---- calculate the inner products ----
    value_type __tmp_s, __ii, __jj;
    std::complex<value_type> __tmp_a;
    
    for(index_type __i = __new_low; __i < __new_up; __i++)
      for(index_type __j = __i + 1; __j <= __new_up; __j++) {
        __ii = (__i == -1 || __i == 0 ? -1.0 : 1.0);
        __jj = (__j == -1 || __j == 0 ? -1.0 : 1.0);
        
        __tmp_s = 2.0*__ii*__jj*dot(__p[__i], __p[__j]);
        __tmp_a = scalar_mp(__ii*__p[__i], __jj*__p[__j]);
        
        _M_s.ptr[__j][__i] = _M_s.ptr[__i][__j] = __tmp_s;
	_M_a.ptr[__j][__i] = -(_M_a.ptr[__i][__j] = __tmp_a);
        _M_b.ptr[__j][__i] = -(_M_b.ptr[__i][__j] = -__tmp_s/__tmp_a);
      }
  }

  template<class _Lv> 
  void innerprod<_Lv>::calculate_s(const event_type& __p) 
  {
    //----- resize the objects -----
    index_type __new_low = __p.lower();
    index_type __new_up  = __p.upper();
    
    _M_s.resize(__new_low, __new_up);
    
    //---- calculate the inner products ----
    value_type __ii, __jj;
    
    for(index_type __i = __new_low; __i < __new_up; __i++)
      for(index_type __j = __i + 1; __j <= __new_up; __j++) {
        __ii = (__i == -1 || __i == 0 ? -1.0 : 1.0);
        __jj = (__j == -1 || __j == 0 ? -1.0 : 1.0);
        
	_M_s.ptr[__j][__i] = _M_s.ptr[__i][__j] = 
	  2.0*__ii*__jj*dot(__p[__i], __p[__j]);
      }
  }

  template<class _Lv>
  void innerprod<_Lv>::calculate_ab(const event_type& __p) 
  {
    //----- resize the objects -----
    index_type __new_low = __p.lower();
    index_type __new_up  = __p.upper();
    
    _M_a.resize(__new_low, __new_up);
    _M_b.resize(__new_low, __new_up);
    
    //---- calculate the inner products ----
    value_type __tmp_s, __ii, __jj;
    std::complex<value_type> __tmp_a;
    
    for(index_type __i = __new_low; __i < __new_up; __i++)
      for(index_type __j = __i + 1; __j <= __new_up; __j++) {
        __ii = (__i == -1 || __i == 0 ? -1.0 : 1.0);
        __jj = (__j == -1 || __j == 0 ? -1.0 : 1.0);
        
        __tmp_s = 2.0*__ii*__jj*dot(__p[__i], __p[__j]);
        __tmp_a = scalar_mp(__ii*__p[__i], __jj*__p[__j]);
        
	_M_a.ptr[__j][__i] = -(_M_a.ptr[__i][__j] = __tmp_a);
        _M_b.ptr[__j][__i] = -(_M_b.ptr[__i][__j] = -__tmp_s/__tmp_a);
      }
  }
}   //  namespace nlo

#endif

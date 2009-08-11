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
#ifndef __NLO_NLO_PROCESS_I0F0_H__
#define __NLO_NLO_PROCESS_I0F0_H__ 1


//   nlojet++ includes
#include <bits/nlo-dipole_mom.h>
#include <bits/nlo-process.h>


namespace nlo {
  

  //
  //  Declaraton of the abstract class process in the 
  //  case of e+e- annihilation. 
  template<typename _Weight, class _Event, class _EvenTraits>
  class process<_Weight, _Event, _EvenTraits, 0U, 0U>
  {
  public:
    //   types 
    typedef _Weight weight_type;
    typedef _Event event_type;

    //   member access
    unsigned int npar() const { return _M_npar;}
    unsigned int npow() const { return _M_npow;}
    double alpha() const { return _M_alpha;}

    //  set the alpha parameter
    void alpha(double __new_alpha) { _M_alpha = __new_alpha;}
    
    //   destructor
    virtual ~process() {} 
    
    //   born contributions
    virtual void born_term(const _Event&, _Weight&) = 0;
    
    //   real contributions
    virtual void real_term(const _Event&, _Weight&) = 0;
    
    //   finite contributions (1-loop, ...)
    virtual void fini_term(const _Event&, _Weight&) = 0;
    
    //   dipole contributions
    virtual void dipole_term(const _Event&, const _Event&, int, int, int, _Weight&) = 0;
    
    //   generate the dipole momenta
    bool dipole_mom(const _Event& p, int i, int j, int k, _Event& q) {
      return dipole_mom_fff<_Event>(_M_alpha, p, i, j, k, q);
    }

  protected:
    //  constructors
    process(unsigned int np, unsigned int nw, double al = 1.0)
      : _M_npar(np), _M_npow(nw), _M_alpha(al) {}
    
  private:
    //   data members
    unsigned int _M_npar, _M_npow;
    double _M_alpha;
  };
  
 
  template<typename _Weight, class _Event, class _EvenTraits>
  class amplitude<process<_Weight, _Event, _EvenTraits, 0U, 0U> >
  {
  public:
    //   public types
    typedef process<_Weight, _Event, _EvenTraits, 0U, 0U> process_type;
    typedef typename process_type::event_type event_type;
    typedef typename process_type::weight_type weight_type;
    
    //   type of the contributions
    enum contrib_type { born = 0, real, sub, fini};
    enum integral_type { lo = 0, nlo};

    //   constructor
    explicit amplitude(process_type *proc, integral_type itype)
      : _M_proc(proc), _M_itype(itype) {}
    
    //   leading  order contribution
    void leading_order(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = born;  
    }
    
    //    next-to-leading order contributions
    void next_to_leading_order_real(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = real;  
    }
    
    void next_to_leading_order_fini(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_calc = false; _M_contr = fini;  
    }
    
    void next_to_leading_order_sub(double w, const _Event& p, const _Event& dp,
				   int i, int j, int k) 
    { 
      _M_p = &p; _M_dp = &dp; _M_weight = w;
      _M_i = i; _M_j = j; _M_k = k;
      _M_calc = false; _M_contr = sub;  
    }
    
    //    get the value of the current amplitude
    const _Weight& operator()() const;
    
    //    get number of the jets
    unsigned int npar() const { return _M_proc -> npar();}
    unsigned int npow() const { return _M_proc -> npow();}
    
    //    get the contributionb type of the current amplitude
    contrib_type contrib() const { return _M_contr;}
    integral_type integral() const { return _M_itype;}

  private:
    //   pointer to the process
    process_type *_M_proc;
    
    //   store the current process type and function
    contrib_type _M_contr;
    integral_type _M_itype;

    //   store the original and dipole event
    double _M_weight;
    const _Event *_M_p, *_M_dp;
    int _M_i, _M_j, _M_k;
    
    //   value of the current amp
    mutable bool _M_calc;
    mutable _Weight _M_amp;
  };    

  template<typename _Weight, class _Event, class _EvenTraits> const _Weight& 
  amplitude<process<_Weight, _Event, _EvenTraits, 0U, 0U> >::operator()() const
  {
    if(!_M_calc) {
      switch(_M_contr) {
      case born: _M_proc -> born_term(*_M_p, _M_amp); break;
      case real: _M_proc -> real_term(*_M_p, _M_amp); break;
      case fini: _M_proc -> fini_term(*_M_p, _M_amp); break;
      case sub:  
	_M_proc -> dipole_term(*_M_p, *_M_dp, _M_i, _M_j, _M_k, _M_amp); 
	break;
      }
      
      _M_calc = true;
      _M_amp *= _M_weight;
    }
    
    return _M_amp;
  }
}    //   namespace nlo

#endif

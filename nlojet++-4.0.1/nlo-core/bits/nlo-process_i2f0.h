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
#ifndef __NLO_PROCESS_I2F0_H__
#define __NLO_PROCESS_I2F0_H__ 1


//   nlojet++ includes
#include <bits/nlo-dipole_mom.h>
#include <bits/nlo-process.h>



namespace nlo {
  
 
  //
  //  Declaraton of the abstract class process in the 
  //  case of hadron-hadron collision. 
  template<typename _Weight, class _Event, class _EvenTraits>
  class process<_Weight, _Event, _EvenTraits, 2U, 0U>
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
    
    //  number of quark flavours
    unsigned int nu() const { return _M_nu;}
    unsigned int nd() const { return _M_nd;}
    unsigned int nf() const { return _M_nu + _M_nd;}
    
    //   destructor
    virtual ~process() {} 
    
    //   born contributions
    virtual void born_term(const _Event&, _Weight&) = 0;
    
    //   real contributions
    virtual void real_term(const _Event&, _Weight&) = 0;
    
    //  finite contributions (1-loop,...)
    virtual void fini_term(double, double, double, double, const _Event&, _Weight *) = 0;
    
    //   dipole contributions
    virtual void dipole_term(const _Event&, const _Event&, int, int, int, _Weight&) = 0;
    
    //   select the dipoles (extra contraints for the dipole indices 
    //     e.g: in photoproduction k > -1 no photon spectator)
    virtual bool dipole_index(int, int, int) { 
      return true;
    }

    //   generate the dipole momenta
    bool dipole_mom(const _Event&, int, int, int, _Event&);

  protected:
    //  constructors
    explicit process(unsigned int np, unsigned int nw, unsigned int nu, 
		     unsigned int nd, double al = 1.0)
      : _M_npar(np), _M_npow(nw), _M_nu(nu), _M_nd(nd), _M_alpha(al) {}
    
  private:
    //   data members
    unsigned int _M_npar, _M_npow, _M_nu, _M_nd;
    double _M_alpha;
  };
  
  template<typename _Weight, class _Event, class _EvenTraits>
  bool process<_Weight, _Event, _EvenTraits, 2U, 0U>::
  dipole_mom(const _Event& p, int i, int j, int k, _Event& q)
  {
    q[hadron(-1)] = p[hadron(-1)];
    q[hadron( 0)] = p[hadron( 0)];
    if(i > 0 && j > 0 && k > 0)
      return dipole_mom_fff<_Event>(_M_alpha, p, i, j, k, q);
    else if(i > 0 && j > 0 && k <= 0)  
      return dipole_mom_ffi<_Event>(_M_alpha, p, i, j, k, q);
    else if(i <= 0 && j > 0 && k > 0) 
      return dipole_mom_iff<_Event>(_M_alpha, p, i, j, k, q);
    else if(i <= 0 && j > 0 && k <= 0) 
      return dipole_mom_ifi<_Event>(_M_alpha, p, i, j, k, q);
    else return false;
  }

  template<class _Weight>
  struct pdf_and_coupling<_Weight, 2U, 0U> 
  {
    //   Destructor
    virtual ~pdf_and_coupling() {}
    
    //   The QCD coupling
    virtual double alpha_qcd(unsigned int, double) = 0;
    
    //   the parton distribution function
    virtual _Weight pdf(double, double, double, unsigned int=2U, unsigned int=3U) = 0;
  };
  
  template<typename _Weight, class _Event, class _EvenTraits>
  class amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >
  {
  public:
    //   public types
    typedef process<_Weight, _Event, _EvenTraits, 2U, 0U> process_type;
    typedef typename process_type::event_type event_type;
    typedef typename process_type::weight_type weight_type;
    typedef pdf_and_coupling<_Weight, 2U, 0U> pdf_type;
    
    //   type of the contributions
    enum contrib_type { notdef = -1, born = 0, real, sub, fini};
    enum integral_type { lo = 0, nlo};
    
    //   constructor
    explicit amplitude(process_type *proc, integral_type itype)
      : _M_proc(proc), _M_itype(itype) 
    { _M_eta1 = _M_eta2 = _M_mr2 = _M_mf2 = -1.0;}
    
    //   leading  order contribution
    void leading_order(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_pdf_als = 0;
      _M_calc = false; _M_contr = born; _M_coef = 1.0; 
    }
    
    //    next-to-leading order contributions
    void next_to_leading_order_real(double w, const _Event& p) {
      _M_p = &p; _M_weight = w; _M_pdf_als = 0;
      _M_calc = false; _M_contr = real; _M_coef = 1.0; 
    }
    
    void next_to_leading_order_fini(double x1, double xjac1, double x2,
				    double xjac2, double w, const _Event& p) 
    {
      _M_p = &p; _M_weight = w; _M_pdf_als = 0;
      _M_calc = false; _M_contr = fini; _M_coef = 1.0; 

      //   fini specific options
      _M_fini.x1 = x1; _M_fini.xjac1 = xjac1; 
      _M_fini.x2 = x2; _M_fini.xjac2 = xjac2;
      _M_fini.mode = 0;
    }
    
    void next_to_leading_order_finix1() { _M_fini.mode = 0;}
    void next_to_leading_order_finix2() { _M_fini.mode = 1;}
    void next_to_leading_order_fini1()  { _M_fini.mode = 2;}
        
    void next_to_leading_order_sub(double w, const _Event& p, 
				   const _Event& dp, int i, int j, int k) 
    { 
      _M_p = &p; _M_weight = w; _M_pdf_als = 0;
      _M_calc = false; _M_contr = sub; _M_coef = 1.0; 

      //   dipole specific options
      _M_dipole.dp = &dp; _M_dipole.i = i; 
      _M_dipole.j = j; _M_dipole.k = k; 
    }
    
    //    get the value of the current amplitude
    const _Weight& operator()(double, double) const;
    
    //    get number of the jets
    unsigned int npar() const { return _M_proc -> npar();}
    unsigned int npow() const { return _M_proc -> npow();}
    
    //  number of quark flavours
    unsigned int nu() const { return _M_proc -> nu();}
    unsigned int nd() const { return _M_proc -> nd();}
    unsigned int nf() const { return _M_proc -> nf();}

    //    get the contributionb type of the current amplitude
    contrib_type contrib() const { return _M_contr;}
    integral_type integral() const { return _M_itype;}

    void pdf_and_qcd_coupling(pdf_type *pa, double cf=1.0) const {
      _M_coef = cf; _M_pdf_als = pa;
      _M_eta1 = _M_eta2 = _M_mr2 = _M_mf2 = -1.0;
    }

    void pdf_and_qcd_coupling(pdf_type& pa, double cf=1.0) const {
      _M_coef = cf; _M_pdf_als = &pa;
      _M_eta1 = _M_eta2 = _M_mr2 = _M_mf2 = -1.0;
    }

  private:
    //   pointer to the process
    process_type *_M_proc;
    
    //   store the current process type and function
    contrib_type _M_contr;
    integral_type _M_itype;

    //   store the event
    const _Event *_M_p;
    double _M_weight;
    
    //   value of the current amp
    mutable bool _M_calc;
    mutable _Weight _M_amp[2]; 
    
    //   scales 
    mutable double _M_mr2, _M_mf2;
    
    //   the pdf and the couplings
    mutable pdf_type *_M_pdf_als;
    mutable double  _M_als, _M_eta1, _M_eta2, _M_coef;
    mutable _Weight _M_pdf;
    
    //   contribution specific variables
    //   dipole contribution
    struct {
      //   dipole event
      const _Event *dp;
      
      //   dipole indices
      int i, j, k;
      
      void operator()(process_type *proc, const _Event *p, _Weight& res) const {
	proc -> dipole_term(*p, *dp, i, j, k, res);
      }
    } _M_dipole;
    
    //   finite contributions
    struct {
      //  finite type
      unsigned int mode;
      
      //  the x1, x2 integral
      double x1, xjac1, x2, xjac2;
      
      //   log of the scales 
      mutable double lxr, lxf;
      
      //  scale independent decomposation
      mutable _Weight amp[7];
      
      void operator()(process_type *proc, const _Event *p) const {
	proc -> fini_term(x1, xjac1, x2, xjac2, *p, amp);
      }
      
      void finix1(_Weight& res) const { 
	res = amp[0] + amp[3]*lxf;
      }
      
      void finix2(_Weight& res) const { 
	res = amp[1] + amp[4]*lxf;
      }
      
      void fini1(_Weight& res) const { 
	res = amp[2] + amp[5]*lxf + amp[6]*lxr;
      }
    } _M_fini;
    
    //   private member functions
    void update_scales(double, double) const; 
    
    //   contributions
    void amp_born() const;
    void amp_real() const;
    void amp_dipole() const;
    void amp_fini() const;
  };    
 
  
  template<typename _Weight, class _Event, class _EvenTraits> 
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::
  update_scales(double mr2, double mf2) const
  {
    if(_M_contr == fini) {
      if(_M_mr2 != mr2) _M_fini.lxr = std::log(mr2); 
      if(_M_mf2 != mf2) _M_fini.lxf = std::log(mf2);
    }
    
    if(_M_pdf_als) {
      const _Event& p = *_M_p;
      double s = 2.0*(p[hadron(-1)]*p[hadron(0)]);
      double eta1 = 2.0*(p[-1]*p[hadron(0)])/s;
      double eta2 = 2.0*(p[0]*p[hadron(-1)])/s;
      
      if(_M_contr == fini && _M_fini.mode == 0) eta1 /= _M_fini.x1;
      if(_M_contr == fini && _M_fini.mode == 1) eta2 /= _M_fini.x2;
      
      if(_M_mr2 != mr2) {
	_M_mr2 = mr2;
	_M_als = _M_pdf_als -> alpha_qcd(this -> nf(), _M_mr2);
	
	if(_M_contr == born) _M_als = std::pow(_M_als, (int) (_M_proc->npow()));
	else _M_als = std::pow(_M_als, (int) (_M_proc->npow() + 1));
	_M_als *= _M_coef;
      }
      
      if(_M_mf2 != mf2 || _M_eta1 != eta1 || _M_eta2 != eta2) {
	_M_mf2 = mf2; _M_eta1 = eta1; _M_eta2 = eta2;
	_M_pdf = _M_pdf_als -> pdf(eta1, eta2, mf2, this->nu(), this->nd());
      }
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits> 
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_born() const
  {
    if(!_M_calc) { 
      //----- calculate the matrix element squares -----
      _M_proc -> born_term(*_M_p, _M_amp[0]);
      
      //----- phase space weight -----
      _M_amp[0] *= _M_weight;
      _M_calc = true;
    }

    if(_M_pdf_als) _M_amp[1] = _M_amp[0]*_M_pdf*_M_als;
  }
  
  template<typename _Weight, class _Event, class _EvenTraits> 
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_real() const
  {
    if(!_M_calc) { 
      //----- calculate the matrix element squares -----
      _M_proc -> real_term(*_M_p, _M_amp[0]);
      
      //----- phase space weight -----
      _M_amp[0] *= _M_weight;
      _M_calc = true;
    }

    if(_M_pdf_als) _M_amp[1] = _M_amp[0]*_M_pdf*_M_als;
  }
  
  template<typename _Weight, class _Event, class _EvenTraits> 
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_dipole() const
  {
    if(!_M_calc) { 
      //----- dipole contributions -----
      _M_dipole(_M_proc, _M_p, _M_amp[0]);
      
      //----- phase space weight -----
      _M_amp[0] *= _M_weight;
      _M_calc = true;
    }
    
    if(_M_pdf_als) _M_amp[1] = _M_amp[0]*_M_pdf*_M_als;
  }

  template<typename _Weight, class _Event, class _EvenTraits> 
  void amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::amp_fini() const
  {
    if(!_M_calc) { 
      //----- calculate the matrix element squares -----
      _M_fini(_M_proc, _M_p);
      
      //----- phase space weight -----
      for(unsigned int i = 0; i < 7; i++)
	_M_fini.amp[i] *= _M_weight;
      _M_calc = true;
    }
    
    switch(_M_fini.mode) {
    case 0: _M_fini.finix1(_M_amp[0]); break;
    case 1: _M_fini.finix2(_M_amp[0]); break;
    case 2: _M_fini.fini1(_M_amp[0]);  break;
    }
    
    if(_M_pdf_als) _M_amp[1] = _M_amp[0]*_M_pdf*_M_als;
  }
  
  template<typename _Weight, class _Event, class _EvenTraits> 
  const _Weight& amplitude<process<_Weight, _Event, _EvenTraits, 2U, 0U> >::
  operator()(double mr2, double mf2) const
  {
    update_scales(mr2, mf2);
    
    switch(_M_contr) {
    case born: this -> amp_born();   break;
    case real: this -> amp_real();   break;
    case sub:  this -> amp_dipole(); break;  
    case fini: this -> amp_fini();   break;
    default: throw "amplitude<..., 2, 0> : no contribution defined"; break;
    }
    
    return (_M_pdf_als ? _M_amp[1] : _M_amp[0]);
  }
}    //   namespace nlo

#endif

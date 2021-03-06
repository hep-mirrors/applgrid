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
#ifndef __NLO_NLO_INTEGRAL_I1F0_H__
#define __NLO_NLO_INTEGRAL_I1F0_H__ 1

//   standard c++ includes
#include <cmath>
#include <iostream>

//   nlojet++ includes
#include <bits/nlo-phasespace.h>
#include <bits/nlo-jetfunc.h>
#include <bits/nlo-integral.h>
#include <bits/nlo-process_i1f0.h>


namespace nlo {

  template<typename _Weight, class _Event, class _EvenTraits>
  class integral<jetfunc<amplitude<process<_Weight, _Event, _EvenTraits, 1U, 0U> > > >
  {
  public:
    //    types
    typedef process<_Weight, _Event, _EvenTraits, 1U, 0U> process_type;
    typedef amplitude<process_type> amplitude_type;
    typedef jetfunc<amplitude_type> jetfunc_type;
    typedef phasespace<_Event> phasespace_type;

    //   constructors
    integral(phasespace_type *ps, process_type *proc, jetfunc_type *jet)
      : _M_ps(ps), _M_proc(proc), _M_jet(jet) {}
        
    //   calculate the full, the born and the nlo contributions 
    void calculate_born(unsigned long);
    void calculate(unsigned long, unsigned int = 1U);
    void calculate_nlo(unsigned long, unsigned int = 1U);
    
  protected:
    //   data members
    phasespace_type *_M_ps;
    process_type    *_M_proc;
    jetfunc_type    *_M_jet;
  };


  template<typename _Weight, class _Event, class _EvenTraits>
  void integral<jetfunc<amplitude<process<_Weight, _Event, _EvenTraits, 1U, 0U> > > >::
  calculate_born(unsigned long mxne) 
  {
    unsigned short nj = _M_proc -> npar();
    double weight = 0.0;
    bool again;
    
    _Event p(nj);
    amplitude_type amp(_M_proc, amplitude_type::lo);
    
    //----- initialize the jet function -----
    _M_jet -> initfunc(nj);
    
    std::cout<<"\nStarting the calculation"<<std::endl;
    for(unsigned long int iter = 0; iter < mxne; iter++) {
      //----- generate the event-----
      do{
	again = false;
	try {
	  weight = _M_ps -> operator()(p);
	} catch(numeric_error) {
	  again = true;
	  _M_jet -> end_of_event();
	} catch(fp_exception) {
	  again = true;
	  _M_jet -> end_of_event();
	} catch(...) {
	  std::cerr<<"integral<..., 1, 0> : It was an unexpected error\n"
		   <<"while the phase space was being generated."<<std::endl;
	  throw;
	}
      } while(again);
      
      //----- flux factor & average by the initial state helicity states -----
      weight /= 16.0*(p[-1]*p[0]);
      
      //----- born term -----
      amp.leading_order(weight, p);
      _M_jet -> userfunc(p, amp);
      _M_jet -> end_of_event();
    }
  }

  template<typename _Weight, class _Event, class _EvenTraits>
  void integral<jetfunc<amplitude<process<_Weight, _Event, _EvenTraits, 1U, 0U> > > >::
  calculate_nlo(unsigned long int mxne, unsigned int time_rate_rf) 
  {
    unsigned int nj = _M_proc -> npar();
    unsigned int np = nj+1, time_rf = 0U;
    int ii, i, j, k;

    double weight1 = 0.0, weight2 = 0.0;
    double eta, x, xjac;
    bool again;

    _Event p1(nj), p2(np), p3(np);
    amplitude_type amp(_M_proc, amplitude_type::nlo);
    const random_generator *rng = _M_ps -> random();

    //----- initialize the jet function -----
    _M_jet -> initfunc(nj);

    //  the time rate parameter must be greater than zero 
    //  and should be less than 10
    if(time_rate_rf == 0U) time_rate_rf = 1U;
    if(time_rate_rf > 10U) time_rate_rf = 10U;
    
    std::cout<<"\nStarting the calculation"<<std::endl;
    for(unsigned long int iter = 1; iter <= mxne; iter++) {
      //----- generate the phase space -----
      do{
	again = false;
	try {
	  weight1 = _M_ps -> operator()(p1);
	  weight2 = _M_ps -> operator()(p1, p2);
	} catch(numeric_error) {
	  again = true;
	  _M_jet -> end_of_event();
	} catch(fp_exception) {
	  again = true;
	  _M_jet -> end_of_event();
	} catch(...) {
	  std::cerr<<"integral<..., 1, 0> : It was an unexpected error\n"
		   <<"while the phase space was being generated."<<std::endl;
	  throw;
	}
      } while(again);
      
      weight2 *= weight1;
      weight1 /= 16.0*(p1[-1]*p1[0]);
      weight2 /= 16.0*(p2[-1]*p2[0]);

      if(++time_rf == time_rate_rf) {
	time_rf = 0U;
	//----- finite terms -----
	eta = (p1[-1]*p1[0])/(p1[-1]*p1[hadron(0)]);
	xjac = -std::log(eta);
	x = eta*std::exp((rng->operator()())*xjac);
	
	p3[hadron(0)] = p1[hadron(0)];
	for(ii = p1.lower(); ii < (int) np; ii++) p3[ii] = p1[ii];
	p3[0] /= x; p3[np] = (1.0-x)*p3[0];
	
	amp.next_to_leading_order_fini(x, x*xjac, time_rate_rf*weight1, p1);
	
	amp.next_to_leading_order_fini1();
	_M_jet -> userfunc(p1, amp);
	
	amp.next_to_leading_order_finix();
	_M_jet -> userfunc(p3, amp);
      }
      
      //----- real term -----
      amp.next_to_leading_order_real(weight2, p2);
      _M_jet -> userfunc(p2, amp);
      
      //----- substraction terms -----
      for(i = 0; i < (int) np; i++)
	for(j = i + 1; j <= (int) np; j++)
	  for(k = 0; k <= (int) np; k++)
	    if(k != j && k != i)
	      if(_M_proc -> dipole_mom(p2, i, j, k, p1)) {
		amp.next_to_leading_order_sub(-weight2, p2, p1, i, j, k);
		
		if(k == 0 || i == 0) {
		  p3[hadron(0)] = p1[hadron(0)];
		  for(ii = p1.lower(); ii < (int) np; ii++) p3[ii] = p1[ii];
		  p3[0] = p2[0]; p3[np] = p2[0] - p1[0];
		  
		  _M_jet -> userfunc(p3, amp);
		} else _M_jet -> userfunc(p1, amp);
	      }
      
      //----- end of event ----- 
      _M_jet -> end_of_event();  
    }
  }


  template<typename _Weight, class _Event, class _EvenTraits>
  void integral<jetfunc<amplitude<process<_Weight, _Event, _EvenTraits, 1U, 0U> > > >::
  calculate(unsigned long int mxne, unsigned int time_rate_rf) 
  {
    unsigned int nj = _M_proc -> npar();
    unsigned int np = nj+1, time_rf = 0U;
    int ii, i, j, k;

    double weight1 = 0.0, weight2 = 0.0;
    double eta, x, xjac;
    bool again;

    _Event p1(nj), p2(np), p3(np);
    amplitude_type amp(_M_proc, amplitude_type::nlo);
    const random_generator *rng = _M_ps -> random();

    //----- initialize the jet function -----
    _M_jet -> initfunc(nj);

    //  the time rate parameter must be greater than zero 
    //  and should be less than 10
    if(time_rate_rf == 0U) time_rate_rf = 1U;
    if(time_rate_rf > 10U) time_rate_rf = 10U;
    
    std::cout<<"\nStarting the calculation"<<std::endl;
    for(unsigned long int iter = 1; iter <= mxne; iter++) {
      //----- generate the phase space -----
      do{
	again = false;
	try {
	  weight1 = _M_ps -> operator()(p1);
	  weight2 = _M_ps -> operator()(p1, p2);
	} catch(numeric_error) {
	  again = true;
	  _M_jet -> end_of_event();
	} catch(fp_exception) {
	  again = true;
	  _M_jet -> end_of_event();
	} catch(...) {
	  std::cerr<<"integral<..., 1, 0> : It was an unexpected error\n"
		   <<"while the phase space was being generated."<<std::endl;
	  throw;
	}
      } while(again);
      
      weight2 *= weight1;
      weight1 /= 16.0*(p1[-1]*p1[0]);
      weight2 /= 16.0*(p2[-1]*p2[0]);

      if(++time_rf == time_rate_rf) {
	time_rf = 0U;

	//----- born term -----
	amp.leading_order(time_rate_rf*weight1, p1);
	_M_jet -> userfunc(p1, amp);

	//----- finite terms -----
	eta = (p1[-1]*p1[0])/(p1[-1]*p1[hadron(0)]);
	xjac = -std::log(eta);
	x = eta*std::exp((rng->operator()())*xjac);
	
	p3[hadron(0)] = p1[hadron(0)];
	for(ii = p1.lower(); ii < (int) np; ii++) p3[ii] = p1[ii];
	p3[0] /= x; p3[np] = (1.0-x)*p3[0];
	
	amp.next_to_leading_order_fini(x, x*xjac, time_rate_rf*weight1, p1);
	
	amp.next_to_leading_order_fini1();
	_M_jet -> userfunc(p1, amp);
	
	amp.next_to_leading_order_finix();
	_M_jet -> userfunc(p3, amp);
      }
      
      //----- real term -----
      amp.next_to_leading_order_real(weight2, p2);
      _M_jet -> userfunc(p2, amp);
      
      //----- substraction terms -----
      for(i = 0; i < (int) np; i++)
	for(j = i + 1; j <= (int) np; j++)
	  for(k = 0; k <= (int) np; k++)
	    if(k != j && k != i)
	      if(_M_proc -> dipole_mom(p2, i, j, k, p1)) {
		amp.next_to_leading_order_sub(-weight2, p2, p1, i, j, k);
		
		if(k == 0 || i == 0) {
		  p3[hadron(0)] = p1[hadron(0)];
		  for(ii = p1.lower(); ii < (int) np; ii++) p3[ii] = p1[ii];
		  p3[0] = p2[0]; p3[np] = p2[0] - p1[0];
		  
		  _M_jet -> userfunc(p3, amp);
		} else _M_jet -> userfunc(p1, amp);
	      }
      
      //----- end of event ----- 
      _M_jet -> end_of_event();  
    }
  }

  
}    //   namespace nlo

#endif

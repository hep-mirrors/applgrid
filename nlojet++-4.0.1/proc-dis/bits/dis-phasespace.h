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
#ifndef __NLO_DIS_PHASESPACE_H__
#define __NLO_DIS_PHASESPACE_H__ 1

#include <bits/dis-event.h>
#include <bits/nlo-phasespace.h>
#include <bits/psg-phasespace_n1i1f0.h>
 

namespace nlo {

  //   Shorthand notations
  typedef phasespace<event_dis> phasespace_dis;

  
  class basic_phasespace_dis
	: public basic_phasespace<event_dis>
  {
	//   private types
	typedef phasespace_generator<lorentzvector<double>, event_dis::iterator> _PhaseSpaceGen;

  public:
	typedef event_dis event_type;
	
    //  constructor
    basic_phasespace_dis(const random_generator *rng, double el, double eh)
      : basic_phasespace<event_dis>(rng, el, eh) {}
    
    basic_phasespace_dis(const random_generator *rng, _PhaseSpaceGen *psg, double el, double eh)
      : basic_phasespace<event_dis>(rng, psg, el, eh) {}
    	
	//   set the phase space cut parameters
	void phasespace_cuts(double, double, double, double, double, double);
	
    double operator()(event_type&);
    	
  private:
	//    data members
	//   phase space cut parameters
	double _M_q2min, _M_q2max, _M_xmin, _M_xmax, _M_ymin, _M_ymax;
  };
}

#endif

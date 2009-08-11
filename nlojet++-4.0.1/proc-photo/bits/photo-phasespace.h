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
#ifndef __NLO_PHOTO_PHASESPACE_H__
#define __NLO_PHOTO_PHASESPACE_H__ 1

#include <bits/photo-event.h>
#include <bits/nlo-phasespace.h>
#include <bits/psg-phasespace_n0i2f0.h>


namespace nlo {

  //   Shorthand notations
  typedef phasespace<event_photo> phasespace_photo;

  class basic_phasespace_photo 
    : public basic_phasespace<event_photo>
  {
    //  private types
    typedef basic_phasespace<event_photo> _Base;  
    typedef event_photo::iterator _Iterator;
    typedef phasespace_generator<lorentzvector<double>, _Iterator> _PhaseSpaceGen;
    
  public:
    typedef _Base::event_type event_type;
    
    //  constructor
    basic_phasespace_photo(const random_generator *rng, _PhaseSpaceGen *psg, 
			   double el, double eh, double beta = 2.0, double eps = 1.0) 
      : _Base(rng, psg, 4.0*el*eh, beta, eps), _M_el(el), _M_eh(eh) {}
     

    basic_phasespace_photo(const random_generator *rng, double el, double eh,
			   double beta = 2.0, double eps = 1.0)
      : _Base(rng, 4.0*el*eh, beta, eps), _M_el(el), _M_eh(eh) {}
     
    //  generate the phase space
    double operator()(event_type&);
     
  private:
    double _M_el, _M_eh;
  };
}

#endif

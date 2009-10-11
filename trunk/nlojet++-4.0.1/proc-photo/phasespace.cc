//  Copyright (C) 2003 Zoltan Nagy
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


//   nlo includes
#include "bits/photo-phasespace.h"


namespace nlo {
  
  double basic_phasespace_photo::operator()(event_type& p) 
  {
    double ret_val = this -> _Base::operator()(p);
    
    //----- incoming hadron -----
    p[hadron(-1)] = lorentzvector<double>(0.0, 0.0,  _M_el, _M_el);  
    p[hadron( 0)] = lorentzvector<double>(0.0, 0.0, -_M_eh, _M_eh);
    
    //----- boost the event into the labor frame -----
    double bz = (_M_el-_M_eh)/(_M_eh+_M_el);
    for(event_type::iterator pi = p.begin(); pi < p.end(); pi++) 
      pi -> boost(0, 0, bz);

    return ret_val;
  }
}

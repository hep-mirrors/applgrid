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


//  nlojet++ includes
#include "statistic.h"



namespace nlo {
  
  template<> std::ostream& 
  print(std::ostream& os, const sample_statistic<double, sample_traits<double> >& ss) 
  {
    return os<<ss.mean()<<" +/- "<<ss.error()<<std::endl;
  }
  
  template<> std::ostream& 
  print(std::ostream& os, const sample_statistic<float, sample_traits<float> >& ss) 
  {
    return os<<ss.mean()<<" +/- "<<ss.error()<<std::endl;
  }
  
  template<> std::ostream& 
  print(std::ostream& os, const sample_statistic<long double, sample_traits<long double> >& ss) 
  {
    return os<<ss.mean()<<" +/- "<<ss.error()<<std::endl;
  }


}

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
#ifndef __STATISTIC_H__
#define __STATISTIC_H__ 1

#include <bits/hep-sample_statistic.h>

//   Some specialization
#include <bits/nlo-weight.h>

namespace nlo {
  
  template<> std::ostream& 
  print(std::ostream&, const sample_statistic<double, sample_traits<double> >&); 

  template<> std::ostream& 
  print(std::ostream&, const sample_statistic<float, sample_traits<float> >&); 
  
  template<> std::ostream& 
  print(std::ostream&, const sample_statistic<long double, sample_traits<long double> >&);

  
  template<unsigned int _Size, const char **_Label> std::ostream& 
  print(std::ostream& os, const sample_statistic<weight<_Size, _Label>, sample_traits<weight<_Size, _Label> > >& ss) 
  {
    weight<_Size, _Label> mean(ss.mean()), var(ss.error());   
    for(unsigned int i = 0; i < _Size; i++)
      os<<mean.label(i)<<" : "<<mean[i]<<" +/- "<<var[i]<<"\n";
    return os<<std::endl;
  }
}
#endif

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
#ifndef __DISTRIBUTION_H__
#define __DISTRIBUTION_H__ 1

//   Standard C++ headers
#include <iostream>


//  nlojet++ headers
#include <statistic.h>
#include <bits/hep-dist.h>
#include <bits/hep-distbook.h>




namespace nlo {

  
  template<unsigned int _Size, const char **_Label> std::ostream& 
  print(std::ostream& os, const distribution<weight<_Size, _Label>, void, sample_traits<weight<_Size, _Label> > >& d) 
  {    
    typedef distribution<weight<_Size, _Label>, void, sample_traits<weight<_Size, _Label> > > _Dist;
	return print(os, dynamic_cast<const typename _Dist::statistic_type&>(d));
  }
  

  template<unsigned int _Size, typename _Point, const char **_Label> std::ostream& 
  print(std::ostream& os, const distribution<weight<_Size, _Label>, _Point, sample_traits<weight<_Size, _Label> > >& d) 
  {    
    typedef distribution<weight<_Size, _Label>, _Point, sample_traits<weight<_Size, _Label> > > _Dist;
	
    os<<"#";
    for(unsigned int l = 0; l < _Size; l++) 
      os<<"   "<<_Label[l]; 
    os<<std::endl;
    
    for(typename _Dist::const_iterator i = d.begin(); i < d.end(); i++) {
      _Dist::point_traits::print(os, i->first); os<<"    ";
      _Dist::sample_traits::print(os, i->second.mean(), i->second.error());
      os<<"\n";
    }
    
    return os;
  }



  
  //    a simple point type for the 1D histograms
  struct histpoint1d {
    double xmin, xmid, xmax;
  };
  

  template<>
  struct distpoint_traits<histpoint1d>
  {
    typedef histpoint1d distpoint_type;
    
    static void write_txt(std::ostream& os, const distpoint_type& p) {
      os<<p.xmin<<"  "<<p.xmid<<"  "<<p.xmax;
    }
    
    static void write_bin(std::ostream& os, const distpoint_type& p) {
      os.write((const char *) &p, sizeof(distpoint_type));
    }
    
    static void read_txt(std::istream& is, distpoint_type& p) {
      is>>p.xmin>>p.xmid>>p.xmax;
    }
     
    static void read_bin(std::istream& is, distpoint_type& p) {
      is.read((char *) &p, sizeof(distpoint_type));
    }
   
    static void print(std::ostream& os,  const distpoint_type& p) {
      if(p.xmin>= 0) os<<" "; os<<p.xmin<<"  ";
      if(p.xmid>= 0) os<<" "; os<<p.xmid<<"  ";
      if(p.xmax>= 0) os<<" "; os<<p.xmax;
    }
  };



  //    a simple point type for the 2D histograms
  struct histpoint2d {
    double xmin, xmid, xmax;
    double ymin, ymid, ymax;
  };

  template<>
  struct distpoint_traits<histpoint2d>
  {
    typedef histpoint2d distpoint_type;
    
    static void write_txt(std::ostream& os, const distpoint_type& p) {
      os<<p.xmin<<"  "<<p.xmid<<"  "<<p.xmax<<"  "
	<<p.ymin<<"  "<<p.ymid<<"  "<<p.ymax;
    }
     
    static void write_bin(std::ostream& os, const distpoint_type& p) {
      os.write((const char *) &p, sizeof(distpoint_type));
    }
   
    static void read_txt(std::istream& is, distpoint_type& p) {
      is>>p.xmin>>p.xmid>>p.xmax>>p.ymin>>p.ymid>>p.ymax;
    }
     
    static void read_bin(std::istream& is, distpoint_type& p) {
      is.read((char *) &p, sizeof(distpoint_type));
    }
   
    static void print(std::ostream& os,  const distpoint_type& p) {
      if(p.xmin>= 0) os<<" "; os<<p.xmin<<"  ";
      if(p.xmid>= 0) os<<" "; os<<p.xmid<<"  ";
      if(p.xmax>= 0) os<<" "; os<<p.xmax<<"  ";
      if(p.ymin>= 0) os<<" "; os<<p.ymin<<"  ";
      if(p.ymid>= 0) os<<" "; os<<p.ymid<<"  ";
      if(p.ymax>= 0) os<<" "; os<<p.ymax;
    }
  };


  
}


#endif

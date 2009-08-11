//
// //  Copyright (C) 2002 Zoltan Nagy
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
#ifndef __NLO_HEP_DISTPOINT_TRAITS_H__
#define __NLO_HEP_DISTPOINT_TRAITS_H__ 1

#include <iostream>


namespace nlo {

  template<class _Point>
  struct distpoint_traits
  {
    typedef _Point distpoint_type;
    
    static void write_txt(std::ostream& os, const distpoint_type& p) {
      os<<p;
    }
     
    static void write_bin(std::ostream& os, const distpoint_type& p) {
      os.write((const char *) &p, sizeof(distpoint_type));
    }
   
    static void read_txt(std::istream& is, distpoint_type& p) {
      is>>p;
    }
     
    static void read_bin(std::istream& is, distpoint_type& p) {
      is.read((char *) &p, sizeof(distpoint_type));
    }
   
    static void print(std::ostream& os,  const distpoint_type& p) {
      os<<p;    
    }
  };

	template<>
	struct distpoint_traits<void>
	{
		typedef void distpoint_type;
	};
		
	
  template<>
  struct distpoint_traits<float>
  {
    typedef float distpoint_type;
    
    static void write_txt(std::ostream& os, const distpoint_type& p) {
      os<<p;
    }
     
    static void write_bin(std::ostream& os, const distpoint_type& p) {
      os.write((const char *) &p, sizeof(distpoint_type));
    }
   
    static void read_txt(std::istream& is, distpoint_type& p) {
      is>>p;
    }
     
    static void read_bin(std::istream& is, distpoint_type& p) {
      is.read((char *) &p, sizeof(distpoint_type));
    }
   
    static void print(std::ostream& os,  const distpoint_type& p) {
      if(p <= 0.0f) os<<" "; os<<p;    
    }
  };

  template<>
  struct distpoint_traits<double>
  {
    typedef double distpoint_type;
    
    static void write_txt(std::ostream& os, const distpoint_type& p) {
      os<<p;
    }
     
    static void write_bin(std::ostream& os, const distpoint_type& p) {
      os.write((const char *) &p, sizeof(distpoint_type));
    }
   
    static void read_txt(std::istream& is, distpoint_type& p) {
      is>>p;
    }
     
    static void read_bin(std::istream& is, distpoint_type& p) {
      is.read((char *) &p, sizeof(distpoint_type));
    }
   
    static void print(std::ostream& os,  const distpoint_type& p) {
      if(p <= 0.0) os<<" "; os<<p;    
    }
  };

  template<>
  struct distpoint_traits<long double>
  {
    typedef long double distpoint_type;
    
    static void write_txt(std::ostream& os, const distpoint_type& p) {
      os<<p;
    }
     
    static void write_bin(std::ostream& os, const distpoint_type& p) {
      os.write((const char *) &p, sizeof(distpoint_type));
    }
   
    static void read_txt(std::istream& is, distpoint_type& p) {
      is>>p;
    }
     
    static void read_bin(std::istream& is, distpoint_type& p) {
      is.read((char *) &p, sizeof(distpoint_type));
    }
   
    static void print(std::ostream& os,  const distpoint_type& p) {
      if(p <= 0.0l) os<<" "; os<<p;    
    }
  };

}   //  namespace nlo


#endif

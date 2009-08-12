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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#include <fstream>
#include <string>
#include <list>

#include "distribution.h"
#include <bits/nlo-basic_user_result.h>


namespace nlo {


#define PRINTBOOK(I)							\
  if(r.book##I.size() > 0) {						\
    std::fstream f((dir + "/result-"#I).c_str(), std::ios::out);	\
    f.setf(std::ios::scientific, std::ios::floatfield);			\
    print(f, r.book##I);						\
    f.close();								\
  }

  template<class _Db0>
  void print(const basic_user_result<_Db0>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
  }

  template<class _Db0, class _Db1>
  void print(const basic_user_result<_Db0, _Db1>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
  }
 
  template<class _Db0, class _Db1, class _Db2>
  void print(const basic_user_result<_Db0, _Db1, _Db2>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
    PRINTBOOK(2)
  }

  template<class _Db0, class _Db1, class _Db2, class _Db3>
  void print(const basic_user_result<_Db0, _Db1, _Db2, _Db3>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
    PRINTBOOK(2)
    PRINTBOOK(3)
  }

  template<class _Db0, class _Db1, class _Db2, class _Db3, class _Db4>
  void print(const basic_user_result<_Db0, _Db1, _Db2, _Db3, _Db4>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
    PRINTBOOK(2)
    PRINTBOOK(3)
    PRINTBOOK(4)
  }

  template<class _Db0, class _Db1, class _Db2, class _Db3, class _Db4,
	   class _Db5>
  void print(const basic_user_result<_Db0, _Db1, _Db2, _Db3, _Db4, _Db5>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
    PRINTBOOK(2)
    PRINTBOOK(3)
    PRINTBOOK(4)
    PRINTBOOK(5)
  }

  template<class _Db0, class _Db1, class _Db2, class _Db3, class _Db4,
	   class _Db5, class _Db6>
  void print(const basic_user_result<_Db0, _Db1, _Db2, _Db3, _Db4, _Db5, _Db6>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
    PRINTBOOK(2)
    PRINTBOOK(3)
    PRINTBOOK(4)
    PRINTBOOK(5)
    PRINTBOOK(6)
  }

  template<class _Db0, class _Db1, class _Db2, class _Db3, class _Db4,
	   class _Db5, class _Db6, class _Db7>
  void print(const basic_user_result<_Db0, _Db1, _Db2, _Db3, _Db4, _Db5, _Db6, _Db7>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
    PRINTBOOK(2)
    PRINTBOOK(3)
    PRINTBOOK(4)
    PRINTBOOK(5)
    PRINTBOOK(6)
    PRINTBOOK(7)
  }

  template<class _Db0, class _Db1, class _Db2, class _Db3, class _Db4,
	   class _Db5, class _Db6, class _Db7, class _Db8>
  void print(const basic_user_result<_Db0, _Db1, _Db2, _Db3, _Db4, _Db5, _Db6, _Db7, _Db8>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
    PRINTBOOK(2)
    PRINTBOOK(3)
    PRINTBOOK(4)
    PRINTBOOK(5)
    PRINTBOOK(6)
    PRINTBOOK(7)
    PRINTBOOK(8)
  }

  template<class _Db0, class _Db1, class _Db2, class _Db3, class _Db4,
	   class _Db5, class _Db6, class _Db7, class _Db8, class _Db9>
  void print(const basic_user_result<_Db0, _Db1, _Db2, _Db3, _Db4, _Db5, _Db6, _Db7, _Db8, _Db9>& r, 
	     const std::basic_string<char>& dir) 
  {
    PRINTBOOK(0)
    PRINTBOOK(1)
    PRINTBOOK(2)
    PRINTBOOK(3)
    PRINTBOOK(4)
    PRINTBOOK(5)
    PRINTBOOK(6)
    PRINTBOOK(7)
    PRINTBOOK(8)
    PRINTBOOK(9)
  }


  
  template<class _Result> unsigned long int 
  main_module_add(bool weighted, const std::list<std::basic_string<char> >& names, 
				  const std::basic_string<char>& dir)
  {
    typedef std::list<std::basic_string<char> >::const_iterator iterator;
    
    std::fstream file;
    unsigned long int saved, nsaved = 0UL;
    _Result tmp, *res = 0;
    iterator in = names.begin(), end = names.end();
    
    while(in != end) {
      //----- open the file ----
      file.open(in -> c_str(), std::ios::in);
      
      //----- read the file if it is exist -----
      if(file.is_open()) { 
		nsaved += (saved = read_basic_user_output(file, tmp)); 
		file.flush();
	
		if(saved != 0) {
		  if(!res) res = new _Result(tmp);
		  else {
			if(weighted) res-> assadd_weighted(tmp);
			else res -> operator+=(tmp);
		  }
		}
      }
      
      //----- close the file ----
      file.close(); ++in;
    }
    
    if(!res) return 0; 

    //----- save the results into the directory dir 
    print(*res, dir);
    
    return nsaved;
  } 

    

}






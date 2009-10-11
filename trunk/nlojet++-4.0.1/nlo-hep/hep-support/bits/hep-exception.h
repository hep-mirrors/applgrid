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
#ifndef __NLO_NLO_EXCEPTION_H__
#define __NLO_NLO_EXCEPTION_H__ 1

#include <exception>
#include <csignal>


namespace nlo {
 
  class monte_carlo_bad_event : public std::exception
  {
  public:
    virtual const char* what () const throw() { 
      return "Generated Monte Carlo event arised an exception";
    }
  };
}


namespace nlo {


  class numeric_error : public std::exception
  {
  public:
    virtual const char* what () const throw() { 
      return "Numeric error";
    }
  };


  class fp_exception : public std::exception
  {
  public:
    enum exception_type {inexact = 0, divbyzero, overflow, underflow, invalid};

    //   constructor
    explicit fp_exception(exception_type except) 
      : _M_excpt(except) {}
    
    // type of the exception
    exception_type type() const { return _M_excpt;}
    
    virtual const char* what () const throw() { 
      return "Floating pont exception";
    }
    
  private:
    exception_type _M_excpt;
  };

  class fp_inexact : public fp_exception
  {
  public:
    fp_inexact() : fp_exception(inexact) {}
    
    virtual const char* what () const throw() { 
      return "Floating pont exception : Inexact";
    }
  };

  class fp_divbyzero : public fp_exception
  {
  public:
    fp_divbyzero() : fp_exception(divbyzero) {}
    
    virtual const char* what () const throw() { 
      return "Floating pont exception : Divbyzero";
    }
  };

  class fp_underflow : public fp_exception
  {
  public:
    fp_underflow() : fp_exception(underflow) {}
    
    virtual const char* what () const throw() { 
      return "Floating pont exception : Underflow";
    }
  };

  class fp_overflow : public fp_exception
  {
  public:
    fp_overflow() : fp_exception(overflow) {}
    
    virtual const char* what () const throw() { 
      return "Floating pont exception : Overflow";
    }
  };

  class fp_invalid : public fp_exception
  {
  public:
    fp_invalid() : fp_exception(invalid) {}
    
    virtual const char* what () const throw() { 
      return "Floating pont exception : Invalid";
    }
  };
}


//   Floating point exception handling (GNU libc only)
#include <fenv.h>

namespace nlo {

  inline void __fe_clear_exception(int excpt = FE_ALL_EXCEPT) {
    feclearexcept(FE_ALL_EXCEPT);
  }
  
  inline void __fe_throw_exception()
  {
    int raised = fetestexcept(FE_OVERFLOW | FE_UNDERFLOW | FE_DIVBYZERO | FE_INVALID);
    if(raised & FE_OVERFLOW)  throw fp_overflow();
    if(raised & FE_UNDERFLOW) throw fp_underflow();
    if(raised & FE_DIVBYZERO) throw fp_divbyzero();
    if(raised & FE_INVALID)   throw fp_invalid();
  }
  
  inline void __fe_throw_exception(int except)
  {
    int raised = fetestexcept(except);
    if(raised & FE_OVERFLOW)  throw fp_overflow();
    if(raised & FE_UNDERFLOW) throw fp_underflow();
    if(raised & FE_DIVBYZERO) throw fp_divbyzero();
    if(raised & FE_INVALID)   throw fp_invalid();
    if(raised & FE_INEXACT)   throw fp_inexact();
  }
}

#endif

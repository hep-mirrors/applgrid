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
#ifndef __NLO_PSG_RAMBO__
#define __NLO_PSG_RAMBO__ 1

// nlo++ includes
#include <bits/psg-psgen.h>
#include <bits/hep-lorentzvector.h>


namespace nlo {
  
  /// \class rambo
  ///  Phase Space algorithm : RAMBO
  ///    R. Kleiss, W.J. Stirling and S.D. Ellis
  ///     Comp. Phys. Comm. 40 (1986) 359-373
  template<typename _OutputIterator>
  class rambo
    : public phasespace_generator<lorentzvector<double>, _OutputIterator>
  {
    //  private types
    typedef phasespace_generator<lorentzvector<double>, _OutputIterator> _Base;

   public:
    /// This is an alias to the Lorentz vector type.
    typedef typename _Base::lorentzvector_type lorentzvector_type;
    
    /// This is an alias to the output iterator type. 
    typedef typename _Base::iterator iterator;
    
    /// \fn double operator()(double s, iterator first, iterator last);
    /// \brief This function generates a phase space at the given center of mass energy.
    ///  
    /// \param s is the center of mass energy square.
    /// \param first is an iterator which is point to the first element of 
    /// that array (or any container) which stores the momenta.
    /// \param last is an iterator which is the past-the-end value for
    /// the container.
    /// \return The pahse space weight of the event.
    double operator()(double, iterator, iterator);
    
    /// This constructor initialize the object with the random nuber generator.
    rambo(const random_generator *__rng) : _Base(__rng) {}
  };

  
  template<typename _OutputIterator> double rambo<_OutputIterator>::
  operator()(double s, iterator first, iterator last) 
  {
    extern lorentzvector<double> __rambo_helper_random_momentum(const random_generator *);
    extern double __rambo_helper_weight(unsigned int, double);
    
    lorentzvector<double> psum;
    unsigned int n = 0U;
        
    iterator iter = first;
    while(iter != last) {
      psum += (*iter = __rambo_helper_random_momentum(_Base::_M_rng));
      n++; iter++;
    }
    
    //----- parameters of the conform transformation -----
    double x = std::sqrt(s)/psum.mag();
    threevector<double> bVec = -psum.boostVector();
    
    //----- do the conform transformation -----
    iter = first;
    while(iter != last) {
      iter -> boost(bVec);
      iter -> operator*=(x);
      iter++;
    }
    
    return __rambo_helper_weight(n, s);
  }
 
}

#endif

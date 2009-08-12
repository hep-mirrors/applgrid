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
#ifndef __NLO_PSG_PSGEN__
#define __NLO_PSG_PSGEN__ 1

// nlo++ includes
#include <bits/hep-rng.h>


namespace nlo {
  
  /// \brief Abstract class for defining the generally used phase 
  /// space algorithm.
  /// 
  /// The general idea is that building the process dependent phase 
  /// space generators from simple algorithm like RAMBO, SARGE. These 
  /// algorithms usaualy generate a phase space configuration with a given 
  /// number of outgoing monenta at fix center-of-mass energy. This 
  /// abstract type provides a possible protocol to define these algorithm 
  /// and use them in a same way.
  /// \param _Lv is a Lorentz vector type template parameter.
  /// \param _OutputIterator is an output iterator type template parameter 
  /// which can iterate (point to)  the \b _Lv Lorent vector type variables.
  /// To get more information about the iterators see the Standard C++ 
  /// reference manual.
  /// \see International Standard, Programming Language -- C++, ISO/IEC 14882 \n
  /// lorentzvector
  template<class _Lv, typename _OutputIterator> 
  class phasespace_generator
  {
   public:
	/// Desctructor
	virtual ~phasespace_generator() {}
	
    /// This is an alias to the Lorentz vector type parameter.
    typedef _Lv lorentzvector_type;
    
    /// This is an alias to the output iterator type parameter. 
    typedef _OutputIterator iterator;
    
    /// \fn virtual double operator()(double s, iterator first, iterator last) = 0;
    /// \brief This function generates a phase space at the given center of mass energy.
    ///  
    /// \param s is the center of mass energy square.
    /// \param first is an iterator which is point to the first element of 
    /// that array (or any container) which stores the momenta.
    /// \param last is an iterator which is the past-the-end value for
    /// the container.
    /// \return The pahse space weight of the event.
    virtual double operator()(double, iterator, iterator) = 0;
    
  protected:
    /// This constructor initialize the object with the random nuber generator.
    phasespace_generator(const random_generator *__rng = 0) 
      : _M_rng(__rng) {}
    
    ///  const pointer to the random number generator
    const random_generator *_M_rng;
  };
}

#endif

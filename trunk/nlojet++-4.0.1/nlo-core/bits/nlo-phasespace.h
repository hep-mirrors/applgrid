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
#ifndef __NLO_NLO_PHASESPACE_H__
#define __NLO_NLO_PHASESPACE_H__ 1

//   nlo includes
#include <bits/hep-rng.h>


namespace nlo {
  
  
  /// \briefAbstract class for generating the process dependent events.
  ///  
  /// This class tries to implement generaly the \f$m\f$ and 
  /// \f$m+1\f$ particle phase space integral. Every phase space generator should be inhereted from this class.
  /// \param _Event is the event type template parameter
  template<class _Event> 
  struct phasespace
  {
    /// Alias for the template parameter.
    /// This parameter must be an event type parameter.
    /// \see partonic_event, hadronic_event
    typedef _Event event_type;
    
    ///   destructor
    virtual ~phasespace() {}
    
    ///  \fn virtual double phasespace::operator()(event_type& p) = 0;
    ///  \brief Abstract function for generating the phase space.
    ///
    ///  This function should generate random event for that process 
    ///  which is defined by the corrent event_type variable \b p. To 
    ///  generate the random momenta it should use the \b _M_rng random number generator.
    ///  \param p is a reference to event_type. This is a return parameter.
    ///  \return Returns with a double. It gives the weight of the phase space point.
    virtual double operator()(event_type& p) = 0;

    ///  \fn virtual double phasespace::operator()(const event_type& q, event_type& p) = 0;
    ///  \brief Abstract function for generating an \f$m+n\f$ particle event 
    ///  from an \f$m\f$ particle event.
    ///
    ///  This function should implement any algorithm which can generate an 
    ///  \f$m+n\f$ particle event from an \f$m\f$ particle event. In this 
    ///  program we will use the dipole emission method to do this.  
    ///  Usually it has only two parameter (the input event and the output event) 
    ///  but sometimes the implemented algorithm requires some external parameters 
    ///  and it is usefull if we can pass them through the virtual function calling. 
    ///  \param q is a const reference to the \f$m\f$ particle event.
    ///  \param p is a  reference to the \f$m+n\f$ particle event.
    ///  \return Returns with a double. It gives the weight of this one 
    ///  particle subphase space point.
    virtual double operator()(const event_type&, event_type&) = 0;
        
    /// This function returns the const reference to the random number generator.
    virtual const random_generator * random() const { return _M_rng;}
    
  protected:
    ///  \brief Protected constructor
    ///
    ///  This constructor initialize the class with a given random number generator.
    ///  \param __rng is a const reference to the random number generator
    ///  (\ref random_generator)
    phasespace(const random_generator *__rng = 0) 
      : _M_rng(__rng) {}
    
    ///  constatnt ponter to random number generator
    const random_generator *_M_rng;
  };

  

}   // namespace nlo

#endif

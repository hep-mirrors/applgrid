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
#ifndef __NLO_NLO_JETFUNC_H__
#define __NLO_NLO_JETFUNC_H__ 1


namespace nlo {

  //   abstarct class : jetfunc
  //
  //     calculates the jet function
  //     and fills the histograms
  //
  template<class _Amp>
  class jetfunc
  {
  public:
    //   types
    typedef _Amp amp_type;
    typedef typename _Amp::event_type event_type;
    
    //   destructor
    virtual ~jetfunc() {}
    
    //   initialize the histograms and other storage variables
    virtual void initfunc(unsigned int) = 0;
    
    //   user function
    virtual void userfunc(const event_type&, const _Amp&) = 0;
    
    //   end of event
    virtual void end_of_event() = 0;
  };
}   //  namespace nlo

#endif

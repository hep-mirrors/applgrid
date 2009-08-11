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
#ifndef __NLO_PROCESS_H__
#define __NLO_PROCESS_H__ 1


namespace nlo {

  //  Empty declaration of the template class process.
  template<typename _Weight, class _Event, 
	   class _EventTraits = typename _Event::traits_type,
	   unsigned int _Inh = _EventTraits::incoming_hadrons,
	   unsigned int _Idh = _EventTraits::identified_hadrons> 
  class process;
  
  //  Empty declaration of the template class amplitude.
  template<class _Process> class amplitude;

  //  Empty declaration of the pdf and coupling functions
  template<class _Weight, unsigned int _Inh, unsigned int _Idh>
  struct pdf_and_coupling;
}  //  namespace nlo


#endif

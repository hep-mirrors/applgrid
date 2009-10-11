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
#ifndef __NLO_HHC2PH_PROCESS_H__
#define __NLO_HHC2PH_PROCESS_H__ 1

#include <bits/hhc2ph-event.h>
#include <bits/hhc2ph-weight.h>
#include <bits/nlo-process_i2f0.h>
 

namespace nlo {

   //  Short name declarations (hadron-hadron with 2 photons)
  typedef process<weight_hhc2ph, event_hhc2ph> process_hhc2ph;
  typedef amplitude<process_hhc2ph> amplitude_hhc2ph;

  class pdf_and_coupling_hhc2ph
    : public pdf_and_coupling<weight_hhc2ph,2U,0U>
  {
    virtual void hadronA(double, double, unsigned int, unsigned int, double *) = 0;
    virtual void hadronB(double, double, unsigned int, unsigned int, double *) = 0;
    weight_hhc2ph pdf(double, double, double, unsigned int=2U, unsigned int=3U);
  };
}

#endif

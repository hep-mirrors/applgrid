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
#ifndef __NLO_PHOTO_PROCESS_H__
#define __NLO_PHOTO_PROCESS_H__ 1

#include <bits/photo-event.h>
#include <bits/photo-weight.h>
#include <bits/nlo-process_i2f0.h>
 

namespace nlo {

  //   Shorthand notations
  typedef process<weight_photo, event_photo> process_photo;
  typedef amplitude<process_photo> amplitude_photo;
 
  class pdf_and_coupling_photo
    : public pdf_and_coupling<weight_photo,2U,0U>
  {
    virtual double photon(double) = 0;
    virtual void hadron(double, double, unsigned int, unsigned int, double *) = 0;
    weight_photo pdf(double, double, double, unsigned int=2U, unsigned int=3U);
  };
}

#endif

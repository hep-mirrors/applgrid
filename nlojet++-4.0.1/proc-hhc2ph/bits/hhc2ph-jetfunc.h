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
#ifndef __NLO_HHC2PH_JETFUNC_H__
#define __NLO_HHC2PH_JETFUNC_H__ 1

#include <bits/hhc2ph-process.h>
#include <bits/nlo-jetfunc.h>
#include <bits/nlo-basic_user.h>
#include <bits/nlo-basic_user_set.h>


namespace nlo {

  //   Shorthand notations
  typedef jetfunc<amplitude_hhc2ph> jetfunc_hhc2ph;
  typedef basic_user_base<jetfunc_hhc2ph> user_base_hhc2ph;

  typedef basic_user<jetfunc_hhc2ph, void,        double, weight_conversion> user0d_hhc2ph;
  typedef basic_user<jetfunc_hhc2ph, double,      double, weight_conversion> user1d_hhc2ph;
  typedef basic_user<jetfunc_hhc2ph, histpoint1d, double, weight_conversion> user1h_hhc2ph;
  typedef basic_user<jetfunc_hhc2ph, histpoint2d, double, weight_conversion> user2h_hhc2ph;

  typedef basic_user_set<user0d_hhc2ph, user1d_hhc2ph, user1h_hhc2ph, user2h_hhc2ph> user_set_hhc2ph;
}

#endif

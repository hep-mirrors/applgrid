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
#ifndef __NLO_DIS_JETFUNC_H__
#define __NLO_DIS_JETFUNC_H__ 1

#include <bits/dis-process.h>
#include <bits/nlo-jetfunc.h>
#include <bits/nlo-basic_user.h>
#include <bits/nlo-basic_user_set.h>


namespace nlo {

  //   Shorthand notations
  typedef jetfunc<amplitude_dis> jetfunc_dis;
  typedef basic_user_base<jetfunc_dis> user_base_dis;

  typedef basic_user<jetfunc_dis, void,        double, weight_conversion> user0d_dis;
  typedef basic_user<jetfunc_dis, double,      double, weight_conversion> user1d_dis;
  typedef basic_user<jetfunc_dis, histpoint1d, double, weight_conversion> user1h_dis;
  typedef basic_user<jetfunc_dis, histpoint2d, double, weight_conversion> user2h_dis;

  typedef basic_user_set<user0d_dis, user1d_dis, user1h_dis, user2h_dis> user_set_dis;
}

#endif

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
#ifndef __NLO_EPA_JETFUNC_H__
#define __NLO_EPA_JETFUNC_H__ 1

#include <bits/epa-process.h>
#include <bits/nlo-jetfunc.h>
#include <bits/nlo-basic_user.h>
#include <bits/nlo-basic_user_set.h>


namespace nlo {

  //   Shorthand notations
  typedef jetfunc<amplitude_epa>       jetfunc_epa;
  typedef basic_user_base<jetfunc_epa> user_base_epa;

  typedef basic_user<jetfunc_epa, void>        user0d_epa;
  typedef basic_user<jetfunc_epa, double>      user1d_epa;
  typedef basic_user<jetfunc_epa, histpoint1d> user1h_epa;
  typedef basic_user<jetfunc_epa, histpoint2d> user2h_epa;

  typedef basic_user_set<user0d_epa, user1d_epa, user1h_epa, user2h_epa> user_set_epa;
}

#endif

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
#ifndef __NLO_RANDOM_H__
#define __NLO_RANDOM_H__ 1

//  abstract class for random generators
#include <bits/hep-rng.h>

//  available generators 
#include <bits/rng-cmrg.h>
#include <bits/rng-ranlux.h>  
#include <bits/rng-ranlxs.h>
#include <bits/rng-cmgr.h>
#include <bits/rng-mt19937.h>
#include <bits/rng-ranlxd.h>


//  random distributions
#include <bits/ran-gamma.h>
#include <bits/ran-beta.h> 
#include <bits/ran-gaussian.h> 
#include <bits/ran-sphere.h>


#endif

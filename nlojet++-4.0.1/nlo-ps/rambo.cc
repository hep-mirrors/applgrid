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

// Standard C inludes
#include <cmath>

// nlo++ includes
#include "bits/hep-rng.h"
#include "bits/hep-lorentzvector.h"


namespace nlo {

  extern lorentzvector<double> __rambo_helper_random_momentum(const random_generator *);
  extern double __rambo_helper_weight(unsigned int, double);
  
  
  lorentzvector<double> __rambo_helper_random_momentum(const random_generator *x)
  {
    double E   = -std::log((x->operator()())*(x->operator()())); 
    double pz  = E*(2.0*(x->operator()()) - 1.0);
    double pt  = std::sqrt(E*E - pz*pz);
    double phi = 6.28318530717958647692*(x->operator()());
    
    return lorentzvector<double>(pt*std::cos(phi), pt*std::sin(phi), pz, E);
  }

  double __rambo_helper_weight(unsigned int n, double s)
  {
#define __DEB_8PI   25.13274122871834590768
#define __DEB_16PI2 157.91367041742973790108
    
    static const double fact[] = {
      1.000000000000000000000, 1.000000000000000000000,
      1.000000000000000000000, 2.000000000000000000000,
      3.464101615137754587050, 5.241482788417793214280,
      7.325683002969412729490, 9.711867496604076595340,
      12.39635026120588892703, 15.37624008321082805550,
      
      18.64921022211253121271, 22.21334732883035904420,
      26.06704940665039454817, 30.20895473052848875326,
      34.63789090520187258023, 39.35283742689448285515,
      44.35289755618292354658, 49.63727677424801864857,
      55.20526599889455311234, 61.05622831083046240725,
      
      67.189588314900381553350, 73.604823510775700955150,
      80.301457218028968401580, 87.279052719104399855340,
      94.537208367699191038680, 102.07555347055500771508,
      109.89374479486319958781, 117.99146358623666580202,
      126.36841300677002588944, 135.02431592135521673949,
      
      143.95891297472152603605, 153.17196091274811705250,
      162.66323111025862674191, 172.43250827433702043364,
      182.47958929763210591531, 192.80428224046477961708,
      203.40640542405709375678, 214.28578662004717056856,
      225.44226232377745474056, 236.87567710075244487879,
    };
    
    return std::pow(s/(__DEB_16PI2*fact[n]), (int) (n-2))/__DEB_8PI;
  }
}

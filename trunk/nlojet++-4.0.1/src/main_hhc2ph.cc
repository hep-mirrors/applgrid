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
#include <cstdlib>
#include "random.h"

#include "bits/hhc2ph-event.h"
#include "bits/hhc2ph-phasespace.h"
#include "bits/hhc2ph-process.h"
#include "bits/hhc2ph-jetfunc.h"
#include "bits/hhc2ph-integral.h"

#include "hhc2ph1jet.h"
#include "hhc2ph2jet.h"

#include "nlojet++.h"


//----- the used namespaces -----
using namespace nlo;



void main_module_hhc2ph(const main_input& in) 
{
  typedef void (*ifun_type)(unsigned int&, unsigned int&, unsigned int&);
  typedef void (*psfun_type)(phasespace_hhc2ph *, double&, double&, double&);

  //--- set up the input parameters ---
  unsigned int nu = 2U, nd = 3U, njet = 2U;
  (*((ifun_type) in.find_symbol("inputfunc")))(njet, nu, nd);

  //--- initialize the random number generator ---
  rng_mt19937 rng((unsigned long int) std::rand());
    
  //--- create phase space ---
  double s = 14000*14000, q2min = 118.5*118.5, q2max = 121.5*121.5;  // GeV^2
  phasespace_hhc2ph *ps = 0;
  bool psdealloc = false;
  
  (*((psfun_type) in.find_symbol("psinput")))(ps, s, q2min, q2max);
  
  if(ps == 0) {
    basic_phasespace_hhc2ph *__ps = new basic_phasespace_hhc2ph(&rng, s);
	__ps -> important_sampling(2.0, in.alpha);
	__ps -> phasespace_cuts(q2min, q2max);
	ps = __ps;
    psdealloc = true;
  }
  
  //--- create the process ---
  process_hhc2ph *amp = 0;
  
  switch(njet) {
  case 1: amp = new hhc2ph1jet(rng, in.mchel, nu, nd, in.alpha); break;
  case 2: amp = new hhc2ph2jet(rng, in.mchel, nu, nd, in.alpha); break;
  default: throw "njet > 2 jet calculations aren't implemented"; break;
  }

  //--- create the user function ---
  switch(njet) {
  case 1: in.outfile.replace(in.outfile.find("@NJET@"), 6, "1"); break;
  case 2: in.outfile.replace(in.outfile.find("@NJET@"), 6, "2"); break;
  }

  user_base_hhc2ph *jet = (*((user_base_hhc2ph * (*)()) in.find_symbol("userfunc")))();
  jet -> phys_output(in.outfile, in.nsave, in.txtout);

  //--- create integral ---
  integral_hhc2ph integral(ps, amp, jet);
  
  //--- start the calculation ---
  switch (in.contr) {
  case 0: integral.calculate_born(in.nevent);          break;
  case 1: integral.calculate_nlo (in.nevent, in.time); break;
  case 2: integral.calculate     (in.nevent, in.time); break;
  }
  
  if(amp) delete amp;
  if(psdealloc && ps) delete ps;  
}

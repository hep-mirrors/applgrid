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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#include <cstdlib>
#include "random.h"

#include "bits/epa-event.h"
#include "bits/epa-phasespace.h"
#include "bits/epa-process.h"
#include "bits/epa-jetfunc.h"
#include "bits/epa-integral.h"

#include "epa3jet.h"
#include "epa4jet.h"
#include "epa5jet.h"

#include "nlojet++.h"


//----- the used namespaces -----
using namespace nlo;



void main_module_epa(const main_input& in) 
{
  typedef void (*ifun_type)(unsigned int&, unsigned int&, unsigned int&);
  typedef void (*psfun_type)(phasespace_epa *, double&);
  
  //--- set up the input parameters ---
  unsigned int nu = 2U, nd = 3U, njet = 3U;
  (*((ifun_type) in.find_symbol("inputfunc")))(njet, nu, nd);
    
  //--- initialize the random number generator ---
  rng_mt19937 rng((unsigned long int) std::rand());
    
  //--- create phase space ---
  double s = 8315.068969;
  phasespace_epa *ps = 0;
  bool psdealloc = false;
  
  (*((psfun_type) in.find_symbol("psinput")))(ps, s);
  
  if(ps == 0) {
    ps = new basic_phasespace_epa(&rng, s, 2.0, in.alpha);
    psdealloc = true;
  }
  
  //--- create the process ---
  process_epa *amp = 0;
    
  switch(njet) {
  case 3: amp = new epa3jet(rng, in.mchel, nu, nd, in.alpha); break;
  case 4: amp = new epa4jet(rng, in.mchel, nu, nd, in.alpha); break;
  case 5: amp = new epa5jet(rng, in.mchel, nu, nd, in.alpha); break;
  default: throw "n < 3 or n > 5 jet calculations aren't implemented"; break;
  }
    
  //--- create the user function ---
  switch(njet) {
  case 3: in.outfile.replace(in.outfile.find("@NJET@"), 6, "3"); break;
  case 4: in.outfile.replace(in.outfile.find("@NJET@"), 6, "4"); break;
  case 5: in.outfile.replace(in.outfile.find("@NJET@"), 6, "5"); break;
  }
    
  user_base_epa *jet = (*((user_base_epa * (*)()) in.find_symbol("userfunc")))();
  jet -> phys_output(in.outfile, in.nsave, in.txtout);
    
  //--- create integral ---
  integral_epa integral(ps, amp, jet);
    
  //--- start the calculation ---
  switch (in.contr) {
  case 0: integral.calculate_born(in.nevent);          break;
  case 1: integral.calculate_nlo (in.nevent, in.time); break;
  case 2: integral.calculate     (in.nevent, in.time); break;
  }
    
  if(amp) delete amp;
  if(psdealloc && ps) delete ps;  
}

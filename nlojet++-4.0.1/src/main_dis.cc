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

#include "bits/dis-event.h"
#include "bits/dis-phasespace.h"
#include "bits/dis-process.h"
#include "bits/dis-jetfunc.h"
#include "bits/dis-integral.h"

#include "dis2jet.h"
#include "dis3jet.h"
#include "dis4jet.h"

#include "nlojet++.h"


//----- the used namespaces -----
using namespace nlo;




void main_module_dis(const main_input& in) 
{
  typedef void (*ifun_type)(unsigned int&, unsigned int&, unsigned int&);
  typedef void (*psfun_type)(phasespace_dis *, double&, double&, double&, double&, double&, double&, double&, double&);
  
  //--- set up the input parameters ---
  unsigned int nu = 2U, nd = 3U, njet = 2U;
  (*((ifun_type) in.find_symbol("inputfunc")))(njet, nu, nd);

  //--- initialize the random number generator ---
  rng_mt19937 rng((unsigned long int) std::rand());
    
  //--- create phase space ---
  double el = 27.439, eh = 820.0;
  double q2min = 5.0, q2max = 5000.0, xmin = 0.0, xmax = 1.0, ymin = 0.2, ymax = 0.6;
  
  phasespace_dis *ps = 0;
  bool psdealloc = false;
  
  (*((psfun_type) in.find_symbol("psinput")))(ps, el, eh, q2min, q2max, xmin, xmax, ymin, ymax);
  
  if(ps == 0) {
    basic_phasespace_dis *__ps = new basic_phasespace_dis(&rng, el, eh);
	__ps -> important_sampling(2.0, in.alpha);
	__ps -> phasespace_cuts(q2min, q2max, xmin, xmax, ymin, ymax);
	ps = __ps;
    psdealloc = true;
  }
  
  //--- create the process ---
  process_dis *amp = 0;
  
  switch(njet) {
  case 2: amp = new dis2jet(rng, in.mchel, nu, nd, in.alpha); break;
  case 3: amp = new dis3jet(rng, in.mchel, nu, nd, in.alpha); break;
  case 4: amp = new dis4jet(rng, in.mchel, nu, nd, in.alpha); break;
  default: throw "n < 2 or n > 4 jet calculations aren't implemented"; break;
  }

  //--- create the user function ---
  switch(njet) {
  case 2: in.outfile.replace(in.outfile.find("@NJET@"), 6, "2"); break;
  case 3: in.outfile.replace(in.outfile.find("@NJET@"), 6, "3"); break;
  case 4: in.outfile.replace(in.outfile.find("@NJET@"), 6, "4"); break;
  }
  
  user_base_dis *jet = (*((user_base_dis * (*)()) in.find_symbol("userfunc")))();
  jet -> phys_output(in.outfile, in.nsave, in.txtout);
  
  //--- create integral ---
  integral_dis integral(ps, amp, jet);
  
  //--- start the calculation ---
  switch (in.contr) {
  case 0: integral.calculate_born(in.nevent);          break;
  case 1: integral.calculate_nlo (in.nevent, in.time); break;
  case 2: integral.calculate     (in.nevent, in.time); break;
  }
  
  if(amp) delete amp;
  if(psdealloc && ps) delete ps;  
}

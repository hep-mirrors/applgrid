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

#include "bits/photo-event.h"
#include "bits/photo-phasespace.h"
#include "bits/photo-process.h"
#include "bits/photo-jetfunc.h"
#include "bits/photo-integral.h"

#include "photo2jet.h"
#include "photo3jet.h"
#include "photo4jet.h"

#include "nlojet++.h"


//----- the used namespaces -----
using namespace nlo;
using namespace std;


void main_module_photo(const main_input& in) 
{
  typedef void (*ifun_type)(unsigned int&, unsigned int&, unsigned int&);
  typedef void (*psfun_type)(phasespace_photo *, double&, double&);
  
  //--- set up the input parameters ---
  unsigned int nu = 2U, nd = 3U, njet = 2U;
  (*((ifun_type) in.find_symbol("inputfunc")))(njet, nu, nd);
  
  //--- initialize the random number generator ---
  rng_mt19937 rng((unsigned long int) std::rand());
  
  //--- create phase space ---
  phasespace_photo *ps = 0;
  bool psdealloc = false;
  double el = 27.439, eh = 820.0;

  (*((psfun_type) in.find_symbol("psinput")))(ps, el, eh);
  
  if(ps == 0) {
    ps = new basic_phasespace_photo(&rng, el, eh, 2.0, in.alpha);
    psdealloc = true;
  }
  
  //--- create the process ---
  process_photo *amp = 0;
  
  switch(njet) {
  case 1:
  case 2: amp = new photo2jet(rng, in.mchel, nu, nd, in.alpha); break;
  case 3: amp = new photo3jet(rng, in.mchel, nu, nd, in.alpha); break;
  case 4: amp = new photo4jet(rng, in.mchel, nu, nd, in.alpha); break;
  default: throw "n > 4 jet calculations aren't implemented"; break;
  }

  //--- create the user function ---
  switch(njet) {
  case 1: in.outfile.replace(in.outfile.find("@NJET@"), 6, "1"); break;
  case 2: in.outfile.replace(in.outfile.find("@NJET@"), 6, "2"); break;
  case 3: in.outfile.replace(in.outfile.find("@NJET@"), 6, "3"); break;
  case 4: in.outfile.replace(in.outfile.find("@NJET@"), 6, "4"); break;
  }

  user_base_photo *jet = (*((user_base_photo * (*)()) in.find_symbol("userfunc")))();
  jet -> phys_output(in.outfile, in.nsave, in.txtout);

  //--- create integral ---
  integral_photo integral(ps, amp, jet);
  

  //--- start the calculation ---
  switch (in.contr) {
  case 0: integral.calculate_born(in.nevent);          break;
  case 1: integral.calculate_nlo (in.nevent, in.time); break;
  case 2: integral.calculate     (in.nevent, in.time); break;
  }

  
  if(amp) delete amp;
  if(psdealloc && ps) delete ps;  
}

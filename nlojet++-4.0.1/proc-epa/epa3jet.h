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
#ifndef __NLO_EPA3JET_H__
#define __NLO_EPA3JET_H__ 1

#include <bits/hep-rng.h>
#include <bits/nlo-innerprod.h>
#include <bits/nlo-split.h>
#include <bits/proc-epa.h>
#include <bits/epa-process.h>




namespace nlo {
  

  class epa3jet 
    : public process_epa, private _epa_jet_base 
  {
  public:
    //   constructor
    explicit epa3jet(const random_generator&, bool=false, unsigned int=2, unsigned int=3, double=1.0);
    
    //   destructor
    ~epa3jet();
    
    //  born contributions
    void born_term(const event_type&, weight_type&);
    
    //  real contributions
    void real_term(const event_type&, weight_type&);
    
    //  finite contributions (1-loop, ...)
    void fini_term(const event_type&, weight_type&);
    
    //  dipole contributions
    void dipole_term(const event_type&, const event_type&, int, int, int, weight_type&);
    
  private:
    //   inner products
    innerprod<lorentzvector<double> > _M_ip;
    
    //   amplitudes
    ampq2g1l2 *_M_q2g1;
    ampq2g2l2 *_M_q2g2;
    ampq4l2   *_M_q4;
    
    //  splitting function
    splitfff  _M_split;
    
    //  Monte Carlo helicity sum
    bool _M_mchel;

    //       dipole contributions
    typedef double (epa3jet::*dipole_func)(int, int);
    static dipole_func _S_dipole[6];
    
    //       dipole functions
    double _M_d12(int, int);
    double _M_d13(int, int);
    double _M_d14(int, int);
    double _M_d23(int, int);
    double _M_d34(int, int);
  };
}     //  namespace nlo  

#endif


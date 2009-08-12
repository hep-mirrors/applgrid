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
#ifndef __NLO_HHC2PH1JET_H__
#define __NLo_HHC2PH1JET_H__ 1


#include <bits/hep-rng.h>
#include <bits/nlo-innerprod.h>
#include <bits/nlo-split.h>
#include <bits/proc-hhc2ph.h>
#include <bits/hhc2ph-process.h>




namespace nlo {
  

  class hhc2ph1jet
    : public process_hhc2ph, private _hhc2ph_jet_base 
  {
  public:
    //   constructor
    explicit hhc2ph1jet(const random_generator&, bool=false, unsigned int=2, unsigned int=3, double=1.0);
    
    //   destructor
    ~hhc2ph1jet();
    
    //      born contribution
    void born_term(const event_type&, weight_type&);
    
    //      real contribution
    void real_term(const event_type&, weight_type&);
    
    //  finite contributions (1-loop,...)
    void fini_term(double, double, double, double, const event_type&, weight_type *);

    //         dipole contributions
    void dipole_term(const event_type&, const event_type&, int, int, int, weight_type&);
    
  private:
    //   inner products
    innerprod<lorentzvector<double> > _M_ip;
    
    //   amplitudes
    ampq2g1p2  *_M_q2g1p2;
    ampq2g2p2  *_M_q2g2p2;
    ampq4p2    *_M_q4p2;
    
    //  splitting functions
    splitffi _M_sffi;
    splitiff _M_siff;
    splitifi _M_sifi;
    
    split_ini<lorentzvector<double> > *_M_sini;

    //  Monte Carlo helicity sum
    bool _M_mchel;

    //   Dipole contributions
    typedef void (hhc2ph1jet::*dipole_func)(int, int, weight_type&);
    static dipole_func _S_dipole[5];

    //  Dipole functions
    void _M_di1(int, int, weight_type&); 
    void _M_di2(int, int, weight_type&); 
    void _M_d01(int, int, weight_type&); 
    void _M_d02(int, int, weight_type&); 
    void _M_d12(int, int, weight_type&); 
  };
}     //  namespace nlo   

#endif


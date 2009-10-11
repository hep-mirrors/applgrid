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
#ifndef __NLO_PHOTO3JET_H__
#define __NLO_PHOTO3JET_H__ 1


#include <bits/hep-rng.h>
#include <bits/nlo-innerprod.h>
#include <bits/nlo-split.h>
#include <bits/proc-photo.h>
#include <bits/photo-process.h>



namespace nlo {
  

  class photo3jet
    : public process_photo, private _photo_jet_base 
  {
  public:
    //   constructor
    explicit photo3jet(const random_generator&, bool=false, unsigned int=2, unsigned int=3, double=1.0);
    
    //   destructor
    ~photo3jet();
    
    //      born contribution
    void born_term(const event_type&, weight_type&);
    
    //      real contribution
    void real_term(const event_type&, weight_type&);
    
    //  finite contributions (1-loop,...)
    void fini_term(double, double, double, double, const event_type&, weight_type *);

    //         dipole contributions
    void dipole_term(const event_type&, const event_type&, int, int, int, weight_type&);
    
    //   select the dipoles 
    bool dipole_index(int, int, int);

  private:
    //   inner products
    innerprod<lorentzvector<double> > _M_ip;
    
    //   amplitudes
    ampq2g2p1  *_M_q2g2p1;
    ampq4p1    *_M_q4p1;
    ampq2g3p1  *_M_q2g3p1;
    ampq4g1p1  *_M_q4g1p1;
    ampq2g3    *_M_q2g3;
    ampq4g1    *_M_q4g1;
    
    //  splitting functions
    splitffi _M_sffi;
    splitfff _M_sfff;
    splitiff _M_siff;
    splitifi _M_sifi;
    
    split_fin<lorentzvector<double> > *_M_sfin;
    split_ini<lorentzvector<double> > *_M_sini;

    //  Monte Carlo helicity sum
    bool _M_mchel;

    //   Dipole contributions
    typedef void (photo3jet::*dipole_func)(int, int, weight_type&);
    static dipole_func _S_dipole[14];

    //  Dipole functions
    void _M_di1(int, int, weight_type&); 
    void _M_di2(int, int, weight_type&); 
    void _M_di3(int, int, weight_type&); 
    void _M_di4(int, int, weight_type&); 
    void _M_d01(int, int, weight_type&); 
    void _M_d02(int, int, weight_type&); 
    void _M_d03(int, int, weight_type&); 
    void _M_d04(int, int, weight_type&); 
    void _M_d12(int, int, weight_type&); 
    void _M_d13(int, int, weight_type&); 
    void _M_d14(int, int, weight_type&); 
    void _M_d23(int, int, weight_type&); 
    void _M_d24(int, int, weight_type&); 
    void _M_d34(int, int, weight_type&); 
  };
}     //  namespace nlo   

#endif


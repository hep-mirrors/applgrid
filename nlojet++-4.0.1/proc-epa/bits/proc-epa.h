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
#ifndef __NLO_PROC_EPA_H__
#define __NLO_PROC_EPA_H__ 1

#include <complex>
#include <utility>


namespace nlo {

  //    The amplitude classes
  class ampq2g1l2;
  class ampq2g2l2;
  class ampq2g3l2;
  class ampq4l2;
  class ampq4g1l2;

  
  //    Helper class for implementing the e+e- n-jets processes
  class _epa_jet_base 
  {
  protected:
    //   constructor
    explicit _epa_jet_base(unsigned int = 2U, unsigned int = 3U);
    
    //  types
    typedef std::pair<double, std::complex<double> > _Pair;
    
    //    Three parton contributions
    double amp_tree (ampq2g1l2*);
    double amp_tree_mch(ampq2g1l2*);
    
    double amp_1loop(ampq2g1l2*, double);
    double amp_1loop_mch(ampq2g1l2*, double);

    _Pair amp_cc(ampq2g1l2*, int, int, int, int, int);
    
    //    Four parton contributions
    double amp_tree(ampq2g2l2*, ampq4l2*);
    double amp_tree_mch(ampq2g2l2*, ampq4l2*);
    
    double amp_1loop(ampq2g2l2*, ampq4l2*, double);
    double amp_1loop_mch(ampq2g2l2*, ampq4l2*, double);

    _Pair amp_cc(ampq2g2l2*, int, int, int, int, int, int);
    _Pair amp_cc(ampq4l2*,   int, int, int, int, int, int);
    
    //    Five parton contributions
    double amp_tree(ampq2g3l2*, ampq4g1l2*);
    double amp_tree_mch(ampq2g3l2*, ampq4g1l2*);
    
    //  static members
    static double _M_Dijk(const _Pair& __v, const _Pair& __cc)
      { return __v.first*__cc.first + 2.0*real(__v.second*__cc.second);}
    
    //   data members
    unsigned int _M_nf;
    
  private:
    //   private data members
    double _M_xcharge;
  };
}   //   namespace nlo

#endif

